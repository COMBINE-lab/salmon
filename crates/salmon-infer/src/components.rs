//! Connected-component EM/VBEM.
//!
//! The class ↔ transcript incidence graph falls apart into many independent
//! connected components: a transcript's update reads only the abundances of
//! transcripts it shares a class with, so no arithmetic ever crosses a
//! component boundary. This module lays the eq-classes out component by
//! component and iterates each component on its own.
//!
//! **Layout.** Every transcript some class (with a nonzero count) references is
//! given a dense id; components are contiguous runs of dense ids, ordered by
//! their smallest original transcript id, and within a component the
//! transcripts keep ascending original order. A component's classes keep their
//! original relative order and store component-local labels (`u32`) beside
//! their combined weights (`f64`): 12 bytes per incidence. Single-transcript
//! classes are not iterated at all: their counts are summed (in class order)
//! into a per-transcript constant that starts each M-step's output.
//!
//! **Arithmetic.** Each transcript's new abundance is its singleton constant
//! plus its multi-transcript contributions in class order — a data-defined
//! summation order, so the result never depends on the thread count, and
//! scheduling components onto threads cannot change it.
//!
//! **VBEM.** The VBEM expectation is `exp(ψ(α_t + prior_t) − ψ(Σ_j α_j + prior_j))`.
//! The second term is common to every transcript and so cancels in each
//! class's ratio `e_t w_t / Σ_j e_j w_j`; it only moved the rounding. It is
//! omitted here (`e_t = exp(ψ(α_t + prior_t))`), which is what makes VBEM
//! component-separable at all — with it, every component would have to wait on
//! a global sum each iteration. Over the range VBEM visits (`α + prior ≥ prior
//! = 0.01` for any referenced transcript under the default prior), `exp(ψ(·))`
//! neither overflows nor underflows. It is evaluated by [`exp_digamma`], a
//! branch-free form whose loop vectorizes (and is dispatched to AVX2/AVX-512
//! at run time on x86-64, with bit-identical results).
//!
//! **Schedule.** The stopping rule is salmon's: after `min_iter`, stop at the
//! first iteration where the relative-change criterion holds over all
//! transcripts at once. Since no arithmetic crosses components, that loop is
//! run component-major in tiles (see [`run_tiled`]): each component advances
//! up to [`TILE`] iterations while its data is in cache, and the stopping
//! iteration is reproduced exactly. A component large enough to be worth
//! splitting is stepped in parallel by gathering ([`Pull`]), which adds the
//! same terms in the same order as the class-major loop.
//!
//! VBEM by default instead stops each component on its own ([`run_projected`]):
//! once salmon's criterion holds over the component and its projected
//! remaining movement is small. The global rule keeps every component
//! iterating until the slowest settles, which is most of VBEM's work on a real
//! sample; the criterion alone, applied per component, stops slowly mixing
//! isoform families (*TTN*, *NEB*) far from their fixpoint. Plain EM keeps the
//! global rule.

use rayon::prelude::*;

use crate::packed::{PackedEqClasses, DIGAMMA_MIN, MIN_EQ_CLASS_WEIGHT};
use crate::EmStop;

/// The component-ordered eq-classes.
pub(crate) struct ComponentPlan {
    /// dense id → original transcript id
    txps: Vec<u32>,
    /// component boundaries in dense ids, `ncomp + 1` long
    comp_txp: Vec<usize>,
    /// component boundaries in (multi-transcript) classes, `ncomp + 1` long
    comp_cls: Vec<usize>,
    /// class boundaries in incidences, `nclasses + 1` long
    cls_start: Vec<usize>,
    /// component-local dense label of every incidence
    labels: Vec<u32>,
    /// combined weight of every incidence
    weights: Vec<f64>,
    /// count of every class
    counts: Vec<f64>,
    /// summed single-transcript class counts, per dense id
    singleton: Vec<f64>,
    /// per component, the transposed incidences used to step it in parallel
    /// (built only for large components of a parallel run)
    pull: Vec<Option<Pull>>,
}

/// A component's incidences transposed: for each of its transcripts, the
/// classes containing it in ascending class order, with the weight there.
///
/// Stepping a component this way — every class's scale in parallel, then
/// every transcript's sum in parallel — adds exactly the same terms in exactly
/// the same order as the class-major loop, so the two are bit-identical and
/// choosing between them is a scheduling decision.
struct Pull {
    start: Vec<usize>,
    cls: Vec<u32>,
    w: Vec<f64>,
}

/// Components with at least this many incidences are stepped in parallel
/// (when the run is parallel).
const PULL_MIN_INCIDENCES: usize = 1 << 16;

impl ComponentPlan {
    /// Lay out `p` (with per-class `counts`) by component. `parallel` asks for
    /// the transposed incidences that let a large component be stepped in
    /// parallel; it changes scheduling only.
    pub(crate) fn new(p: &PackedEqClasses, counts: &[u64], parallel: bool) -> Self {
        let pull_min = if parallel && rayon::current_num_threads() > 1 {
            PULL_MIN_INCIDENCES
        } else {
            usize::MAX
        };
        Self::with_pull_min(p, counts, pull_min)
    }

    fn with_pull_min(p: &PackedEqClasses, counts: &[u64], pull_min: usize) -> Self {
        let n = p.num_txps;
        let nclasses = p.num_classes();
        // Union-find over transcripts, joining each multi class's members.
        let mut parent: Vec<u32> = (0..n as u32).collect();
        fn find(parent: &mut [u32], mut x: u32) -> u32 {
            while parent[x as usize] != x {
                let gp = parent[parent[x as usize] as usize];
                parent[x as usize] = gp;
                x = gp;
            }
            x
        }
        let mut active = vec![false; n];
        let mut singleton_full = vec![0.0f64; n];
        for ci in 0..nclasses {
            if counts[ci] == 0 {
                continue;
            }
            let (tids, _) = p.class(ci);
            match tids.len() {
                0 => {}
                1 => {
                    active[tids[0] as usize] = true;
                    singleton_full[tids[0] as usize] += counts[ci] as f64;
                }
                _ => {
                    let mut r0 = find(&mut parent, tids[0]);
                    active[tids[0] as usize] = true;
                    for &t in &tids[1..] {
                        active[t as usize] = true;
                        let r = find(&mut parent, t);
                        if r != r0 {
                            // Union toward the smaller id; any rule works, the
                            // component order below does not depend on it.
                            let (lo, hi) = if r < r0 { (r, r0) } else { (r0, r) };
                            parent[hi as usize] = lo;
                            r0 = lo;
                        }
                    }
                }
            }
        }
        // Components numbered by their smallest transcript id.
        const NONE: u32 = u32::MAX;
        let mut comp_of_root = vec![NONE; n];
        let mut comp_of_txp = vec![NONE; n];
        let mut ncomp = 0u32;
        for t in 0..n {
            if !active[t] {
                continue;
            }
            let r = find(&mut parent, t as u32) as usize;
            if comp_of_root[r] == NONE {
                comp_of_root[r] = ncomp;
                ncomp += 1;
            }
            comp_of_txp[t] = comp_of_root[r];
        }
        let ncomp = ncomp as usize;
        // Counting sort of transcripts by component (stable: ascending tid).
        let mut comp_txp = vec![0usize; ncomp + 1];
        for &c in &comp_of_txp {
            if c != NONE {
                comp_txp[c as usize + 1] += 1;
            }
        }
        for c in 0..ncomp {
            comp_txp[c + 1] += comp_txp[c];
        }
        let nactive = comp_txp[ncomp];
        let mut txps = vec![0u32; nactive];
        let mut dense_of = vec![0u32; n];
        let mut cursor = comp_txp.clone();
        for t in 0..n {
            let c = comp_of_txp[t];
            if c != NONE {
                let d = cursor[c as usize];
                cursor[c as usize] += 1;
                txps[d] = t as u32;
                dense_of[t] = d as u32;
            }
        }
        let singleton = txps.iter().map(|&t| singleton_full[t as usize]).collect();
        // Counting sort of multi classes by component (stable: class order).
        let mut comp_cls = vec![0usize; ncomp + 1];
        let mut ninc = 0usize;
        for ci in 0..nclasses {
            let (tids, _) = p.class(ci);
            if counts[ci] != 0 && tids.len() > 1 {
                comp_cls[comp_of_txp[tids[0] as usize] as usize + 1] += 1;
                ninc += tids.len();
            }
        }
        for c in 0..ncomp {
            comp_cls[c + 1] += comp_cls[c];
        }
        let nmulti = comp_cls[ncomp];
        let mut order = vec![0usize; nmulti];
        let mut cursor = comp_cls.clone();
        for ci in 0..nclasses {
            let (tids, _) = p.class(ci);
            if counts[ci] != 0 && tids.len() > 1 {
                let c = comp_of_txp[tids[0] as usize] as usize;
                order[cursor[c]] = ci;
                cursor[c] += 1;
            }
        }
        let mut cls_start = Vec::with_capacity(nmulti + 1);
        cls_start.push(0);
        let mut labels = Vec::with_capacity(ninc);
        let mut weights = Vec::with_capacity(ninc);
        let mut ccounts = Vec::with_capacity(nmulti);
        for c in 0..ncomp {
            let base = comp_txp[c] as u32;
            for &ci in &order[comp_cls[c]..comp_cls[c + 1]] {
                let (tids, ws) = p.class(ci);
                labels.extend(tids.iter().map(|&t| dense_of[t as usize] - base));
                weights.extend_from_slice(ws);
                ccounts.push(counts[ci] as f64);
                cls_start.push(labels.len());
            }
        }
        let mut plan = Self {
            txps,
            comp_txp,
            comp_cls,
            cls_start,
            labels,
            weights,
            counts: ccounts,
            singleton,
            pull: Vec::new(),
        };
        plan.pull = (0..ncomp)
            .map(|c| {
                // Push and pull are bit-identical, so this may depend on the
                // pool: one thread gains nothing from the transposed copy.
                (plan.comp_incidences(c) >= pull_min).then(|| plan.transpose(c))
            })
            .collect();
        plan
    }

    fn transpose(&self, c: usize) -> Pull {
        let v = self.view(c);
        let nt = v.singleton.len();
        let mut start = vec![0usize; nt + 1];
        for &t in v.labels {
            start[t as usize + 1] += 1;
        }
        for t in 0..nt {
            start[t + 1] += start[t];
        }
        let mut cursor = start.clone();
        let mut cls = vec![0u32; v.labels.len()];
        let mut w = vec![0.0f64; v.labels.len()];
        for k in 0..v.counts.len() {
            let (lo, hi) = (v.cls_start[k] - v.base, v.cls_start[k + 1] - v.base);
            for i in lo..hi {
                let t = v.labels[i] as usize;
                cls[cursor[t]] = k as u32;
                w[cursor[t]] = v.weights[i];
                cursor[t] += 1;
            }
        }
        Pull { start, cls, w }
    }

    pub(crate) fn num_components(&self) -> usize {
        self.comp_txp.len() - 1
    }

    /// Incidences of component `c`.
    fn comp_incidences(&self, c: usize) -> usize {
        self.cls_start[self.comp_cls[c + 1]] - self.cls_start[self.comp_cls[c]]
    }

    fn view(&self, c: usize) -> CompView<'_> {
        let (c0, c1) = (self.comp_cls[c], self.comp_cls[c + 1]);
        let (t0, t1) = (self.comp_txp[c], self.comp_txp[c + 1]);
        let i0 = self.cls_start[c0];
        CompView {
            cls_start: &self.cls_start[c0..=c1],
            base: i0,
            labels: &self.labels[i0..self.cls_start[c1]],
            weights: &self.weights[i0..self.cls_start[c1]],
            counts: &self.counts[c0..c1],
            singleton: &self.singleton[t0..t1],
        }
    }
}

/// One component's slice of the plan.
struct CompView<'a> {
    cls_start: &'a [usize],
    base: usize,
    labels: &'a [u32],
    weights: &'a [f64],
    counts: &'a [f64],
    singleton: &'a [f64],
}

/// Per-component mutable state: abundances, the next M-step's output, and
/// (VBEM) the expectations.
struct CompState<'a> {
    alpha: &'a mut [f64],
    next: &'a mut [f64],
    prior: &'a [f64],
    et: &'a mut [f64],
    /// whether the criterion over this component last failed (tiling's
    /// witness choice; scheduling only)
    changed: bool,
    /// per-class scale scratch for a parallel (pull) step
    scale: Vec<f64>,
}

/// Smallest prior for which every VBEM expectation is a positive normal:
/// `α ≥ 0`, so `α + prior ≥ 0.002` and [`exp_digamma`] is at least
/// `exp(ψ(0.002)) ≈ 1e-217` (it flushes to `0` only below about `0.0014`).
const POSITIVE_PRIOR_MIN: f64 = 0.002;

/// `e^z` for `z ≤ 0`, branch-free (so a loop over it vectorizes): Cody–Waite
/// reduction by `ln 2`, a degree-13 Taylor polynomial on `|r| ≤ ln2/2`
/// (truncation below 4e-18), and the power of two set in the exponent bits.
/// Results below the smallest normal double (`z < −708`) are `0`.
#[inline(always)]
fn exp_nonpos(z: f64) -> f64 {
    // ln 2 split as fdlibm does: the high part has 32 trailing zero bits, so
    // `k·LN2_HI` is exact for every `k` used here.
    const LN2_HI: f64 = f64::from_bits(0x3fe6_2e42_fee0_0000);
    const LN2_LO: f64 = f64::from_bits(0x3dea_39ef_3579_3c76);
    // 1.5·2⁵²: adding it rounds to an integer, held in the low mantissa bits.
    const SHIFTER: f64 = 6_755_399_441_055_744.0;
    let zc = if z < -708.0 { -708.0 } else { z };
    let t = zc * std::f64::consts::LOG2_E + SHIFTER;
    let k = t - SHIFTER;
    let r = (zc - k * LN2_HI) - k * LN2_LO;
    let (r2, r4) = (r * r, r * r * r * r);
    let pair = |a: f64, b: f64| a + r * b;
    let q0 = pair(1.0, 1.0) + r2 * pair(1.0 / 2.0, 1.0 / 6.0);
    let q1 = pair(1.0 / 24.0, 1.0 / 120.0) + r2 * pair(1.0 / 720.0, 1.0 / 5040.0);
    let q2 =
        pair(1.0 / 40_320.0, 1.0 / 362_880.0) + r2 * pair(1.0 / 3_628_800.0, 1.0 / 39_916_800.0);
    let q3 = pair(1.0 / 479_001_600.0, 1.0 / 6_227_020_800.0);
    let poly = (q0 + r4 * q1) + (r4 * r4) * (q2 + r4 * q3);
    let ki = t.to_bits() as i64 - SHIFTER.to_bits() as i64;
    let v = poly * f64::from_bits(((ki + 1023) as u64) << 52);
    if z < -708.0 {
        0.0
    } else {
        v
    }
}

/// `exp(ψ(x))`, the VBEM expectation, for `x > 0`; branch-free, without `ln`
/// and with one inline exponential.
///
/// `ψ(x) = ln y + u(y) − S`, where `y = x + 10` and `S = Σ_{j<10} 1/(x + j)`
/// when `x < 10` (the recurrence; `y = x`, `S = 0` otherwise), and `u` is the
/// asymptotic series of `ψ(y) − ln y` with Bernoulli terms through `y⁻¹²`
/// (truncation below 1e-15 at `y = 10`). So `exp(ψ(x)) = y·e^{u − S}`, with
/// `S` formed as one ratio of pairwise terms. Against 30-digit mpmath on
/// `[1e-2, 1e8]` the relative error is at most 4.5e-14 (below 0.1), 6.5e-15
/// (0.1–1), 1.7e-15 (1–10) and 9.3e-16 (≥ 10) — no worse than `statrs`'s
/// `digamma(x).exp()` in each range and within ~2× the conditioning floor
/// `|ψ(x)|·2⁻⁵²`. Where `exp(ψ(x))` is subnormal (`x` below about `0.0014`,
/// reachable only with `--vbPrior` < 0.01) it returns `0`.
#[inline(always)]
fn exp_digamma(x: f64) -> f64 {
    let small = x < 10.0;
    let (mut num, mut den, mut j) = (0.0f64, 1.0f64, 0.0f64);
    while j < 10.0 {
        let a = x + j;
        let b = a + 1.0;
        num = num * (a * b) + (a + b) * den;
        den *= a * b;
        j += 2.0;
    }
    let y = if small { x + 10.0 } else { x };
    let s = if small { num / den } else { 0.0 };
    let r = 1.0 / y;
    let r2 = r * r;
    let u = -0.5 * r
        - r2 * (1.0 / 12.0
            - r2 * (1.0 / 120.0
                - r2 * (1.0 / 252.0
                    - r2 * (1.0 / 240.0 - r2 * (1.0 / 132.0 - r2 * (691.0 / 32_760.0))))));
    y * exp_nonpos(u - s)
}

/// VBEM expectations `exp(ψ(α + prior))`, `0` where `α + prior` is not above
/// [`DIGAMMA_MIN`]. Branch-free, so the loop vectorizes.
#[inline(always)]
fn fill_et_body(et: &mut [f64], alpha: &[f64], prior: &[f64]) {
    for ((et, &a), &pr) in et.iter_mut().zip(alpha).zip(prior) {
        let ap = a + pr;
        let e = exp_digamma(ap);
        *et = if ap > DIGAMMA_MIN { e } else { 0.0 };
    }
}

/// [`fill_et_body`] compiled for AVX-512F. The same IEEE operations in the
/// same order (Rust never contracts into FMA), so the same bits.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f")]
fn fill_et_avx512(et: &mut [f64], alpha: &[f64], prior: &[f64]) {
    fill_et_body(et, alpha, prior)
}

/// [`fill_et_body`] compiled for AVX2; see [`fill_et_avx512`].
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
fn fill_et_avx2(et: &mut [f64], alpha: &[f64], prior: &[f64]) {
    fill_et_body(et, alpha, prior)
}

/// The expectations for one component (or one chunk of it), on the widest
/// vector unit available.
#[allow(unsafe_code)]
fn fill_et(et: &mut [f64], alpha: &[f64], prior: &[f64]) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::arch::is_x86_feature_detected!("avx512f") {
            // SAFETY: the CPU supports AVX-512F (checked just above).
            return unsafe { fill_et_avx512(et, alpha, prior) };
        }
        if std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: the CPU supports AVX2 (checked just above).
            return unsafe { fill_et_avx2(et, alpha, prior) };
        }
    }
    fill_et_body(et, alpha, prior)
}

/// A class's denominator `Σ basis[t]·w_t` (VBEM: over positive expectations).
#[inline]
fn class_denom<const VB: bool, const G: bool>(tids: &[u32], ws: &[f64], basis: &[f64]) -> f64 {
    let mut denom = 0.0;
    if VB {
        for (&t, &w) in tids.iter().zip(ws) {
            let e = basis[t as usize];
            denom += if !G || e > 0.0 { e * w } else { 0.0 };
        }
    } else {
        for (&t, &w) in tids.iter().zip(ws) {
            denom += basis[t as usize] * w;
        }
    }
    denom
}

/// `max_d` folded with transcript `t`'s relative change, as salmon's criterion.
#[inline]
fn fold_d(max_d: f64, a_in: f64, a_out: f64, cutoff: f64) -> f64 {
    // Filter on `a_out`, as salmon (and `crate::max_rel_diff`) does; a positive
    // cutoff also keeps the denominator non-zero.
    if a_out > cutoff {
        let d = (a_out - a_in).abs() / a_out;
        if d > max_d {
            return d;
        }
    }
    max_d
}

impl CompView<'_> {
    /// The component's fragments: its multi-transcript class counts, then its
    /// single-transcript totals (each summed in layout order).
    fn reads(&self) -> f64 {
        self.counts.iter().sum::<f64>() + self.singleton.iter().sum::<f64>()
    }

    /// One M-step `next = F(alpha)` of this component; returns the convergence
    /// criterion over its transcripts (`-inf` if none is above the cutoff).
    #[inline]
    fn step<const VB: bool, const G: bool>(
        &self,
        pull: Option<&Pull>,
        s: &mut CompState,
        cutoff: f64,
    ) -> (f64, f64) {
        match pull {
            Some(pl) => self.step_pull::<VB, G>(pl, s, cutoff),
            None => self.step_push::<VB, G>(s, cutoff),
        }
    }

    /// [`Self::step_push`], parallel: the same terms in the same order.
    fn step_pull<const VB: bool, const G: bool>(
        &self,
        pl: &Pull,
        s: &mut CompState,
        cutoff: f64,
    ) -> (f64, f64) {
        const CH: usize = 2048;
        if VB {
            s.et.par_chunks_mut(CH)
                .zip(s.alpha.par_chunks(CH))
                .zip(s.prior.par_chunks(CH))
                .for_each(|((et, al), pr)| fill_et(et, al, pr));
        }
        let basis: &[f64] = if VB { s.et } else { s.alpha };
        s.scale
            .par_chunks_mut(CH)
            .enumerate()
            .for_each(|(chunk, sc)| {
                for (j, out) in sc.iter_mut().enumerate() {
                    let k = chunk * CH + j;
                    let (lo, hi) = (
                        self.cls_start[k] - self.base,
                        self.cls_start[k + 1] - self.base,
                    );
                    let denom =
                        class_denom::<VB, G>(&self.labels[lo..hi], &self.weights[lo..hi], basis);
                    // A class the push step skips is marked NaN, which a valid
                    // scale (a count over a positive denominator) never is.
                    *out = if denom > MIN_EQ_CLASS_WEIGHT {
                        self.counts[k] / denom
                    } else {
                        f64::NAN
                    };
                }
            });
        let scale = &s.scale;
        s.next
            .par_chunks_mut(CH)
            .zip(s.alpha.par_chunks(CH))
            .enumerate()
            .map(|(chunk, (nx, al))| {
                let mut max_d = f64::NEG_INFINITY;
                let mut l1 = 0.0f64;
                for (j, (o, &a_in)) in nx.iter_mut().zip(al).enumerate() {
                    let t = chunk * CH + j;
                    let b = basis[t];
                    let mut acc = self.singleton[t];
                    if !(VB && G) || b > 0.0 {
                        for e in pl.start[t]..pl.start[t + 1] {
                            let inv = scale[pl.cls[e] as usize];
                            if inv.is_nan() {
                                continue;
                            }
                            let v = b * pl.w[e];
                            if VB || !v.is_nan() {
                                acc += v * inv;
                            }
                        }
                    }
                    *o = acc;
                    max_d = fold_d(max_d, a_in, acc, cutoff);
                    l1 += (acc - a_in).abs();
                }
                (max_d, l1)
            })
            // Collected in chunk order and folded sequentially, so the L1 sum
            // (unlike a rayon reduce tree) does not depend on the thread count.
            .collect::<Vec<_>>()
            .into_iter()
            .fold((f64::NEG_INFINITY, 0.0), |a, b| {
                (if b.0 > a.0 { b.0 } else { a.0 }, a.1 + b.1)
            })
    }

    #[inline]
    fn step_push<const VB: bool, const G: bool>(
        &self,
        s: &mut CompState,
        cutoff: f64,
    ) -> (f64, f64) {
        if VB {
            fill_et(s.et, s.alpha, s.prior);
        }
        let basis: &[f64] = if VB { s.et } else { s.alpha };
        let next = &mut *s.next;
        next.copy_from_slice(self.singleton);
        let mut lo = 0usize;
        for (k, &count) in self.counts.iter().enumerate() {
            let hi = self.cls_start[k + 1] - self.base;
            let tids = &self.labels[lo..hi];
            let ws = &self.weights[lo..hi];
            lo = hi;
            let denom = class_denom::<VB, G>(tids, ws, basis);
            if denom > MIN_EQ_CLASS_WEIGHT {
                let inv = count / denom;
                for (&t, &w) in tids.iter().zip(ws) {
                    let v = basis[t as usize] * w;
                    if VB {
                        if !G || basis[t as usize] > 0.0 {
                            next[t as usize] += v * inv;
                        }
                    } else if !v.is_nan() {
                        next[t as usize] += v * inv;
                    }
                }
            }
        }
        let mut max_d = f64::NEG_INFINITY;
        let mut l1 = 0.0f64;
        for (&a_in, &a_out) in s.alpha.iter().zip(next.iter()) {
            max_d = fold_d(max_d, a_in, a_out, cutoff);
            l1 += (a_out - a_in).abs();
        }
        (max_d, l1)
    }
}

/// Split `v` into consecutive mutable pieces at `bounds` (`bounds[0] == 0`).
fn split_at_bounds<'a, T>(mut v: &'a mut [T], bounds: &[usize]) -> Vec<&'a mut [T]> {
    let mut out = Vec::with_capacity(bounds.len().saturating_sub(1));
    for w in bounds.windows(2) {
        let (head, tail) = std::mem::take(&mut v).split_at_mut(w[1] - w[0]);
        out.push(head);
        v = tail;
    }
    out
}

/// Components grouped into consecutive batches of at least this many
/// incidences for parallel stepping, so a task is never a handful of tiny
/// components. Scheduling only; it does not touch the arithmetic.
const BATCH_INCIDENCES: usize = 32 * 1024;

/// Iterations a tile may run before every component is brought level again.
const TILE: u32 = 64;

/// The iteration's settings, as [`crate::EmOptions`] gives them.
#[derive(Clone, Copy)]
pub(crate) struct Settings {
    pub(crate) vb: bool,
    pub(crate) cutoff: f64,
    pub(crate) tol: f64,
    pub(crate) min_iter: u32,
    pub(crate) max_iter: u32,
    pub(crate) parallel: bool,
    /// VBEM: whether an expectation can be `0` (some prior below
    /// [`POSITIVE_PRIOR_MIN`]), so that the M-step must skip zero terms. Set
    /// by [`run`]; when `false` every expectation is a positive normal and the
    /// skips, always taken, are compiled out (bit-identical either way).
    pub(crate) guard_zero_et: bool,
    /// When a run stops: [`EmStop::Global`] runs [`run_tiled`],
    /// [`EmStop::Projected`] [`run_projected`].
    pub(crate) stop: EmStop,
}

/// How a component-engine run ended, per component.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ComponentSummary {
    /// the stopping rule the run used
    pub stop: EmStop,
    /// connected components iterated
    pub num_components: usize,
    /// components that never met salmon's relative-change criterion
    pub num_unconverged: usize,
    /// components stopped by `max_iter` rather than by the stopping rule
    pub num_at_max_iter: usize,
    /// fragments (class counts) in the components stopped by `max_iter`
    pub reads_at_max_iter: f64,
    /// fragments in all components
    pub total_reads: f64,
    /// median, over components, of the iterations each ran
    pub median_iters: u32,
}

impl ComponentSummary {
    /// The `meta_info.json` record of this run.
    pub fn to_meta(&self) -> salmon_core::EmComponentsMeta {
        salmon_core::EmComponentsMeta {
            stop_rule: match self.stop {
                EmStop::Projected => "projected",
                EmStop::Global => "global",
            }
            .to_string(),
            num_components: self.num_components,
            num_unconverged: self.num_unconverged,
            num_at_max_iter: self.num_at_max_iter,
            fraction_reads_at_max_iter: self.fraction_at_max_iter(),
            median_iterations: self.median_iters,
        }
    }

    /// Fraction of the fragments held by components cut off by `max_iter`.
    pub fn fraction_at_max_iter(&self) -> f64 {
        if self.total_reads > 0.0 {
            self.reads_at_max_iter / self.total_reads
        } else {
            0.0
        }
    }

    /// One line for the run log, given the most iterations any component ran;
    /// `true` when it deserves a warning (some component never met salmon's
    /// criterion).
    pub fn describe(&self, iters: u32) -> (String, bool) {
        let pct = |n: f64| 100.0 * n / self.total_reads.max(1.0);
        let mut line = match self.stop {
            EmStop::Global => format!(
                "EM: {} iterations over {} connected components (stopping rule: global)",
                iters, self.num_components
            ),
            EmStop::Projected => format!(
                "EM: {} connected components, each iterated until it converged \
                 (median {} iterations, most {})",
                self.num_components, self.median_iters, iters
            ),
        };
        if self.num_unconverged > 0 {
            line += &format!(
                "; {} component(s) never met the convergence criterion",
                self.num_unconverged
            );
        } else if self.num_at_max_iter > 0 {
            line += &format!(
                "; {} slowly mixing component(s) ({:.2}% of fragments) met the \
                 convergence criterion but were still shifting mass between \
                 near-identical transcripts at the {}-iteration cap, and stop there",
                self.num_at_max_iter,
                pct(self.reads_at_max_iter),
                iters
            );
        }
        (line, self.num_unconverged > 0)
    }
}

/// The dense working vectors of a run, split per component.
struct Buffers {
    alpha: Vec<f64>,
    next: Vec<f64>,
    prior: Vec<f64>,
    et: Vec<f64>,
}

impl Buffers {
    fn new(plan: &ComponentPlan, init: &[f64], prior_full: &[f64], vb: bool) -> Self {
        let nd = plan.txps.len();
        let nv = if vb { nd } else { 0 };
        Self {
            alpha: plan.txps.iter().map(|&t| init[t as usize]).collect(),
            next: vec![0.0; nd],
            prior: plan.txps.iter().map(|&t| prior_full[t as usize]).collect(),
            et: vec![0.0; nv],
        }
    }

    /// One [`CompState`] per component, in component order.
    fn states<'a>(&'a mut self, plan: &ComponentPlan) -> Vec<(usize, CompState<'a>)> {
        let bounds = &plan.comp_txp;
        let empty = vec![0usize; bounds.len()];
        let vb_bounds = if self.et.is_empty() { &empty } else { bounds };
        let prior = &self.prior;
        split_at_bounds(&mut self.alpha, bounds)
            .into_iter()
            .zip(split_at_bounds(&mut self.next, bounds))
            .zip(split_at_bounds(&mut self.et, vb_bounds))
            .enumerate()
            .map(|(c, ((alpha, next), et))| {
                let scale = if plan.pull[c].is_some() {
                    vec![0.0; plan.comp_cls[c + 1] - plan.comp_cls[c]]
                } else {
                    Vec::new()
                };
                let state = CompState {
                    alpha,
                    next,
                    prior: &prior[bounds[c]..bounds[c + 1]],
                    et,
                    changed: true,
                    scale,
                };
                (c, state)
            })
            .collect()
    }
}

/// One M-step of component `c`, swapping its buffers so `alpha` holds the
/// new iterate; returns the convergence criterion over the component.
#[inline]
fn step_component(plan: &ComponentPlan, c: usize, s: &mut CompState, set: &Settings) -> f64 {
    step_component_l1(plan, c, s, set).0
}

/// [`step_component`], also returning the L1 size of the step.
#[inline]
fn step_component_l1(
    plan: &ComponentPlan,
    c: usize,
    s: &mut CompState,
    set: &Settings,
) -> (f64, f64) {
    let v = plan.view(c);
    let pull = plan.pull[c].as_ref();
    let d = match (set.vb, set.guard_zero_et) {
        (true, true) => v.step::<true, true>(pull, s, set.cutoff),
        (true, false) => v.step::<true, false>(pull, s, set.cutoff),
        (false, _) => v.step::<false, false>(pull, s, set.cutoff),
    };
    std::mem::swap(&mut s.alpha, &mut s.next);
    d
}

/// The larger of two criterion values (neither is ever NaN).
#[inline]
fn max_of(a: f64, b: f64) -> f64 {
    if b > a {
        b
    } else {
        a
    }
}

/// Run EM or VBEM to convergence over the component plan, starting from
/// `init` (full length, indexed by original transcript id). Returns
/// `(alphas, iterations, converged, summary)`, with every transcript no class
/// references at `0.0`, as the M-step leaves it. `iterations` is the most any
/// component ran; `converged` says whether salmon's criterion was met (over
/// all transcripts under [`EmStop::Global`], over every component under
/// [`EmStop::Projected`]).
pub(crate) fn run(
    plan: &ComponentPlan,
    init: &[f64],
    prior_full: &[f64],
    set: &Settings,
) -> (Vec<f64>, u32, bool, ComponentSummary) {
    let mut buffers = Buffers::new(plan, init, prior_full, set.vb);
    let set = &Settings {
        guard_zero_et: !buffers.prior.iter().all(|&p| p >= POSITIVE_PRIOR_MIN),
        ..*set
    };
    let mut out = vec![0.0f64; init.len()];
    let mut states = buffers.states(plan);
    let reads: Vec<f64> = (0..plan.num_components())
        .map(|c| plan.view(c).reads())
        .collect();
    let total_reads = reads.iter().sum();
    let (iters, converged, summary) = match set.stop {
        EmStop::Global => {
            let (it, ok) = run_tiled(plan, &mut states, set);
            let summary = ComponentSummary {
                stop: set.stop,
                num_components: plan.num_components(),
                num_unconverged: if ok { 0 } else { plan.num_components() },
                num_at_max_iter: if ok { 0 } else { plan.num_components() },
                reads_at_max_iter: if ok { 0.0 } else { total_reads },
                total_reads,
                median_iters: it,
            };
            (it, ok, summary)
        }
        EmStop::Projected => {
            let ends = run_projected(plan, &mut states, set);
            let mut iters: Vec<u32> = ends.iter().map(|e| e.iters).collect();
            iters.sort_unstable();
            let capped = ends.iter().filter(|e| !e.stopped);
            let summary = ComponentSummary {
                stop: set.stop,
                num_components: ends.len(),
                num_unconverged: ends.iter().filter(|e| !e.met).count(),
                num_at_max_iter: capped.clone().count(),
                reads_at_max_iter: capped.map(|e| reads[e.comp]).sum(),
                total_reads,
                median_iters: iters.get(iters.len() / 2).copied().unwrap_or(0),
            };
            let max = iters.last().copied().unwrap_or(0);
            (max, summary.num_unconverged == 0, summary)
        }
    };
    write_out(plan, &states, &mut out);
    (out, iters, converged, summary)
}

/// The global loop, scheduled component-major in tiles of up to [`TILE`]
/// iterations so that a component's data stays in cache across iterations
/// (the schedule kallistors uses).
///
/// Exactly the plain loop's result: each component's arithmetic depends only
/// on its own state, so only the stopping iteration has to be reproduced. A
/// *witness* component — the largest whose criterion last failed — runs first,
/// until its own criterion holds at an eligible iteration (or the tile cap).
/// The global criterion cannot have held at any earlier iteration of the tile,
/// since the witness's alone failed there; so every other component is
/// advanced to exactly that iteration, and the global criterion is evaluated
/// there as the plain loop would.
fn run_tiled(
    plan: &ComponentPlan,
    states: &mut [(usize, CompState)],
    set: &Settings,
) -> (u32, bool) {
    if states.is_empty() {
        return (0, false);
    }
    // Batch boundaries over the (component-ordered) states.
    let mut bounds = vec![0usize];
    let mut size = 0;
    for (i, st) in states.iter().enumerate() {
        size += plan.comp_incidences(st.0) + 1;
        if size >= BATCH_INCIDENCES {
            bounds.push(i + 1);
            size = 0;
        }
    }
    if *bounds.last().expect("starts with 0") != states.len() {
        bounds.push(states.len());
    }
    let mut it = 0u32;
    while it < set.max_iter {
        let cap = TILE.min(set.max_iter - it);
        let w = (0..states.len())
            .max_by_key(|&i| {
                let (c, s) = &states[i];
                (s.changed, plan.comp_incidences(*c), std::cmp::Reverse(i))
            })
            .expect("states is not empty");
        let mut taken = 0;
        let mut dw = f64::NEG_INFINITY;
        {
            let (c, s) = &mut states[w];
            for off in 0..cap {
                dw = step_component(plan, *c, s, set);
                taken = off + 1;
                if it + taken >= set.min_iter && dw < set.tol {
                    break;
                }
            }
            s.changed = dw >= set.tol;
        }
        let run_batch = |batch: &mut [(usize, CompState)], first: usize| -> f64 {
            let mut m = f64::NEG_INFINITY;
            for (k, (c, s)) in batch.iter_mut().enumerate() {
                if first + k == w {
                    continue;
                }
                let mut d = f64::NEG_INFINITY;
                for _ in 0..taken {
                    d = step_component(plan, *c, s, set);
                }
                s.changed = d >= set.tol;
                m = max_of(m, d);
            }
            m
        };
        let batches = split_at_bounds(states, &bounds);
        let others = if set.parallel {
            batches
                .into_par_iter()
                .zip(bounds.par_iter())
                .map(|(b, &first)| run_batch(b, first))
                .reduce(|| f64::NEG_INFINITY, max_of)
        } else {
            batches
                .into_iter()
                .zip(bounds.iter())
                .map(|(b, &first)| run_batch(b, first))
                .fold(f64::NEG_INFINITY, max_of)
        };
        it += taken;
        let d = max_of(dw, others);
        if it >= set.min_iter && d < set.tol {
            return (it, true);
        }
    }
    (it, false)
}

/// Each component iterated on its own until salmon's criterion holds over its
/// own transcripts (after `min_iter`), or `max_iter`. Returns the most
/// iterations any component ran, and whether every component converged.
///
/// Successive step-size ratios the contraction-rate estimate is taken over.
const RATE_WINDOW: usize = 5;

/// Cap on the estimated contraction rate, so `1 / (1 − ρ)` stays finite.
const RATE_CAP: f64 = 1.0 - 1e-4;

/// Largest projected remaining movement at which a component may stop, as a
/// fraction of its fragments.
const PROJECTED_EPS: f64 = 1e-3;

/// How one component's iteration ended under [`EmStop::Projected`].
struct CompEnd {
    comp: usize,
    iters: u32,
    /// salmon's criterion held over the component at some eligible iteration
    met: bool,
    /// the stopping rule (not `max_iter`) ended it
    stopped: bool,
}

/// Each component iterated on its own until salmon's criterion holds over its
/// transcripts (after `min_iter`) **and** its projected remaining movement is
/// below [`PROJECTED_EPS`] of its fragments, or until `max_iter`.
///
/// Under the global rule every component runs as many iterations as the
/// slowest one needs; most of that work is on components that met the
/// criterion long before. Stopping each component on the criterion alone is
/// not enough, though: a family of near-identical isoforms (*TTN*, *NEB*)
/// trades mass so slowly that each step is tiny — the criterion holds — while
/// the family is still far from its fixpoint. For a linearly convergent
/// iteration with contraction rate `ρ`, the movement still to come after a
/// step of L1 size `l` is about `l·ρ/(1 − ρ)`; `ρ` is estimated as the largest
/// ratio of successive step sizes over the last [`RATE_WINDOW`] steps (capped
/// at [`RATE_CAP`]), and `l/(1 − ρ)` must be small next to the component's
/// fragment count. Fast-mixing components stop where the criterion alone
/// would; slow ones keep iterating.
///
/// A component's iterate depends only on its own data and its stopping test
/// only on its own step sizes (each summed in a data-defined order), so
/// scheduling is free and the result is independent of the thread count.
fn run_projected(
    plan: &ComponentPlan,
    states: &mut [(usize, CompState)],
    set: &Settings,
) -> Vec<CompEnd> {
    let one = |(c, s): &mut (usize, CompState)| -> CompEnd {
        let mass = plan.view(*c).reads();
        let mut ratios = [f64::INFINITY; RATE_WINDOW];
        let mut prev_l1 = f64::NAN;
        let mut met = false;
        let mut it = 0;
        while it < set.max_iter {
            let (d, l1) = step_component_l1(plan, *c, s, set);
            ratios[it as usize % RATE_WINDOW] = if prev_l1 > 0.0 {
                l1 / prev_l1
            } else if l1 == 0.0 {
                0.0
            } else {
                f64::INFINITY
            };
            prev_l1 = l1;
            it += 1;
            if it < set.min_iter {
                continue;
            }
            // No transcript above the cutoff: nothing the criterion measures.
            if d == f64::NEG_INFINITY {
                return CompEnd {
                    comp: *c,
                    iters: it,
                    met: true,
                    stopped: true,
                };
            }
            if d < set.tol {
                met = true;
                let rho = ratios.iter().copied().fold(0.0f64, f64::max).min(RATE_CAP);
                if l1 / (1.0 - rho) < PROJECTED_EPS * mass {
                    return CompEnd {
                        comp: *c,
                        iters: it,
                        met,
                        stopped: true,
                    };
                }
            }
        }
        CompEnd {
            comp: *c,
            iters: it,
            met,
            stopped: false,
        }
    };
    if set.parallel {
        // Largest first, so the long tail of small components fills in behind
        // the big ones. Scheduling only.
        let mut order: Vec<&mut (usize, CompState)> = states.iter_mut().collect();
        order.sort_by_key(|(c, _)| std::cmp::Reverse(plan.comp_incidences(*c)));
        order.into_par_iter().with_max_len(1).map(&one).collect()
    } else {
        states.iter_mut().map(one).collect()
    }
}

/// Scatter each component's latest iterate to original transcript ids. The
/// M-steps swap a component's two buffer references, so its latest iterate
/// is wherever its `alpha` reference now points.
fn write_out(plan: &ComponentPlan, states: &[(usize, CompState)], out: &mut [f64]) {
    for (c, s) in states {
        let t0 = plan.comp_txp[*c];
        for (&t, &a) in plan.txps[t0..t0 + s.alpha.len()].iter().zip(s.alpha.iter()) {
            out[t as usize] = a;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::packed::{em_step_seq, vbem_step_seq};
    use salmon_eqclass::{EquivalenceClassBuilder, TranscriptGroup};

    /// Classes drawn within `blocks` blocks of `per_block` transcripts (so the
    /// graph has many components), plus single-transcript classes; weights and
    /// effective lengths chosen so the sums are not exactly representable.
    fn fixture(blocks: usize, per_block: usize, classes: usize) -> PackedEqClasses {
        let b = EquivalenceClassBuilder::new();
        let mut seed = 0x9E37_79B9_7F4A_7C15u64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let num_txps = blocks * per_block + 7;
        for c in 0..classes {
            let base = (c % blocks) * per_block;
            let n = if c % 9 == 0 { 1 } else { 2 + (c % 5) };
            let mut tids: Vec<u32> = (0..n)
                .map(|_| (base + (rnd() % per_block as u64) as usize) as u32)
                .collect();
            tids.sort_unstable();
            tids.dedup();
            let ws: Vec<f64> = (0..tids.len())
                .map(|i| 0.1 + ((rnd() % 1000) as f64) / 997.0 + i as f64 * 1e-3)
                .collect();
            b.add_group(TranscriptGroup::new(tids), ws, 1 + (rnd() % 50));
        }
        let mut eq = b.finish();
        let eff: Vec<f64> = (0..num_txps)
            .map(|i| 500.0 + (i % 997) as f64 * 1.7)
            .collect();
        eq.update_eff_lengths(&eff);
        PackedEqClasses::from_collapsed(&eq, num_txps)
    }

    fn settings(vb: bool, parallel: bool) -> Settings {
        Settings {
            vb,
            cutoff: 1e-2,
            tol: 1e-2,
            min_iter: 50,
            max_iter: 10_000,
            parallel,
            // Guarded here; `run` drops the guard for these priors, so the
            // tests below also check that dropping it changes no bit.
            guard_zero_et: true,
            stop: EmStop::Global,
        }
    }

    fn start(p: &PackedEqClasses) -> (Vec<f64>, Vec<f64>) {
        let total: u64 = p.counts.iter().sum();
        let init = vec![total as f64 / p.num_txps as f64; p.num_txps];
        (init, vec![1e-2; p.num_txps])
    }

    /// The plain loop: every component one M-step per iteration.
    fn untiled(
        plan: &ComponentPlan,
        init: &[f64],
        prior: &[f64],
        set: &Settings,
    ) -> (Vec<f64>, u32, bool) {
        let mut buffers = Buffers::new(plan, init, prior, set.vb);
        let mut states = buffers.states(plan);
        let mut it = 0;
        let mut converged = false;
        while it < set.max_iter {
            let mut d = f64::NEG_INFINITY;
            for (c, s) in states.iter_mut() {
                d = max_of(d, step_component(plan, *c, s, set));
            }
            it += 1;
            if it >= set.min_iter && d < set.tol {
                converged = true;
                break;
            }
        }
        let mut out = vec![0.0; init.len()];
        write_out(plan, &states, &mut out);
        (out, it, converged)
    }

    fn bits(v: &[f64]) -> Vec<u64> {
        v.iter().map(|x| x.to_bits()).collect()
    }

    /// Tiling only reschedules: same iterates, same stopping iteration.
    #[test]
    fn tiled_schedule_matches_the_plain_loop_bitwise() {
        let p = fixture(150, 20, 6000);
        let (init, prior) = start(&p);
        let plan = ComponentPlan::with_pull_min(&p, &p.counts, usize::MAX);
        assert!(plan.num_components() > 100);
        for vb in [false, true] {
            let set = settings(vb, false);
            let (a, ia, ca, _) = run(&plan, &init, &prior, &set);
            let (b, ib, cb) = untiled(&plan, &init, &prior, &set);
            assert_eq!((ia, ca), (ib, cb), "vb={vb}");
            assert!(ia > TILE, "the run must span several tiles");
            assert_eq!(bits(&a), bits(&b), "vb={vb}");
        }
    }

    /// Gathering adds the same terms in the same order as the class-major
    /// loop, so neither the stepping strategy nor the pool size changes a bit.
    #[test]
    fn parallel_gather_matches_sequential_push_bitwise() {
        let p = fixture(3, 400, 9000);
        let (init, prior) = start(&p);
        let push = ComponentPlan::with_pull_min(&p, &p.counts, usize::MAX);
        let pull = ComponentPlan::with_pull_min(&p, &p.counts, 0);
        assert!(pull.pull.iter().all(Option::is_some));
        for vb in [false, true] {
            let reference = run(&push, &init, &prior, &settings(vb, false));
            for threads in [1, 4] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let got = pool.install(|| run(&pull, &init, &prior, &settings(vb, true)));
                assert_eq!((got.1, got.2), (reference.1, reference.2));
                assert_eq!(bits(&got.0), bits(&reference.0), "vb={vb} t{threads}");
            }
        }
    }

    /// Against the global kernels, only rounding moves: single-transcript
    /// counts are added first instead of in class order, and VBEM drops the
    /// common `ψ(Σα)` term that cancels in every class's ratio.
    #[test]
    fn matches_the_global_kernels_to_rounding() {
        let p = fixture(40, 30, 5000);
        let (init, prior) = start(&p);
        let plan = ComponentPlan::with_pull_min(&p, &p.counts, usize::MAX);
        for vb in [false, true] {
            let set = settings(vb, false);
            let (got, iters, converged, _) = run(&plan, &init, &prior, &set);
            let mut a = init.clone();
            let mut b = vec![0.0; p.num_txps];
            let mut et = vec![0.0; p.num_txps];
            let mut scratch = Vec::new();
            let mut it = 0;
            loop {
                if vb {
                    vbem_step_seq::<false>(
                        &p,
                        &p.counts,
                        &prior,
                        &a,
                        &mut b,
                        &mut et,
                        &mut scratch,
                    );
                } else {
                    em_step_seq::<false>(&p, &p.counts, &a, &mut b, &mut scratch);
                }
                it += 1;
                let d = crate::max_rel_diff(&a, &b, set.cutoff);
                std::mem::swap(&mut a, &mut b);
                if it >= set.min_iter && d < set.tol {
                    break;
                }
            }
            assert!(converged);
            assert_eq!(iters, it, "vb={vb}");
            for (x, y) in got.iter().zip(&a) {
                assert!(
                    (x - y).abs() <= 1e-9 * y.abs().max(1.0),
                    "vb={vb}: {x} vs {y}"
                );
            }
        }
    }

    /// The projected rule: neither the pool size nor stepping by gathering
    /// changes a bit or the summary.
    #[test]
    fn projected_stop_is_independent_of_threads_and_stepping() {
        let p = fixture(60, 15, 3000);
        let (init, prior) = start(&p);
        let seq = ComponentPlan::with_pull_min(&p, &p.counts, usize::MAX);
        let pull = ComponentPlan::with_pull_min(&p, &p.counts, 0);
        for vb in [false, true] {
            let set = Settings {
                stop: EmStop::Projected,
                ..settings(vb, false)
            };
            let (got, iters, converged, summary) = run(&seq, &init, &prior, &set);
            assert!(converged);
            assert_eq!(summary.num_components, seq.num_components());
            assert_eq!(summary.num_at_max_iter, 0);
            for threads in [1, 4] {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .unwrap();
                let par = pool.install(|| {
                    run(
                        &pull,
                        &init,
                        &prior,
                        &Settings {
                            parallel: true,
                            ..set
                        },
                    )
                });
                assert_eq!(bits(&got), bits(&par.0), "vb={vb} t{threads}");
                assert_eq!((iters, converged, &summary), (par.1, par.2, &par.3));
            }
        }
    }

    /// Two near-identical transcripts that differ by a handful of unique
    /// reads: each step moves little, so salmon's criterion holds long before
    /// the split is near its fixpoint. The projected rule keeps iterating and
    /// lands far closer to it.
    #[test]
    fn projected_stop_keeps_a_slowly_mixing_pair_iterating() {
        let b = EquivalenceClassBuilder::new();
        b.add_group(TranscriptGroup::new(vec![0, 1]), vec![1.0, 1.0], 100_000);
        b.add_group(TranscriptGroup::new(vec![0]), vec![1.0], 60);
        b.add_group(TranscriptGroup::new(vec![1]), vec![1.0], 20);
        let mut eq = b.finish();
        eq.update_eff_lengths(&[1000.0, 1000.0]);
        let p = PackedEqClasses::from_collapsed(&eq, 2);
        let (init, prior) = start(&p);
        let plan = ComponentPlan::with_pull_min(&p, &p.counts, usize::MAX);
        for vb in [false, true] {
            let global = run(&plan, &init, &prior, &settings(vb, false));
            let fixpoint = run(
                &plan,
                &init,
                &prior,
                &Settings {
                    tol: 1e-12,
                    max_iter: 1_000_000,
                    ..settings(vb, false)
                },
            );
            let set = Settings {
                stop: EmStop::Projected,
                ..settings(vb, false)
            };
            let (got, iters, converged, summary) = run(&plan, &init, &prior, &set);
            assert!(converged && summary.num_at_max_iter == 0, "vb={vb}");
            assert!(iters > global.1, "vb={vb}: {iters} vs {}", global.1);
            let err = |a: &[f64]| (a[0] - fixpoint.0[0]).abs();
            assert!(err(&got) * 10.0 < err(&global.0), "vb={vb}");
            assert!(err(&got) < PROJECTED_EPS * 100_080.0, "vb={vb}");
        }
    }

    /// A component cut off by `max_iter` is counted with its fragments, and
    /// still counts as converged once salmon's criterion has held over it.
    #[test]
    fn projected_stop_reports_components_cut_off_by_max_iter() {
        let b = EquivalenceClassBuilder::new();
        b.add_group(TranscriptGroup::new(vec![0, 1]), vec![1.0, 1.0], 100_000);
        b.add_group(TranscriptGroup::new(vec![0]), vec![1.0], 60);
        b.add_group(TranscriptGroup::new(vec![1]), vec![1.0], 20);
        b.add_group(TranscriptGroup::new(vec![2, 3]), vec![1.0, 1.0], 50);
        b.add_group(TranscriptGroup::new(vec![2]), vec![1.0], 40);
        let mut eq = b.finish();
        eq.update_eff_lengths(&[1000.0; 4]);
        let p = PackedEqClasses::from_collapsed(&eq, 4);
        let (init, prior) = start(&p);
        let plan = ComponentPlan::with_pull_min(&p, &p.counts, usize::MAX);
        let set = Settings {
            stop: EmStop::Projected,
            max_iter: 200,
            ..settings(false, false)
        };
        let (_, iters, converged, summary) = run(&plan, &init, &prior, &set);
        assert_eq!(iters, 200);
        assert!(converged, "the criterion held over both components");
        assert_eq!(summary.num_components, 2);
        assert_eq!(summary.num_unconverged, 0);
        assert_eq!(summary.num_at_max_iter, 1);
        assert_eq!(summary.reads_at_max_iter, 100_080.0);
        assert_eq!(summary.total_reads, 100_170.0);
        // Cut off before the criterion is even checked: not converged.
        let short = Settings {
            max_iter: 10,
            ..set
        };
        let (_, _, converged, summary) = run(&plan, &init, &prior, &short);
        assert!(!converged && summary.num_unconverged >= 1);
    }

    /// `exp_digamma` agrees with `statrs`'s `digamma(x).exp()` over the range
    /// VBEM visits under any prior ≥ 0.01, to within the conditioning of
    /// `exp(ψ(x))` (relative `|ψ(x)|·ε`; both sides carry about twice that).
    #[test]
    fn exp_digamma_matches_statrs() {
        use statrs::function::gamma::digamma;
        let n = 200_000;
        let mut worst = 0.0f64;
        for k in 0..n {
            let x = 10f64.powf(-2.0 + 10.0 * (k as f64 + 0.5) / n as f64);
            let (got, want) = (exp_digamma(x), digamma(x).exp());
            let rel = (got - want).abs() / want;
            let bound = 8.0 * f64::EPSILON * (digamma(x).abs() + 4.0);
            assert!(
                rel <= bound,
                "x={x:e}: {got:e} vs {want:e} (rel {rel:e} > {bound:e})"
            );
            worst = worst.max(rel / bound);
        }
        assert!(worst > 0.0);
    }

    /// The vectorized dispatch computes the same bits as the portable loop.
    #[test]
    fn fill_et_dispatch_is_bit_identical() {
        let n = 10_007;
        let alpha: Vec<f64> = (0..n)
            .map(|i| match i % 4 {
                0 => 0.0,
                1 => 1e-9 * i as f64,
                2 => (i as f64 * 0.618_033_988_7).fract() * 20.0,
                _ => 1e6 * (i as f64).sqrt(),
            })
            .collect();
        for pr in [0.0, 1e-12, 1e-3, 1e-2, 7.5] {
            let prior = vec![pr; n];
            let (mut a, mut b) = (vec![0.0; n], vec![0.0; n]);
            fill_et(&mut a, &alpha, &prior);
            fill_et_body(&mut b, &alpha, &prior);
            assert_eq!(bits(&a), bits(&b), "prior {pr}");
        }
    }
}
