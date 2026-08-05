//! Sequence-specific bias model (`SBModel`).
//!
//! # The effect being modelled
//!
//! RNA-seq protocols fragment and prime with enzymes and random hexamers that
//! are not indifferent to sequence. Certain short motifs are cut or primed far
//! more readily than others, so fragments *start* at some sequences much more
//! often than at others. A transcript rich in favoured motifs yields more
//! fragments at the same abundance — again chemistry, not biology.
//!
//! # How it is measured
//!
//! A faithful port of salmon's `SBModel` (`src/model/SBModel.cpp`): a
//! variable-order Markov model over the 9-base sequence context surrounding a
//! fragment's start position (3 bases before the start, the start, and 5 after).
//! Per-position Markov orders are `{0,1,2,2,2,2,2,2,2}`.
//!
//! **What "variable-order Markov" means here.** Rather than one probability for
//! each of the 4^9 possible 9-mers — which would need far more data than exists —
//! the model predicts each position from a few preceding ones: position 0 from
//! nothing (order 0), position 1 from its predecessor (order 1), and the rest
//! from their two predecessors (order 2). The context's probability is then the
//! product of nine small conditional probabilities. That is enough structure to
//! capture real motif preference while keeping the parameter count at 64 per
//! position, which a run easily estimates.
//!
//! Counts are accumulated from observed fragment-start contexts (the *observed*
//! model) and from the transcriptome (the *expected* model); the ratio of the two
//! scores the sequence bias at any position, which is used to correct effective
//! lengths.
//!
//! The 2-bit base encoding only needs to be self-consistent (the bias is a
//! ratio of two models built with the same encoding), so we use A=0, C=1,
//! G=2, T=3.

use realfft::RealFftPlanner;
use std::cell::RefCell;

thread_local! {
    /// Per-thread FFT planner (the per-transcript correction runs in a parallel
    /// sweep). Padding to a power of two keeps the set of distinct plan sizes
    /// tiny, so plans are reused across transcripts.
    static FFT_PLANNER: RefCell<RealFftPlanner<f64>> = RefCell::new(RealFftPlanner::<f64>::new());
}

/// Cross-correlation `xc[Δ] = Σ_k a[k]·b[k+Δ]` for `Δ in [0, max_lag]`, via a
/// real FFT (zero-padded so `b` beyond its length contributes 0 — i.e. linear,
/// not circular, correlation). Correlation theorem: `corr(a,b) =
/// IFFT(conj(FFT(a))·FFT(b))`. rustfft is unnormalized, so divide by `n`.
///
/// **Why this is here.** The effective-length sweep asks, for every fragment
/// length Δ, "what is the total of (start factor × end factor) over all
/// positions?". Done directly that is O(L) work per length and O(L²) overall —
/// prohibitive for a 100 kb transcript. A cross-correlation computes the answer
/// for *every* Δ at once, and the correlation theorem turns it into three FFTs,
/// giving O(L log L).
///
/// The zero padding matters: without it the FFT would wrap around, so fragments
/// running off the right end would fold back onto the left end (circular rather
/// than linear correlation) and produce fragments that do not exist.
fn xcorr_fft(fw: &[f64], rc: &[f64], max_lag: usize) -> Vec<f64> {
    let l = fw.len();
    debug_assert_eq!(l, rc.len());
    let n = (l + max_lag + 1).next_power_of_two().max(2);
    FFT_PLANNER.with(|p| {
        let mut planner = p.borrow_mut();
        let r2c = planner.plan_fft_forward(n);
        let c2r = planner.plan_fft_inverse(n);
        let mut a = r2c.make_input_vec();
        let mut b = r2c.make_input_vec();
        a[..l].copy_from_slice(fw);
        b[..l].copy_from_slice(rc);
        let mut fa = r2c.make_output_vec();
        let mut fb = r2c.make_output_vec();
        r2c.process(&mut a, &mut fa).expect("rfft fw");
        r2c.process(&mut b, &mut fb).expect("rfft rc");
        for (x, y) in fa.iter_mut().zip(&fb) {
            *x = x.conj() * *y;
        }
        let mut out = c2r.make_output_vec();
        c2r.process(&mut fa, &mut out).expect("irfft");
        let scale = 1.0 / n as f64;
        out[..=max_lag].iter().map(|v| v * scale).collect()
    })
}

/// Per-position Markov orders (salmon's "simple" model). Length is the context.
///
/// The first positions necessarily have lower order — there is nothing before
/// position 0 to condition on, and only one base before position 1.
const ORDER: [u32; 9] = [0, 1, 2, 2, 2, 2, 2, 2, 2];
/// Context length (= ORDER.len()): 3 left + start + 5 right.
pub const CONTEXT_LENGTH: usize = 9;
/// Bases before the fragment-start position.
pub const CONTEXT_LEFT: usize = 3;
/// Bases at/after the fragment-start position.
pub const CONTEXT_RIGHT: usize = 5;
/// Rows in the probability table: 4^(maxOrder+1) = 4^3.
///
/// One row per (two-base context, predicted base) combination; lower-order
/// positions use only the first few rows.
const ROWS: usize = 64;
/// Pseudocount prior.
///
/// Keeps an unobserved transition from being exactly zero, which would make its
/// log negative infinity and poison every context containing it.
const PRIOR: f64 = 1e-10;
/// Floor used when taking the log of a zero probability.
///
/// A finite floor rather than -inf, so a single impossible transition cannot
/// annihilate an otherwise plausible context.
const LOG_SMALL: f64 = -11.512_925_464_970_229; // ln(1e-5)

/// 2-bit encode an ASCII base (non-ACGT -> 0).
#[inline]
fn base2bit(b: u8) -> u32 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

#[inline]
fn complement_bit(x: u32) -> u32 {
    3 - x // A<->T (0<->3), C<->G (1<->2)
}

/// The sequence-specific bias Markov model.
#[derive(Debug, Clone)]
pub struct SBModel {
    /// log (after [`normalize`](Self::normalize)) or linear (before) transition
    /// probabilities, laid out position-major: `probs[pos * ROWS + idx]`. Before
    /// `normalize` this is materialized from the integer `probs_fp`.
    probs: Vec<f64>,
    /// Fixed-point integer accumulator for the transition counts (`weight *
    /// BIAS_WEIGHT_SCALE`, truncated), summed order-independently across worker
    /// threads and materialized into `probs` (plus the `PRIOR`) at `normalize`.
    probs_fp: Vec<u64>,
    /// per-position base marginals: `marginals[pos * 4 + base]`
    marginals: Vec<f64>,
    shifts: [u32; CONTEXT_LENGTH],
    masks: [u32; CONTEXT_LENGTH],
    trained: bool,
}

impl Default for SBModel {
    fn default() -> Self {
        Self::new()
    }
}

impl SBModel {
    pub fn new() -> Self {
        let mut shifts = [0u32; CONTEXT_LENGTH];
        let mut masks = [0u32; CONTEXT_LENGTH];
        for i in 0..CONTEXT_LENGTH {
            // base i occupies the high bits; isolate the (order+1)-mer ending at i
            shifts[i] = (2 * CONTEXT_LENGTH as u32) - 2 * (i as u32 + 1);
            let width = 2 * (ORDER[i] + 1);
            masks[i] = (1u32 << width) - 1;
        }
        Self {
            probs: vec![PRIOR; ROWS * CONTEXT_LENGTH],
            probs_fp: vec![0u64; ROWS * CONTEXT_LENGTH],
            marginals: vec![PRIOR; 4 * CONTEXT_LENGTH],
            shifts,
            masks,
            trained: false,
        }
    }

    /// Encode a 9-base context (`CONTEXT_LENGTH` bytes) into a 2-bit-per-base
    /// integer with base 0 in the high bits. `rev_comp` reverse-complements it.
    ///
    /// Packing the whole context into one `u32` means every per-position lookup
    /// below is a shift and a mask rather than a slice access.
    fn encode(context: &[u8], rev_comp: bool) -> u32 {
        debug_assert_eq!(context.len(), CONTEXT_LENGTH);
        let mut mer = 0u32;
        if rev_comp {
            // reverse complement: last base becomes first
            for &b in context.iter().rev() {
                mer = (mer << 2) | complement_bit(base2bit(b));
            }
        } else {
            for &b in context {
                mer = (mer << 2) | base2bit(b);
            }
        }
        mer
    }

    #[inline]
    fn index_at(&self, mer: u32, pos: usize) -> usize {
        ((mer >> self.shifts[pos]) & self.masks[pos]) as usize
    }

    /// The flattened transition table (`probs[pos * ROWS + idx]`), for dumping to
    /// the aux bias files. Linear counts before [`normalize`](Self::normalize),
    /// conditional log-probabilities after.
    pub fn dump(&self) -> &[f64] {
        &self.probs
    }

    /// Accumulate one observed context with the given weight.
    pub fn add_context(&mut self, context: &[u8], rev_comp: bool, weight: f64) {
        debug_assert!(!self.trained, "cannot add to a normalized model");
        let mer = Self::encode(context, rev_comp);
        let w = crate::bias_mass_to_fp(weight);
        for pos in 0..CONTEXT_LENGTH {
            let idx = self.index_at(mer, pos);
            self.probs_fp[pos * ROWS + idx] += w;
        }
    }

    /// Convert accumulated counts into conditional log-probabilities. Idempotent
    /// guard: a model can only be normalized once.
    ///
    /// Two steps: divide each group of four counts (the four possible bases given
    /// one preceding context) by their total, turning counts into conditional
    /// probabilities; then take logs so `evaluate_log` can add instead of
    /// multiply.
    pub fn normalize(&mut self) {
        if self.trained {
            return;
        }
        // Materialize the integer counts into `probs`, reintroducing the `PRIOR`
        // pseudocount in f64 (matching the pre-fixed-point `probs = PRIOR + Σw`).
        for (p, &fp) in self.probs.iter_mut().zip(&self.probs_fp) {
            *p = PRIOR + fp as f64 / crate::BIAS_WEIGHT_SCALE;
        }
        for pos in 0..CONTEXT_LENGTH {
            let num_states = 4usize.pow(ORDER[pos]);
            for s in 0..num_states {
                let node = s * 4;
                let base = pos * ROWS + node;
                let tot: f64 = self.probs[base..base + 4].iter().sum();
                if tot > 0.0 {
                    for j in 0..4 {
                        self.probs[base + j] /= tot;
                        self.marginals[pos * 4 + j] += self.probs[base + j];
                    }
                }
            }
            for j in 0..4 {
                self.marginals[pos * 4 + j] /= num_states as f64;
            }
        }
        for p in &mut self.probs {
            *p = if *p > 0.0 { p.ln() } else { LOG_SMALL };
        }
        self.trained = true;
    }

    /// Log-probability the (normalized) model assigns to a context.
    ///
    /// A sum of nine per-position log-probabilities, which is the product of the
    /// nine conditional probabilities — the chain rule for a Markov model.
    pub fn evaluate_log(&self, context: &[u8], rev_comp: bool) -> f64 {
        debug_assert!(self.trained, "evaluate_log requires a normalized model");
        let mer = Self::encode(context, rev_comp);
        let mut lp = 0.0;
        for pos in 0..CONTEXT_LENGTH {
            let idx = self.index_at(mer, pos);
            lp += self.probs[pos * ROWS + idx];
        }
        lp
    }

    pub fn is_trained(&self) -> bool {
        self.trained
    }

    /// Add another (un-normalized) model's counts into this one. Both must be
    /// pre-normalization; used to merge per-thread observed models.
    pub fn combine_counts(&mut self, other: &SBModel) {
        debug_assert!(!self.trained && !other.trained, "combine before normalize");
        // Integer sum of the raw counts (the PRIOR is reintroduced once, in f64,
        // at `normalize`), so the merge is order/thread-count independent.
        for (a, b) in self.probs_fp.iter_mut().zip(&other.probs_fp) {
            *a += *b;
        }
    }
}

/// Reverse-complement a DNA byte slice (ACGT; other bases map to `A`).
///
/// DNA is double-stranded: reading the other strand means walking backwards and
/// swapping each base for its pair. The 3' end of a fragment is sequenced from
/// that strand, so its bias is scored against the reverse complement.
pub(crate) fn revcomp_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'A',
        })
        .collect()
}

/// Minimum transcript abundance to contribute to / be corrected by the bias
/// background (salmon's `minAlpha`).
pub(crate) const MIN_ALPHA: f64 = 1e-8;
/// Minimum reliable CDF mass for a transcript (salmon's `minCDFMass`).
pub(crate) const MIN_CDF_MASS: f64 = 1e-10;
/// Fragment-length sampling stride in the effective-length convolution
/// (salmon's `pdfSampFactor` = `biasSpeedSamp` default).
pub const FLD_SAMP_STRIDE: usize = 5;

/// Linear cumulative fragment-length distribution plus the `[low, high]`
/// fragment-length quantile bounds (0.5% / 99.5%), mirroring the `cdf`,
/// `fldLow`, `fldHigh` salmon computes in `updateEffectiveLengths`.
///
/// The bounds bracket the lengths worth iterating over: the outer 1% of the
/// distribution contributes almost nothing to the convolution but would extend
/// the sweep over a long tail of near-zero-probability lengths.
pub fn fld_cdf_and_bounds(pmf_lin: &[f64]) -> (Vec<f64>, usize, usize) {
    let mut cdf = vec![0.0f64; pmf_lin.len()];
    let mut acc = 0.0;
    let (mut lo, mut hi) = (0usize, 1usize);
    let (mut lb, mut ub) = (false, false);
    for i in 0..pmf_lin.len() {
        acc += pmf_lin[i];
        cdf[i] = acc;
        if !lb && acc >= 0.005 {
            lb = true;
            lo = i;
        }
        if !ub && acc >= 0.995 {
            ub = true;
            hi = i;
        }
    }
    (cdf, lo, hi)
}

/// Per-transcript conditional fragment-length CDF: salmon's
/// `conditionalCDF(x) = (x > cdfMaxArg) ? 1.0 : cdf[x] / cdfMaxVal`, where
/// `cdfMaxArg = min(cdf.len()-1, refLen)` normalizes the FLD to the fragment
/// lengths that fit in this transcript.
///
/// A 500-base transcript cannot host a 700-base fragment, so its length
/// distribution is the global one restricted to what fits, renormalized to sum to
/// 1 — otherwise short transcripts would appear to be missing probability mass.
#[inline]
pub(crate) fn conditional_cdf(cdf: &[f64], cdf_max_arg: usize, cdf_max_val: f64, x: i32) -> f64 {
    if x > cdf_max_arg as i32 {
        1.0
    } else if x <= 0 {
        cdf[0] / cdf_max_val
    } else {
        cdf[x as usize] / cdf_max_val
    }
}

/// Build the expected forward/RC sequence-bias models by sliding the context
/// window over each expressed transcript. Each context is weighted by the
/// transcript's abundance density (`alpha / effLen`) times the conditional FLD
/// mass that can start there (`conditionalCDF(maxFragLen)`), matching salmon's
/// expected-model construction in `updateEffectiveLengths`.
///
/// The "expected" model answers: if fragmentation were indifferent to sequence,
/// which 9-mers would we see at fragment starts, given these transcripts and
/// these abundances? Dividing the observed model by it leaves the protocol's
/// sequence preference alone. Weighting by abundance is essential — a motif that
/// is common in a highly expressed transcript should be expected to be common.
pub fn build_expected<'a, F>(
    num_targets: usize,
    seq_of: F,
    alphas: &[f64],
    eff_lens: &[f64],
    cdf: &[f64],
) -> (SBModel, SBModel)
where
    F: Fn(usize) -> &'a [u8] + Sync,
{
    use rayon::prelude::*;
    let k = CONTEXT_LENGTH;
    let cu = CONTEXT_LEFT as i32;
    // Each expressed transcript contributes independently to the expected
    // forward/RC context counts, an O(refLen) sweep per transcript. salmon
    // parallelizes this over transcripts; do the same with rayon (per-thread
    // `SBModel` partials reduced via `combine_counts`). `seq_of` must be `Sync`
    // to share across threads (it is: a closure over the index). `num_targets`
    // excludes decoys (the contiguous tail): decoys are never expressed and so
    // contribute nothing, but skipping them outright guarantees no O(refLen)
    // decoy sweep can ever run.
    let per_tid = |tid: usize| -> Option<(SBModel, SBModel)> {
        if alphas[tid] < MIN_ALPHA || eff_lens[tid] <= 0.0 {
            return None;
        }
        let seq = seq_of(tid);
        let ref_len = seq.len();
        if ref_len < k {
            return None;
        }
        let cdf_max_arg = (cdf.len() - 1).min(ref_len);
        let cdf_max_val = cdf[cdf_max_arg];
        if cdf_max_val < MIN_CDF_MASS {
            return None;
        }
        let weight = alphas[tid] / eff_lens[tid];
        let rc = revcomp_bytes(seq);
        let mut fw = SBModel::new();
        let mut rc_m = SBModel::new();
        // fragStartPos in 0..(refLen - K) (salmon's loop bound)
        for frag_start in 0..(ref_len - k) {
            let max_frag_len = ref_len as i32 - (frag_start as i32 + cu);
            if max_frag_len >= 0 && (max_frag_len as usize) < ref_len {
                let cdensity = conditional_cdf(cdf, cdf_max_arg, cdf_max_val, max_frag_len);
                let w = weight * cdensity;
                fw.add_context(&seq[frag_start..frag_start + k], false, w);
                rc_m.add_context(&rc[frag_start..frag_start + k], false, w);
            }
        }
        Some((fw, rc_m))
    };
    let (mut exp_fw, mut exp_rc) = (0..num_targets)
        .into_par_iter()
        .fold(
            || (SBModel::new(), SBModel::new()),
            |mut acc, tid| {
                if let Some((fw, rc_m)) = per_tid(tid) {
                    acc.0.combine_counts(&fw);
                    acc.1.combine_counts(&rc_m);
                }
                acc
            },
        )
        .reduce(
            || (SBModel::new(), SBModel::new()),
            |mut a, b| {
                a.0.combine_counts(&b.0);
                a.1.combine_counts(&b.1);
                a
            },
        );
    exp_fw.normalize();
    exp_rc.normalize();
    (exp_fw, exp_rc)
}

/// Bias-corrected effective length of one transcript, matching salmon's
/// `updateEffectiveLengths` (`src/util/SalmonUtils.cpp`).
///
/// `cdf` is the linear cumulative FLD; `fld_low`/`fld_high` the 0.5%/99.5%
/// fragment-length quantiles (from [`fld_cdf_and_bounds`]). `elen` is the
/// transcript's *unbiased* effective length (used for the lower barrier and the
/// `unprocessedLen` guard). `stride` subsamples fragment lengths
/// ([`FLD_SAMP_STRIDE`] matches salmon).
///
/// Per-position 5'/3' bias factors `exp(obsLog − expLog)` are placed at the
/// fragment *read-start* (`fragStart + contextBefore`), the 3' factors reversed
/// to forward fragment-end coordinates, then convolved with the conditional FLD:
/// `effLen = Σ_l flWeight(l) · Σ_s fw[s]·rc[s+l−1]`. The result is floored at
/// `min(elen, max(1, unprocessedLen))` (salmon's lower "barrier"; there is no
/// upper cap, so a strongly-biased transcript's effLen can exceed its length).
#[allow(clippy::too_many_arguments)]
pub fn corrected_effective_length(
    seq: &[u8],
    cdf: &[f64],
    fld_low: usize,
    fld_high: usize,
    obs_fw: &SBModel,
    exp_fw: &SBModel,
    obs_rc: &SBModel,
    exp_rc: &SBModel,
    elen: f64,
    stride: usize,
) -> f64 {
    let k = CONTEXT_LENGTH;
    let cu = CONTEXT_LEFT; // contextBefore(false)
    let ref_len = seq.len();
    let unprocessed = (ref_len as i32 - elen as i32).max(0);
    let cdf_max_arg = (cdf.len() - 1).min(ref_len);
    let cdf_max_val = cdf[cdf_max_arg];
    if ref_len < k || unprocessed <= 0 || cdf_max_val < MIN_CDF_MASS {
        return elen;
    }
    let cond = |x: i32| conditional_cdf(cdf, cdf_max_arg, cdf_max_val, x);

    // Per-position 5' and 3' sequence-bias factors, placed at the read-start.
    let rc_seq = revcomp_bytes(seq);
    let mut fw = vec![1.0f64; ref_len];
    let mut rc = vec![1.0f64; ref_len];
    for frag_start in 0..(ref_len - k) {
        let read_start = frag_start + cu;
        if read_start < ref_len {
            fw[read_start] =
                log_bias(obs_fw, exp_fw, &seq[frag_start..frag_start + k], false).exp();
            rc[read_start] =
                log_bias(obs_rc, exp_rc, &rc_seq[frag_start..frag_start + k], false).exp();
        }
    }
    rc.reverse(); // align RC factors with forward fragment-end coordinates

    // Convolve the bias factors with the conditional FLD over [fld_low, fld_high].
    let stride = stride.max(1) as i32;
    let max_len = (ref_len as i32).min(fld_high as i32 + 1);
    let mut fl = fld_low as i32;
    let mut done = fl >= max_len;
    let sp = if fl > 0 { fl - 1 } else { 0 };
    let mut prev_mass = cond(sp);
    let mut eff = 0.0f64;
    while !done {
        if fl >= max_len {
            done = true;
            fl = max_len - 1;
        }
        let fl_weight = cond(fl) - prev_mass;
        prev_mass = cond(fl);
        let mut mass = 0.0f64;
        let mut kstart = 0i32;
        while kstart < ref_len as i32 - fl {
            let frag_start = kstart as usize;
            let frag_end = (kstart + fl - 1) as usize;
            if frag_end < ref_len {
                mass += fw[frag_start] * rc[frag_end];
            } else {
                break;
            }
            kstart += 1;
        }
        eff += fl_weight * mass;
        fl += stride;
    }

    // Lower barrier (salmon default; no upper cap).
    let offset = (unprocessed as f64).max(1.0);
    eff.max(elen.min(offset))
}

/// Bias-corrected effective length when the per-fragment factor is **separable**
/// as `a[start]·b[end]`. This holds for any combination of sequence and
/// positional bias (each contributes an independent 5′ start factor and 3′ end
/// factor); GC bias is *not* separable (its windowed-GC binning couples start
/// and length) and stays on the scalar convolution.
///
/// "Separable" means a fragment's weight depends on where it starts and where it
/// ends, but not on the two jointly. That is exactly the condition under which
/// the length sweep becomes a cross-correlation and the FFT applies.
///
/// The length sweep `mass(fl) = Σ_k a[k]·b[k+fl-1]` is then the cross-correlation
/// of `a` and `b` evaluated at every lag, computed once via a real FFT in
/// `O(L log L)` instead of `O(L · n_len)`. `a`/`b` must both have length
/// `ref_len`; `cond` is the conditional fragment-length CMF; `unprocessed` is
/// `max(0, ref_len − elen)`. Mirrors the scalar loop's `stride` (biasSpeedSamp)
/// sampling and boundary exclusion exactly, so it is a drop-in replacement.
#[allow(clippy::too_many_arguments)]
pub fn eff_len_from_xcorr(
    a: &[f64],
    b: &[f64],
    cond: impl Fn(i32) -> f64,
    fld_low: usize,
    fld_high: usize,
    elen: f64,
    unprocessed: i32,
    stride: usize,
    no_length_threshold: bool,
) -> f64 {
    let ref_len = a.len();
    debug_assert_eq!(ref_len, b.len());
    let max_len = (ref_len as i32).min(fld_high as i32 + 1);
    if (fld_low as i32) >= max_len {
        let offset = (unprocessed as f64).max(1.0);
        return elen.max(elen.min(offset));
    }
    // Lags needed: Δ = fl-1 for fl in [fld_low, max_len). The scalar inner loop
    // stops at kstart < ref_len-fl (so frag_end ≤ ref_len-2), i.e. it excludes
    // the single fragment ending at ref_len-1; the zero-padded xcorr includes it
    // (term a[ref_len-fl]·b[ref_len-1]), so subtract that for exact parity.
    let max_lag = (max_len - 2).max(0) as usize;
    let xc = xcorr_fft(a, b, max_lag);
    let b_last = b[ref_len - 1];

    // Mirror the scalar loop's `stride` (biasSpeedSamp) sampling EXACTLY so this
    // is a drop-in faster replacement, not an accuracy change: same fragment
    // lengths, same FLD weights. `xc[fl-1] - boundary` is the scalar's inner
    // position sum (the boundary term is the one fragment ending at ref_len-1
    // that the scalar's `kstart < ref_len-fl` bound excludes).
    let stride = stride.max(1) as i32;
    let mut eff = 0.0f64;
    let mut fl = fld_low as i32;
    let mut done = fl >= max_len;
    let sp = if fl > 0 { fl - 1 } else { 0 };
    let mut prev_mass = cond(sp);
    while !done {
        if fl >= max_len {
            done = true;
            fl = max_len - 1;
        }
        let fl_weight = cond(fl) - prev_mass;
        prev_mass = cond(fl);
        if fl >= 1 {
            let delta = (fl - 1) as usize;
            let boundary = a[(ref_len as i32 - fl) as usize] * b_last;
            eff += fl_weight * (xc[delta] - boundary);
        }
        fl += stride;
    }
    if no_length_threshold {
        if eff > 1.0 {
            eff
        } else {
            elen
        }
    } else {
        let offset = (unprocessed as f64).max(1.0);
        eff.max(elen.min(offset))
    }
}

/// FFT form of [`corrected_effective_length`] (sequence-only): builds the 5′/3′
/// sequence factor arrays, then evaluates the length sweep as their
/// cross-correlation via [`eff_len_from_xcorr`]. Kept as the validated
/// seq-only reference (the combined no-GC path in [`crate::bias`] builds the
/// same factors, optionally fused with positional bias, and calls the same
/// core). Numerically identical to [`corrected_effective_length`] up to FFT
/// round-off.
#[allow(clippy::too_many_arguments)]
pub fn corrected_effective_length_fft(
    seq: &[u8],
    cdf: &[f64],
    fld_low: usize,
    fld_high: usize,
    obs_fw: &SBModel,
    exp_fw: &SBModel,
    obs_rc: &SBModel,
    exp_rc: &SBModel,
    elen: f64,
    stride: usize,
    no_length_threshold: bool,
) -> f64 {
    let k = CONTEXT_LENGTH;
    let cu = CONTEXT_LEFT;
    let ref_len = seq.len();
    let unprocessed = (ref_len as i32 - elen as i32).max(0);
    let cdf_max_arg = (cdf.len() - 1).min(ref_len);
    let cdf_max_val = cdf[cdf_max_arg];
    if ref_len < k || unprocessed <= 0 || cdf_max_val < MIN_CDF_MASS {
        return elen;
    }
    let cond = |x: i32| conditional_cdf(cdf, cdf_max_arg, cdf_max_val, x);

    let rc_seq = revcomp_bytes(seq);
    let mut fw = vec![1.0f64; ref_len];
    let mut rc = vec![1.0f64; ref_len];
    for frag_start in 0..(ref_len - k) {
        let read_start = frag_start + cu;
        if read_start < ref_len {
            fw[read_start] =
                log_bias(obs_fw, exp_fw, &seq[frag_start..frag_start + k], false).exp();
            rc[read_start] =
                log_bias(obs_rc, exp_rc, &rc_seq[frag_start..frag_start + k], false).exp();
        }
    }
    rc.reverse();

    eff_len_from_xcorr(
        &fw,
        &rc,
        cond,
        fld_low,
        fld_high,
        elen,
        unprocessed,
        stride,
        no_length_threshold,
    )
}

/// Log bias of `observed` relative to `expected` for a context:
/// `log P_obs(context) - log P_exp(context)`. The fragment-level bias weight is
/// `exp` of this.
///
/// A difference of logs is a ratio of probabilities: how much more often this
/// context was seen at a fragment start than chance predicts. Above 1 means the
/// protocol favours it.
pub fn log_bias(observed: &SBModel, expected: &SBModel, context: &[u8], rev_comp: bool) -> f64 {
    observed.evaluate_log(context, rev_comp) - expected.evaluate_log(context, rev_comp)
}

/// Precomputed `observed − expected` log-transition table for a fixed model
/// pair. The effective-length correction evaluates `log_bias` for a context at
/// every transcript position; `log_bias` evaluates BOTH models (each
/// re-encoding the context and sweeping all `CONTEXT_LENGTH` positions), so a
/// context costs two encodes + two table sweeps. Folding the pair into a single
/// difference table `diff[pos·ROWS+idx] = obs − exp` (built once per quant run,
/// the models being fixed during correction) collapses that to **one** encode +
/// **one** sweep — ~1.4× on the per-position factor build, the dominant cost of
/// the seqBias sweep.
///
/// `eval` equals `log_bias(obs, exp, ctx, rc)` up to floating-point
/// reassociation: it sums `Σ(obs−exp)` rather than `(Σobs) − (Σexp)`, a
/// difference of ~1e-15 per context (machine epsilon), far below quant-output
/// resolution.
pub struct LogBiasTable {
    diff: Vec<f64>,
    shifts: [u32; CONTEXT_LENGTH],
    masks: [u32; CONTEXT_LENGTH],
}

impl LogBiasTable {
    /// Build the difference table from a normalized observed/expected pair.
    pub fn new(observed: &SBModel, expected: &SBModel) -> Self {
        debug_assert!(observed.trained && expected.trained);
        let diff = observed
            .probs
            .iter()
            .zip(&expected.probs)
            .map(|(&o, &e)| o - e)
            .collect();
        Self {
            diff,
            shifts: observed.shifts,
            masks: observed.masks,
        }
    }

    /// `log_bias` for a context (one encode, one table sweep).
    #[inline]
    pub fn eval(&self, context: &[u8], rev_comp: bool) -> f64 {
        let mer = SBModel::encode(context, rev_comp);
        let mut lp = 0.0;
        for pos in 0..CONTEXT_LENGTH {
            let idx = ((mer >> self.shifts[pos]) & self.masks[pos]) as usize;
            lp += self.diff[pos * ROWS + idx];
        }
        lp
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Build a realistic-ish trained obs/exp pair: expected from a sweep of the
    // transcript, observed = expected with a few enriched contexts.
    //
    // Starting the observed model as a copy of the expected one means the tests
    // measure exactly the injected enrichment, with no incidental background
    // difference between the two.
    fn trained_pair(seq: &[u8]) -> (SBModel, SBModel) {
        let rc = revcomp_bytes(seq);
        let mut exp = SBModel::new();
        for p in 0..=(seq.len() - CONTEXT_LENGTH) {
            exp.add_context(&seq[p..p + CONTEXT_LENGTH], false, 1.0);
            exp.add_context(&rc[p..p + CONTEXT_LENGTH], false, 1.0);
        }
        let mut obs = exp.clone();
        for _ in 0..500 {
            obs.add_context(b"AAACCCGGG", false, 1.0);
            obs.add_context(b"TTTGGGCCC", true, 1.0);
        }
        obs.normalize();
        exp.normalize();
        (obs, exp)
    }

    /// A profiling harness rather than a correctness test (hence `#[ignore]`):
    /// it times the successive optimizations of the per-position factor build,
    /// which is the dominant cost of the seqBias sweep.
    #[test]
    #[ignore = "profiling bench; run with --ignored --nocapture"]
    fn bench_factor_build() {
        use std::time::Instant;
        let bases = b"ACGTACGTAGGCCTTAACCGGTTACGTACGTTTAGCGATCG";
        let seq: Vec<u8> = (0..2000)
            .map(|i| bases[(i * 7 + 3) % bases.len()])
            .collect();
        let (obs, exp) = trained_pair(&seq);
        let rc_seq = revcomp_bytes(&seq);
        let k = CONTEXT_LENGTH;
        let n = seq.len() - k;
        let iters = 8000usize; // ~16M positions, real-workload scale

        // V0: current path — log_bias (two evaluate_log, each re-encodes) + exp.
        let t = Instant::now();
        let mut acc = 0.0f64;
        for _ in 0..iters {
            for fs in 0..n {
                acc += log_bias(&obs, &exp, &seq[fs..fs + k], false).exp();
                acc += log_bias(&obs, &exp, &rc_seq[fs..fs + k], false).exp();
            }
        }
        let v0 = t.elapsed().as_secs_f64();

        // No-exp: V0 minus the exp() to isolate exp cost.
        let t = Instant::now();
        let mut acc1 = 0.0f64;
        for _ in 0..iters {
            for fs in 0..n {
                acc1 += log_bias(&obs, &exp, &seq[fs..fs + k], false);
                acc1 += log_bias(&obs, &exp, &rc_seq[fs..fs + k], false);
            }
        }
        let v_noexp = t.elapsed().as_secs_f64();

        // Encode-only: isolate the encode cost (2 encodes per position as today).
        let t = Instant::now();
        let mut enc = 0u64;
        for _ in 0..iters {
            for fs in 0..n {
                enc ^= SBModel::encode(&seq[fs..fs + k], false) as u64;
                enc ^= SBModel::encode(&rc_seq[fs..fs + k], false) as u64;
            }
        }
        let v_enc = t.elapsed().as_secs_f64();

        // V1: encode ONCE per context, evaluate obs and exp from the shared mer
        // (byte-identical: encode is deterministic, sum order unchanged).
        let eval_mer = |m: &SBModel, mer: u32| -> f64 {
            let mut lp = 0.0;
            for pos in 0..CONTEXT_LENGTH {
                lp += m.probs[pos * ROWS + m.index_at(mer, pos)];
            }
            lp
        };
        let t = Instant::now();
        let mut acc_v1 = 0.0f64;
        let mut max_d1 = 0.0f64;
        for it in 0..iters {
            for fs in 0..n {
                let mf = SBModel::encode(&seq[fs..fs + k], false);
                let mr = SBModel::encode(&rc_seq[fs..fs + k], false);
                let bf = (eval_mer(&obs, mf) - eval_mer(&exp, mf)).exp();
                let br = (eval_mer(&obs, mr) - eval_mer(&exp, mr)).exp();
                acc_v1 += bf + br;
                if it == 0 {
                    let rf = log_bias(&obs, &exp, &seq[fs..fs + k], false).exp();
                    let rr = log_bias(&obs, &exp, &rc_seq[fs..fs + k], false).exp();
                    max_d1 = max_d1.max((bf - rf).abs()).max((br - rr).abs());
                }
            }
        }
        let v1 = t.elapsed().as_secs_f64();

        // V2: precomputed diff table d[pos*ROWS+idx] = obs - exp (one eval per
        // direction). NON-byte-identical (Σ(a-b) vs Σa-Σb reassociation).
        let mut diff = vec![0.0f64; ROWS * CONTEXT_LENGTH];
        for (i, d) in diff.iter_mut().enumerate() {
            *d = obs.probs[i] - exp.probs[i];
        }
        let eval_diff = |mer: u32| -> f64 {
            let mut lp = 0.0;
            for pos in 0..CONTEXT_LENGTH {
                lp += diff[pos * ROWS + obs.index_at(mer, pos)];
            }
            lp
        };
        let t = Instant::now();
        let mut acc_v2 = 0.0f64;
        let mut max_d2 = 0.0f64;
        for it in 0..iters {
            for fs in 0..n {
                let bf = eval_diff(SBModel::encode(&seq[fs..fs + k], false)).exp();
                let br = eval_diff(SBModel::encode(&rc_seq[fs..fs + k], false)).exp();
                acc_v2 += bf + br;
                if it == 0 {
                    let rf = log_bias(&obs, &exp, &seq[fs..fs + k], false).exp();
                    let rr = log_bias(&obs, &exp, &rc_seq[fs..fs + k], false).exp();
                    max_d2 = max_d2.max((bf - rf).abs()).max((br - rr).abs());
                }
            }
        }
        let v2 = t.elapsed().as_secs_f64();

        eprintln!("--- factor-build bench ({iters} iters x {n} pos x2) ---");
        eprintln!("V0 current (log_bias+exp)   : {v0:.3}s   acc={acc:.3}");
        eprintln!(
            "  no-exp (log_bias only)    : {v_noexp:.3}s  acc={acc1:.3}  => exp cost ~{:.3}s",
            v0 - v_noexp
        );
        eprintln!("  encode-only (2x/pos)      : {v_enc:.3}s   enc={enc}");
        eprintln!("V1 encode-once (byte-ident) : {v1:.3}s   acc={acc_v1:.3}  max|Δ|={max_d1:.3e}  speedup={:.2}x", v0 / v1);
        eprintln!("V2 diff-table (reassoc)     : {v2:.3}s   acc={acc_v2:.3}  max|Δ|={max_d2:.3e}  speedup={:.2}x", v0 / v2);
    }

    /// No enrichment in, no correction out: identical observed and expected
    /// models must score every context at log-bias ~0 (factor ~1).
    #[test]
    fn uniform_contexts_give_near_zero_bias() {
        // Both models trained on the same uniform set of contexts -> bias ~ 0.
        let ctxs: Vec<Vec<u8>> = (0..256)
            .map(|i| {
                let bases = b"ACGT";
                (0..CONTEXT_LENGTH)
                    .map(|p| bases[((i >> (p * 2)) & 3) as usize])
                    .collect()
            })
            .collect();
        let mut obs = SBModel::new();
        let mut exp = SBModel::new();
        for c in &ctxs {
            obs.add_context(c, false, 1.0);
            exp.add_context(c, false, 1.0);
        }
        obs.normalize();
        exp.normalize();
        for c in &ctxs {
            assert!(log_bias(&obs, &exp, c, false).abs() < 1e-9);
        }
    }

    /// And the converse: a context deliberately over-represented in the observed
    /// model must score positively, so the correction actually responds to bias.
    #[test]
    fn enriched_context_has_positive_bias() {
        // observed enriched for a specific context vs a uniform expected model
        let target: Vec<u8> = b"ACGTACGTA".to_vec();
        let bases = b"ACGT";
        let uniform: Vec<Vec<u8>> = (0..4096)
            .map(|i| {
                (0..CONTEXT_LENGTH)
                    .map(|p| bases[((i >> (p * 2)) & 3) as usize])
                    .collect()
            })
            .collect();

        let mut exp = SBModel::new();
        for c in &uniform {
            exp.add_context(c, false, 1.0);
        }
        exp.normalize();

        let mut obs = SBModel::new();
        for c in &uniform {
            obs.add_context(c, false, 1.0);
        }
        for _ in 0..5000 {
            obs.add_context(&target, false, 1.0); // enrich
        }
        obs.normalize();

        assert!(
            log_bias(&obs, &exp, &target, false) > 0.5,
            "enriched context should have positive log-bias"
        );
    }

    /// With no bias to correct, the corrected effective length must collapse to
    /// the ordinary one — the correction may not shift results merely by being
    /// switched on.
    #[test]
    fn unbiased_correction_reduces_to_standard_eff_len() {
        // obs == exp -> all bias factors 1 -> corrected effLen == standard
        // effLen = sum_l pmf(l)*(refLen - l). Point-mass FLD at l=100.
        let bases = b"ACGTACGTAGGCCTTAACCGGTTACGTACGT";
        let seq: Vec<u8> = (0..400).map(|i| bases[i % bases.len()]).collect();
        let mut m = SBModel::new();
        let rc = revcomp_bytes(&seq);
        for p in 0..=(seq.len() - CONTEXT_LENGTH) {
            m.add_context(&seq[p..p + CONTEXT_LENGTH], false, 1.0);
            m.add_context(&rc[p..p + CONTEXT_LENGTH], false, 1.0);
        }
        let mut obs = m.clone();
        let mut exp = m.clone();
        obs.normalize();
        exp.normalize();

        let mut pmf = vec![0.0; 200];
        pmf[100] = 1.0;
        let (cdf, lo, hi) = fld_cdf_and_bounds(&pmf);
        // unbiased effLen at point-mass 100 on a 400nt transcript = 400 - 100 = 300
        let eff = corrected_effective_length(&seq, &cdf, lo, hi, &obs, &exp, &obs, &exp, 300.0, 1);
        assert!((eff - 300.0).abs() < 1e-6, "got {eff}");
    }

    /// The FFT path exists purely to be faster, so it has to agree with the
    /// direct scalar convolution to within floating-point round-off — including
    /// the fiddly boundary term the scalar loop excludes and the FFT includes.
    #[test]
    fn fft_matches_exact_scalar_corrected_eff_len() {
        // Build a genuinely biased obs/exp pair (so per-position factors != 1),
        // a spread FLD, and check the FFT cross-correlation form equals the exact
        // (stride=1) scalar convolution up to FFT round-off.
        let bases = b"ACGTACGTAGGCCTTAACCGGTTACGTACGTTTAGCGATCG";
        let seq: Vec<u8> = (0..1500)
            .map(|i| bases[(i * 7 + 3) % bases.len()])
            .collect();
        let rc = revcomp_bytes(&seq);
        let mut exp = SBModel::new();
        for p in 0..=(seq.len() - CONTEXT_LENGTH) {
            exp.add_context(&seq[p..p + CONTEXT_LENGTH], false, 1.0);
            exp.add_context(&rc[p..p + CONTEXT_LENGTH], false, 1.0);
        }
        let mut obs = exp.clone();
        // enrich a couple of contexts so obs != exp
        let t1 = b"AAACCCGGG";
        let t2 = b"TTTGGGCCC";
        for _ in 0..500 {
            obs.add_context(t1, false, 1.0);
            obs.add_context(t2, true, 1.0);
        }
        obs.normalize();
        exp.normalize();

        // spread FLD (Gaussian-ish around 250)
        let mut pmf = vec![0.0f64; 600];
        for (l, v) in pmf.iter_mut().enumerate() {
            let d = l as f64 - 250.0;
            *v = (-d * d / (2.0 * 40.0 * 40.0)).exp();
        }
        let (cdf, lo, hi) = fld_cdf_and_bounds(&pmf);

        // FFT must match the scalar at the SAME stride (drop-in, not an accuracy
        // change) — check both the exact (stride=1) and strided (stride=5) cases.
        for stride in [1usize, 5] {
            let scalar = corrected_effective_length(
                &seq, &cdf, lo, hi, &obs, &exp, &obs, &exp, 1200.0, stride,
            );
            let fft = corrected_effective_length_fft(
                &seq, &cdf, lo, hi, &obs, &exp, &obs, &exp, 1200.0, stride, false,
            );
            let rel = (scalar - fft).abs() / scalar.abs();
            assert!(
                rel < 1e-9,
                "FFT vs scalar mismatch at stride={stride}: scalar={scalar} fft={fft} rel={rel:.3e}"
            );
        }
    }

    /// Same obligation for the shared cross-correlation core when sequence and
    /// positional factors are fused into one start/end array pair.
    #[test]
    fn eff_len_from_xcorr_matches_scalar_combined_factors() {
        // The generic core handles ANY separable per-fragment factor
        // a[start]·b[end] — i.e. seq-only, pos-only, or seq+pos (GC is not
        // separable). Validate it against an explicit scalar double-loop on
        // arbitrary positive factor arrays (a stand-in for seqFW·posFW etc.),
        // at both stride 1 and 5, so the pos / seq+pos dispatch in `bias.rs` is
        // covered independently of how the factors were built.
        let ref_len = 1300usize;
        let a: Vec<f64> = (0..ref_len)
            .map(|i| 0.5 + 1.5 * ((i as f64 * 0.013).sin() * 0.5 + 0.5))
            .collect();
        let b: Vec<f64> = (0..ref_len)
            .map(|i| 0.4 + 1.8 * ((i as f64 * 0.021 + 1.0).cos() * 0.5 + 0.5))
            .collect();

        let mut pmf = vec![0.0f64; 600];
        for (l, v) in pmf.iter_mut().enumerate() {
            let d = l as f64 - 250.0;
            *v = (-d * d / (2.0 * 40.0 * 40.0)).exp();
        }
        let (cdf, lo, hi) = fld_cdf_and_bounds(&pmf);
        let elen = 1100.0f64;
        let unprocessed = (ref_len as i32 - elen as i32).max(0);
        let cdf_max_arg = (cdf.len() - 1).min(ref_len);
        let cdf_max_val = cdf[cdf_max_arg];
        let cond = |x: i32| conditional_cdf(&cdf, cdf_max_arg, cdf_max_val, x);

        for stride in [1usize, 5] {
            // Explicit scalar reference: same fragment lengths/weights and the
            // same `kstart < ref_len - fl` bound as the combined scalar loop.
            let max_len = (ref_len as i32).min(hi as i32 + 1);
            let st = stride.max(1) as i32;
            let mut fl = lo as i32;
            let mut done = fl >= max_len;
            let sp = if fl > 0 { fl - 1 } else { 0 };
            let mut prev = cond(sp);
            let mut eff = 0.0f64;
            while !done {
                if fl >= max_len {
                    done = true;
                    fl = max_len - 1;
                }
                let w = cond(fl) - prev;
                prev = cond(fl);
                let kmax = ref_len as i32 - fl;
                let mut mass = 0.0f64;
                let mut k = 0i32;
                while k < kmax {
                    mass += a[k as usize] * b[(k + fl - 1) as usize];
                    k += 1;
                }
                eff += w * mass;
                fl += st;
            }
            let offset = (unprocessed as f64).max(1.0);
            let scalar = eff.max(elen.min(offset));

            let fft = eff_len_from_xcorr(&a, &b, cond, lo, hi, elen, unprocessed, stride, false);
            let rel = (scalar - fft).abs() / scalar.abs();
            assert!(
                rel < 1e-9,
                "combined-factor FFT vs scalar mismatch at stride={stride}: scalar={scalar} fft={fft} rel={rel:.3e}"
            );
        }
    }

    /// The 3' factor is scored against the reverse complement, so encoding a
    /// context reverse-complemented must equal encoding its reverse complement
    /// directly; otherwise the two ends would be scored against different models.
    #[test]
    fn revcomp_encoding_is_consistent() {
        // RC of a context evaluated forward equals the context evaluated as RC.
        let ctx: Vec<u8> = b"ACGTACGTA".to_vec();
        let rc: Vec<u8> = ctx
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => x,
            })
            .collect();
        assert_eq!(SBModel::encode(&ctx, true), SBModel::encode(&rc, false));
    }

    /// Decoys sit past `num_targets` and must never enter the expected model: a
    /// genome decoy swept as a transcript would dominate the background.
    #[test]
    fn build_expected_respects_num_targets_bound() {
        // Five real transcripts plus a sixth "decoy" with a very distinctive
        // composition (poly-AC). The decoy must influence the expected model only
        // when `num_targets` includes it, and must be skipped (cheaply) when its
        // alpha is zero — the two ways decoys are kept out of the bias models.
        let bases = b"ACGTACGTAGGCCTTAACCGGTTACGTACGT";
        let mut refs: Vec<Vec<u8>> = (0..5)
            .map(|s| (0..200).map(|i| bases[(i + s) % bases.len()]).collect())
            .collect();
        refs.push(
            (0..400)
                .map(|i| if i % 2 == 0 { b'A' } else { b'C' })
                .collect(),
        );
        let num_refs = refs.len();
        let alphas = vec![1.0; num_refs];
        let eff_lens = vec![150.0; num_refs];
        let mut pmf = vec![0.0; 200];
        pmf[100] = 1.0;
        let (cdf, _lo, _hi) = fld_cdf_and_bounds(&pmf);

        // Exclude the decoy (num_targets = 5) vs include it (num_targets = 6).
        let (a_fw, _) = build_expected(5, |t| refs[t].as_slice(), &alphas, &eff_lens, &cdf);
        let (b_fw, _) = build_expected(6, |t| refs[t].as_slice(), &alphas, &eff_lens, &cdf);
        assert!(a_fw.is_trained() && b_fw.is_trained());
        assert!(a_fw.dump().iter().all(|v| v.is_finite()));
        let diff: f64 = a_fw
            .dump()
            .iter()
            .zip(b_fw.dump())
            .map(|(x, y)| (x - y).abs())
            .sum();
        assert!(
            diff > 1e-6,
            "a target beyond num_targets must not contribute (diff={diff})"
        );

        // With the decoy's alpha zeroed the MIN_ALPHA guard skips it, so including
        // it (num_targets = 6) must match excluding it (num_targets = 5).
        let mut alphas0 = alphas.clone();
        alphas0[5] = 0.0;
        let (c_fw, _) = build_expected(6, |t| refs[t].as_slice(), &alphas0, &eff_lens, &cdf);
        let diff2: f64 = a_fw
            .dump()
            .iter()
            .zip(c_fw.dump())
            .map(|(x, y)| (x - y).abs())
            .sum();
        assert!(
            diff2 < 1e-9,
            "zero-alpha target must not contribute (diff={diff2})"
        );
    }
}
