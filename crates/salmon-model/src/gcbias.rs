//! Fragment-GC bias model (`GCFragModel`).
//!
//! A port of salmon's `GCFragModel` (`include/.../model/GCFragModel.hpp`): a
//! `condBins × numGCBins` table of fragment counts, where a fragment is keyed by
//! its **GC fraction** (0–100, → `numGCBins` bins) and a **conditioning context**
//! (a coarse sequence/context fraction, → `condBins` bins). An *observed* model
//! (from mapped fragments) and an *expected* model (from the transcriptome) are
//! each normalized per conditioning row into a distribution, then divided —
//! [`gc_ratio`] — to give the per-`(context, GC)` bias factor used to weight
//! fragments in the bias-corrected effective length.
//!
//! salmon accumulates the observed model in log space and the expected in linear
//! space, but both are normalized to linear distributions before the ratio is
//! taken, so the normalized result is independent of the accumulation space.
//! We therefore accumulate in **linear** space (summing weights) throughout —
//! numerically equivalent after [`GcFragModel::normalize`], and simpler.

/// salmon's default number of conditioning (context) bins (`numConditionalGCBins`).
pub const DEFAULT_COND_BINS: usize = 3;
/// salmon's default number of fragment-GC bins for accumulation
/// (`salmon::defaults::numFragGCBins`).
pub const DEFAULT_GC_BINS: usize = 25;
/// salmon's `ratio()` clamp (`gcBias = gcCounts.ratio(transcriptGCDist, 1000.0)`).
pub const GC_MAX_RATIO: f64 = 1000.0;
/// salmon's per-row normalization prior.
const NORM_PRIOR: f64 = 0.1;

/// Bin a 0–100 fraction into `n` bins (salmon's `GCDesc::fragBin(n)`/`contextBin(n)`:
/// `min(n-1, floor(frac / (100/n)))`). With `n == 101` this is the identity.
#[inline]
pub fn bin_frac(frac: i32, n: usize) -> usize {
    if n == 101 {
        return (frac.clamp(0, 100)) as usize;
    }
    let w = 100.0 / n as f64;
    let b = (frac as f64 / w) as i32;
    b.clamp(0, n as i32 - 1) as usize
}

/// `(condBins × numGCBins)` table of fragment-GC counts (then a normalized
/// distribution / ratio after [`normalize`](Self::normalize) / [`gc_ratio`]).
#[derive(Debug, Clone)]
pub struct GcFragModel {
    cond_bins: usize,
    gc_bins: usize,
    /// row-major `counts[ctx * gc_bins + gc]`
    counts: Vec<f64>,
    normalized: bool,
}

impl GcFragModel {
    pub fn new(cond_bins: usize, gc_bins: usize) -> Self {
        Self {
            cond_bins,
            gc_bins,
            counts: vec![0.0; cond_bins * gc_bins],
            normalized: false,
        }
    }

    /// salmon's default 3 × 101 model.
    pub fn default_model() -> Self {
        Self::new(DEFAULT_COND_BINS, DEFAULT_GC_BINS)
    }

    #[inline]
    fn idx(&self, ctx_frac: i32, gc_frac: i32) -> usize {
        let ctx = if self.cond_bins > 1 {
            bin_frac(ctx_frac, self.cond_bins)
        } else {
            0
        };
        let gc = bin_frac(gc_frac, self.gc_bins);
        ctx * self.gc_bins + gc
    }

    /// The flattened `cond_bins × gc_bins` table (`counts[ctx * gc_bins + gc]`,
    /// row-major), for dumping to the aux bias files.
    pub fn dump(&self) -> &[f64] {
        &self.counts
    }

    /// Accumulate `weight` for a fragment with the given GC and context fractions.
    pub fn inc(&mut self, gc_frac: i32, ctx_frac: i32, weight: f64) {
        debug_assert!(!self.normalized, "cannot inc a normalized model");
        let i = self.idx(ctx_frac, gc_frac);
        self.counts[i] += weight;
    }

    /// Value at a `(context, GC)` cell (a normalized density after `normalize`,
    /// or a clamped ratio for a model produced by [`gc_ratio`]).
    #[inline]
    pub fn get(&self, gc_frac: i32, ctx_frac: i32) -> f64 {
        self.counts[self.idx(ctx_frac, gc_frac)]
    }

    /// Merge another (compatible) model's counts. Both must be pre-normalization.
    pub fn combine_counts(&mut self, other: &GcFragModel) {
        debug_assert!(!self.normalized && !other.normalized);
        for (a, b) in self.counts.iter_mut().zip(&other.counts) {
            *a += *b;
        }
    }

    /// Normalize each conditioning row into a distribution over GC bins, with a
    /// pseudocount `prior` (salmon's default 0.1). Idempotent.
    pub fn normalize(&mut self) {
        if self.normalized {
            return;
        }
        for r in 0..self.cond_bins {
            let base = r * self.gc_bins;
            let row = &mut self.counts[base..base + self.gc_bins];
            let row_mass: f64 = row.iter().map(|c| NORM_PRIOR + *c).sum();
            if row_mass > 0.0 {
                let norm = 1.0 / row_mass;
                for c in row.iter_mut() {
                    *c = (NORM_PRIOR + *c) * norm;
                }
            }
        }
        self.normalized = true;
    }
}

/// Per-cell bias ratio `observed / expected`, clamped to `[1/max_ratio, max_ratio]`
/// (salmon's `GCFragModel::ratio`). Both models are normalized first.
pub fn gc_ratio(observed: &mut GcFragModel, expected: &mut GcFragModel, max_ratio: f64) -> GcFragModel {
    observed.normalize();
    expected.normalize();
    let min_ratio = 1.0 / max_ratio;
    let mut out = GcFragModel::new(observed.cond_bins, observed.gc_bins);
    for i in 0..out.counts.len() {
        let e = expected.counts[i];
        let rat = if e != 0.0 { observed.counts[i] / e } else { max_ratio };
        out.counts[i] = rat.clamp(min_ratio, max_ratio);
    }
    out.normalized = true; // a ratio model is used directly, not re-normalized
    out
}

// ===========================================================================
// Fragment-GC content (salmon's `Transcript::GCCount_` / `gcFrac` / `gcDesc`)
// ===========================================================================

use crate::seqbias::{
    conditional_cdf, log_bias, revcomp_bytes, SBModel, CONTEXT_LEFT, CONTEXT_LENGTH, MIN_ALPHA,
    MIN_CDF_MASS,
};

/// salmon's fragment-length sampling stride for the GC convolution
/// (`pdfSampFactor` = `biasSpeedSamp`, default 5).
pub const GC_SAMP_STRIDE: usize = 5;

/// `gcDesc` context-window geometry (salmon `Transcript::gcDesc`):
/// `outsideContext = 3`, `insideContext = 2`.
const OUTSIDE_5P: i32 = 4; // outsideContext + 1
const OUTSIDE_3P: i32 = 3; // outsideContext
const INSIDE_5P: i32 = 1; // insideContext - 1
const INSIDE_3P: i32 = 2; // insideContext

/// Round half-to-even, matching C++ `std::lrint` under the default rounding mode.
#[inline]
fn lrint(x: f64) -> i32 {
    x.round_ties_even() as i32
}

/// Cumulative G+C counts (salmon's `Transcript::GCCount_`): `prefix[p]` is the
/// number of `G`/`C` bases in `seq[0..=p]`.
pub fn gc_prefix(seq: &[u8]) -> Vec<u32> {
    let mut prefix = Vec::with_capacity(seq.len());
    let mut acc = 0u32;
    for &b in seq {
        if matches!(b, b'G' | b'g' | b'C' | b'c') {
            acc += 1;
        }
        prefix.push(acc);
    }
    prefix
}

/// GC count in the closed interval `[a, b]` from a cumulative-count prefix.
#[inline]
fn gc_in(prefix: &[u32], a: i32, b: i32) -> i64 {
    let lo = if a > 0 { prefix[(a - 1) as usize] as i64 } else { 0 };
    prefix[b as usize] as i64 - lo
}

/// Fragment GC fraction (0–100) over the closed interval `[s, e]`
/// (salmon's `Transcript::gcFrac`): `round(100·GC[s,e] / (e−s+1))`.
#[inline]
pub fn gc_frac(prefix: &[u32], s: i32, e: i32) -> i32 {
    lrint(100.0 * gc_in(prefix, s, e) as f64 / (e - s + 1) as f64)
}

/// Fragment GC descriptor `(fragFrac, contextFrac)` for the closed fragment
/// `[s, e]` (salmon's `Transcript::gcDesc`). The context fraction is the GC
/// content over a 5-base 5' window `[s−3, s+1]` and a 5-base 3' window
/// `[e−1, e+3]` (edge-clamped). Returns `None` when the context window is empty
/// (matching salmon's `valid = false`).
#[inline]
pub fn gc_desc(prefix: &[u32], ref_len: i32, s: i32, e: i32) -> Option<(i32, i32)> {
    let last = ref_len - 1;
    let cs = if s > 0 { prefix[(s - 1) as usize] as i64 } else { 0 };
    let ce = prefix[e as usize] as i64;

    let fs = s - OUTSIDE_5P;
    let fe = s + INSIDE_5P;
    let ts = e - INSIDE_3P;
    let te = e + OUTSIDE_3P;

    let fp_left = fs >= 0;
    let fp_right = fe <= last;
    let tp_left = ts >= 0;
    let tp_right = te <= last;

    let fps = if fp_left { prefix[fs as usize] as i64 } else { 0 };
    let fpe = if fp_right { prefix[fe as usize] as i64 } else { ce };
    let tps = if tp_left { prefix[ts as usize] as i64 } else { 0 };
    let tpe = if tp_right { prefix[te as usize] as i64 } else { ce };

    let fs_c = fs.max(0);
    let fe_c = fe.min(last);
    let ts_c = ts.max(0);
    let te_c = te.min(last);
    let fp_context_size = if !fp_left { fe_c + 1 } else { fe_c - fs_c };
    let tp_context_size = if !tp_left { te_c + 1 } else { te_c - ts_c };
    let context_size = (fp_context_size + tp_context_size) as f64;
    if context_size == 0.0 {
        return None;
    }

    let frag_frac = lrint(100.0 * (ce - cs) as f64 / (e - s + 1) as f64);
    let context_frac = lrint(100.0 * ((fpe - fps) + (tpe - tps)) as f64 / context_size);
    Some((frag_frac, context_frac))
}

/// Build the *expected* fragment-GC model (salmon's `transcriptGCDist`): slide
/// every fragment `[fragStart, fragStart+fl−1]` over each expressed transcript,
/// weighting each by `(alpha/effLen)·(conditionalCDF(fl) − prevMass)` — the same
/// abundance-density × conditional-FLD weighting used for the expected
/// sequence-bias model. `k` is the leading offset salmon excludes from the
/// fragment-start loop (`9` with `--seqBias`, `1` otherwise); `stride` subsamples
/// fragment lengths ([`GC_SAMP_STRIDE`]).
#[allow(clippy::too_many_arguments)]
pub fn build_expected_gc<'a, FS, FP>(
    num_refs: usize,
    seq_of: FS,
    prefix_of: FP,
    alphas: &[f64],
    eff_lens: &[f64],
    cdf: &[f64],
    fld_low: usize,
    fld_high: usize,
    cond_bins: usize,
    gc_bins: usize,
    k: usize,
    stride: usize,
) -> GcFragModel
where
    FS: Fn(usize) -> &'a [u8] + Sync,
    FP: Fn(usize) -> &'a [u32] + Sync,
{
    let stride = stride.max(1) as i32;
    // The expected-GC distribution is a sum of independent per-transcript
    // contributions, each an O(refLen · fldRange/stride) double loop — billions
    // of `gc_desc` calls in total. salmon parallelizes this over transcripts;
    // do the same with rayon (per-thread `GcFragModel` partials, reduced via
    // `combine_counts`). `seq_of`/`prefix_of` must be `Sync` to share across
    // threads (they are: closures over `&[Vec<_>]`).
    use rayon::prelude::*;
    let per_tid = |tid: usize| -> Option<GcFragModel> {
        if alphas[tid] < MIN_ALPHA || eff_lens[tid] <= 0.0 {
            return None;
        }
        let seq = seq_of(tid);
        let ref_len = seq.len();
        if ref_len <= k {
            return None;
        }
        let cdf_max_arg = (cdf.len() - 1).min(ref_len);
        let cdf_max_val = cdf[cdf_max_arg];
        if cdf_max_val < MIN_CDF_MASS {
            return None;
        }
        let prefix = prefix_of(tid);
        let weight = alphas[tid] / eff_lens[tid];
        let cond = |x: i32| conditional_cdf(cdf, cdf_max_arg, cdf_max_val, x);
        let sp = if fld_low > 0 { fld_low as i32 - 1 } else { 0 };
        let mut model = GcFragModel::new(cond_bins, gc_bins);
        for frag_start in 0..(ref_len - k) {
            let mut prev = cond(sp);
            let mut fl = fld_low as i32;
            while fl <= fld_high as i32 {
                let frag_end = frag_start as i32 + fl - 1;
                if (frag_end as usize) < ref_len {
                    if let Some((ff, cf)) = gc_desc(prefix, ref_len as i32, frag_start as i32, frag_end)
                    {
                        model.inc(ff, cf, weight * (cond(fl) - prev));
                    }
                    prev = cond(fl);
                } else {
                    break;
                }
                fl += stride;
            }
        }
        Some(model)
    };
    (0..num_refs)
        .into_par_iter()
        .fold(
            || GcFragModel::new(cond_bins, gc_bins),
            |mut acc, tid| {
                if let Some(m) = per_tid(tid) {
                    acc.combine_counts(&m);
                }
                acc
            },
        )
        .reduce(
            || GcFragModel::new(cond_bins, gc_bins),
            |mut a, b| {
                a.combine_counts(&b);
                a
            },
        )
}

/// Bias-corrected effective length including fragment-GC bias (and, when
/// `seq_models` is provided, sequence-specific bias too). Mirrors salmon's
/// combined `updateEffectiveLengths` convolution: per-position 5'/3' sequence
/// factors are multiplied by the per-fragment `gcBias.get({fragFrac, contextFrac})`
/// ratio, convolved with the conditional FLD, and floored at the lower barrier.
///
/// `gc_bias` is the normalized observed/expected ratio model ([`gc_ratio`]).
/// `seq_models` is `(obs_fw, exp_fw, obs_rc, exp_rc)` when `--seqBias` is also on.
#[allow(clippy::too_many_arguments)]
pub fn gc_corrected_effective_length(
    seq: &[u8],
    prefix: &[u32],
    cdf: &[f64],
    fld_low: usize,
    fld_high: usize,
    gc_bias: &GcFragModel,
    seq_models: Option<(&SBModel, &SBModel, &SBModel, &SBModel)>,
    elen: f64,
    stride: usize,
) -> f64 {
    let k = if seq_models.is_some() { CONTEXT_LENGTH } else { 1 };
    let ref_len = seq.len();
    let unprocessed = (ref_len as i32 - elen as i32).max(0);
    let cdf_max_arg = (cdf.len() - 1).min(ref_len);
    let cdf_max_val = cdf[cdf_max_arg];
    if ref_len < k || unprocessed <= 0 || cdf_max_val < MIN_CDF_MASS {
        return elen;
    }
    let cond = |x: i32| conditional_cdf(cdf, cdf_max_arg, cdf_max_val, x);

    // Per-position sequence-bias factors (1.0 when not seq-correcting).
    let mut fw = vec![1.0f64; ref_len];
    let mut rc = vec![1.0f64; ref_len];
    if let Some((obs_fw, exp_fw, obs_rc, exp_rc)) = seq_models {
        let cu = CONTEXT_LEFT;
        let rc_seq = revcomp_bytes(seq);
        for frag_start in 0..(ref_len - CONTEXT_LENGTH) {
            let read_start = frag_start + cu;
            if read_start < ref_len {
                fw[read_start] =
                    log_bias(obs_fw, exp_fw, &seq[frag_start..frag_start + CONTEXT_LENGTH], false)
                        .exp();
                rc[read_start] = log_bias(
                    obs_rc,
                    exp_rc,
                    &rc_seq[frag_start..frag_start + CONTEXT_LENGTH],
                    false,
                )
                .exp();
            }
        }
        rc.reverse();
    }

    // Convolve seq×gc fragment factors with the conditional FLD.
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
            let frag_start = kstart;
            let frag_end = kstart + fl - 1;
            if (frag_end as usize) < ref_len {
                let mut frag_factor = fw[frag_start as usize] * rc[frag_end as usize];
                if let Some((ff, cf)) = gc_desc(prefix, ref_len as i32, frag_start, frag_end) {
                    frag_factor *= gc_bias.get(ff, cf);
                }
                mass += frag_factor;
            } else {
                break;
            }
            kstart += 1;
        }
        eff += fl_weight * mass;
        fl += stride;
    }

    let offset = (unprocessed as f64).max(1.0);
    eff.max(elen.min(offset))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bin_frac_identity_at_101() {
        assert_eq!(bin_frac(0, 101), 0);
        assert_eq!(bin_frac(57, 101), 57);
        assert_eq!(bin_frac(100, 101), 100);
    }

    #[test]
    fn bin_frac_coarse() {
        // 3 bins: width 33.33 -> [0,33),[33,66),[66,100]
        assert_eq!(bin_frac(0, 3), 0);
        assert_eq!(bin_frac(30, 3), 0);
        assert_eq!(bin_frac(40, 3), 1);
        assert_eq!(bin_frac(70, 3), 2);
        assert_eq!(bin_frac(100, 3), 2); // clamped to n-1
    }

    #[test]
    fn normalize_makes_rows_sum_to_one() {
        let mut m = GcFragModel::new(2, 5);
        for gc in 0..5 {
            m.inc(gc * 25, 10, (gc + 1) as f64); // row 0 (ctx 10 -> bin 0)
        }
        m.normalize();
        let row0: f64 = (0..5).map(|c| m.counts[c]).sum();
        assert!((row0 - 1.0).abs() < 1e-9, "row0 sum = {row0}");
    }

    #[test]
    fn ratio_of_identical_models_is_one() {
        let mut obs = GcFragModel::new(1, 101);
        let mut exp = GcFragModel::new(1, 101);
        for gc in 0..=100 {
            obs.inc(gc, 0, (gc + 1) as f64);
            exp.inc(gc, 0, (gc + 1) as f64);
        }
        let r = gc_ratio(&mut obs, &mut exp, GC_MAX_RATIO);
        for gc in 0..=100 {
            assert!((r.get(gc, 0) - 1.0).abs() < 1e-9, "gc {gc}: {}", r.get(gc, 0));
        }
    }

    #[test]
    fn gc_frac_basic() {
        // 10-base sequence, 5 G/C -> 50%
        let seq = b"ACGTACGTAC"; // A C G T A C G T A C : GC at 1,2,5,6,9 = 5
        let p = gc_prefix(seq);
        assert_eq!(gc_frac(&p, 0, 9), 50);
        // all-GC window
        let seq2 = b"GGGGCCCC";
        let p2 = gc_prefix(seq2);
        assert_eq!(gc_frac(&p2, 0, 7), 100);
        // single base
        assert_eq!(gc_frac(&p, 2, 2), 100); // 'G'
        assert_eq!(gc_frac(&p, 0, 0), 0); // 'A'
    }

    #[test]
    fn gc_desc_matches_frag_and_context() {
        // 40-base sequence; check an interior fragment
        let seq: Vec<u8> = b"ACGT".iter().cycle().take(40).copied().collect();
        let p = gc_prefix(&seq);
        let (ff, cf) = gc_desc(&p, 40, 10, 29).unwrap();
        // ACGT repeat is exactly 50% GC everywhere
        assert_eq!(ff, 50, "fragFrac");
        assert!((0..=100).contains(&cf), "contextFrac in range: {cf}");
        // fragFrac must equal gc_frac over the same interval
        assert_eq!(ff, gc_frac(&p, 10, 29));
    }

    #[test]
    fn gc_desc_context_window_geometry() {
        // Make the 5' context window [s-3, s+1] all GC and the rest AT, so the
        // contextFrac is dominated by the GC island near the fragment ends.
        let mut seq = vec![b'A'; 60];
        for i in 17..=21 {
            seq[i] = b'G'; // 5' window for s=20 is [17,21]
        }
        for i in 39..=43 {
            seq[i] = b'C'; // 3' window for e=40 is [39,43]
        }
        let p = gc_prefix(&seq);
        let (_ff, cf) = gc_desc(&p, 60, 20, 40).unwrap();
        // both 5-base context windows are fully G/C -> 100%
        assert_eq!(cf, 100, "contextFrac");
    }

    #[test]
    fn ratio_clamped() {
        let mut obs = GcFragModel::new(1, 2);
        let mut exp = GcFragModel::new(1, 2);
        obs.inc(0, 0, 1e6); // bin 0 hugely enriched in obs
        exp.inc(75, 0, 1e6); // bin 1 enriched in exp
        let r = gc_ratio(&mut obs, &mut exp, 1000.0);
        for v in [r.get(0, 0), r.get(75, 0)] {
            assert!(v <= 1000.0 + 1e-6 && v >= 1.0 / 1000.0 - 1e-12, "ratio {v} out of clamp");
        }
    }
}
