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
/// salmon's default number of fragment-GC bins for accumulation.
pub const DEFAULT_GC_BINS: usize = 101;
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
