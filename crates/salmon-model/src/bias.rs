//! Unified bias-corrected effective length: composes sequence-specific, GC, and
//! positional bias exactly as salmon's `updateEffectiveLengths` does — a single
//! conditional-FLD convolution whose per-fragment factor is the product of the
//! enabled bias terms.
//!
//! `fragFactor = seqFW[fragStart]·seqRC[fragEnd] · gcBias({fragFrac, ctxFrac}) ·
//! posFW[fragStart]·posRC[fragEnd]`, summed over fragment starts and convolved
//! with the conditional fragment-length distribution.

use crate::gcbias::{GcContext, GcFragModel};
use crate::posbias::{length_class_index, SimplePosBias, NUM_LENGTH_CLASSES};
use crate::seqbias::{
    conditional_cdf, log_bias, revcomp_bytes, SBModel, CONTEXT_LEFT, CONTEXT_LENGTH, MIN_ALPHA,
    MIN_CDF_MASS,
};

/// salmon's `EPSILON` (mass cutoff for adding positional expected mass).
const EPSILON: f64 = 0.375e-10;

/// Additive (Laplace) smoothing fraction for the positional-bias factor: the
/// smoothing constant is `POS_SMOOTH_FRAC · mean(expected density)`.
const POS_SMOOTH_FRAC: f64 = 0.1;

/// Per-position positional-bias factors `obs/exp`, additively smoothed toward 1.
///
/// Instead of dividing two projected densities that the old `0.001` floor drove
/// to a tiny, noisy denominator in the tails (amplifying small model
/// differences), we add a smoothing constant `c = POS_SMOOTH_FRAC·mean(exp)` to
/// both: `factor = (obs + c) / (exp + c)`. Where the expected density is
/// substantial the ratio is preserved; where it vanishes (uninformative tails)
/// the factor shrinks to 1 (no bias), which is the correct default.
pub fn positional_factor(obs: &[f64], exp: &[f64]) -> Vec<f64> {
    let n = exp.len();
    if n == 0 {
        return Vec::new();
    }
    let mean_exp: f64 = exp.iter().sum::<f64>() / n as f64;
    let c = (POS_SMOOTH_FRAC * mean_exp).max(f64::MIN_POSITIVE);
    obs.iter()
        .zip(exp)
        .map(|(&o, &e)| (o + c) / (e + c))
        .collect()
}

/// The enabled bias terms for one transcript's effective-length correction.
#[derive(Clone, Copy, Default)]
pub struct BiasInputs<'a> {
    /// `(obs_fw, exp_fw, obs_rc, exp_rc)` sequence-bias models (`--seqBias`)
    pub seq: Option<(&'a SBModel, &'a SBModel, &'a SBModel, &'a SBModel)>,
    /// GC ratio model + this transcript's cumulative G+C prefix (`--gcBias`)
    pub gc: Option<(&'a GcFragModel, &'a [u32])>,
    /// per-position 5'/3' positional-bias factors (`--posBias`), transcript-length sized
    pub pos: Option<(&'a [f64], &'a [f64])>,
}

impl BiasInputs<'_> {
    fn any(&self) -> bool {
        self.seq.is_some() || self.gc.is_some() || self.pos.is_some()
    }
}

/// Bias-corrected effective length composing every enabled bias term, matching
/// salmon's combined `updateEffectiveLengths` convolution. Floors at the lower
/// barrier `min(elen, max(1, unprocessedLen))` (no upper cap).
pub fn corrected_effective_length_full(
    seq: &[u8],
    cdf: &[f64],
    fld_low: usize,
    fld_high: usize,
    bias: &BiasInputs,
    elen: f64,
    stride: usize,
) -> f64 {
    if !bias.any() {
        return elen;
    }
    let k = if bias.seq.is_some() { CONTEXT_LENGTH } else { 1 };
    let ref_len = seq.len();
    let unprocessed = (ref_len as i32 - elen as i32).max(0);
    let cdf_max_arg = (cdf.len() - 1).min(ref_len);
    let cdf_max_val = cdf[cdf_max_arg];
    if ref_len < k || unprocessed <= 0 || cdf_max_val < MIN_CDF_MASS {
        return elen;
    }
    let cond = |x: i32| conditional_cdf(cdf, cdf_max_arg, cdf_max_val, x);

    // Per-position sequence-bias factors. Only built (and applied in the inner
    // loop) when `--seqBias` is on; otherwise the factors are all 1.0, so we
    // skip both the per-transcript allocation and the per-fragment multiply
    // entirely (the common `--gcBias`-only case).
    let have_seq = bias.seq.is_some();
    let mut fw: Vec<f64> = Vec::new();
    let mut rc: Vec<f64> = Vec::new();
    if let Some((obs_fw, exp_fw, obs_rc, exp_rc)) = bias.seq {
        fw = vec![1.0f64; ref_len];
        rc = vec![1.0f64; ref_len];
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

    let (gc_model, gc_prefix) = match bias.gc {
        Some((m, p)) => (Some(m), p),
        None => (None, &[][..]),
    };
    // Precompute the per-position 5'/3' context-GC arrays once per transcript
    // (salmon's `populateContextCounts`) so the inner convolution does cheap
    // array lookups instead of re-deriving the context geometry per fragment.
    let gc_ctx = gc_model.map(|_| GcContext::build(gc_prefix));
    let (pos_fw, pos_rc) = match bias.pos {
        Some((a, b)) => (Some(a), Some(b)),
        None => (None, None),
    };

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
        // Hoist the bound: for kstart in [0, kmax) we have
        // frag_end = kstart+fl-1 <= ref_len-2 < ref_len, so the old per-iteration
        // `frag_end < ref_len` guard is always true and is dropped (it kept
        // `ref_len` live in the inner loop, forcing a spill/reload — the single
        // hottest source line in `perf annotate`).
        let kmax = ref_len as i32 - fl;
        // Dispatch the bias combination ONCE per fragment length rather than
        // re-testing every bias model's presence per fragment: the common
        // `--gcBias`-only case gets a tight loop with no per-fragment branches.
        match (have_seq, gc_model.zip(gc_ctx.as_ref()), pos_fw.zip(pos_rc)) {
            (false, Some((gc, ctx)), None) => {
                let mut kstart = 0i32;
                while kstart < kmax {
                    let frag_end = kstart + fl - 1;
                    if let Some((ff, cf)) = ctx.desc(gc_prefix, kstart, frag_end) {
                        mass += gc.get(ff, cf);
                    } else {
                        mass += 1.0;
                    }
                    kstart += 1;
                }
            }
            _ => {
                let mut kstart = 0i32;
                while kstart < kmax {
                    let frag_start = kstart;
                    let frag_end = kstart + fl - 1;
                    let mut frag_factor = if have_seq {
                        fw[frag_start as usize] * rc[frag_end as usize]
                    } else {
                        1.0
                    };
                    if let (Some(gc), Some(ctx)) = (gc_model, gc_ctx.as_ref()) {
                        if let Some((ff, cf)) = ctx.desc(gc_prefix, frag_start, frag_end) {
                            frag_factor *= gc.get(ff, cf);
                        }
                    }
                    if let (Some(pf), Some(pr)) = (pos_fw, pos_rc) {
                        frag_factor *= pf[frag_start as usize] * pr[frag_end as usize];
                    }
                    mass += frag_factor;
                    kstart += 1;
                }
            }
        }
        eff += fl_weight * mass;
        fl += stride;
    }

    let offset = (unprocessed as f64).max(1.0);
    eff.max(elen.min(offset))
}

/// Build the *expected* positional-bias models (5'/3', one per length class),
/// mirroring salmon's expected-pos accumulation in `updateEffectiveLengths`:
/// for each expressed transcript and fragment start, add `log(weight·density)`
/// to the length-class bin (forward density = fragments that can start here,
/// reverse density = fragments that can end here). Models are finalized.
#[allow(clippy::too_many_arguments)]
pub fn build_expected_pos<'a, FL>(
    num_refs: usize,
    ref_len_of: FL,
    alphas: &[f64],
    eff_lens: &[f64],
    cdf: &[f64],
    quantiles: &[u32],
    k: usize,
) -> (Vec<SimplePosBias>, Vec<SimplePosBias>)
where
    FL: Fn(usize) -> usize,
{
    let mut pos5: Vec<SimplePosBias> = (0..NUM_LENGTH_CLASSES).map(|_| SimplePosBias::default()).collect();
    let mut pos3: Vec<SimplePosBias> = (0..NUM_LENGTH_CLASSES).map(|_| SimplePosBias::default()).collect();
    for tid in 0..num_refs {
        if alphas[tid] < MIN_ALPHA || eff_lens[tid] <= 0.0 {
            continue;
        }
        let ref_len = ref_len_of(tid) as i32;
        if (ref_len as usize) <= k {
            continue;
        }
        let unprocessed = ref_len - eff_lens[tid] as i32;
        if unprocessed <= 0 {
            continue;
        }
        let cdf_max_arg = (cdf.len() - 1).min(ref_len as usize);
        let cdf_max_val = cdf[cdf_max_arg];
        if cdf_max_val < MIN_CDF_MASS {
            continue;
        }
        let weight = alphas[tid] / eff_lens[tid];
        let lc = length_class_index(quantiles, ref_len as u32);
        let cond = |x: i32| conditional_cdf(cdf, cdf_max_arg, cdf_max_val, x);
        for frag_start in 0..(ref_len - k as i32) {
            let max_fw = ref_len - frag_start + 1;
            let max_rc = frag_start;
            let density_fw = cond(max_fw);
            let density_rc = cond(max_rc);
            if weight * density_fw > EPSILON {
                pos5[lc].add_mass(frag_start, ref_len, (weight * density_fw).ln());
            }
            if weight * density_rc > EPSILON {
                pos3[lc].add_mass(frag_start, ref_len, (weight * density_rc).ln());
            }
        }
    }
    for m in pos5.iter_mut().chain(pos3.iter_mut()) {
        m.finalize();
    }
    (pos5, pos3)
}
