//! Unified bias-corrected effective length: composes sequence-specific, GC, and
//! positional bias exactly as salmon's `updateEffectiveLengths` does — a single
//! conditional-FLD convolution whose per-fragment factor is the product of the
//! enabled bias terms.
//!
//! `fragFactor = seqFW[fragStart]·seqRC[fragEnd] · gcBias({fragFrac, ctxFrac}) ·
//! posFW[fragStart]·posRC[fragEnd]`, summed over fragment starts and convolved
//! with the conditional fragment-length distribution.

use crate::gcbias::{gc_desc, GcFragModel};
use crate::posbias::{length_class_index, SimplePosBias, NUM_LENGTH_CLASSES};
use crate::seqbias::{
    conditional_cdf, log_bias, revcomp_bytes, SBModel, CONTEXT_LEFT, CONTEXT_LENGTH, MIN_ALPHA,
    MIN_CDF_MASS,
};

/// salmon's `EPSILON` (mass cutoff for adding positional expected mass).
const EPSILON: f64 = 0.375e-10;

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

    // Per-position sequence-bias factors (1.0 when not seq-correcting).
    let mut fw = vec![1.0f64; ref_len];
    let mut rc = vec![1.0f64; ref_len];
    if let Some((obs_fw, exp_fw, obs_rc, exp_rc)) = bias.seq {
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
        let mut kstart = 0i32;
        while kstart < ref_len as i32 - fl {
            let frag_start = kstart;
            let frag_end = kstart + fl - 1;
            if (frag_end as usize) < ref_len {
                let mut frag_factor = fw[frag_start as usize] * rc[frag_end as usize];
                if let Some(gc) = gc_model {
                    if let Some((ff, cf)) = gc_desc(gc_prefix, ref_len as i32, frag_start, frag_end)
                    {
                        frag_factor *= gc.get(ff, cf);
                    }
                }
                if let (Some(pf), Some(pr)) = (pos_fw, pos_rc) {
                    frag_factor *= pf[frag_start as usize] * pr[frag_end as usize];
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
