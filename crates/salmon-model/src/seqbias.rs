//! Sequence-specific bias model (`SBModel`).
//!
//! A faithful port of salmon's `SBModel` (`src/model/SBModel.cpp`): a
//! variable-order Markov model over the 9-base sequence context surrounding a
//! fragment's start position (3 bases before the start, the start, and 5 after).
//! Per-position Markov orders are `{0,1,2,2,2,2,2,2,2}`. Counts are accumulated
//! from observed fragment-start contexts (the *observed* model) and from the
//! transcriptome (the *expected* model); the ratio of the two scores the
//! sequence bias at any position, which is used to correct effective lengths.
//!
//! The 2-bit base encoding only needs to be self-consistent (the bias is a
//! ratio of two models built with the same encoding), so we use A=0, C=1,
//! G=2, T=3.

/// Per-position Markov orders (salmon's "simple" model). Length is the context.
const ORDER: [u32; 9] = [0, 1, 2, 2, 2, 2, 2, 2, 2];
/// Context length (= ORDER.len()): 3 left + start + 5 right.
pub const CONTEXT_LENGTH: usize = 9;
/// Bases before the fragment-start position.
pub const CONTEXT_LEFT: usize = 3;
/// Bases at/after the fragment-start position.
pub const CONTEXT_RIGHT: usize = 5;
/// Rows in the probability table: 4^(maxOrder+1) = 4^3.
const ROWS: usize = 64;
/// Pseudocount prior.
const PRIOR: f64 = 1e-10;
/// Floor used when taking the log of a zero probability.
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
    /// probabilities, laid out position-major: `probs[pos * ROWS + idx]`
    probs: Vec<f64>,
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
            marginals: vec![PRIOR; 4 * CONTEXT_LENGTH],
            shifts,
            masks,
            trained: false,
        }
    }

    /// Encode a 9-base context (`CONTEXT_LENGTH` bytes) into a 2-bit-per-base
    /// integer with base 0 in the high bits. `rev_comp` reverse-complements it.
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

    /// Accumulate one observed context with the given weight.
    pub fn add_context(&mut self, context: &[u8], rev_comp: bool, weight: f64) {
        debug_assert!(!self.trained, "cannot add to a normalized model");
        let mer = Self::encode(context, rev_comp);
        for pos in 0..CONTEXT_LENGTH {
            let idx = self.index_at(mer, pos);
            self.probs[pos * ROWS + idx] += weight;
        }
    }

    /// Convert accumulated counts into conditional log-probabilities. Idempotent
    /// guard: a model can only be normalized once.
    pub fn normalize(&mut self) {
        if self.trained {
            return;
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
        for (a, b) in self.probs.iter_mut().zip(&other.probs) {
            *a += *b - PRIOR; // avoid double-counting the prior
        }
    }
}

/// Reverse-complement a DNA byte slice (ACGT; other bases map to `A`).
fn revcomp_bytes(seq: &[u8]) -> Vec<u8> {
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
const MIN_ALPHA: f64 = 1e-8;
/// Minimum reliable CDF mass for a transcript (salmon's `minCDFMass`).
const MIN_CDF_MASS: f64 = 1e-10;
/// Fragment-length sampling stride in the effective-length convolution
/// (salmon's `pdfSampFactor` = `biasSpeedSamp` default).
pub const FLD_SAMP_STRIDE: usize = 5;

/// Linear cumulative fragment-length distribution plus the `[low, high]`
/// fragment-length quantile bounds (0.5% / 99.5%), mirroring the `cdf`,
/// `fldLow`, `fldHigh` salmon computes in `updateEffectiveLengths`.
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
#[inline]
fn conditional_cdf(cdf: &[f64], cdf_max_arg: usize, cdf_max_val: f64, x: i32) -> f64 {
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
pub fn build_expected<'a, F>(
    num_refs: usize,
    seq_of: F,
    alphas: &[f64],
    eff_lens: &[f64],
    cdf: &[f64],
) -> (SBModel, SBModel)
where
    F: Fn(usize) -> &'a [u8],
{
    let k = CONTEXT_LENGTH;
    let cu = CONTEXT_LEFT as i32;
    let mut exp_fw = SBModel::new();
    let mut exp_rc = SBModel::new();
    for tid in 0..num_refs {
        if alphas[tid] < MIN_ALPHA || eff_lens[tid] <= 0.0 {
            continue;
        }
        let seq = seq_of(tid);
        let ref_len = seq.len();
        if ref_len < k {
            continue;
        }
        let cdf_max_arg = (cdf.len() - 1).min(ref_len);
        let cdf_max_val = cdf[cdf_max_arg];
        if cdf_max_val < MIN_CDF_MASS {
            continue;
        }
        let weight = alphas[tid] / eff_lens[tid];
        let rc = revcomp_bytes(seq);
        // fragStartPos in 0..(refLen - K) (salmon's loop bound)
        for frag_start in 0..(ref_len - k) {
            let max_frag_len = ref_len as i32 - (frag_start as i32 + cu);
            if max_frag_len >= 0 && (max_frag_len as usize) < ref_len {
                let cdensity = conditional_cdf(cdf, cdf_max_arg, cdf_max_val, max_frag_len);
                let w = weight * cdensity;
                exp_fw.add_context(&seq[frag_start..frag_start + k], false, w);
                exp_rc.add_context(&rc[frag_start..frag_start + k], false, w);
            }
        }
    }
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
            fw[read_start] = log_bias(obs_fw, exp_fw, &seq[frag_start..frag_start + k], false).exp();
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

/// Log bias of `observed` relative to `expected` for a context:
/// `log P_obs(context) - log P_exp(context)`. The fragment-level bias weight is
/// `exp` of this.
pub fn log_bias(observed: &SBModel, expected: &SBModel, context: &[u8], rev_comp: bool) -> f64 {
    observed.evaluate_log(context, rev_comp) - expected.evaluate_log(context, rev_comp)
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
