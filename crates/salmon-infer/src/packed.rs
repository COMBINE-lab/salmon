//! Packed, parallel-friendly equivalence-class representation for inference.
//!
//! Mirrors the flat CSR layout used by `piscem-infer` (`PackedEqMap`): instead
//! of a `Vec<(TranscriptGroup, TGValue)>` of small per-class allocations, all
//! class labels and weights are concatenated into flat arrays indexed by a CSR
//! `starts` offset vector. This is cache-friendly and trivially parallelizable
//! (a class is just a pair of slices), which matters for the main EM and even
//! more for bootstrap/Gibbs, which run the optimizer many times.
//!
//! For class `i` the targets are `labels[starts[i]..starts[i+1]]`, with aligned
//! `combined` (= `weight / effLen`, what EM multiplies by `alpha`) and raw
//! `weights` (what Gibbs multiplies by the sampled fraction `mu`). `counts[i]`
//! is the class's fragment count (overridable per run so bootstrap can resample).

use rayon::prelude::*;
use salmon_eqclass::CollapsedEqClasses;
use statrs::function::gamma::digamma;

/// Minimum `alpha + prior` for which VBEM evaluates `digamma`.
const DIGAMMA_MIN: f64 = 1e-10;

/// Flat CSR equivalence classes (only `valid` groups are retained).
#[derive(Debug, Clone)]
pub struct PackedEqClasses {
    /// flat transcript ids; class `i` spans `labels[starts[i]..starts[i+1]]`
    pub labels: Vec<u32>,
    /// CSR offsets into `labels`/`combined`/`weights`; length `num_classes + 1`
    pub starts: Vec<u32>,
    /// flat `combined_weights` (`weight/effLen`), aligned to `labels`; used by EM
    pub combined: Vec<f64>,
    /// flat raw conditional weights, aligned to `labels`; used by Gibbs
    pub weights: Vec<f64>,
    /// per-class fragment counts
    pub counts: Vec<u64>,
    /// total transcript count (length of an abundance vector)
    pub num_txps: usize,
    /// total fragment count across classes
    pub total_count: u64,
}

impl PackedEqClasses {
    /// Build the packed layout from a finalized [`CollapsedEqClasses`] whose
    /// `combined_weights` are already populated
    /// ([`update_eff_lengths`](salmon_eqclass::CollapsedEqClasses::update_eff_lengths)).
    pub fn from_collapsed(eq: &CollapsedEqClasses, num_txps: usize) -> Self {
        let n = eq.classes.len();
        let mut labels = Vec::new();
        let mut starts = Vec::with_capacity(n + 1);
        let mut combined = Vec::new();
        let mut weights = Vec::new();
        let mut counts = Vec::with_capacity(n);
        starts.push(0u32);
        let mut total = 0u64;
        for (group, value) in &eq.classes {
            if !group.valid {
                continue;
            }
            labels.extend_from_slice(&group.txps);
            combined.extend_from_slice(&value.combined_weights);
            weights.extend_from_slice(&value.weights);
            counts.push(value.count);
            total += value.count;
            starts.push(labels.len() as u32);
        }
        Self {
            labels,
            starts,
            combined,
            weights,
            counts,
            num_txps,
            total_count: total,
        }
    }

    /// Number of (valid) classes.
    #[inline]
    pub fn num_classes(&self) -> usize {
        self.counts.len()
    }

    /// Targets and combined-weights slices for class `i`.
    #[inline]
    pub fn class(&self, i: usize) -> (&[u32], &[f64]) {
        let s = self.starts[i] as usize;
        let e = self.starts[i + 1] as usize;
        (&self.labels[s..e], &self.combined[s..e])
    }
}

/// Smallest denominator weight below which a class is treated as degenerate.
const MIN_EQ_CLASS_WEIGHT: f64 = f64::MIN_POSITIVE;

/// One sequential EM M-step: `alpha_out[t] += count·(alpha_in[t]·w_t)/Σ_j(alpha_in[j]·w_j)`,
/// with single-transcript classes assigned their full count. `counts` overrides
/// the per-class counts (so bootstrap can pass resampled counts).
pub(crate) fn em_step_seq(
    p: &PackedEqClasses,
    counts: &[u64],
    alpha_in: &[f64],
    alpha_out: &mut [f64],
    scratch: &mut Vec<f64>,
) {
    alpha_out.iter_mut().for_each(|a| *a = 0.0);
    for ci in 0..p.num_classes() {
        let count = counts[ci] as f64;
        let (tids, ws) = p.class(ci);
        if tids.len() > 1 {
            scratch.clear();
            let mut denom = 0.0;
            for (&tid, &w) in tids.iter().zip(ws) {
                let v = alpha_in[tid as usize] * w;
                scratch.push(v);
                denom += v;
            }
            if denom > MIN_EQ_CLASS_WEIGHT {
                let inv = count / denom;
                for (&tid, &v) in tids.iter().zip(scratch.iter()) {
                    if !v.is_nan() {
                        alpha_out[tid as usize] += v * inv;
                    }
                }
            }
        } else {
            alpha_out[tids[0] as usize] += count;
        }
    }
}

/// Parallel EM M-step (rayon fold/reduce over classes into per-thread
/// accumulators). Equivalent to [`em_step_seq`] but uses all cores; the
/// reduction is a deterministic sum tree.
pub(crate) fn em_step_par(
    p: &PackedEqClasses,
    counts: &[u64],
    alpha_in: &[f64],
    alpha_out: &mut [f64],
) {
    let num_txps = p.num_txps;
    let acc = (0..p.num_classes())
        .into_par_iter()
        .fold(
            || (vec![0.0f64; num_txps], Vec::<f64>::with_capacity(64)),
            |(mut acc, mut scratch), ci| {
                let count = counts[ci] as f64;
                let (tids, ws) = p.class(ci);
                if tids.len() > 1 {
                    scratch.clear();
                    let mut denom = 0.0;
                    for (&tid, &w) in tids.iter().zip(ws) {
                        let v = alpha_in[tid as usize] * w;
                        scratch.push(v);
                        denom += v;
                    }
                    if denom > MIN_EQ_CLASS_WEIGHT {
                        let inv = count / denom;
                        for (&tid, &v) in tids.iter().zip(scratch.iter()) {
                            if !v.is_nan() {
                                acc[tid as usize] += v * inv;
                            }
                        }
                    }
                } else {
                    acc[tids[0] as usize] += count;
                }
                (acc, scratch)
            },
        )
        .map(|(acc, _)| acc)
        .reduce(
            || vec![0.0f64; num_txps],
            |mut a, b| {
                for (x, y) in a.iter_mut().zip(b) {
                    *x += y;
                }
                a
            },
        );
    alpha_out.copy_from_slice(&acc);
}

/// `exp_theta[i] = exp(digamma(alpha_in[i]+prior_i) - digamma(Σ_j alpha_in[j]+prior_j))`,
/// the VBEM mean-field expectation substituted for `alpha` in the M-step.
fn fill_exp_theta(alpha_in: &[f64], prior_alphas: &[f64], exp_theta: &mut [f64]) {
    let alpha_sum: f64 = alpha_in.iter().zip(prior_alphas).map(|(a, p)| a + p).sum();
    let log_norm = digamma(alpha_sum);
    for i in 0..alpha_in.len() {
        let ap = alpha_in[i] + prior_alphas[i];
        exp_theta[i] = if ap > DIGAMMA_MIN {
            (digamma(ap) - log_norm).exp()
        } else {
            0.0
        };
    }
}

/// One sequential VBEM M-step (uses `exp_theta` in place of `alpha`).
pub(crate) fn vbem_step_seq(
    p: &PackedEqClasses,
    counts: &[u64],
    prior_alphas: &[f64],
    alpha_in: &[f64],
    alpha_out: &mut [f64],
    exp_theta: &mut [f64],
    scratch: &mut Vec<f64>,
) {
    fill_exp_theta(alpha_in, prior_alphas, exp_theta);
    alpha_out.iter_mut().for_each(|a| *a = 0.0);
    for ci in 0..p.num_classes() {
        let count = counts[ci] as f64;
        let (tids, ws) = p.class(ci);
        if tids.len() > 1 {
            scratch.clear();
            let mut denom = 0.0;
            for (&tid, &w) in tids.iter().zip(ws) {
                let et = exp_theta[tid as usize];
                let v = if et > 0.0 { et * w } else { 0.0 };
                scratch.push(v);
                denom += v;
            }
            if denom > MIN_EQ_CLASS_WEIGHT {
                let inv = count / denom;
                for (&tid, &v) in tids.iter().zip(scratch.iter()) {
                    if v > 0.0 {
                        alpha_out[tid as usize] += v * inv;
                    }
                }
            }
        } else {
            alpha_out[tids[0] as usize] += count;
        }
    }
}

/// Parallel VBEM M-step.
pub(crate) fn vbem_step_par(
    p: &PackedEqClasses,
    counts: &[u64],
    prior_alphas: &[f64],
    alpha_in: &[f64],
    alpha_out: &mut [f64],
    exp_theta: &mut [f64],
) {
    fill_exp_theta(alpha_in, prior_alphas, exp_theta);
    let num_txps = p.num_txps;
    let acc = (0..p.num_classes())
        .into_par_iter()
        .fold(
            || (vec![0.0f64; num_txps], Vec::<f64>::with_capacity(64)),
            |(mut acc, mut scratch), ci| {
                let count = counts[ci] as f64;
                let (tids, ws) = p.class(ci);
                if tids.len() > 1 {
                    scratch.clear();
                    let mut denom = 0.0;
                    for (&tid, &w) in tids.iter().zip(ws) {
                        let et = exp_theta[tid as usize];
                        let v = if et > 0.0 { et * w } else { 0.0 };
                        scratch.push(v);
                        denom += v;
                    }
                    if denom > MIN_EQ_CLASS_WEIGHT {
                        let inv = count / denom;
                        for (&tid, &v) in tids.iter().zip(scratch.iter()) {
                            if v > 0.0 {
                                acc[tid as usize] += v * inv;
                            }
                        }
                    }
                } else {
                    acc[tids[0] as usize] += count;
                }
                (acc, scratch)
            },
        )
        .map(|(acc, _)| acc)
        .reduce(
            || vec![0.0f64; num_txps],
            |mut a, b| {
                for (x, y) in a.iter_mut().zip(b) {
                    *x += y;
                }
                a
            },
        );
    alpha_out.copy_from_slice(&acc);
}
