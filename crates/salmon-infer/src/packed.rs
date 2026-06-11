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

/// Reduce per-shard dense accumulators into `alpha_out` (one writer per `tid`,
/// no contention). Parallelized over transcripts.
fn reduce_shards(shards: &[Vec<f64>], alpha_out: &mut [f64]) {
    alpha_out.par_iter_mut().enumerate().for_each(|(tid, out)| {
        let mut s = 0.0;
        for buf in shards {
            s += buf[tid];
        }
        *out = s;
    });
}

/// Parallel EM M-step. Each shard owns a private dense `num_txps` buffer and
/// processes a contiguous slice of the classes with plain (non-atomic) adds;
/// the shards are then summed into `alpha_out`. This avoids both the per-task
/// allocation of a naive fold/reduce and the cross-thread CAS contention of a
/// single shared `AtomicF64` array (which, on hot transcripts, dominated the
/// M-step). The buffers are allocated once in [`run_em_counts`] and reused.
pub(crate) fn em_step_par(
    p: &PackedEqClasses,
    counts: &[u64],
    alpha_in: &[f64],
    alpha_out: &mut [f64],
    shards: &mut [Vec<f64>],
) {
    let nclasses = p.num_classes();
    let chunk = nclasses.div_ceil(shards.len().max(1));
    shards.par_iter_mut().enumerate().for_each(|(s, buf)| {
        buf.iter_mut().for_each(|x| *x = 0.0);
        let start = s * chunk;
        let end = ((s + 1) * chunk).min(nclasses);
        for ci in start..end {
            let count = counts[ci] as f64;
            let (tids, ws) = p.class(ci);
            if tids.len() > 1 {
                let mut denom = 0.0;
                for (&tid, &w) in tids.iter().zip(ws) {
                    denom += alpha_in[tid as usize] * w;
                }
                if denom > MIN_EQ_CLASS_WEIGHT {
                    let inv = count / denom;
                    for (&tid, &w) in tids.iter().zip(ws) {
                        let v = alpha_in[tid as usize] * w;
                        if !v.is_nan() {
                            buf[tid as usize] += v * inv;
                        }
                    }
                }
            } else {
                buf[tids[0] as usize] += count;
            }
        }
    });
    reduce_shards(shards, alpha_out);
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

/// Parallel VBEM M-step. Sharded private buffers + reduce (see [`em_step_par`]).
pub(crate) fn vbem_step_par(
    p: &PackedEqClasses,
    counts: &[u64],
    prior_alphas: &[f64],
    alpha_in: &[f64],
    alpha_out: &mut [f64],
    exp_theta: &mut [f64],
    shards: &mut [Vec<f64>],
) {
    fill_exp_theta(alpha_in, prior_alphas, exp_theta);
    let nclasses = p.num_classes();
    let chunk = nclasses.div_ceil(shards.len().max(1));
    let exp_theta: &[f64] = exp_theta;
    shards.par_iter_mut().enumerate().for_each(|(s, buf)| {
        buf.iter_mut().for_each(|x| *x = 0.0);
        let start = s * chunk;
        let end = ((s + 1) * chunk).min(nclasses);
        for ci in start..end {
            let count = counts[ci] as f64;
            let (tids, ws) = p.class(ci);
            if tids.len() > 1 {
                let mut denom = 0.0;
                for (&tid, &w) in tids.iter().zip(ws) {
                    let et = exp_theta[tid as usize];
                    if et > 0.0 {
                        denom += et * w;
                    }
                }
                if denom > MIN_EQ_CLASS_WEIGHT {
                    let inv = count / denom;
                    for (&tid, &w) in tids.iter().zip(ws) {
                        let et = exp_theta[tid as usize];
                        if et > 0.0 {
                            buf[tid as usize] += et * w * inv;
                        }
                    }
                }
            } else {
                buf[tids[0] as usize] += count;
            }
        }
    });
    reduce_shards(shards, alpha_out);
}
