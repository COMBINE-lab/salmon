//! Packed, parallel-friendly equivalence-class representation for inference.
//!
//! # What "packed" means and why it matters
//!
//! The builder in `salmon-eqclass` produces a `Vec<(TranscriptGroup, TGValue)>`:
//! a vector of structs, each owning its own little vectors of transcript ids and
//! weights. That shape is convenient to build concurrently, but it is the wrong
//! shape to iterate a thousand times: each class costs several pointer hops, the
//! data is scattered across the heap, and the CPU's cache prefetcher cannot
//! help.
//!
//! This module flattens it into the classic **CSR** ("compressed sparse row")
//! layout: one long array of transcript ids, one long array of weights, and an
//! offsets array saying where each class begins. Everything a class needs is then
//! contiguous, iterating classes walks memory in a straight line, and splitting
//! the work across threads is just splitting a range of indices.
//!
//! Mirrors the flat CSR layout used by `piscem-infer` (`PackedEqMap`). This is
//! cache-friendly and trivially parallelizable (a class is just a pair of
//! slices), which matters for the main EM and even more for bootstrap/Gibbs,
//! which run the optimizer many times.
//!
//! For class `i` the targets are `labels[starts[i]..starts[i+1]]`, with aligned
//! `combined` (= `weight / effLen`, what EM multiplies by `alpha`) and raw
//! `weights` (what Gibbs multiplies by the sampled fraction `mu`). `counts[i]`
//! is the class's fragment count (overridable per run so bootstrap can resample).
//!
//! # The M-step, in one sentence
//!
//! Every routine below is a variation on the same operation: for each class,
//! split its fragment count across its member transcripts in proportion to
//! `abundance × weight`, and accumulate the result. That is the "M-step" of EM.

use rayon::prelude::*;
use salmon_eqclass::CollapsedEqClasses;
use statrs::function::gamma::digamma;

/// Minimum `alpha + prior` for which VBEM evaluates `digamma`.
/// Below this the function heads to negative infinity and the result is
/// indistinguishable from zero anyway, so it is short-circuited.
const DIGAMMA_MIN: f64 = 1e-10;

/// Flat CSR equivalence classes (only `valid` groups are retained).
#[derive(Debug, Clone)]
pub struct PackedEqClasses {
    /// flat transcript ids; class `i` spans `labels[starts[i]..starts[i+1]]`
    pub labels: Vec<u32>,
    /// CSR offsets into `labels`/`combined`/`weights`; length `num_classes + 1`
    ///
    /// The extra trailing entry is what lets class `i`'s span be written without
    /// a special case for the last class.
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

/// Which posterior-sampling method a run will use, if any.
///
/// Bootstrap and Gibbs are mutually exclusive — the CLI rejects both together —
/// so this is the single place that fact is written down. Deriving it once and
/// passing it around replaces a `bool` at each call site whose meaning ("does
/// this run need the raw per-mapping weights?") was only recoverable by knowing
/// which sampler reads which array.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub enum PosteriorMethod {
    /// Point estimate only; no posterior replicates.
    #[default]
    None,
    /// Non-parametric bootstrap over the equivalence classes.
    Bootstrap { replicates: u32 },
    /// Collapsed Gibbs sampling.
    Gibbs { samples: u32 },
}

impl PosteriorMethod {
    /// Resolve from the raw option values.
    ///
    /// `skip_quant` suppresses both. Bootstrap wins if both counts are somehow
    /// non-zero, matching the historical dispatch order — though the CLI now
    /// rejects that combination, so it should be unreachable from a real run.
    pub fn resolve(skip_quant: bool, num_bootstraps: u32, num_gibbs_samples: u32) -> Self {
        if skip_quant {
            Self::None
        } else if num_bootstraps > 0 {
            Self::Bootstrap {
                replicates: num_bootstraps,
            }
        } else if num_gibbs_samples > 0 {
            Self::Gibbs {
                samples: num_gibbs_samples,
            }
        } else {
            Self::None
        }
    }

    /// Whether the run needs [`PackedEqClasses::weights`], the raw per-mapping
    /// conditional weights.
    ///
    /// Only the Gibbs sampler reads them; the EM and the bootstrap use
    /// `combined`. On a human-scale run that array is ~8 bytes per CSR
    /// incidence, so skipping it saves on the order of 160 MB that would be
    /// allocated, filled and never read.
    pub fn needs_raw_weights(self) -> bool {
        matches!(self, Self::Gibbs { .. })
    }
}

impl PackedEqClasses {
    /// Build the packed layout from a finalized [`CollapsedEqClasses`] whose
    /// `combined_weights` are already populated
    /// ([`update_eff_lengths`](salmon_eqclass::CollapsedEqClasses::update_eff_lengths)).
    pub fn from_collapsed(eq: &CollapsedEqClasses, num_txps: usize) -> Self {
        Self::from_collapsed_for(eq, num_txps, PosteriorMethod::Gibbs { samples: 1 })
    }

    /// As [`from_collapsed`](Self::from_collapsed), but skips the raw weights
    /// leaves the raw `weights` array empty.
    ///
    /// `weights` is read only by the Gibbs sampler; the EM path uses `combined`.
    /// On a run without `--numGibbsSamples` it is therefore 8 bytes per CSR
    /// incidence that is allocated, filled and never read — on the order of
    /// 160 MB for a human-scale run. Callers that know Gibbs will not run pass
    /// `false`.
    pub fn from_collapsed_for(
        eq: &CollapsedEqClasses,
        num_txps: usize,
        method: PosteriorMethod,
    ) -> Self {
        let n = eq.classes.len();
        // Size the flat arrays up front. Growing them from empty costs a
        // reallocation plus a full copy at every doubling, and leaves a transient
        // peak of ~1.5x the final size while both buffers are live. The exact
        // total is one cheap pass over the class list (it only reads `Vec::len`).
        let total_labels: usize = eq
            .classes
            .iter()
            .filter(|(g, _)| g.valid)
            .map(|(g, _)| g.txps.len())
            .sum();
        // The CSR offsets below are `u32`. Past `u32::MAX` incidences the `as u32`
        // cast would wrap silently and `class()` would hand the optimizer slices
        // belonging to other classes — wrong abundances with no diagnostic. Fail
        // loudly instead.
        assert!(
            total_labels <= u32::MAX as usize,
            "equivalence-class CSR needs {total_labels} incidences, which exceeds \
             the u32 offset limit ({}); please report this input",
            u32::MAX
        );
        let mut labels = Vec::with_capacity(total_labels);
        let mut starts = Vec::with_capacity(n + 1);
        let mut combined = Vec::with_capacity(total_labels);
        let keep_weights = method.needs_raw_weights();
        let mut weights = Vec::with_capacity(if keep_weights { total_labels } else { 0 });
        let mut counts = Vec::with_capacity(n);
        // The first class starts at offset 0; every later entry is appended after
        // that class's data has been copied in.
        starts.push(0u32);
        let mut total = 0u64;
        for (group, value) in &eq.classes {
            // Classes the optimizer marked degenerate are simply not carried over.
            if !group.valid {
                continue;
            }
            labels.extend_from_slice(&group.txps);
            combined.extend_from_slice(&value.combined_weights);
            if keep_weights {
                weights.extend_from_slice(&value.weights);
            }
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

    /// Rewrite `combined` in place from `eq`'s current `combined_weights`,
    /// leaving `labels`, `starts`, `counts` and `weights` untouched.
    ///
    /// Bias correction recomputes effective lengths and calls
    /// [`CollapsedEqClasses::update_eff_lengths`], which changes *only* the
    /// combined weights. Rebuilding the whole packed layout to pick that up
    /// re-copies the labels and counts too — hundreds of MB on a human-scale run
    /// — for arrays bit-identical to the ones already held.
    ///
    /// `eq` must be the class set the layout was built from. Classes are never
    /// invalidated after construction (nothing clears `TranscriptGroup::valid`),
    /// so the filter below matches the original build; the assertion pins that.
    pub fn refresh_combined(&mut self, eq: &CollapsedEqClasses) {
        let mut at = 0usize;
        for (group, value) in &eq.classes {
            if !group.valid {
                continue;
            }
            let n = value.combined_weights.len();
            self.combined[at..at + n].copy_from_slice(&value.combined_weights);
            at += n;
        }
        assert_eq!(
            at,
            self.combined.len(),
            "refresh_combined: class set changed since the packed layout was built"
        );
    }

    /// Number of (valid) classes.
    #[inline]
    pub fn num_classes(&self) -> usize {
        self.counts.len()
    }

    /// Targets and combined-weights slices for class `i`.
    ///
    /// Returning borrowed slices means no allocation and no copying: this is
    /// called once per class per iteration, millions of times per run.
    #[inline]
    pub fn class(&self, i: usize) -> (&[u32], &[f64]) {
        let s = self.starts[i] as usize;
        let e = self.starts[i + 1] as usize;
        (&self.labels[s..e], &self.combined[s..e])
    }
}

/// Smallest denominator weight below which a class is treated as degenerate.
/// Dividing by anything smaller would produce an infinity.
const MIN_EQ_CLASS_WEIGHT: f64 = f64::MIN_POSITIVE;

/// Final masked M-step that redistributes truncated mass instead of rescaling.
///
/// After the EM/VBEM has converged, transcripts whose abundance is below
/// `min_alpha` are negligible noise; salmon zeroes them. Rather than zero them
/// and rescale the survivors back up to the total (which can inflate a transcript
/// above the eq-class mass that supports it), this marks those transcripts
/// **inactive** and runs one more M-step that distributes each class's count only
/// among its remaining active members — so the truncated mass flows to genuine
/// eq-class co-members and the total is preserved. Inactive transcripts are given
/// zero weight (in VBEM this also keeps the Dirichlet prior from reviving them).
///
/// The distinction matters: rescaling would move a truncated transcript's mass to
/// transcripts that share no evidence with it, whereas redistribution moves it
/// only to transcripts the same fragments were also compatible with.
///
/// A class whose every member is inactive has nowhere to put its count; that mass
/// cannot be redistributed and is summed into the returned `dropped` value
/// (reported as `inference_truncated_mass`, not rescaled away). Returns the
/// redistributed alphas (`Σ == total_count − dropped`) and `dropped`.
pub(crate) fn redistribute_truncated(
    p: &PackedEqClasses,
    counts: &[u64],
    alpha_conv: &[f64],
    prior_alphas: &[f64],
    min_alpha: f64,
    use_vbem: bool,
) -> (Vec<f64>, f64) {
    let n = p.num_txps;
    let inactive: Vec<bool> = alpha_conv.iter().map(|&a| a < min_alpha).collect();
    // Per-transcript distribution basis: VBEM uses expTheta(α+prior); EM uses α.
    // Inactive transcripts get 0 either way (the digamma normalization is a common
    // factor that cancels in each class's ratio, so it is omitted here).
    let mut basis = vec![0.0f64; n];
    for i in 0..n {
        if inactive[i] {
            continue;
        }
        basis[i] = if use_vbem {
            let ap = alpha_conv[i] + prior_alphas[i];
            if ap > DIGAMMA_MIN {
                digamma(ap).exp()
            } else {
                0.0
            }
        } else {
            alpha_conv[i]
        };
    }
    let mut alpha_out = vec![0.0f64; n];
    let mut dropped = 0.0f64;
    // Reused across classes so the loop allocates nothing; 64 comfortably covers
    // a typical class size.
    let mut scratch: Vec<f64> = Vec::with_capacity(64);
    for ci in 0..p.num_classes() {
        let count = counts[ci] as f64;
        let (tids, ws) = p.class(ci);
        if tids.len() > 1 {
            // Ambiguous class: split the count in proportion to basis × weight.
            scratch.clear();
            let mut denom = 0.0;
            for (&tid, &w) in tids.iter().zip(ws) {
                let v = basis[tid as usize] * w;
                scratch.push(v);
                denom += v;
            }
            if denom > MIN_EQ_CLASS_WEIGHT {
                // Multiply by the reciprocal once rather than dividing per member.
                let inv = count / denom;
                for (&tid, &v) in tids.iter().zip(scratch.iter()) {
                    if v > 0.0 {
                        alpha_out[tid as usize] += v * inv;
                    }
                }
            } else {
                dropped += count; // every member truncated: cannot redistribute
            }
        } else if inactive[tids[0] as usize] {
            dropped += count; // single-transcript class, its transcript truncated
        } else {
            // Unambiguous class: its whole count belongs to that transcript.
            alpha_out[tids[0] as usize] += count;
        }
    }
    (alpha_out, dropped)
}

/// One sequential EM M-step: `alpha_out[t] += count·(alpha_in[t]·w_t)/Σ_j(alpha_in[j]·w_j)`,
/// with single-transcript classes assigned their full count. `counts` overrides
/// the per-class counts (so bootstrap can pass resampled counts).
///
/// Read the formula as: a fragment that could have come from several transcripts
/// is split between them in proportion to how likely each is — which depends on
/// the current abundance estimate, which is why this has to be iterated.
pub(crate) fn em_step_seq(
    p: &PackedEqClasses,
    counts: &[u64],
    alpha_in: &[f64],
    alpha_out: &mut [f64],
    scratch: &mut Vec<f64>,
) {
    // The output accumulates from zero every step; it is not an update in place.
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
                    // NaN would silently poison the transcript's total for the
                    // rest of the run; skip rather than propagate.
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
///
/// Parallelizing over *transcripts* rather than shards is what removes the
/// contention: each output slot is written by exactly one task, which reads that
/// slot from every shard.
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
///
/// The trade is memory for speed: one dense `num_txps` buffer per shard, in
/// exchange for every add being an ordinary non-atomic one. A fixed shard count
/// and a fixed class partition also keep the summation order stable, so the
/// parallel result reproduces run to run.
pub(crate) fn em_step_par(
    p: &PackedEqClasses,
    counts: &[u64],
    alpha_in: &[f64],
    alpha_out: &mut [f64],
    shards: &mut [Vec<f64>],
) {
    let nclasses = p.num_classes();
    // Contiguous slices: `div_ceil` makes the last shard the short one.
    let chunk = nclasses.div_ceil(shards.len().max(1));
    shards.par_iter_mut().enumerate().for_each(|(s, buf)| {
        // Buffers are reused across iterations, so clear before accumulating.
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
                    // Recomputes `alpha_in[tid] * w` rather than keeping a scratch
                    // vector: the multiply is cheaper than the extra memory
                    // traffic inside a parallel closure.
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
///
/// **What VBEM changes.** Plain EM asks "what single abundance vector best
/// explains the data?". Variational Bayes instead keeps a distribution over
/// abundance vectors (a Dirichlet) and works with the expectation of its log.
/// That expectation is a difference of `digamma` functions, and `exp_theta` is
/// simply its exponential — so every M-step below is the EM one with `alpha`
/// replaced by `exp_theta`. In practice VBEM shrinks weakly-supported
/// transcripts toward zero rather than letting them float on noise.
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
                // A zeroed expectation contributes nothing; guarding here keeps
                // `0 * w` out of the sum entirely.
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
    // Computed once, before the parallel region: it depends on a global sum over
    // all transcripts, so it cannot be sharded.
    fill_exp_theta(alpha_in, prior_alphas, exp_theta);
    let nclasses = p.num_classes();
    let chunk = nclasses.div_ceil(shards.len().max(1));
    // Reborrow as immutable so the closure below can share it across threads.
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

#[cfg(test)]
mod tests {
    use super::*;
    use salmon_eqclass::{EquivalenceClassBuilder, TranscriptGroup};

    /// Build a packed set from `(transcript ids, count)` pairs with uniform
    /// weights and unit effective lengths, so the tests exercise the
    /// redistribution logic and nothing else.
    fn packed(classes: &[(Vec<u32>, u64)], num_txps: usize) -> PackedEqClasses {
        let b = EquivalenceClassBuilder::new();
        for (txps, count) in classes {
            b.add_group(
                TranscriptGroup::new(txps.clone()),
                vec![1.0; txps.len()],
                *count,
            );
        }
        let mut eq = b.finish();
        eq.update_eff_lengths(&vec![1.0; num_txps]);
        PackedEqClasses::from_collapsed(&eq, num_txps)
    }

    /// The core claim of `redistribute_truncated`: truncated mass goes to
    /// transcripts that shared evidence with it, and the total is unchanged.
    #[test]
    fn redistribute_moves_truncated_mass_to_comembers_no_rescale() {
        // Shared class {0,1} with count 100; transcript 1 is truncated. Its share
        // must flow to its co-member (0), not be recovered by rescaling everything.
        let p = packed(&[(vec![0, 1], 100)], 2);
        let alpha_conv = vec![100.0, 1e-12];
        let (out, dropped) =
            redistribute_truncated(&p, &p.counts, &alpha_conv, &[0.0, 0.0], 1e-8, false);
        assert_eq!(dropped, 0.0);
        assert!(
            (out[0] - 100.0).abs() < 1e-9,
            "co-member should get the mass: {out:?}"
        );
        assert_eq!(out[1], 0.0, "truncated transcript stays 0");
        assert!(
            ((out[0] + out[1]) - 100.0).abs() < 1e-9,
            "mass preserved exactly"
        );
    }

    /// Mass with nowhere to go must be reported, not silently absorbed — that is
    /// what makes the `inference_truncated_mass` metric trustworthy.
    #[test]
    fn redistribute_reports_fully_truncated_class_mass() {
        // Class {0} count 5 (active) + class {1} count 3 whose only transcript is
        // truncated -> that 3 cannot be redistributed and is reported as dropped.
        let p = packed(&[(vec![0], 5), (vec![1], 3)], 2);
        let alpha_conv = vec![10.0, 1e-12];
        let (out, dropped) =
            redistribute_truncated(&p, &p.counts, &alpha_conv, &[0.0, 0.0], 1e-8, false);
        assert_eq!(out[0], 5.0);
        assert_eq!(out[1], 0.0);
        assert_eq!(dropped, 3.0, "fully-truncated class mass must be reported");
    }

    /// The subtle VBEM interaction: because the Dirichlet prior keeps `expTheta`
    /// positive even at zero abundance, the explicit inactive mask is what stops a
    /// truncated transcript from coming back to life in the final step.
    #[test]
    fn redistribute_vbem_prior_does_not_revive_truncated() {
        // Under VBEM the Dirichlet prior makes expTheta nonzero even at alpha=0;
        // the inactive mask must keep a truncated transcript at 0 regardless.
        let p = packed(&[(vec![0, 1], 100)], 2);
        let alpha_conv = vec![100.0, 1e-12];
        let (out, dropped) =
            redistribute_truncated(&p, &p.counts, &alpha_conv, &[0.01, 0.01], 1e-8, true);
        assert_eq!(dropped, 0.0);
        assert_eq!(
            out[1], 0.0,
            "VBEM prior must not revive a truncated transcript"
        );
        assert!(
            (out[0] - 100.0).abs() < 1e-9,
            "all mass to the surviving co-member"
        );
    }
}
