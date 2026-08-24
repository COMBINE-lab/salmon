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
use std::sync::atomic::{AtomicBool, AtomicU32, AtomicU64, AtomicUsize, Ordering};
use std::sync::OnceLock;

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
    ///
    /// `u64` rather than `u32` because these are cumulative *incidence* counts,
    /// not transcript ids: they grow with the whole dataset, and past
    /// `u32::MAX` a narrower offset wrapped silently and handed the optimizer
    /// slices belonging to other classes. One entry per class, so the width
    /// costs ~3.5% of the packed structure at any scale — a fixed and small
    /// price for having no ceiling at all (#1097).
    pub starts: Vec<u64>,
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
        let total_labels: usize = eq.classes.iter().map(|(g, _)| g.txps.len()).sum();
        let mut labels = Vec::with_capacity(total_labels);
        let mut starts = Vec::with_capacity(n + 1);
        let mut combined = Vec::with_capacity(total_labels);
        let keep_weights = method.needs_raw_weights();
        let mut weights = Vec::with_capacity(if keep_weights { total_labels } else { 0 });
        let mut counts = Vec::with_capacity(n);
        // The first class starts at offset 0; every later entry is appended after
        // that class's data has been copied in.
        starts.push(0u64);
        let mut total = 0u64;
        for (group, value) in &eq.classes {
            labels.extend_from_slice(&group.txps);
            combined.extend_from_slice(&value.combined_weights);
            if keep_weights {
                weights.extend_from_slice(&value.weights);
            }
            counts.push(value.count);
            total += value.count;
            starts.push(labels.len() as u64);
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
    /// `eq` must be the class set the layout was built from: this walks the
    /// classes in the same order the build did, so the offsets line up. The
    /// assertion pins that.
    pub fn refresh_combined(&mut self, eq: &CollapsedEqClasses) {
        let mut at = 0usize;
        for (_group, value) in &eq.classes {
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
    ///
    /// The `as usize` narrowing cannot lose anything: an offset only exceeds
    /// `usize` on a 32-bit target, where `labels` would have had to exceed a
    /// 16 GB allocation to get there.
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

/// The fixed work partition for the parallel M-step, plus the set of
/// transcripts each shard can touch.
///
/// Two problems live here, and they share one solution.
///
/// **Determinism.** The shard count used to be `rayon::current_num_threads()`,
/// so `-p` changed the class partition and therefore the order of the
/// floating-point sums. Different thread counts produced different `quant.sf`
/// files — on 26M real fragments across 238k transcripts, four rows disagreed
/// between `-p 4` and `-p 32`, while `-p 64` and `-p 80` (both clamped to 64)
/// agreed exactly. That directly contradicts the guarantee 2.6.0 is built on.
/// The partition is now a function of the *data* alone, so the result is
/// identical at any `-p`; rayon schedules these shards onto whatever threads
/// exist, which is a scheduling decision and cannot change the arithmetic.
///
/// **Scaling.** Clearing and reducing dense `num_txps` buffers cost
/// `O(nshards × num_txps)` per iteration no matter how much real work there
/// was — 30M operations per iteration at 64 shards over 238k transcripts,
/// against 6.9M eq-class entries of actual work on the measured workload. That
/// overhead grew with the thread count while the useful work stayed fixed,
/// which is why the EM phase stopped scaling past ~32 threads and then
/// regressed. Since the class labels
/// never change between iterations, each shard's reachable transcript set can
/// be computed once, and the per-iteration cost falls to `O(entries)`.
#[derive(Clone, Copy, Debug, Default)]
struct Contributor {
    /// Logical shard id. `MAX_SHARDS` is deliberately below `u16::MAX`.
    shard: u16,
    /// Position of this transcript in that shard's compressed accumulator.
    local: u32,
}

pub(crate) struct ShardPlan {
    /// Contiguous class boundaries, length `num_shards + 1`.
    boundaries: Vec<usize>,
    /// Number of reachable transcripts in each shard; also its buffer length.
    local_lens: Vec<usize>,
    /// Local accumulator index for every packed incidence, aligned with
    /// `PackedEqClasses::labels` and `combined`.
    entry_local: Vec<u32>,
    /// CSR offsets into `contributors`, one range per transcript. `u64` is
    /// deliberate: the packed layout supports more than 4G incidences and its
    /// auxiliary indices must not silently reintroduce that old ceiling.
    contrib_start: Vec<u64>,
    /// Contributing shards in ascending shard order, with the local buffer index
    /// needed to read each compressed accumulator directly.
    contributors: Vec<Contributor>,
}

impl ShardPlan {
    /// Logical parallelism is derived only from the packed data, never the Rayon
    /// pool. Several logical shards per high-core worker leave enough tasks for
    /// work stealing without making `-p` part of the arithmetic.
    const MAX_SHARDS: usize = 128;
    const MIN_INCIDENCES_PER_SHARD: usize = 4096;

    pub(crate) fn new(p: &PackedEqClasses) -> Self {
        let nclasses = p.num_classes();
        let total_incidences = p.labels.len();
        let max_shards = Self::MAX_SHARDS.min(nclasses.max(1));
        let nshards = total_incidences
            .div_ceil(Self::MIN_INCIDENCES_PER_SHARD)
            .clamp(1, max_shards);

        // Place deterministic boundaries at equal cumulative-incidence targets.
        // Classes stay intact, so a single exceptionally large class may still
        // dominate one shard, but ordinary class-size skew no longer creates a
        // long tail merely because every shard received the same class count.
        let mut boundaries = Vec::with_capacity(nshards + 1);
        boundaries.push(0);
        for s in 1..nshards {
            let target = ((total_incidences as u128 * s as u128).div_ceil(nshards as u128)) as u64;
            let min_ci = boundaries[s - 1] + 1;
            let max_ci = nclasses - (nshards - s);
            let rel = p.starts[min_ci..=max_ci].partition_point(|&offset| offset < target);
            let upper = (min_ci + rel).min(max_ci);
            let lower = upper.saturating_sub(1).max(min_ci);
            let lower_distance = p.starts[lower].abs_diff(target);
            let upper_distance = p.starts[upper].abs_diff(target);
            boundaries.push(if lower_distance <= upper_distance {
                lower
            } else {
                upper
            });
        }
        boundaries.push(nclasses);

        // Build the sorted touched set and the incidence->local-index map
        // together. The latter removes every binary search/hash lookup from the
        // thousands of repeated M-steps; the binary searches happen only here.
        let shard_maps: Vec<(Vec<u32>, Vec<u32>)> = (0..nshards)
            .into_par_iter()
            .map(|s| {
                let entry_start = p.starts[boundaries[s]] as usize;
                let entry_end = p.starts[boundaries[s + 1]] as usize;
                let labels = &p.labels[entry_start..entry_end];
                let mut touched = labels.to_vec();
                touched.sort_unstable();
                touched.dedup();
                let locals = labels
                    .iter()
                    .map(|tid| {
                        u32::try_from(
                            touched
                                .binary_search(tid)
                                .expect("a shard's touched set came from these labels"),
                        )
                        .expect("a shard cannot touch more than u32::MAX transcripts")
                    })
                    .collect();
                (touched, locals)
            })
            .collect();
        let local_lens = shard_maps
            .iter()
            .map(|(touched, _)| touched.len())
            .collect();
        let entry_local = shard_maps
            .iter()
            .flat_map(|(_, locals)| locals.iter().copied())
            .collect();

        // Invert the touched sets into CSR. Counting pass, prefix sum, fill pass —
        // the shards are visited in ascending order in the fill, so each
        // transcript's contributor list is ascending too.
        let num_txps = p.num_txps;
        let mut counts = vec![0u64; num_txps + 1];
        for (touched, _) in &shard_maps {
            for &tid in touched {
                counts[tid as usize + 1] = counts[tid as usize + 1]
                    .checked_add(1)
                    .expect("contributor count exceeds u64");
            }
        }
        for i in 0..num_txps {
            counts[i + 1] = counts[i + 1]
                .checked_add(counts[i])
                .expect("contributor CSR offset exceeds u64");
        }
        let contrib_start = counts;
        let mut cursor = contrib_start.clone();
        let contributor_len = usize::try_from(contrib_start[num_txps])
            .expect("contributor array must fit in address space");
        let mut contributors = vec![Contributor::default(); contributor_len];
        for (s, (touched, _)) in shard_maps.iter().enumerate() {
            let shard = u16::try_from(s).expect("MAX_SHARDS fits in u16");
            for (local, &tid) in touched.iter().enumerate() {
                let at = &mut cursor[tid as usize];
                contributors[*at as usize] = Contributor {
                    shard,
                    local: u32::try_from(local)
                        .expect("a shard cannot touch more than u32::MAX transcripts"),
                };
                *at = at.checked_add(1).expect("contributor cursor exceeds u64");
            }
        }

        Self {
            boundaries,
            local_lens,
            entry_local,
            contrib_start,
            contributors,
        }
    }

    pub(crate) fn num_shards(&self) -> usize {
        self.local_lens.len()
    }

    pub(crate) fn allocate_buffers(&self) -> Vec<Vec<f64>> {
        let mut buffers = Vec::with_capacity(self.num_shards());
        buffers.extend(self.local_lens.iter().map(|&len| vec![0.0; len]));
        buffers
    }

    fn range(&self, s: usize) -> (usize, usize) {
        (self.boundaries[s], self.boundaries[s + 1])
    }

    fn entry_locals(&self, p: &PackedEqClasses, ci: usize) -> &[u32] {
        let lo = p.starts[ci] as usize;
        let hi = p.starts[ci + 1] as usize;
        &self.entry_local[lo..hi]
    }
}

/// Sum the shard accumulators into `alpha_out`, reading only the entries the
/// shards actually wrote.
///
/// Parallel over transcripts and deterministic: each output element is an
/// independent sum over that transcript's contributing shards, taken in
/// ascending shard order from the plan's inverted index. Total work is
/// `Σ|touched|` — on the order of the eq-class entry count — where the previous
/// dense reduction cost `nshards × num_txps` regardless of the data.
fn reduce_shards_sparse(shards: &[Vec<f64>], plan: &ShardPlan, alpha_out: &mut [f64]) {
    alpha_out.par_iter_mut().enumerate().for_each(|(tid, out)| {
        let lo = plan.contrib_start[tid] as usize;
        let hi = plan.contrib_start[tid + 1] as usize;
        let mut s = 0.0;
        for contributor in &plan.contributors[lo..hi] {
            s += shards[contributor.shard as usize][contributor.local as usize];
        }
        *out = s;
    });
}

/// Parallel EM M-step. Each shard owns a private compressed buffer containing
/// only the transcripts reachable from its class range. A precomputed local
/// index aligned with every packed incidence makes accumulation a direct add;
/// there is no lookup in the repeated hot path. Each shard processes a
/// contiguous slice of the classes with plain (non-atomic) adds;
/// the shards are then summed into `alpha_out`. This avoids both the per-task
/// allocation of a naive fold/reduce and the cross-thread CAS contention of a
/// single shared `AtomicF64` array (which, on hot transcripts, dominated the
/// M-step). The buffers are allocated once in [`run_em_counts`] and reused.
///
/// Fixed data-derived shard boundaries and an ascending contributor reduction
/// keep the floating-point order stable regardless of how Rayon schedules work.
pub(crate) fn em_step_par(
    p: &PackedEqClasses,
    counts: &[u64],
    alpha_in: &[f64],
    alpha_out: &mut [f64],
    shards: &mut [Vec<f64>],
    plan: &ShardPlan,
) {
    shards.par_iter_mut().enumerate().for_each(|(s, buf)| {
        // The whole compressed buffer was touched in the previous iteration;
        // clearing it contiguously is both complete and cache-friendly.
        buf.fill(0.0);
        let (start, end) = plan.range(s);
        for ci in start..end {
            let count = counts[ci] as f64;
            let (tids, ws) = p.class(ci);
            let locals = plan.entry_locals(p, ci);
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
                    for i in 0..tids.len() {
                        let v = alpha_in[tids[i] as usize] * ws[i];
                        if !v.is_nan() {
                            buf[locals[i] as usize] += v * inv;
                        }
                    }
                }
            } else {
                buf[locals[0] as usize] += count;
            }
        }
    });
    reduce_shards_sparse(shards, plan, alpha_out);
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
fn fill_exp_theta_seq(alpha_in: &[f64], prior_alphas: &[f64], exp_theta: &mut [f64]) {
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

/// Parallel VBEM expectation fill used by the point estimate. Its sum remains
/// sequential and therefore pool-independent; the element-wise map is parallel.
fn fill_exp_theta_par(alpha_in: &[f64], prior_alphas: &[f64], exp_theta: &mut [f64]) {
    let alpha_sum: f64 = alpha_in.iter().zip(prior_alphas).map(|(a, p)| a + p).sum();
    let log_norm = digamma(alpha_sum);
    exp_theta.par_iter_mut().enumerate().for_each(|(i, et)| {
        let ap = alpha_in[i] + prior_alphas[i];
        *et = if ap > DIGAMMA_MIN {
            (digamma(ap) - log_norm).exp()
        } else {
            0.0
        };
    });
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
    // Bootstrap already parallelizes over replicates. Keeping this entire
    // per-replicate kernel sequential avoids nested Rayon scheduling.
    fill_exp_theta_seq(alpha_in, prior_alphas, exp_theta);
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
    plan: &ShardPlan,
) {
    // Computed once, before the parallel region: it depends on a global sum over
    // all transcripts, so it cannot be sharded.
    fill_exp_theta_par(alpha_in, prior_alphas, exp_theta);
    // Reborrow as immutable so the closure below can share it across threads.
    let exp_theta: &[f64] = exp_theta;
    shards.par_iter_mut().enumerate().for_each(|(s, buf)| {
        buf.fill(0.0);
        let (start, end) = plan.range(s);
        for ci in start..end {
            let count = counts[ci] as f64;
            let (tids, ws) = p.class(ci);
            let locals = plan.entry_locals(p, ci);
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
                    for i in 0..tids.len() {
                        let et = exp_theta[tids[i] as usize];
                        if et > 0.0 {
                            buf[locals[i] as usize] += et * ws[i] * inv;
                        }
                    }
                }
            } else {
                buf[locals[0] as usize] += count;
            }
        }
    });
    reduce_shards_sparse(shards, plan, alpha_out);
}

/// Minimum fixed-point budget for the persistent worker path.
///
/// The one-step benchmark and short bias-seed runs are better served by the
/// ordinary Rayon kernels: their plan/workspace setup is not amortized. The
/// production point estimate normally runs hundreds or thousands of steps.
pub(crate) const PERSISTENT_MIN_STEPS: u32 = 32;

/// Fork-join is still cheaper for smaller pools; the persistent phase queue is
/// selected where repeated scheduler overhead becomes the measured bottleneck.
pub(crate) const PERSISTENT_MIN_WORKERS: usize = 32;

#[inline]
fn load_f64(cell: &AtomicU64) -> f64 {
    f64::from_bits(cell.load(Ordering::Relaxed))
}

#[inline]
fn store_f64(cell: &AtomicU64, value: f64) {
    cell.store(value.to_bits(), Ordering::Relaxed);
}

fn atomic_f64_vec(values: impl IntoIterator<Item = f64>) -> Vec<AtomicU64> {
    values
        .into_iter()
        .map(|value| AtomicU64::new(value.to_bits()))
        .collect()
}

/// Run an unaccelerated fixed-point EM/VBEM loop as one persistent Rayon
/// broadcast rather than constructing four independent parallel task graphs per
/// iteration.
///
/// Rayon already owns persistent OS threads. `broadcast` places one long-lived
/// closure on each of them. Workers claim fixed vector chunks or complete
/// logical shards from the current phase. The last completed job publishes the
/// next phase, so an idle or preempted worker is not itself a barrier participant.
///
/// The numerical contract is unchanged:
///
/// - shard boundaries and within-shard class order still come from `ShardPlan`;
/// - each transcript still reduces contributors in ascending shard order;
/// - convergence maxima still use fixed 2,048-transcript chunks and an
///   ascending final fold;
/// - worker count affects only ownership and completion time, never grouping.
///
/// Atomic cells provide phase-shared storage without atomic arithmetic. Every
/// cell has exactly one writer in a phase, and readers access it only after the
/// release/acquire phase transition. Relaxed loads/stores therefore compile to
/// ordinary scalar accesses on the supported targets while keeping this
/// implementation data-race-free without an unsafe phase-sharing abstraction.
#[allow(clippy::too_many_arguments)]
pub(crate) fn fixed_point_persistent(
    p: &PackedEqClasses,
    counts: &[u64],
    prior_alphas: &[f64],
    initial_alphas: Vec<f64>,
    use_vbem: bool,
    min_iter: u32,
    max_iter: u32,
    rel_diff_tol: f64,
    alpha_check_cutoff: f64,
) -> (Vec<f64>, u32, bool) {
    if max_iter == 0 {
        return (initial_alphas, 0, false);
    }

    let plan = ShardPlan::new(p);
    let num_workers = rayon::current_num_threads().max(1);
    let num_txps = p.num_txps;
    let num_chunks = num_txps.div_ceil(crate::VECTOR_REDUCTION_CHUNK);

    let alpha = [
        atomic_f64_vec(initial_alphas),
        atomic_f64_vec(std::iter::repeat_n(0.0, num_txps)),
    ];
    let exp_theta = atomic_f64_vec(std::iter::repeat_n(0.0, num_txps));
    let rel_diff_partials = atomic_f64_vec(std::iter::repeat_n(f64::NEG_INFINITY, num_chunks));
    let shard_buffers: Vec<OnceLock<Vec<AtomicU64>>> =
        (0..plan.num_shards()).map(|_| OnceLock::new()).collect();

    let current = AtomicUsize::new(0);
    let iterations = AtomicU32::new(0);
    let log_norm = AtomicU64::new(0.0f64.to_bits());
    // High 32 bits are a monotonically increasing phase generation; low 32
    // bits are the next job within that phase. One CAS therefore prevents a
    // worker paused during an old phase from claiming a new phase's job.
    let work_state = AtomicU64::new(0);
    let completed_jobs = AtomicUsize::new(0);
    let stop = AtomicBool::new(false);
    let converged = AtomicBool::new(false);

    rayon::broadcast(|ctx| {
        debug_assert_eq!(ctx.num_threads(), num_workers);

        loop {
            if stop.load(Ordering::Acquire) {
                break;
            }

            let snapshot = work_state.load(Ordering::Acquire);
            let phase_id = (snapshot >> 32) as u32;
            let next_job = snapshot as u32;
            // VBEM: sum, expectation map, accumulate, reduce, control.
            // EM: accumulation, reduction, control.
            let phase_kind = if use_vbem {
                phase_id % 5
            } else {
                2 + phase_id % 3
            };
            let job_count = match phase_kind {
                0 | 4 => 1,
                1 | 3 => num_chunks,
                2 => plan.num_shards(),
                _ => unreachable!(),
            };
            let job_count = u32::try_from(job_count)
                .expect("persistent EM phase cannot contain more than 2^32 jobs");

            // Empty inputs still need to advance through a zero-job phase.
            if job_count == 0 {
                let next_phase = (u64::from(phase_id) + 1) << 32;
                let _ = work_state.compare_exchange(
                    snapshot,
                    next_phase,
                    Ordering::AcqRel,
                    Ordering::Acquire,
                );
                continue;
            }

            if next_job >= job_count {
                // This worker owns no outstanding job. Yield until the last
                // completed job advances the generation; unlike a barrier, this
                // worker's own scheduling cannot delay that transition.
                'wait: loop {
                    for _ in 0..64 {
                        if stop.load(Ordering::Acquire)
                            || (work_state.load(Ordering::Acquire) >> 32) != u64::from(phase_id)
                        {
                            break 'wait;
                        }
                        std::hint::spin_loop();
                    }
                    std::thread::yield_now();
                }
                continue;
            }

            if work_state
                .compare_exchange_weak(snapshot, snapshot + 1, Ordering::AcqRel, Ordering::Acquire)
                .is_err()
            {
                continue;
            }

            let job = next_job as usize;
            let src_index = current.load(Ordering::Relaxed);
            let dst_index = src_index ^ 1;
            let src = &alpha[src_index];
            let dst = &alpha[dst_index];
            let next_iteration = iterations.load(Ordering::Relaxed) + 1;

            match phase_kind {
                0 => {
                    // Preserve the historical left-to-right VBEM sum.
                    let mut alpha_sum = 0.0;
                    for i in 0..num_txps {
                        alpha_sum += load_f64(&src[i]) + prior_alphas[i];
                    }
                    store_f64(&log_norm, digamma(alpha_sum));
                }
                1 => {
                    let lo = job * crate::VECTOR_REDUCTION_CHUNK;
                    let hi = (lo + crate::VECTOR_REDUCTION_CHUNK).min(num_txps);
                    let norm = load_f64(&log_norm);
                    for i in lo..hi {
                        let ap = load_f64(&src[i]) + prior_alphas[i];
                        let value = if ap > DIGAMMA_MIN {
                            (digamma(ap) - norm).exp()
                        } else {
                            0.0
                        };
                        store_f64(&exp_theta[i], value);
                    }
                }
                2 => {
                    let shard = job;
                    let buf = shard_buffers[shard].get_or_init(|| {
                        atomic_f64_vec(std::iter::repeat_n(0.0, plan.local_lens[shard]))
                    });
                    for value in buf {
                        store_f64(value, 0.0);
                    }

                    let (class_lo, class_hi) = plan.range(shard);
                    for ci in class_lo..class_hi {
                        let count = counts[ci] as f64;
                        let (tids, ws) = p.class(ci);
                        let locals = plan.entry_locals(p, ci);
                        if tids.len() > 1 {
                            let mut denom = 0.0;
                            if use_vbem {
                                for (&tid, &weight) in tids.iter().zip(ws) {
                                    let et = load_f64(&exp_theta[tid as usize]);
                                    if et > 0.0 {
                                        denom += et * weight;
                                    }
                                }
                            } else {
                                for (&tid, &weight) in tids.iter().zip(ws) {
                                    denom += load_f64(&src[tid as usize]) * weight;
                                }
                            }

                            if denom > MIN_EQ_CLASS_WEIGHT {
                                let inv = count / denom;
                                for i in 0..tids.len() {
                                    let basis = if use_vbem {
                                        load_f64(&exp_theta[tids[i] as usize])
                                    } else {
                                        load_f64(&src[tids[i] as usize])
                                    };
                                    let value = if use_vbem && basis <= 0.0 {
                                        continue;
                                    } else {
                                        basis * ws[i]
                                    };
                                    if !use_vbem && value.is_nan() {
                                        continue;
                                    }
                                    let cell = &buf[locals[i] as usize];
                                    store_f64(cell, load_f64(cell) + value * inv);
                                }
                            }
                        } else {
                            let cell = &buf[locals[0] as usize];
                            store_f64(cell, load_f64(cell) + count);
                        }
                    }
                }
                3 => {
                    // Fuse sparse reduction with the convergence scan. Output
                    // jobs retain the fixed 2,048-transcript boundaries used by
                    // `max_rel_diff_par`, independent of worker count.
                    let lo = job * crate::VECTOR_REDUCTION_CHUNK;
                    let hi = (lo + crate::VECTOR_REDUCTION_CHUNK).min(num_txps);
                    let mut max_d = f64::NEG_INFINITY;
                    for tid in lo..hi {
                        let contrib_lo = plan.contrib_start[tid] as usize;
                        let contrib_hi = plan.contrib_start[tid + 1] as usize;
                        let mut sum = 0.0;
                        for contributor in &plan.contributors[contrib_lo..contrib_hi] {
                            let shard = shard_buffers[contributor.shard as usize]
                                .get()
                                .expect("the accumulation phase initialized every shard");
                            sum += load_f64(&shard[contributor.local as usize]);
                        }
                        store_f64(&dst[tid], sum);

                        if next_iteration >= min_iter {
                            let before = load_f64(&src[tid]);
                            if before > alpha_check_cutoff && sum > 0.0 {
                                let d = (sum - before).abs() / sum;
                                if d > max_d {
                                    max_d = d;
                                }
                            }
                        }
                    }
                    store_f64(&rel_diff_partials[job], max_d);
                }
                4 => {
                    let mut max_d = f64::NEG_INFINITY;
                    if next_iteration >= min_iter {
                        for partial in &rel_diff_partials {
                            let d = load_f64(partial);
                            if d > max_d {
                                max_d = d;
                            }
                        }
                    }

                    iterations.store(next_iteration, Ordering::Relaxed);
                    current.store(dst_index, Ordering::Relaxed);
                    if next_iteration >= max_iter || (max_d.is_finite() && max_d < rel_diff_tol) {
                        converged
                            .store(max_d.is_finite() && max_d < rel_diff_tol, Ordering::Relaxed);
                        stop.store(true, Ordering::Release);
                    }
                }
                _ => unreachable!(),
            }

            if completed_jobs.fetch_add(1, Ordering::AcqRel) + 1 == job_count as usize {
                completed_jobs.store(0, Ordering::Relaxed);
                work_state.store((u64::from(phase_id) + 1) << 32, Ordering::Release);
            }
        }
    });

    let final_index = current.load(Ordering::Relaxed);
    let result = alpha[final_index].iter().map(load_f64).collect();
    (
        result,
        iterations.load(Ordering::Relaxed),
        converged.load(Ordering::Relaxed),
    )
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

#[cfg(test)]
mod shard_plan_determinism {
    use super::*;
    use crate::{EmOptions, EmResult};
    use salmon_eqclass::{EquivalenceClassBuilder, TranscriptGroup};

    fn fixture(num_txps: usize, num_classes: usize) -> (PackedEqClasses, Vec<f64>) {
        let b = EquivalenceClassBuilder::new();
        // Values chosen so the sums are not exactly representable: with tidy
        // powers of two, every grouping agrees and the test cannot fail.
        let mut seed = 0x9E37_79B9_7F4A_7C15u64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for c in 0..num_classes {
            let n = 2 + (c % 4);
            let mut tids: Vec<u32> = (0..n).map(|_| (rnd() % num_txps as u64) as u32).collect();
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
        (PackedEqClasses::from_collapsed(&eq, num_txps), eff)
    }

    fn assert_same_result(label: &str, threads: usize, base: &EmResult, other: &EmResult) {
        assert_eq!(
            base.iters, other.iters,
            "{label}: iteration count at t{threads}"
        );
        assert_eq!(
            base.converged, other.converged,
            "{label}: convergence flag at t{threads}"
        );
        assert_eq!(
            base.dropped_mass.to_bits(),
            other.dropped_mass.to_bits(),
            "{label}: dropped mass at t{threads}"
        );
        assert_eq!(
            base.alphas.len(),
            other.alphas.len(),
            "{label}: thread count changed the result shape"
        );
        let diffs = base
            .alphas
            .iter()
            .zip(&other.alphas)
            .filter(|(a, b)| a.to_bits() != b.to_bits())
            .count();
        assert_eq!(
            diffs,
            0,
            "{label}: result differs between 1 and {threads} threads in {diffs} of {} \
             transcripts",
            base.alphas.len()
        );
    }

    fn dense_em_reference(
        p: &PackedEqClasses,
        counts: &[u64],
        alpha_in: &[f64],
        plan: &ShardPlan,
    ) -> Vec<f64> {
        let mut shards = vec![vec![0.0; p.num_txps]; plan.num_shards()];
        for (s, buf) in shards.iter_mut().enumerate() {
            let (start, end) = plan.range(s);
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
        }
        (0..p.num_txps)
            .map(|tid| shards.iter().fold(0.0, |sum, buf| sum + buf[tid]))
            .collect()
    }

    fn dense_vbem_reference(
        p: &PackedEqClasses,
        counts: &[u64],
        prior: &[f64],
        alpha_in: &[f64],
        plan: &ShardPlan,
    ) -> Vec<f64> {
        let mut exp_theta = vec![0.0; p.num_txps];
        fill_exp_theta_par(alpha_in, prior, &mut exp_theta);
        let mut shards = vec![vec![0.0; p.num_txps]; plan.num_shards()];
        for (s, buf) in shards.iter_mut().enumerate() {
            let (start, end) = plan.range(s);
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
        }
        (0..p.num_txps)
            .map(|tid| shards.iter().fold(0.0, |sum, buf| sum + buf[tid]))
            .collect()
    }

    /// The EM result must not depend on the thread count.
    ///
    /// It used to: the shard count was `rayon::current_num_threads()`, so `-p`
    /// changed the class partition and therefore the grouping of the
    /// floating-point sums. On 26M real fragments this produced different
    /// `quant.sf` files at `-p 4` and `-p 32` — the guarantee the deterministic
    /// mode exists to provide.
    ///
    /// The existing suite missed it because its fixtures were far too small:
    /// with a handful of classes every partition is the same partition. This
    /// builds enough incidences to cross `MIN_INCIDENCES_PER_SHARD` several
    /// times, so the plan really does shard, and enough transcript overlap that
    /// the per-shard partial sums genuinely differ in grouping.
    #[test]
    fn em_result_is_independent_of_thread_count() {
        const NUM_TXPS: usize = 5_000;
        const NUM_CLASSES: usize = 20_000;
        let (packed, eff) = fixture(NUM_TXPS, NUM_CLASSES);

        // Cover every kernel the shard plan feeds, not just the one
        // `EmOptions::default()` happens to select. The default is plain EM
        // without acceleration, but production runs VBEM (`opt_type: "vb"`),
        // and `--emAccel squarem` drives the same M-steps through a different
        // iterate sequence. A test that exercised only the default would have
        // left the actual shipping path unguarded — which is how the original
        // defect survived.
        let configs = [
            ("plain EM", false, false, crate::EmAccel::None),
            (
                "VBEM (production default)",
                true,
                false,
                crate::EmAccel::None,
            ),
            (
                "VBEM per-nucleotide prior",
                true,
                true,
                crate::EmAccel::None,
            ),
            ("plain EM + SQUAREM", false, false, crate::EmAccel::Squarem),
            ("VBEM + SQUAREM", true, false, crate::EmAccel::Squarem),
            ("plain EM + DAAREM", false, false, crate::EmAccel::Daarem),
            ("VBEM + DAAREM", true, false, crate::EmAccel::Daarem),
        ];

        for (label, use_vbem, per_nucleotide_prior, accel) in configs {
            let opts = EmOptions {
                max_iter: 200,
                use_vbem,
                per_nucleotide_prior,
                accel,
                ..EmOptions::default()
            };

            // Run the same problem inside rayon pools of different sizes. Only
            // the scheduling differs; the arithmetic must not.
            let run = |threads: usize| -> EmResult {
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(threads)
                    .build()
                    .expect("thread pool");
                pool.install(|| {
                    crate::optimize_packed_with_init(
                        &packed,
                        &opts,
                        true,
                        crate::InitAlphas::NONE,
                        crate::EffLens::new(&eff),
                    )
                })
            };

            let base = run(1);
            for threads in [2usize, 3, 8, 16, 32] {
                let other = run(threads);
                assert_same_result(label, threads, &base, &other);
            }
        }
    }

    #[test]
    fn compressed_plan_indices_and_ranges_are_consistent() {
        let (packed, _) = fixture(3_000, 12_000);
        let plan = ShardPlan::new(&packed);
        assert_eq!(plan.boundaries[0], 0);
        assert_eq!(*plan.boundaries.last().unwrap(), packed.num_classes());
        assert!(plan.boundaries.windows(2).all(|w| w[0] < w[1]));
        assert_eq!(plan.entry_local.len(), packed.labels.len());

        // Invert contributor metadata back to (shard, local)->tid, then verify
        // every incidence's precomputed local index reaches its original tid.
        let mut local_tid: Vec<Vec<Option<u32>>> =
            plan.local_lens.iter().map(|&len| vec![None; len]).collect();
        for tid in 0..packed.num_txps {
            let lo = plan.contrib_start[tid] as usize;
            let hi = plan.contrib_start[tid + 1] as usize;
            let contributors = &plan.contributors[lo..hi];
            assert!(contributors.windows(2).all(|w| w[0].shard < w[1].shard));
            for contributor in contributors {
                local_tid[contributor.shard as usize][contributor.local as usize] =
                    Some(tid as u32);
            }
        }
        for s in 0..plan.num_shards() {
            let (start, end) = plan.range(s);
            for ci in start..end {
                let (tids, _) = packed.class(ci);
                let locals = plan.entry_locals(&packed, ci);
                for (&tid, &local) in tids.iter().zip(locals) {
                    assert_eq!(local_tid[s][local as usize], Some(tid));
                }
            }
        }
    }

    #[test]
    fn compressed_steps_match_dense_reference_for_the_same_plan() {
        let (packed, _) = fixture(3_000, 12_000);
        let plan = ShardPlan::new(&packed);
        let alpha: Vec<f64> = (0..packed.num_txps)
            .map(|i| 0.25 + (i % 31) as f64 / 17.0)
            .collect();
        let prior: Vec<f64> = (0..packed.num_txps)
            .map(|i| 0.01 + (i % 7) as f64 * 1e-4)
            .collect();

        let em_reference = dense_em_reference(&packed, &packed.counts, &alpha, &plan);
        let mut em_out = vec![0.0; packed.num_txps];
        let mut em_shards = plan.allocate_buffers();
        em_step_par(
            &packed,
            &packed.counts,
            &alpha,
            &mut em_out,
            &mut em_shards,
            &plan,
        );
        assert!(em_reference
            .iter()
            .zip(&em_out)
            .all(|(a, b)| a.to_bits() == b.to_bits()));

        let vb_reference = dense_vbem_reference(&packed, &packed.counts, &prior, &alpha, &plan);
        let mut vb_out = vec![0.0; packed.num_txps];
        let mut exp_theta = vec![0.0; packed.num_txps];
        let mut vb_shards = plan.allocate_buffers();
        vbem_step_par(
            &packed,
            &packed.counts,
            &prior,
            &alpha,
            &mut vb_out,
            &mut exp_theta,
            &mut vb_shards,
            &plan,
        );
        assert!(vb_reference
            .iter()
            .zip(&vb_out)
            .all(|(a, b)| a.to_bits() == b.to_bits()));
    }

    #[test]
    fn persistent_loop_matches_repeated_rayon_steps_bit_for_bit() {
        const STEPS: u32 = 50;
        let (packed, _) = fixture(3_000, 12_000);
        let initial: Vec<f64> = (0..packed.num_txps)
            .map(|i| 0.25 + (i % 31) as f64 / 17.0)
            .collect();
        let prior: Vec<f64> = (0..packed.num_txps)
            .map(|i| 0.01 + (i % 7) as f64 * 1e-4)
            .collect();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build()
            .expect("thread pool");

        for use_vbem in [false, true] {
            let reference = pool.install(|| {
                let plan = ShardPlan::new(&packed);
                let mut shards = plan.allocate_buffers();
                let mut exp_theta = vec![0.0; packed.num_txps];
                let mut src = initial.clone();
                let mut dst = vec![0.0; packed.num_txps];
                for _ in 0..STEPS {
                    if use_vbem {
                        vbem_step_par(
                            &packed,
                            &packed.counts,
                            &prior,
                            &src,
                            &mut dst,
                            &mut exp_theta,
                            &mut shards,
                            &plan,
                        );
                    } else {
                        em_step_par(&packed, &packed.counts, &src, &mut dst, &mut shards, &plan);
                    }
                    std::mem::swap(&mut src, &mut dst);
                }
                src
            });

            let (persistent, iterations, converged) = pool.install(|| {
                fixed_point_persistent(
                    &packed,
                    &packed.counts,
                    &prior,
                    initial.clone(),
                    use_vbem,
                    STEPS,
                    STEPS,
                    -1.0,
                    0.01,
                )
            });
            assert_eq!(iterations, STEPS);
            assert!(!converged);
            assert!(reference
                .iter()
                .zip(&persistent)
                .all(|(a, b)| a.to_bits() == b.to_bits()));
        }
    }

    #[test]
    fn empty_plan_has_one_empty_shard() {
        let b = EquivalenceClassBuilder::new();
        let mut eq = b.finish();
        eq.update_eff_lengths(&[1.0; 4]);
        let packed = PackedEqClasses::from_collapsed(&eq, 4);
        let plan = ShardPlan::new(&packed);
        assert_eq!(plan.boundaries, [0, 0]);
        assert_eq!(plan.local_lens, [0]);
        assert!(plan.entry_local.is_empty());
        assert!(plan.contributors.is_empty());
        assert_eq!(plan.contrib_start, [0; 5]);
    }

    #[test]
    fn persistent_loop_handles_empty_classes_and_zero_job_phases() {
        let b = EquivalenceClassBuilder::new();
        let mut eq = b.finish();
        eq.update_eff_lengths(&[1.0; 4]);
        let packed = PackedEqClasses::from_collapsed(&eq, 4);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("thread pool");

        for use_vbem in [false, true] {
            let (result, iterations, converged) = pool.install(|| {
                fixed_point_persistent(
                    &packed,
                    &packed.counts,
                    &[0.01; 4],
                    vec![1.0; 4],
                    use_vbem,
                    1,
                    3,
                    0.01,
                    0.01,
                )
            });
            assert_eq!(result, [0.0; 4]);
            assert_eq!(iterations, 3);
            assert!(!converged);
        }

        let b = EquivalenceClassBuilder::new();
        let mut eq = b.finish();
        eq.update_eff_lengths(&[]);
        let packed = PackedEqClasses::from_collapsed(&eq, 0);
        let (result, iterations, converged) = pool.install(|| {
            fixed_point_persistent(
                &packed,
                &packed.counts,
                &[],
                Vec::new(),
                false,
                1,
                2,
                0.01,
                0.01,
            )
        });
        assert!(result.is_empty());
        assert_eq!(iterations, 2);
        assert!(!converged);
    }

    /// The plan's shard count must come from the data, never from the pool.
    #[test]
    fn shard_count_ignores_the_thread_pool() {
        let b = EquivalenceClassBuilder::new();
        for c in 0..12_000u32 {
            b.add_group(
                TranscriptGroup::new(vec![c % 500, 500 + c % 300]),
                vec![1.0, 2.0],
                3,
            );
        }
        let mut eq = b.finish();
        eq.update_eff_lengths(&vec![1000.0; 900]);
        let packed = PackedEqClasses::from_collapsed(&eq, 900);

        let n = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| ShardPlan::new(&packed).num_shards())
        };
        let base = n(1);
        for t in [2usize, 7, 32] {
            assert_eq!(n(t), base, "shard count changed with a {t}-thread pool");
        }
    }
}
