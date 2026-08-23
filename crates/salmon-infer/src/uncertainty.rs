//! Posterior uncertainty: multinomial **bootstrap** (`CollapsedEMOptimizer::gatherBootstraps`)
//! and the non-collapsed **Gibbs sampler** (`CollapsedGibbsSampler`).
//!
//! # Why uncertainty at all
//!
//! The EM returns a single number per transcript, but that number is not equally
//! trustworthy for every transcript. A transcript with plenty of unique evidence
//! is pinned down precisely; one that shares all its sequence with a paralogue is
//! barely constrained, and its estimate could shift a long way given slightly
//! different reads. Downstream differential-expression tools need to know which
//! is which, so salmon can produce a *distribution* of plausible estimates
//! instead of only the point estimate.
//!
//! Two methods are offered, answering slightly different questions:
//!
//! * **Bootstrap** — "if I had sequenced this sample again, how much would the
//!   answer move?". Resample the observed fragment counts with replacement and
//!   re-run the EM. Spread across replicates estimates sampling variability.
//! * **Gibbs sampling** — "given this data, what is the posterior distribution of
//!   abundances?". Repeatedly reassign each ambiguous fragment at random in
//!   proportion to current abundances, and redraw abundances from the resulting
//!   counts. After enough rounds the chain wanders the posterior, and snapshots
//!   of it are posterior samples.
//!
//! Both operate on the flat [`PackedEqClasses`] layout and produce one abundance
//! vector per replicate/sample. Bootstrap parallelizes across replicates (each
//! runs a sequential EM on resampled class counts); Gibbs parallelizes across
//! independent chains (each runs sequential thinned rounds) — in both cases the
//! unit of parallelism is an entirely independent computation, which is the
//! cleanest kind. RNG is PCG (`rand_pcg`), seeded per replicate/chain for
//! reproducibility — results are statistically equivalent to salmon's but not
//! bit-identical (different RNG).

use rand::{Rng, SeedableRng};
use rand_distr::{Binomial, Distribution, Gamma};
use rand_pcg::Pcg64Mcg;
use rayon::prelude::*;

use crate::packed::PackedEqClasses;
use crate::{run_em_counts, EmOptions};

/// Smallest class denominator below which mass is redistributed evenly (Gibbs).
const MIN_EQ_CLASS_WEIGHT: f64 = f64::MIN_POSITIVE;

/// Draw a multinomial count vector: `total` draws over categories with the given
/// (non-normalized) `weights`, via the conditional-binomial method — `O(k)` in
/// the number of categories rather than `O(total)` individual draws.
///
/// **The conditional-binomial trick.** Throwing `total` balls into `k` buckets
/// one at a time would cost `O(total)`, and `total` is the fragment count —
/// hundreds of millions. Instead, ask "how many of the remaining balls land in
/// bucket 0?", which is a single binomial draw, then recurse on the rest with
/// bucket 0 removed and the probabilities renormalized. That is `k` draws total,
/// and the result has exactly the right joint distribution.
fn multinomial(total: u64, weights: &[f64], rng: &mut impl Rng) -> Vec<u64> {
    let n = weights.len();
    let mut out = vec![0u64; n];
    if total == 0 || n == 0 {
        return out;
    }
    let mut remaining = total;
    let mut remaining_w: f64 = weights.iter().sum();
    for i in 0..n {
        // All balls placed; the rest of the buckets stay empty.
        if remaining == 0 {
            break;
        }
        // Last bucket (or no weight left anywhere): it takes the remainder by
        // construction, which also guarantees the counts sum to `total` exactly.
        if i == n - 1 || remaining_w <= 0.0 {
            out[i] = remaining;
            break;
        }
        // Probability this bucket takes a ball, among those not yet allocated.
        // Clamped because the running subtraction of `remaining_w` accumulates
        // rounding error.
        let p = (weights[i] / remaining_w).clamp(0.0, 1.0);
        let k = if p >= 1.0 {
            remaining
        } else if p <= 0.0 {
            0
        } else {
            Binomial::new(remaining, p).unwrap().sample(rng)
        };
        out[i] = k;
        remaining -= k;
        remaining_w -= weights[i];
    }
    out
}

/// Run `num_bootstraps` multinomial bootstrap replicates. Each resamples the
/// per-class counts (multinomial over the original counts, `total_count` draws),
/// runs EM/VBEM to convergence (min 50 iters), and rescales the abundances to sum
/// to the input total `p.total_count` (the summed equivalence-class counts that
/// were resampled). The rescale restores the mass the `min_alpha` truncation
/// below removes, so each replicate sums to exactly what the point estimate does
/// (`Σ alphas == p.total_count`). We rescale to this *intrinsic* total — not an
/// externally-computed mapped-fragment counter — so the replicates stay
/// self-consistent with the resampled input regardless of how (or whether)
/// `num_mapped` is computed elsewhere. (This is also unconditional, unlike C++'s
/// mode-gated `useScaledCounts`, which corrects prior-inflated VBEM alphas we do
/// not have.) Returns one abundance vector per replicate.
///
/// Resampling *classes* rather than individual fragments is what makes this
/// affordable: a class is a group of identical fragments, so one multinomial draw
/// over classes is statistically the same as drawing every fragment.
pub fn bootstrap(
    p: &PackedEqClasses,
    opts: &EmOptions,
    eff_lens: Option<&[f64]>,
    num_bootstraps: u32,
    seed: u64,
) -> Vec<Vec<f64>> {
    let sample_weights: Vec<f64> = p.counts.iter().map(|&c| c as f64).collect();
    let total = p.total_count;
    (0..num_bootstraps)
        .into_par_iter()
        .map(|bs| {
            // Seed derived from (run seed, replicate index) so results reproduce
            // exactly and do not depend on which thread ran which replicate. The
            // odd multiplier scatters adjacent indices to distant seeds.
            let mut rng =
                Pcg64Mcg::seed_from_u64(seed ^ (bs as u64).wrapping_mul(0x9E3779B97F4A7C15));
            let resampled = multinomial(total, &sample_weights, &mut rng);
            // Effective lengths reach the replicate EM so `--perNucleotidePrior`
            // shapes the prior here exactly as it does the point estimate;
            // passing `None` silently fell back to the flat prior, which made
            // the flag a no-op inside every replicate (#1140, audit D13).
            //
            // NOTE the argument order: `init_alphas` precedes `eff_lens` and
            // both are `Option<&[f64]>`, so swapping them type-checks. The
            // first attempt at this fix did exactly that — the prior stayed
            // flat AND every replicate warm-started from the effective-length
            // vector. Keep `None` explicit for the uniform start.
            let (alphas, _, _) = run_em_counts(p, &resampled, opts, false, 50, None, eff_lens);
            // Finalize like the point estimate: truncate the negligible
            // abundances, then redistribute that mass to eq-class co-members via a
            // masked final M-step over the *resampled* counts (no rescale-up). The
            // replicate's `dropped` mass is negligible and not reported (only the
            // point estimate surfaces `inference_truncated_mass`).
            //
            // Finalizing identically to the point estimate is what makes the
            // spread interpretable: any systematic difference in post-processing
            // would show up as bias, not variance.
            let (alphas, _dropped) =
                crate::finalize_truncate_redistribute(p, &resampled, alphas, opts, eff_lens);
            alphas
        })
        .collect()
}

/// Gibbs sampling parameters (salmon defaults).
#[derive(Debug, Clone)]
pub struct GibbsOptions {
    /// number of posterior samples to draw
    pub num_samples: u32,
    /// internal thinning rounds between recorded samples (salmon default 16)
    ///
    /// Consecutive states of a Markov chain are highly correlated, so recording
    /// every one would give many samples carrying little independent information.
    /// Keeping only every 16th decorrelates them.
    pub thinning: u32,
    /// base prior value (per-transcript, or per-nucleotide when `!per_transcript_prior`)
    pub prior: f64,
    /// whether the prior is per-transcript (else scaled by effective length)
    pub per_transcript_prior: bool,
}

impl Default for GibbsOptions {
    fn default() -> Self {
        Self {
            num_samples: 0,
            thinning: 16,
            prior: 1e-3,
            per_transcript_prior: true,
        }
    }
}

/// salmon's Gibbs rate parameter `beta`.
/// The rate of the Gamma prior on each transcript's expression rate.
const GIBBS_BETA: f64 = 0.1;

/// One Gibbs round (salmon's `sampleRoundNonCollapsedMultithreaded_`): draw the
/// transcript fractions `mu` from their Gamma posterior, then resample each
/// equivalence class's count multinomially across its transcripts.
///
/// These are the two conditional distributions the sampler alternates between:
/// abundances given a fragment assignment, and a fragment assignment given
/// abundances. Alternating between the full conditionals of a model is exactly
/// what makes this a Gibbs sampler.
///
/// The Gamma appears because it is the conjugate prior for a Poisson count: with
/// a Gamma prior and Poisson-distributed observed counts, the posterior is again
/// a Gamma, with the observed count folded into its shape parameter.
#[allow(clippy::too_many_arguments)]
fn gibbs_round(
    p: &PackedEqClasses,
    active: &[u32],
    eff_lens: &[f64],
    prior_alphas: &[f64],
    txp_count: &mut [f64],
    mu: &mut [f64],
    rng: &mut impl Rng,
) {
    // Sample mu[i] ~ Gamma(txpCount[i] + prior[i], 1/(beta + effLen[i])); reset count.
    //
    // Shape = observed count + prior (more reads ⇒ tighter, higher draw); scale
    // shrinks with effective length, converting counts back to a rate per
    // available position. The counts are zeroed as we go because the second half
    // of the round refills them from scratch.
    for &i in active {
        let i = i as usize;
        let ci = txp_count[i] + prior_alphas[i];
        let scale = 1.0 / (GIBBS_BETA + eff_lens[i]);
        mu[i] = if ci > 0.0 {
            Gamma::new(ci, scale).unwrap().sample(rng)
        } else {
            0.0
        };
        txp_count[i] = 0.0;
    }
    // Resample each class's reads across its transcripts.
    let mut probs: Vec<f64> = Vec::with_capacity(64);
    for ci in 0..p.num_classes() {
        let class_count = p.counts[ci];
        let s = p.starts[ci] as usize;
        let e = p.starts[ci + 1] as usize;
        let tids = &p.labels[s..e];
        // Note: the *raw* weights here, not the effective-length-divided
        // `combined` ones — the length correction already lives in `mu`'s scale.
        let weights = &p.weights[s..e];
        if tids.len() > 1 {
            probs.clear();
            let mut denom = 0.0;
            for (&tid, &w) in tids.iter().zip(weights) {
                // The 1000× factor only keeps the products away from denormal
                // territory; the multinomial normalizes them anyway.
                let v = 1000.0 * mu[tid as usize] * w;
                probs.push(v);
                denom += v;
            }
            if denom <= MIN_EQ_CLASS_WEIGHT {
                // fall back to uniform over the class
                //
                // Every member drew a vanishing rate; assigning uniformly keeps
                // the chain moving instead of stalling on a zero denominator.
                for v in probs.iter_mut() {
                    *v = 1.0;
                }
            }
            let draws = multinomial(class_count, &probs, rng);
            for (&tid, &k) in tids.iter().zip(&draws) {
                txp_count[tid as usize] += k as f64;
            }
        } else {
            // Unambiguous class: nothing random about where its reads go.
            txp_count[tids[0] as usize] += class_count as f64;
        }
    }
}

/// Draw `opts.num_samples` Gibbs posterior samples. `init_alphas` is the point
/// estimate (EM result) each chain restarts from; `eff_lens` the effective
/// lengths. Chains (1/2/4/8 by sample count, like salmon) run in parallel.
/// Returns one abundance vector per sample (scaled to the input total
/// `p.total_count`, matching the point estimate — not an external mapped-fragment
/// counter).
///
/// Starting every chain from the converged EM estimate is what allows the burn-in
/// to be skipped: the chain begins in a high-probability region rather than
/// having to walk there.
pub fn gibbs_sample(
    p: &PackedEqClasses,
    eff_lens: &[f64],
    init_alphas: &[f64],
    opts: &GibbsOptions,
    seed: u64,
) -> Vec<Vec<f64>> {
    let num_txps = p.num_txps;
    let num_samples = opts.num_samples as usize;
    if num_samples == 0 {
        return Vec::new();
    }

    // Active transcripts = those appearing in some class.
    //
    // A transcript no fragment is compatible with has no posterior to sample;
    // restricting the per-round loop to active ones skips them entirely, which
    // matters because most of a transcriptome is typically unexpressed.
    let mut active_flag = vec![false; num_txps];
    for &t in &p.labels {
        active_flag[t as usize] = true;
    }
    let active: Vec<u32> = (0..num_txps as u32)
        .filter(|&t| active_flag[t as usize])
        .collect();

    // Per-transcript prior (per-txp = constant; per-nucleotide = prior·max(1,effLen)).
    //
    // The per-nucleotide form makes the prior's strength proportional to how much
    // sequence a transcript offers, so it does not dominate short transcripts.
    let prior_alphas: Vec<f64> = (0..num_txps)
        .map(|i| {
            if opts.per_transcript_prior {
                opts.prior
            } else {
                opts.prior * eff_lens[i].max(1.0)
            }
        })
        .collect();

    // Initial counts (0 for inactive transcripts).
    let mut init = init_alphas.to_vec();
    for i in 0..num_txps {
        if !active_flag[i] {
            init[i] = 0.0;
        }
    }

    // Chain layout: salmon uses 1/2/4/8 chains by sample count.
    //
    // More chains means more parallelism, but each still pays its own thinning
    // cost per sample, so splitting a small request across many chains would be
    // pure overhead.
    let nchains: usize = if num_samples >= 200 {
        8
    } else if num_samples >= 100 {
        4
    } else if num_samples >= 50 {
        2
    } else {
        1
    };
    let step = num_samples / nchains;
    // chain c produces samples [c*step .. c*step+len_c)
    // The last chain absorbs the remainder of an uneven division.
    let bounds: Vec<(usize, usize)> = (0..nchains)
        .map(|c| {
            let start = c * step;
            let end = if c == nchains - 1 {
                num_samples
            } else {
                (c + 1) * step
            };
            (start, end)
        })
        .collect();

    let mut all: Vec<Vec<f64>> = vec![Vec::new(); num_samples];
    // Run chains in parallel; each writes its contiguous block.
    //
    // Each chain returns its block plus the offset to place it at, so the final
    // ordering is by sample index and not by which chain finished first.
    let blocks: Vec<(usize, Vec<Vec<f64>>)> = bounds
        .par_iter()
        .enumerate()
        .map(|(c, &(start, end))| {
            let mut rng =
                Pcg64Mcg::seed_from_u64(seed ^ (c as u64).wrapping_mul(0xD1B54A32D192ED03));
            let mut txp_count = init.clone();
            let mut mu = vec![0.0f64; num_txps];
            let mut out: Vec<Vec<f64>> = Vec::with_capacity(end - start);
            for _ in start..end {
                // Advance the chain `thinning` steps, then record one state.
                for _ in 0..opts.thinning {
                    gibbs_round(
                        p,
                        &active,
                        eff_lens,
                        &prior_alphas,
                        &mut txp_count,
                        &mut mu,
                        &mut rng,
                    );
                }
                // Extrapolate counts from the final fractions mu, then normalize
                // to total_count. We TRUNCATE the negligible rate values FIRST and
                // normalize the survivors, so each sample sums to *exactly*
                // total_count with no mass lost to a post-normalization truncation.
                // The normalization is intrinsic to converting the μ rate back to
                // counts (and is paradox-free: μ is anchored to the conserving
                // discrete assignment and the scale factor is ≤ 1).
                let mut sample = vec![0.0f64; num_txps];
                let mut denom = 0.0;
                for t in 0..num_txps {
                    // rate × available positions = expected fragment count.
                    let ext = mu[t] * eff_lens[t];
                    if ext > 1e-8 {
                        sample[t] = ext;
                        denom += ext;
                    }
                }
                if denom > 0.0 {
                    let scale = p.total_count as f64 / denom;
                    for s in &mut sample {
                        *s *= scale;
                    }
                }
                out.push(sample);
            }
            (start, out)
        })
        .collect();
    for (start, out) in blocks {
        for (j, s) in out.into_iter().enumerate() {
            all[start + j] = s;
        }
    }
    all
}

/// Per-transcript unique / ambiguous fragment counts (salmon's `ambig_info.tsv`):
/// `unique[t]` sums counts of single-transcript classes for `t`; `ambig[t]` sums
/// counts of every multi-transcript class containing `t`.
///
/// A cheap, deterministic proxy for confidence: a transcript whose evidence is
/// mostly unique is well-determined, one whose evidence is mostly ambiguous is
/// not. Note `ambig` deliberately double-counts across transcripts — each member
/// of a shared class is credited the whole class count, because the question
/// being answered is per transcript.
pub fn ambiguity_counts(p: &PackedEqClasses) -> (Vec<u32>, Vec<u32>) {
    let mut unique = vec![0u32; p.num_txps];
    let mut ambig = vec![0u32; p.num_txps];
    for ci in 0..p.num_classes() {
        let s = p.starts[ci] as usize;
        let e = p.starts[ci + 1] as usize;
        let tids = &p.labels[s..e];
        let count = p.counts[ci] as u32;
        if tids.len() > 1 {
            for &t in tids {
                ambig[t as usize] += count;
            }
        } else {
            unique[tids[0] as usize] += count;
        }
    }
    (unique, ambig)
}

#[cfg(test)]
mod tests {
    use super::*;
    use salmon_eqclass::{EquivalenceClassBuilder, TranscriptGroup};

    /// Uniform weights and unit effective lengths, so the tests measure the
    /// resampling behaviour and nothing else.
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

    /// With only unique evidence the bootstrap must be centred on the point
    /// estimate and must conserve the total in every replicate.
    #[test]
    fn bootstrap_mean_near_point_estimate() {
        // unique evidence -> every bootstrap recovers ~the same counts
        let p = packed(&[(vec![0], 300), (vec![1], 700)], 2);
        let bs = bootstrap(&p, &EmOptions::default(), None, 50, 12345);
        assert_eq!(bs.len(), 50);
        let m0: f64 = bs.iter().map(|b| b[0]).sum::<f64>() / 50.0;
        let m1: f64 = bs.iter().map(|b| b[1]).sum::<f64>() / 50.0;
        // means within a few % of the point estimate, totals conserved
        assert!((m0 - 300.0).abs() < 30.0, "m0={m0}");
        assert!((m1 - 700.0).abs() < 30.0, "m1={m1}");
        for b in &bs {
            assert!(((b[0] + b[1]) - 1000.0).abs() < 1e-6);
        }
    }

    /// The whole point of the bootstrap: ambiguity must show up as spread.
    #[test]
    fn bootstrap_variance_grows_with_ambiguity() {
        // a fully shared class has higher per-transcript bootstrap variance
        let p = packed(&[(vec![0], 10), (vec![1], 10), (vec![0, 1], 980)], 2);
        let bs = bootstrap(&p, &EmOptions::default(), None, 100, 7);
        let m0: f64 = bs.iter().map(|b| b[0]).sum::<f64>() / 100.0;
        let var0: f64 = bs.iter().map(|b| (b[0] - m0).powi(2)).sum::<f64>() / 100.0;
        assert!(
            var0 > 0.0,
            "ambiguous transcript should have nonzero bootstrap variance"
        );
    }

    /// Gibbs samples are stochastic, so this checks the invariant that must hold
    /// for every sample — the scale — rather than the values themselves.
    #[test]
    fn gibbs_runs_and_conserves_scale() {
        let p = packed(&[(vec![0], 300), (vec![1], 700)], 2);
        let opts = GibbsOptions {
            num_samples: 20,
            thinning: 8,
            ..Default::default()
        };
        let samples = gibbs_sample(&p, &[1.0, 1.0], &[300.0, 700.0], &opts, 99);
        assert_eq!(samples.len(), 20);
        for s in &samples {
            let tot = s[0] + s[1];
            assert!(
                (tot - 1000.0).abs() < 50.0,
                "gibbs total {tot} not near 1000"
            );
        }
    }

    /// Unique counts are per transcript; ambiguous counts credit the full class
    /// count to every member, hence 100 twice rather than 50 each.
    #[test]
    fn ambiguity_counts_split() {
        let p = packed(&[(vec![0], 30), (vec![1], 70), (vec![0, 1], 100)], 2);
        let (uniq, amb) = ambiguity_counts(&p);
        assert_eq!(uniq, vec![30, 70]);
        assert_eq!(amb, vec![100, 100]);
    }
}

#[cfg(test)]
mod per_nucleotide_prior_reaches_replicates {
    use super::*;
    use salmon_eqclass::{EquivalenceClassBuilder, TranscriptGroup};

    /// `--perNucleotidePrior` must shape the prior inside every bootstrap
    /// replicate, not only the point estimate.
    ///
    /// The bug this pins was invisible to the type system: `run_em_counts`
    /// takes `init_alphas: Option<&[f64]>` immediately before
    /// `eff_lens: Option<&[f64]>`, so passing the effective lengths one slot
    /// early compiled, silently left the prior flat, AND warm-started every
    /// replicate from the effective-length vector. The reported posterior
    /// spread was then drawn under a different prior than the estimate it was
    /// supposed to quantify.
    ///
    /// The transcripts differ in effective length, which is the only thing that
    /// makes a per-nucleotide prior (`prior * effLen`) distinguishable from a
    /// per-transcript one (`prior`), and the class is contested so the prior
    /// actually moves the split.
    #[test]
    fn bootstrap_replicates_follow_the_per_nucleotide_prior() {
        let eff_lens = [1000.0f64, 50.0];
        let b = EquivalenceClassBuilder::new();
        b.add_group(TranscriptGroup::new(vec![0, 1]), vec![1.0, 1.0], 800);
        b.add_group(TranscriptGroup::new(vec![0]), vec![1.0], 200);
        let mut eq = b.finish();
        eq.update_eff_lengths(&eff_lens);
        let packed = PackedEqClasses::from_collapsed(&eq, 2);

        let mut opts = EmOptions {
            use_vbem: true,
            vb_prior: 10.0,
            ..EmOptions::default()
        };

        opts.per_nucleotide_prior = false;
        let flat = bootstrap(&packed, &opts, Some(&eff_lens), 16, 0xC0FFEE);
        opts.per_nucleotide_prior = true;
        let per_nt = bootstrap(&packed, &opts, Some(&eff_lens), 16, 0xC0FFEE);

        assert_eq!(flat.len(), 16);
        assert_eq!(per_nt.len(), 16);
        let mean = |r: &[Vec<f64>], t: usize| -> f64 {
            r.iter().map(|v| v[t]).sum::<f64>() / r.len() as f64
        };
        let (a, b) = (mean(&flat, 0), mean(&per_nt, 0));
        // Same seed, same data, same resampling — only the prior differs, so
        // any difference at all is the prior reaching the replicates.
        assert!(
            (a - b).abs() > 1e-6,
            "--perNucleotidePrior must change the replicates: flat mean {a}, \
             per-nucleotide mean {b}"
        );
    }
}
