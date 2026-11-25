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

use crate::packed::{PackedEqClasses, MIN_EQ_CLASS_WEIGHT};
use crate::{run_em_counts, BootstrapMode, EffLens, EmAccel, EmOptions, EmStop, InitAlphas};

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
    let mut out = Vec::new();
    multinomial_into(total, weights, rng, &mut out);
    out
}

/// [`multinomial`] writing into `out` (resized to `weights.len()`), so hot loops
/// such as the Gibbs rounds can reuse one buffer instead of allocating per class.
fn multinomial_into(total: u64, weights: &[f64], rng: &mut impl Rng, out: &mut Vec<u64>) {
    let n = weights.len();
    out.clear();
    out.resize(n, 0);
    if total == 0 || n == 0 {
        return;
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
}

/// Run `num_bootstraps` multinomial bootstrap replicates. Each resamples the
/// per-class counts (multinomial over the original counts, `total_count` draws),
/// optimizes the abundances on the resampled counts as `opts.bootstrap` says
/// ([`BootstrapMode`]), and finalizes them like the point estimate. Each
/// replicate sums to the input total `p.total_count` (the summed
/// equivalence-class counts that were resampled): the finalization
/// redistributes the mass its `min_alpha` truncation removes. We hold to this
/// *intrinsic* total — not an externally-computed mapped-fragment counter — so
/// the replicates stay self-consistent with the resampled input regardless of
/// how (or whether) `num_mapped` is computed elsewhere. Returns one abundance
/// vector per replicate.
///
/// Under the default [`BootstrapMode::RandomStartEm`], every replicate runs
/// plain EM from a uniform-random split of the fragments (Dirichlet(1) over the
/// transcripts), with each component stopping on its own
/// ([`EmStop::Componentwise`]), whichever
/// optimizer produced the point estimate. Under [`BootstrapMode::Legacy`] it
/// runs the point estimate's own optimizer (prior, stopping rule, accelerator)
/// from a uniform start, as salmon ≤ 2.12 did.
///
/// Resampling *classes* rather than individual fragments is what makes this
/// affordable: a class is a group of identical fragments, so one multinomial draw
/// over classes is statistically the same as drawing every fragment.
pub fn bootstrap(
    p: &PackedEqClasses,
    opts: &EmOptions,
    eff_lens: EffLens<'_>,
    num_bootstraps: u32,
    seed: u64,
) -> Vec<Vec<f64>> {
    let sample_weights: Vec<f64> = p.counts.iter().map(|&c| c as f64).collect();
    let total = p.total_count;
    let random_start = opts.bootstrap == BootstrapMode::RandomStartEm;
    // The replicates' optimizer. The random-start mode swaps in plain EM with
    // per-component stopping; the rest (tolerances, caps, truncation) is the
    // point estimate's.
    let rep_opts = if random_start {
        EmOptions {
            use_vbem: false,
            accel: EmAccel::None,
            em_stop: EmStop::Componentwise,
            ..opts.clone()
        }
    } else {
        opts.clone()
    };
    let unit = Gamma::new(1.0, 1.0).expect("Gamma(1, 1)");
    (0..num_bootstraps)
        .into_par_iter()
        .map(|bs| {
            // Seed derived from (run seed, replicate index) so results reproduce
            // exactly and do not depend on which thread ran which replicate. The
            // odd multiplier scatters adjacent indices to distant seeds.
            let mut rng =
                Pcg64Mcg::seed_from_u64(seed ^ (bs as u64).wrapping_mul(0x9E3779B97F4A7C15));
            let resampled = multinomial(total, &sample_weights, &mut rng);
            // The random start: normalized Gamma(1, 1) draws are a Dirichlet(1)
            // split, scaled to the resampled total. It is drawn after the
            // resample, from the same stream, so a replicate is a function of
            // (seed, index) alone. Its transcripts with no class get no mass at
            // the first M-step, whatever they start with.
            let init_vec: Option<Vec<f64>> = random_start.then(|| {
                let w: Vec<f64> = (0..p.num_txps).map(|_| unit.sample(&mut rng)).collect();
                let sw: f64 = w.iter().sum();
                w.iter().map(|x| x / sw * total as f64).collect()
            });
            let init = match &init_vec {
                Some(v) => InitAlphas::new(v),
                None => InitAlphas::NONE,
            };
            // Effective lengths reach the replicate EM so `--perNucleotidePrior`
            // shapes the prior here exactly as it does the point estimate (when
            // the replicates run VBEM, i.e. legacy mode); passing `None` silently
            // fell back to the flat prior, which made the flag a no-op inside
            // every replicate (#1140, audit D13).
            //
            // NOTE the argument order: `init_alphas` precedes `eff_lens` and the
            // two newtypes exist because swapping them used to type-check. The
            // first attempt at this fix did exactly that — the prior stayed
            // flat AND every replicate warm-started from the effective-length
            // vector.
            let (alphas, _, _, _, _) =
                run_em_counts(p, &resampled, &rep_opts, false, 50, init, eff_lens);
            // Finalize like the point estimate: truncate the negligible
            // abundances, then redistribute that mass to eq-class co-members via a
            // masked final M-step over the *resampled* counts (no rescale-up). The
            // replicate's `dropped` mass is negligible and not reported (only the
            // point estimate surfaces `inference_truncated_mass`).
            //
            // Finalizing the same way as the point estimate keeps post-processing
            // from showing up as spread.
            let (alphas, _dropped) =
                crate::finalize_truncate_redistribute(p, &resampled, alphas, &rep_opts, eff_lens);
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

impl GibbsOptions {
    /// The prior that follows the main optimizer (salmon): with VBEM a
    /// per-transcript `max(1.0, vbPrior)`, or `vbPrior` per nucleotide under
    /// `--perNucleotidePrior` (a VBEM option) to match the point estimate; with
    /// plain EM a flat 1e-3 per transcript. Every driver builds its options here,
    /// so the rule cannot drift between them again.
    pub fn for_em(em: &crate::EmOptions, num_samples: u32, thinning: u32) -> Self {
        let per_nucleotide = em.use_vbem && em.per_nucleotide_prior;
        let prior = if per_nucleotide {
            em.vb_prior
        } else if em.use_vbem {
            em.vb_prior.max(1.0)
        } else {
            1e-3
        };
        Self {
            num_samples,
            thinning,
            prior,
            per_transcript_prior: !per_nucleotide,
        }
    }
}

/// salmon's Gibbs rate parameter `beta`.
/// The rate of the Gamma prior on each transcript's expression rate.
const GIBBS_BETA: f64 = 0.1;

/// Fixed work-block sizes for one Gibbs round. Each block draws from its own
/// stream, seeded from (chain, round, block), so the samples depend on these
/// constants and the seed but never on how many threads run the blocks.
const GIBBS_TXP_BLOCK: usize = 4096;
const GIBBS_CLASS_BLOCK: usize = 2048;

/// The RNG for one block of one round of one chain (splitmix64 finalizer, so
/// neighbouring indices land on unrelated seeds).
fn block_rng(chain_seed: u64, round: u64, block: u64) -> Pcg64Mcg {
    let mut z = chain_seed
        ^ round.wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ block.wrapping_mul(0xC2B2_AE3D_27D4_EB4F);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    Pcg64Mcg::seed_from_u64(z ^ (z >> 31))
}

/// Per-chain buffers reused by every round.
struct GibbsScratch {
    /// `mu` for the active transcripts, in `active` order
    mu_active: Vec<f64>,
    /// each class's multinomial draw, aligned with `p.labels`
    draws: Vec<u64>,
}

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
    chain_seed: u64,
    round: u64,
    sc: &mut GibbsScratch,
    serial: bool,
) {
    if serial {
        gibbs_round_serial(
            p,
            active,
            eff_lens,
            prior_alphas,
            txp_count,
            mu,
            chain_seed,
            round,
        );
        return;
    }
    // Sample mu[i] ~ Gamma(txpCount[i] + prior[i], 1/(beta + effLen[i])); reset count.
    //
    // Shape = observed count + prior (more reads ⇒ tighter, higher draw); scale
    // shrinks with effective length, converting counts back to a rate per
    // available position. The counts are zeroed afterwards because the second
    // half of the round refills them from scratch.
    {
        let tc: &[f64] = txp_count;
        sc.mu_active
            .par_chunks_mut(GIBBS_TXP_BLOCK)
            .zip(active.par_chunks(GIBBS_TXP_BLOCK))
            .enumerate()
            .for_each(|(b, (out, ids))| {
                let mut rng = block_rng(chain_seed, round, 2 * b as u64);
                for (m, &i) in out.iter_mut().zip(ids) {
                    let i = i as usize;
                    let ci = tc[i] + prior_alphas[i];
                    let scale = 1.0 / (GIBBS_BETA + eff_lens[i]);
                    *m = if ci > 0.0 {
                        Gamma::new(ci, scale).unwrap().sample(&mut rng)
                    } else {
                        0.0
                    };
                }
            });
    }
    for (&i, &m) in active.iter().zip(&sc.mu_active) {
        mu[i as usize] = m;
        txp_count[i as usize] = 0.0;
    }

    // Resample each class's reads across its transcripts, in fixed blocks of
    // classes; each block writes the label-aligned slice of `draws` it owns.
    let nclasses = p.num_classes();
    let mu_ro: &[f64] = mu;
    let mut rest: &mut [u64] = &mut sc.draws;
    let mut blocks = Vec::with_capacity(nclasses.div_ceil(GIBBS_CLASS_BLOCK));
    for c0 in (0..nclasses).step_by(GIBBS_CLASS_BLOCK) {
        let c1 = (c0 + GIBBS_CLASS_BLOCK).min(nclasses);
        let len = (p.starts[c1] - p.starts[c0]) as usize;
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(len);
        rest = tail;
        blocks.push((c0, c1, head));
    }
    blocks.into_par_iter().for_each(|(c0, c1, out)| {
        let mut rng = block_rng(chain_seed, round, 2 * (c0 / GIBBS_CLASS_BLOCK) as u64 + 1);
        let base = p.starts[c0] as usize;
        let mut probs: Vec<f64> = Vec::with_capacity(64);
        let mut draws: Vec<u64> = Vec::with_capacity(64);
        for ci in c0..c1 {
            let s = p.starts[ci] as usize;
            let e = p.starts[ci + 1] as usize;
            let dst = &mut out[s - base..e - base];
            draw_class(p, ci, mu_ro, &mut rng, &mut probs, &mut draws);
            dst.copy_from_slice(&draws);
        }
    });
    // Integer-valued sums, so this serial fold is exact in any order.
    for (&tid, &k) in p.labels.iter().zip(&sc.draws) {
        txp_count[tid as usize] += k as f64;
    }
}

/// Draw class `ci`'s reads across its transcripts into `draws` (aligned with the
/// class's labels).
fn draw_class(
    p: &PackedEqClasses,
    ci: usize,
    mu: &[f64],
    rng: &mut Pcg64Mcg,
    probs: &mut Vec<f64>,
    draws: &mut Vec<u64>,
) {
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
        multinomial_into(class_count, probs, rng, draws);
    } else {
        // Unambiguous class: nothing random about where its reads go.
        draws.clear();
        draws.extend(tids.iter().map(|_| class_count));
    }
}

/// [`gibbs_round`] when the chains already use every thread: the same blocks
/// with the same streams, run in order, adding straight into `txp_count`
/// instead of going through the `draws` buffer and a separate fold. The
/// samples are identical to the parallel path's.
#[allow(clippy::too_many_arguments)]
fn gibbs_round_serial(
    p: &PackedEqClasses,
    active: &[u32],
    eff_lens: &[f64],
    prior_alphas: &[f64],
    txp_count: &mut [f64],
    mu: &mut [f64],
    chain_seed: u64,
    round: u64,
) {
    for (b, ids) in active.chunks(GIBBS_TXP_BLOCK).enumerate() {
        let mut rng = block_rng(chain_seed, round, 2 * b as u64);
        for &i in ids {
            let i = i as usize;
            let ci = txp_count[i] + prior_alphas[i];
            let scale = 1.0 / (GIBBS_BETA + eff_lens[i]);
            mu[i] = if ci > 0.0 {
                Gamma::new(ci, scale).unwrap().sample(&mut rng)
            } else {
                0.0
            };
            txp_count[i] = 0.0;
        }
    }
    let nclasses = p.num_classes();
    let mut probs: Vec<f64> = Vec::with_capacity(64);
    let mut draws: Vec<u64> = Vec::with_capacity(64);
    for c0 in (0..nclasses).step_by(GIBBS_CLASS_BLOCK) {
        let c1 = (c0 + GIBBS_CLASS_BLOCK).min(nclasses);
        let mut rng = block_rng(chain_seed, round, 2 * (c0 / GIBBS_CLASS_BLOCK) as u64 + 1);
        for ci in c0..c1 {
            let s = p.starts[ci] as usize;
            let e = p.starts[ci + 1] as usize;
            if e - s == 1 {
                // Unambiguous class: no draw, and no detour through `draws`.
                txp_count[p.labels[s] as usize] += p.counts[ci] as f64;
                continue;
            }
            draw_class(p, ci, mu, &mut rng, &mut probs, &mut draws);
            for (&tid, &k) in p.labels[s..e].iter().zip(&draws) {
                txp_count[tid as usize] += k as f64;
            }
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
    // With no thread to spare beyond one per chain, splitting a round into
    // parallel blocks only adds the fold; run the blocks in order instead.
    let serial = rayon::current_num_threads() <= nchains;
    let blocks: Vec<(usize, Vec<Vec<f64>>)> = bounds
        .par_iter()
        .enumerate()
        .map(|(c, &(start, end))| {
            let chain_seed = seed ^ (c as u64).wrapping_mul(0xD1B54A32D192ED03);
            let mut txp_count = init.clone();
            let mut mu = vec![0.0f64; num_txps];
            // Per-chain scratch reused by every round (no per-class allocation).
            let mut sc = if serial {
                GibbsScratch {
                    mu_active: Vec::new(),
                    draws: Vec::new(),
                }
            } else {
                GibbsScratch {
                    mu_active: vec![0.0; active.len()],
                    draws: vec![0; p.labels.len()],
                }
            };
            let mut round = 0u64;
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
                        chain_seed,
                        round,
                        &mut sc,
                        serial,
                    );
                    round += 1;
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
pub fn ambiguity_counts(p: &PackedEqClasses) -> (Vec<u64>, Vec<u64>) {
    let mut unique = vec![0u64; p.num_txps];
    let mut ambig = vec![0u64; p.num_txps];
    for ci in 0..p.num_classes() {
        let s = p.starts[ci] as usize;
        let e = p.starts[ci + 1] as usize;
        let tids = &p.labels[s..e];
        let count = p.counts[ci];
        if tids.len() > 1 {
            for &t in tids {
                ambig[t as usize] += count;
            }
        } else if tids.len() == 1 {
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
        let bs = bootstrap(&p, &EmOptions::default(), EffLens::NONE, 50, 12345);
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
        let bs = bootstrap(&p, &EmOptions::default(), EffLens::NONE, 100, 7);
        let m0: f64 = bs.iter().map(|b| b[0]).sum::<f64>() / 100.0;
        let var0: f64 = bs.iter().map(|b| (b[0] - m0).powi(2)).sum::<f64>() / 100.0;
        assert!(
            var0 > 0.0,
            "ambiguous transcript should have nonzero bootstrap variance"
        );
    }

    /// Two transcripts no read tells apart: every split of their shared reads
    /// fits equally well. Replicates from one fixed start split them the same
    /// way every time, showing no uncertainty at all; random-start EM
    /// replicates land in different places, showing that the data leave the
    /// split open. Either way the resampled pair total is conserved.
    #[test]
    fn random_start_replicates_show_what_the_data_cannot_split() {
        let p = packed(&[(vec![0, 1], 1000), (vec![2], 200)], 3);
        let spread = |mode| {
            let opts = EmOptions {
                use_vbem: true,
                bootstrap: mode,
                ..EmOptions::default()
            };
            let bs = bootstrap(&p, &opts, EffLens::NONE, 40, 0xB007);
            for b in &bs {
                assert!(((b[0] + b[1] + b[2]) - 1200.0).abs() < 1e-6);
            }
            bs.iter().map(|b| (b[0] - b[1]).abs()).sum::<f64>() / bs.len() as f64
        };
        let legacy = spread(crate::BootstrapMode::Legacy);
        let random = spread(crate::BootstrapMode::RandomStartEm);
        assert!(
            legacy < 1e-6,
            "fixed start: the pair is always split evenly ({legacy})"
        );
        assert!(
            random > 100.0,
            "random starts must spread the split ({random})"
        );
    }

    /// Random starts come from each replicate's own stream, so replicates
    /// reproduce exactly and do not depend on the thread count.
    #[test]
    fn random_start_replicates_are_thread_count_independent() {
        let p = packed(
            &[
                (vec![0, 1], 500),
                (vec![1, 2], 300),
                (vec![0], 40),
                (vec![3], 90),
            ],
            4,
        );
        let run = |n| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
                .install(|| bootstrap(&p, &EmOptions::default(), EffLens::NONE, 12, 5))
        };
        assert_eq!(run(1), run(8));
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

    /// Rounds run in fixed blocks with their own streams, so the samples must
    /// not depend on the thread count (several blocks of each kind here).
    #[test]
    fn gibbs_is_thread_count_independent() {
        let ntx = 2 * GIBBS_TXP_BLOCK + 37;
        let mut classes = Vec::new();
        for c in 0..(2 * GIBBS_CLASS_BLOCK + 11) {
            let a = (c * 7) % ntx;
            let b = (c * 13 + 1) % ntx;
            let mut t = vec![a as u32];
            if a != b {
                t.push(b as u32);
                t.sort();
            }
            classes.push((t, 1 + (c as u64 % 50)));
        }
        let p = packed(&classes, ntx);
        let eff = vec![100.0; ntx];
        let init = vec![1.0; ntx];
        let opts = GibbsOptions {
            num_samples: 4,
            thinning: 2,
            ..Default::default()
        };
        let run = |n| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(n)
                .build()
                .unwrap()
                .install(|| gibbs_sample(&p, &eff, &init, &opts, 7))
        };
        assert_eq!(run(1), run(8));
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

        // The legacy replicates run the point estimate's VBEM, with its prior;
        // the default random-start replicates run plain EM, which has none.
        let mut opts = EmOptions {
            use_vbem: true,
            vb_prior: 10.0,
            bootstrap: crate::BootstrapMode::Legacy,
            ..EmOptions::default()
        };

        opts.per_nucleotide_prior = false;
        let flat = bootstrap(&packed, &opts, EffLens::new(&eff_lens), 16, 0xC0FFEE);
        opts.per_nucleotide_prior = true;
        let per_nt = bootstrap(&packed, &opts, EffLens::new(&eff_lens), 16, 0xC0FFEE);

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

    /// Every driver takes its Gibbs prior from `for_em`. The default path kept
    /// `max(1, vbPrior)` per nucleotide after only the `--online` copy was fixed,
    /// a prior 100x the point estimate's at `--vbPrior 0.01`.
    #[test]
    fn gibbs_prior_follows_the_optimizer() {
        let em = |use_vbem, per_nucleotide_prior| crate::EmOptions {
            use_vbem,
            per_nucleotide_prior,
            vb_prior: 0.01,
            ..Default::default()
        };
        let g = |o: &crate::EmOptions| {
            let g = GibbsOptions::for_em(o, 10, 16);
            (g.prior, g.per_transcript_prior)
        };
        assert_eq!(
            g(&em(true, false)),
            (1.0, true),
            "VBEM: max(1, vbPrior) per transcript"
        );
        assert_eq!(
            g(&em(true, true)),
            (0.01, false),
            "VBEM per nucleotide: vbPrior"
        );
        assert_eq!(
            g(&em(false, true)),
            (1e-3, true),
            "plain EM ignores the VBEM-only flag"
        );
        assert_eq!(g(&em(false, false)), (1e-3, true));
    }
}
