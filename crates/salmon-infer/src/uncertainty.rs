//! Posterior uncertainty: multinomial **bootstrap** (`CollapsedEMOptimizer::gatherBootstraps`)
//! and the non-collapsed **Gibbs sampler** (`CollapsedGibbsSampler`).
//!
//! Both operate on the flat [`PackedEqClasses`] layout and produce one abundance
//! vector per replicate/sample. Bootstrap parallelizes across replicates (each
//! runs a sequential EM on resampled class counts); Gibbs parallelizes across
//! independent chains (each runs sequential thinned rounds). RNG is PCG
//! (`rand_pcg`), seeded per replicate/chain for reproducibility — results are
//! statistically equivalent to salmon's but not bit-identical (different RNG).

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
fn multinomial(total: u64, weights: &[f64], rng: &mut impl Rng) -> Vec<u64> {
    let n = weights.len();
    let mut out = vec![0u64; n];
    if total == 0 || n == 0 {
        return out;
    }
    let mut remaining = total;
    let mut remaining_w: f64 = weights.iter().sum();
    for i in 0..n {
        if remaining == 0 {
            break;
        }
        if i == n - 1 || remaining_w <= 0.0 {
            out[i] = remaining;
            break;
        }
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
/// runs EM/VBEM to convergence (min 50 iters), and — when `scale_counts` (salmon's
/// `useScaledCounts`) — rescales the abundances to sum to `num_mapped_frags`.
/// Returns one abundance vector per replicate.
pub fn bootstrap(
    p: &PackedEqClasses,
    opts: &EmOptions,
    num_bootstraps: u32,
    num_mapped_frags: u64,
    scale_counts: bool,
    seed: u64,
) -> Vec<Vec<f64>> {
    let sample_weights: Vec<f64> = p.counts.iter().map(|&c| c as f64).collect();
    let total = p.total_count;
    (0..num_bootstraps)
        .into_par_iter()
        .map(|bs| {
            let mut rng =
                Pcg64Mcg::seed_from_u64(seed ^ (bs as u64).wrapping_mul(0x9E3779B97F4A7C15));
            let resampled = multinomial(total, &sample_weights, &mut rng);
            let (mut alphas, _, _) = run_em_counts(p, &resampled, opts, false, 50, None, None);
            // truncate tiny values
            for a in &mut alphas {
                if *a < opts.min_alpha {
                    *a = 0.0;
                }
            }
            if scale_counts {
                let sum: f64 = alphas.iter().sum();
                if sum > 0.0 {
                    let scale = num_mapped_frags as f64 / sum;
                    for a in &mut alphas {
                        *a *= scale;
                    }
                }
            }
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
const GIBBS_BETA: f64 = 0.1;

/// One Gibbs round (salmon's `sampleRoundNonCollapsedMultithreaded_`): draw the
/// transcript fractions `mu` from their Gamma posterior, then resample each
/// equivalence class's count multinomially across its transcripts.
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
        let weights = &p.weights[s..e];
        if tids.len() > 1 {
            probs.clear();
            let mut denom = 0.0;
            for (&tid, &w) in tids.iter().zip(weights) {
                let v = 1000.0 * mu[tid as usize] * w;
                probs.push(v);
                denom += v;
            }
            if denom <= MIN_EQ_CLASS_WEIGHT {
                // fall back to uniform over the class
                for v in probs.iter_mut() {
                    *v = 1.0;
                }
            }
            let draws = multinomial(class_count, &probs, rng);
            for (&tid, &k) in tids.iter().zip(&draws) {
                txp_count[tid as usize] += k as f64;
            }
        } else {
            txp_count[tids[0] as usize] += class_count as f64;
        }
    }
}

/// Draw `opts.num_samples` Gibbs posterior samples. `init_alphas` is the point
/// estimate (EM result) each chain restarts from; `eff_lens` the effective
/// lengths. Chains (1/2/4/8 by sample count, like salmon) run in parallel.
/// Returns one abundance vector per sample (scaled to `num_mapped_frags`).
pub fn gibbs_sample(
    p: &PackedEqClasses,
    eff_lens: &[f64],
    init_alphas: &[f64],
    opts: &GibbsOptions,
    num_mapped_frags: u64,
    seed: u64,
) -> Vec<Vec<f64>> {
    let num_txps = p.num_txps;
    let num_samples = opts.num_samples as usize;
    if num_samples == 0 {
        return Vec::new();
    }

    // Active transcripts = those appearing in some class.
    let mut active_flag = vec![false; num_txps];
    for &t in &p.labels {
        active_flag[t as usize] = true;
    }
    let active: Vec<u32> = (0..num_txps as u32)
        .filter(|&t| active_flag[t as usize])
        .collect();

    // Per-transcript prior (per-txp = constant; per-nucleotide = prior·max(1,effLen)).
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
                // Extrapolate scaled counts from the final fractions mu.
                let denom: f64 = (0..num_txps).map(|t| mu[t] * eff_lens[t]).sum();
                let scale = if denom > 0.0 {
                    num_mapped_frags as f64 / denom
                } else {
                    0.0
                };
                let mut sample = vec![0.0f64; num_txps];
                for t in 0..num_txps {
                    let a = mu[t] * eff_lens[t] * scale;
                    sample[t] = if a > 1e-8 { a } else { 0.0 };
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

    #[test]
    fn bootstrap_mean_near_point_estimate() {
        // unique evidence -> every bootstrap recovers ~the same counts
        let p = packed(&[(vec![0], 300), (vec![1], 700)], 2);
        let bs = bootstrap(&p, &EmOptions::default(), 50, 1000, true, 12345);
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

    #[test]
    fn bootstrap_variance_grows_with_ambiguity() {
        // a fully shared class has higher per-transcript bootstrap variance
        let p = packed(&[(vec![0], 10), (vec![1], 10), (vec![0, 1], 980)], 2);
        let bs = bootstrap(&p, &EmOptions::default(), 100, 1000, true, 7);
        let m0: f64 = bs.iter().map(|b| b[0]).sum::<f64>() / 100.0;
        let var0: f64 = bs.iter().map(|b| (b[0] - m0).powi(2)).sum::<f64>() / 100.0;
        assert!(
            var0 > 0.0,
            "ambiguous transcript should have nonzero bootstrap variance"
        );
    }

    #[test]
    fn gibbs_runs_and_conserves_scale() {
        let p = packed(&[(vec![0], 300), (vec![1], 700)], 2);
        let opts = GibbsOptions {
            num_samples: 20,
            thinning: 8,
            ..Default::default()
        };
        let samples = gibbs_sample(&p, &[1.0, 1.0], &[300.0, 700.0], &opts, 1000, 99);
        assert_eq!(samples.len(), 20);
        for s in &samples {
            let tot = s[0] + s[1];
            assert!(
                (tot - 1000.0).abs() < 50.0,
                "gibbs total {tot} not near 1000"
            );
        }
    }

    #[test]
    fn ambiguity_counts_split() {
        let p = packed(&[(vec![0], 30), (vec![1], 70), (vec![0, 1], 100)], 2);
        let (uniq, amb) = ambiguity_counts(&p);
        assert_eq!(uniq, vec![30, 70]);
        assert_eq!(amb, vec![100, 100]);
    }
}
