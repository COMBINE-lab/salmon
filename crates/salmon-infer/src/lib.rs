//! Collapsed EM / VBEM abundance estimation over equivalence classes.
//!
//! Ports salmon's `CollapsedEMOptimizer` (`src/inference/CollapsedEMOptimizer.cpp`):
//! given a finalized set of equivalence classes (each a transcript label, a
//! count, and per-transcript `combined_weights`), iteratively estimate the
//! expected number of fragments originating from each transcript.
//!
//! The update rules match the C++ exactly:
//! - **EM**: `alphaOut[t] += count * (alphaIn[t] * w_t) / sum_j(alphaIn[j] * w_j)`,
//!   with single-transcript classes assigned their full count.
//! - **VBEM**: replaces `alphaIn[t]` with `expTheta[t] = exp(digamma(alphaIn[t] +
//!   prior_t) - digamma(sum_j(alphaIn[j] + prior_j)))`.
//!
//! The M-steps run in parallel with rayon. Optional SQUAREM acceleration
//! ([`EmAccel::Squarem`]) extrapolates over the plain iteration to reach the same
//! fixpoint in far fewer M-steps.

use salmon_eqclass::CollapsedEqClasses;

mod online;
mod packed;
pub mod uncertainty;

pub use online::OnlineInference;
pub use packed::PackedEqClasses;
pub use uncertainty::{ambiguity_counts, bootstrap, gibbs_sample, GibbsOptions};

/// EM/VBEM acceleration scheme applied on top of the fixed-point map.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EmAccel {
    /// Plain fixed-point iteration (one M-step per iteration). Default; the output
    /// is unchanged from historical salmon.
    #[default]
    None,
    /// SQUAREM (SqS3): extrapolate over two M-steps with a stabilizing third,
    /// reaching the same fixpoint in far fewer M-steps. Not byte-identical to
    /// `None` (a different iterate sequence, same fixpoint within `rel_diff_tol`).
    Squarem,
}

/// Optimizer configuration. Defaults mirror salmon's command-line defaults.
#[derive(Debug, Clone)]
pub struct EmOptions {
    pub max_iter: u32,
    pub min_iter: u32,
    /// relative-difference convergence tolerance
    pub rel_diff_tol: f64,
    /// only transcripts with `alpha` above this participate in the convergence check
    pub alpha_check_cutoff: f64,
    /// abundances below this are truncated to zero on output
    pub min_alpha: f64,
    /// use Variational Bayes EM instead of plain EM
    pub use_vbem: bool,
    /// per-transcript Dirichlet prior (VBEM only)
    pub vb_prior: f64,
    /// interpret `vb_prior` as a per-nucleotide prior (`vb_prior * effLen`)
    /// instead of a flat per-transcript prior (salmon's `--perNucleotidePrior`;
    /// salmon's default is the per-transcript interpretation).
    pub per_nucleotide_prior: bool,
    /// convergence acceleration scheme (default [`EmAccel::None`]).
    pub accel: EmAccel,
}

impl Default for EmOptions {
    fn default() -> Self {
        Self {
            max_iter: 10_000,
            min_iter: 50,
            rel_diff_tol: 0.01,
            alpha_check_cutoff: 1e-2,
            min_alpha: 1e-8,
            use_vbem: false,
            vb_prior: 1e-2,
            per_nucleotide_prior: false,
            accel: EmAccel::None,
        }
    }
}

/// Result of an optimization run.
#[derive(Debug, Clone)]
pub struct EmResult {
    /// estimated fragment counts per transcript (indexed by transcript id)
    pub alphas: Vec<f64>,
    /// iterations actually run
    pub iters: u32,
    /// whether the relative-difference criterion was met before `max_iter`
    pub converged: bool,
    /// total count in equivalence classes whose every transcript was truncated
    /// below `min_alpha` in the final redistribution step — mass that could not be
    /// reassigned to any surviving transcript (reported, not rescaled away). 0 in
    /// the normal case.
    pub dropped_mass: f64,
}

/// Relative-difference convergence check, matching salmon: the max over
/// transcripts (with `alpha_in` above the cutoff) of
/// `|alpha_out - alpha_in| / alpha_out`.
fn max_rel_diff(alpha_in: &[f64], alpha_out: &[f64], cutoff: f64) -> f64 {
    let mut max_d = f64::NEG_INFINITY;
    for i in 0..alpha_in.len() {
        if alpha_in[i] > cutoff && alpha_out[i] > 0.0 {
            let d = (alpha_out[i] - alpha_in[i]).abs() / alpha_out[i];
            if d > max_d {
                max_d = d;
            }
        }
    }
    max_d
}

/// Run the optimizer to convergence (parallel EM/VBEM over the packed layout).
///
/// `eq` must already have `combined_weights` populated (call
/// [`CollapsedEqClasses::update_eff_lengths`](salmon_eqclass::CollapsedEqClasses::update_eff_lengths)).
/// `num_txps` is the total transcript count (output length). Abundances are
/// initialized uniformly over the total fragment count. Internally builds a
/// flat CSR [`PackedEqClasses`] and uses rayon-parallel M-steps.
pub fn optimize(
    eq: &CollapsedEqClasses,
    num_txps: usize,
    opts: &EmOptions,
    eff_lens: Option<&[f64]>,
) -> EmResult {
    let packed = PackedEqClasses::from_collapsed(eq, num_txps);
    optimize_packed_with_init(&packed, opts, true, None, eff_lens)
}

/// As [`optimize`], but warm-starts the abundances from `init_alphas` (per
/// transcript id) when its length matches `num_txps` — used to seed the EM with
/// salmon's count-blended initialization (online estimates blended with uniform),
/// which reduces the iteration count to convergence.
pub fn optimize_with_init(
    eq: &CollapsedEqClasses,
    num_txps: usize,
    opts: &EmOptions,
    init_alphas: Option<&[f64]>,
    eff_lens: Option<&[f64]>,
) -> EmResult {
    let packed = PackedEqClasses::from_collapsed(eq, num_txps);
    optimize_packed_with_init(&packed, opts, true, init_alphas, eff_lens)
}

/// Core convergence loop over a [`PackedEqClasses`]. `parallel` selects the
/// rayon M-step (for the single main run) vs. the sequential one (used by
/// bootstrap, which parallelizes across replicates instead). The per-class
/// `counts` are the packed structure's own (bootstrap passes resampled counts
/// through [`run_em_counts`]).
pub fn optimize_packed(p: &PackedEqClasses, opts: &EmOptions, parallel: bool) -> EmResult {
    optimize_packed_with_init(p, opts, parallel, None, None)
}

/// As [`optimize_packed`], but seeds the abundances from `init_alphas` (a warm
/// start, e.g. salmon's online-estimate-blended-with-uniform initialization)
/// when supplied; otherwise starts uniform.
pub fn optimize_packed_with_init(
    p: &PackedEqClasses,
    opts: &EmOptions,
    parallel: bool,
    init_alphas: Option<&[f64]>,
    eff_lens: Option<&[f64]>,
) -> EmResult {
    let (alphas, iters, converged) = run_em_counts(
        p,
        &p.counts,
        opts,
        parallel,
        opts.min_iter,
        init_alphas,
        eff_lens,
    );
    // Truncate negligible abundances (matches salmon's cutoff), but rather than
    // zero-and-rescale, redistribute the truncated mass to eq-class co-members via
    // one masked final M-step (see `redistribute_truncated`). `min_alpha <= 0`
    // (e.g. the bias warm-up) keeps the continuous vector untouched.
    let (alphas, dropped_mass) =
        finalize_truncate_redistribute(p, &p.counts, alphas, opts, eff_lens);
    EmResult {
        alphas,
        iters,
        converged,
        dropped_mass,
    }
}

/// VBEM Dirichlet prior per transcript: flat `vb_prior`, or `vb_prior·max(1,effLen)`
/// under `--perNucleotidePrior`.
pub(crate) fn prior_alphas_vec(
    opts: &EmOptions,
    eff_lens: Option<&[f64]>,
    num_txps: usize,
) -> Vec<f64> {
    match (opts.per_nucleotide_prior, eff_lens) {
        (true, Some(el)) if el.len() == num_txps => {
            el.iter().map(|&l| opts.vb_prior * l.max(1.0)).collect()
        }
        _ => vec![opts.vb_prior; num_txps],
    }
}

/// Apply the post-convergence min-alpha truncation as a mass-preserving
/// redistribution (not a rescale). Returns the finalized alphas and the mass that
/// could not be redistributed (fully-truncated classes). A non-positive
/// `min_alpha` is a no-op (used by the bias warm-up, whose alphas continue into
/// the warm-started EM).
pub(crate) fn finalize_truncate_redistribute(
    p: &PackedEqClasses,
    counts: &[u64],
    alphas: Vec<f64>,
    opts: &EmOptions,
    eff_lens: Option<&[f64]>,
) -> (Vec<f64>, f64) {
    if opts.min_alpha <= 0.0 {
        return (alphas, 0.0);
    }
    let prior = prior_alphas_vec(opts, eff_lens, p.num_txps);
    packed::redistribute_truncated(p, counts, &alphas, &prior, opts.min_alpha, opts.use_vbem)
}

/// Run EM/VBEM to convergence on `p` with explicit per-class `counts`, returning
/// `(alphas, iters, converged)` *without* the final min-alpha truncation (so
/// bootstrap can apply its own scaling first). `min_iter` is the minimum number
/// of iterations before the convergence check engages.
pub(crate) fn run_em_counts(
    p: &PackedEqClasses,
    counts: &[u64],
    opts: &EmOptions,
    parallel: bool,
    min_iter: u32,
    init_alphas: Option<&[f64]>,
    eff_lens: Option<&[f64]>,
) -> (Vec<f64>, u32, bool) {
    let num_txps = p.num_txps;
    let total: u64 = counts.iter().sum();
    let init = if num_txps > 0 {
        total as f64 / num_txps as f64
    } else {
        0.0
    };
    // Warm start from a supplied initialization (e.g. the online-phase abundance
    // estimates blended with uniform, matching salmon's count-blended init) when
    // its length matches; otherwise start uniform over the total fragment count.
    let mut alphas = match init_alphas {
        Some(a) if a.len() == num_txps => a.to_vec(),
        _ => vec![init; num_txps],
    };
    let mut alphas_prime = vec![0.0f64; num_txps];
    // VBEM prior: flat per-transcript `vb_prior` (salmon's default), or — under
    // `--perNucleotidePrior` — `vb_prior * effLen` per transcript.
    let prior_alphas = prior_alphas_vec(opts, eff_lens, num_txps);
    let mut exp_theta = vec![0.0f64; num_txps];
    let mut scratch: Vec<f64> = Vec::with_capacity(64);
    // Per-shard dense accumulators reused across all parallel M-steps (allocated
    // once here, not per-task per-iteration). Each shard processes a contiguous
    // slice of the classes with plain adds, then they are summed into `alpha_out`
    // — avoiding the cross-thread CAS contention of a single shared atomic array.
    // Capped at 64 shards: beyond that, the per-iteration zero/reduce overhead
    // outweighs the extra accumulation parallelism.
    let mut shards: Vec<Vec<f64>> = if parallel {
        let nshards = rayon::current_num_threads().clamp(1, 64);
        vec![vec![0.0f64; num_txps]; nshards]
    } else {
        Vec::new()
    };

    // One fixed-point M-step `dst = F(src)`, dispatched by (VBEM?, parallel?) and
    // reusing the once-allocated shard/scratch/exp_theta buffers. Both the plain
    // and SQUAREM loops drive `F` through this closure, so the accelerated path
    // reuses the exact same (already tuned) kernels.
    let mut f = |src: &[f64], dst: &mut [f64]| match (opts.use_vbem, parallel) {
        (false, true) => packed::em_step_par(p, counts, src, dst, &mut shards),
        (false, false) => packed::em_step_seq(p, counts, src, dst, &mut scratch),
        (true, true) => packed::vbem_step_par(
            p,
            counts,
            &prior_alphas,
            src,
            dst,
            &mut exp_theta,
            &mut shards,
        ),
        (true, false) => packed::vbem_step_seq(
            p,
            counts,
            &prior_alphas,
            src,
            dst,
            &mut exp_theta,
            &mut scratch,
        ),
    };

    let mut converged = false;
    let mut it = 0u32;
    match opts.accel {
        EmAccel::None => {
            while it < opts.max_iter {
                f(&alphas, &mut alphas_prime);
                it += 1;
                if it >= min_iter {
                    let d = max_rel_diff(&alphas, &alphas_prime, opts.alpha_check_cutoff);
                    std::mem::swap(&mut alphas, &mut alphas_prime);
                    if d.is_finite() && d < opts.rel_diff_tol {
                        converged = true;
                        break;
                    }
                } else {
                    std::mem::swap(&mut alphas, &mut alphas_prime);
                }
            }
        }
        EmAccel::Squarem => {
            converged = squarem_loop(
                &mut f,
                &mut alphas,
                &mut alphas_prime,
                num_txps,
                opts,
                min_iter,
                &mut it,
            );
        }
    }
    (alphas, it, converged)
}

/// SQUAREM (SqS3) acceleration of the fixed-point iteration.
///
/// Each cycle takes two M-steps (`x1 = F(x0)`, `x2 = F(x1)`), forms `r = x1 - x0`
/// and `v = (x2 - x1) - r`, and extrapolates
/// `xp = x0 - 2·alpha·r + alpha²·v` with steplength `alpha = -‖r‖/‖v‖` clamped to
/// `alpha ≤ -1` (so it is never shorter than the plain double-step, for which
/// `alpha = -1` gives exactly `xp = x2`). If any component of `xp` goes negative
/// the steplength is halved back toward `-1` (salmon's non-negativity fallback;
/// worst case `xp = x2`). A stabilizing third M-step `x_next = F(xp)` guarantees a
/// valid EM iterate (non-negative, mass-conserving), and convergence is checked on
/// the full cycle `x0 → x_next` with the same [`max_rel_diff`] criterion as the
/// plain loop.
///
/// `x0` aliases the caller's `alphas` (updated in place to the final estimate) and
/// `x_next` reuses the caller's `alphas_prime` scratch; the vector/norm math is
/// done sequentially (over `num_txps`, cheap next to a per-class M-step) so the
/// result is deterministic regardless of thread count. `it` accumulates the total
/// M-step count (comparable to the plain loop's iteration count).
#[allow(clippy::too_many_arguments)]
fn squarem_loop(
    f: &mut impl FnMut(&[f64], &mut [f64]),
    x0: &mut Vec<f64>,
    x_next: &mut Vec<f64>,
    num_txps: usize,
    opts: &EmOptions,
    min_iter: u32,
    it: &mut u32,
) -> bool {
    let n = num_txps;
    let mut x1 = vec![0.0f64; n];
    let mut x2 = vec![0.0f64; n];
    let mut r = vec![0.0f64; n];
    let mut v = vec![0.0f64; n];
    let mut xp = vec![0.0f64; n];

    while *it < opts.max_iter {
        f(x0, &mut x1);
        *it += 1;
        f(&x1, &mut x2);
        *it += 1;

        // r = x1 - x0 ; v = (x2 - x1) - r ; and their squared norms.
        let mut r_sq = 0.0f64;
        let mut v_sq = 0.0f64;
        for i in 0..n {
            let ri = x1[i] - x0[i];
            let vi = (x2[i] - x1[i]) - ri;
            r[i] = ri;
            v[i] = vi;
            r_sq += ri * ri;
            v_sq += vi * vi;
        }

        // Vanishing (or non-finite) second differences ⇒ at the fixpoint: take the
        // plain step x2. `v_sq <= 0.0` is false for NaN, so the finiteness check
        // also guards NaN/inf.
        if v_sq <= 0.0 || !v_sq.is_finite() {
            std::mem::swap(x0, &mut x2);
            return true;
        }

        // SqS3 steplength, clamped so |alpha| ≥ 1 (alpha = -1 ⇒ xp = x2).
        let mut alpha = -(r_sq.sqrt() / v_sq.sqrt());
        if alpha > -1.0 {
            alpha = -1.0;
        }

        // Extrapolate, backing off toward the double-step if a component < 0.
        let mut tries = 0u32;
        loop {
            let a2 = alpha * alpha;
            let mut ok = true;
            for i in 0..n {
                let val = x0[i] - 2.0 * alpha * r[i] + a2 * v[i];
                xp[i] = val;
                if val < 0.0 {
                    ok = false;
                }
            }
            if ok {
                break;
            }
            if alpha >= -1.0 || tries >= 30 {
                // Give up accelerating this cycle: fall back to the plain x2.
                xp.copy_from_slice(&x2);
                break;
            }
            alpha = (alpha - 1.0) / 2.0; // move toward -1
            tries += 1;
        }

        // Stabilizing M-step, then converge-check on the whole cycle.
        f(&xp, x_next);
        *it += 1;
        if *it >= min_iter {
            let d = max_rel_diff(x0, x_next, opts.alpha_check_cutoff);
            std::mem::swap(x0, x_next);
            if d.is_finite() && d < opts.rel_diff_tol {
                return true;
            }
        } else {
            std::mem::swap(x0, x_next);
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use salmon_eqclass::{EquivalenceClassBuilder, TranscriptGroup};

    /// Build a collapsed eq-class set and set unit effective lengths.
    fn build(classes: &[(Vec<u32>, u64)], num_txps: usize) -> CollapsedEqClasses {
        let b = EquivalenceClassBuilder::new();
        for (txps, count) in classes {
            let w = vec![1.0; txps.len()];
            b.add_group(TranscriptGroup::new(txps.clone()), w, *count);
        }
        let mut eq = b.finish();
        eq.update_eff_lengths(&vec![1.0; num_txps]);
        eq
    }

    #[test]
    fn unique_classes_recover_exact_counts() {
        // Two transcripts, only unique evidence: EM must return those counts.
        let eq = build(&[(vec![0], 30), (vec![1], 70)], 2);
        let res = optimize(&eq, 2, &EmOptions::default(), None);
        assert!((res.alphas[0] - 30.0).abs() < 1e-6);
        assert!((res.alphas[1] - 70.0).abs() < 1e-6);
    }

    #[test]
    fn shared_class_splits_by_unique_evidence() {
        // 10 unique to t0, 90 unique to t1, 100 shared between them.
        // The EM fixpoint allocates the shared class proportionally to the
        // current abundances; with equal eff lengths the stable split tracks
        // the unique ratio, so totals converge to 0.1*200=20 and 0.9*200=180.
        let eq = build(&[(vec![0], 10), (vec![1], 90), (vec![0, 1], 100)], 2);
        let res = optimize(&eq, 2, &EmOptions::default(), None);
        let total = res.alphas[0] + res.alphas[1];
        assert!((total - 200.0).abs() < 1e-6, "total={total}");
        assert!((res.alphas[0] - 20.0).abs() < 1e-2, "a0={}", res.alphas[0]);
        assert!((res.alphas[1] - 180.0).abs() < 1e-2, "a1={}", res.alphas[1]);
    }

    #[test]
    fn conserves_total_count() {
        let eq = build(&[(vec![0, 1, 2], 50), (vec![1, 2], 30), (vec![2], 20)], 3);
        let res = optimize(&eq, 3, &EmOptions::default(), None);
        let total: f64 = res.alphas.iter().sum();
        assert!((total - 100.0).abs() < 1e-6, "total={total}");
    }

    #[test]
    fn vbem_runs_and_conserves_approximately() {
        let eq = build(&[(vec![0], 30), (vec![1], 70), (vec![0, 1], 100)], 2);
        let opts = EmOptions {
            use_vbem: true,
            ..Default::default()
        };
        let res = optimize(&eq, 2, &opts, None);
        let total: f64 = res.alphas.iter().sum();
        // VBEM with a tiny prior stays very close to the EM total.
        assert!((total - 200.0).abs() < 1.0, "total={total}");
        assert!(res.alphas[1] > res.alphas[0]);
    }

    /// A moderately ambiguous eq-class set that takes the plain EM many iterations
    /// to converge — enough that SQUAREM's acceleration is exercised and its
    /// same-fixpoint guarantee is meaningfully tested.
    fn ambiguous() -> CollapsedEqClasses {
        build(
            &[
                (vec![0], 40),
                (vec![1], 35),
                (vec![2], 25),
                (vec![0, 1], 300),
                (vec![1, 2], 260),
                (vec![0, 2], 180),
                (vec![0, 1, 2], 500),
            ],
            3,
        )
    }

    #[test]
    fn squarem_matches_plain_em_fixpoint() {
        let eq = ambiguous();
        let plain = optimize(&eq, 3, &EmOptions::default(), None);
        let sq = optimize(
            &eq,
            3,
            &EmOptions {
                accel: EmAccel::Squarem,
                ..Default::default()
            },
            None,
        );
        for t in 0..3 {
            let rel = (sq.alphas[t] - plain.alphas[t]).abs() / plain.alphas[t].max(1.0);
            assert!(
                rel < 1e-3,
                "txp {t}: squarem {} vs plain {} (rel {rel})",
                sq.alphas[t],
                plain.alphas[t]
            );
        }
        let tp: f64 = plain.alphas.iter().sum();
        let ts: f64 = sq.alphas.iter().sum();
        assert!((tp - ts).abs() < 1e-6, "totals differ: {tp} vs {ts}");
        // Same fixpoint, but reached in no more M-steps than plain (usually far fewer).
        assert!(
            sq.iters <= plain.iters,
            "squarem used more M-steps ({}) than plain ({})",
            sq.iters,
            plain.iters
        );
    }

    #[test]
    fn squarem_matches_plain_vbem_fixpoint() {
        let eq = ambiguous();
        let base = EmOptions {
            use_vbem: true,
            ..Default::default()
        };
        let plain = optimize(&eq, 3, &base, None);
        let sq = optimize(
            &eq,
            3,
            &EmOptions {
                accel: EmAccel::Squarem,
                ..base.clone()
            },
            None,
        );
        for t in 0..3 {
            let rel = (sq.alphas[t] - plain.alphas[t]).abs() / plain.alphas[t].max(1.0);
            assert!(
                rel < 1e-3,
                "vbem txp {t}: {} vs {}",
                sq.alphas[t],
                plain.alphas[t]
            );
        }
    }

    #[test]
    fn squarem_conserves_total_count() {
        let eq = ambiguous();
        let res = optimize(
            &eq,
            3,
            &EmOptions {
                accel: EmAccel::Squarem,
                ..Default::default()
            },
            None,
        );
        let total: f64 = res.alphas.iter().sum();
        // total input count = 40+35+25+300+260+180+500 = 1340
        assert!((total - 1340.0).abs() < 1e-6, "total={total}");
    }

    #[test]
    fn squarem_reduces_iters_on_large_ambiguous_case() {
        // Two nearly-indistinguishable transcripts sharing a huge class, with tiny
        // asymmetric unique evidence: the plain EM crawls from the 50/50 uniform
        // start toward the ~25/75 unique ratio over many iterations. A tight
        // tolerance makes that iteration count large, so SQUAREM's fewer-steps
        // property is clearly exercised.
        let ntxp = 2usize;
        let eq = build(&[(vec![0], 1), (vec![1], 3), (vec![0, 1], 100_000)], ntxp);
        let tight = EmOptions {
            rel_diff_tol: 1e-5,
            min_iter: 1,
            alpha_check_cutoff: 1e-8,
            ..Default::default()
        };
        let plain = optimize(&eq, ntxp, &tight, None);
        let sq = optimize(
            &eq,
            ntxp,
            &EmOptions {
                accel: EmAccel::Squarem,
                ..tight.clone()
            },
            None,
        );
        println!(
            "squarem_reduces_iters: plain M-steps={} (converged={}) squarem M-steps={} (converged={})",
            plain.iters, plain.converged, sq.iters, sq.converged
        );
        // Analytic fixpoint: the shared class splits proportionally to abundance, so
        // f = a0/(a0+a1) solves f·100004 = 1 + 100000·f ⇒ f = 1/4, giving
        // a0 = 25001, a1 = 75003. SQUAREM must land there.
        assert!(sq.converged, "squarem did not converge");
        assert!((sq.alphas[0] - 25001.0).abs() < 1.0, "a0={}", sq.alphas[0]);
        assert!((sq.alphas[1] - 75003.0).abs() < 1.0, "a1={}", sq.alphas[1]);
        // On this deliberately slow mixer the plain iteration has not converged even
        // at the 10 000-step cap, while SQUAREM needs only a handful of M-steps.
        assert!(
            sq.iters * 20 < plain.iters,
            "expected SQUAREM ({}) to need far fewer M-steps than plain ({})",
            sq.iters,
            plain.iters
        );
    }

    #[test]
    fn effective_length_shifts_allocation() {
        // One shared class, equal weights, but t0 is 3x longer -> more of the
        // shared mass should go to the shorter t1.
        let b = EquivalenceClassBuilder::new();
        b.add_group(TranscriptGroup::new(vec![0, 1]), vec![1.0, 1.0], 100);
        let mut eq = b.finish();
        eq.update_eff_lengths(&[300.0, 100.0]);
        let res = optimize(&eq, 2, &EmOptions::default(), None);
        assert!(res.alphas[1] > res.alphas[0], "{:?}", res.alphas);
    }
}
