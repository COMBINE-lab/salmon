//! Collapsed EM / VBEM abundance estimation over equivalence classes.
//!
//! # The estimation problem
//!
//! We know, for each group of identical-looking fragments (an *equivalence
//! class*), how many fragments there were and which transcripts they were
//! compatible with. We do not know which transcript each fragment actually came
//! from. If we knew the abundances we could apportion the ambiguous fragments
//! sensibly; if we knew the apportionment we could compute the abundances. That
//! circularity is exactly what **EM** (expectation-maximization) resolves: guess
//! the abundances, apportion accordingly (E-step), recompute the abundances from
//! that apportionment (M-step), repeat. Each round provably improves the fit, and
//! the sequence converges to a fixed point.
//!
//! Ports salmon's `CollapsedEMOptimizer` (`src/inference/CollapsedEMOptimizer.cpp`):
//! given a finalized set of equivalence classes (each a transcript label, a
//! count, and per-transcript `combined_weights`), iteratively estimate the
//! expected number of fragments originating from each transcript.
//!
//! The update rules match the C++ exactly:
//! - **EM**: `alphaOut[t] += count * (alphaIn[t] * w_t) / sum_j(alphaIn[j] * w_j)`,
//!   with single-transcript classes assigned their full count. Read it as:
//!   split each class's count in proportion to `abundance × compatibility`.
//! - **VBEM**: replaces `alphaIn[t]` with `expTheta[t] = exp(digamma(alphaIn[t] +
//!   prior_t) - digamma(sum_j(alphaIn[j] + prior_j)))` — the variational-Bayes
//!   analogue, which adds a Dirichlet prior and pulls weakly-supported
//!   transcripts toward zero instead of letting them float on noise.
//!
//! The M-steps run in parallel with rayon. Optional SQUAREM acceleration
//! ([`EmAccel::Squarem`]) extrapolates over the plain iteration to reach the same
//! fixpoint in far fewer M-steps.

use salmon_eqclass::CollapsedEqClasses;

mod daarem;
mod online;
mod packed;
pub mod uncertainty;

pub use online::OnlineInference;
pub use packed::PackedEqClasses;
pub use uncertainty::{ambiguity_counts, bootstrap, gibbs_sample, GibbsOptions};

/// EM/VBEM acceleration scheme applied on top of the fixed-point map.
///
/// All three reach the *same* answer; they differ only in how many M-steps it
/// takes to get there, and (for the accelerated ones) in the exact sequence of
/// intermediate iterates. That is why the default is the unaccelerated path:
/// identical output to historical salmon, with acceleration strictly opt-in.
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
    /// DAAREM: damped Anderson acceleration with restarts and residual-monotonicity
    /// control, extrapolating over a window of the last several iterates (a
    /// multi-secant quasi-Newton step) rather than just two. Converges faster than
    /// SQUAREM on high-dimensional, ill-conditioned problems. Same opt-in,
    /// same-fixpoint semantics as [`EmAccel::Squarem`]; not byte-identical to
    /// `None`. See the `daarem` module.
    Daarem,
}

/// Optimizer configuration. Defaults mirror salmon's command-line defaults.
#[derive(Debug, Clone)]
pub struct EmOptions {
    /// hard stop, so a pathological input cannot loop forever
    pub max_iter: u32,
    /// floor on iterations before the convergence check engages, so an early
    /// accidental "nothing moved" cannot end the run prematurely
    pub min_iter: u32,
    /// relative-difference convergence tolerance
    pub rel_diff_tol: f64,
    /// only transcripts with `alpha` above this participate in the convergence check
    ///
    /// Without the cutoff, a transcript oscillating between 1e-30 and 2e-30 would
    /// register a 100% relative change forever and block convergence, despite
    /// being numerically irrelevant.
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
    ///
    /// Fractional, because ambiguous fragments are split between transcripts.
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
///
/// The *maximum* rather than an average: convergence means every meaningful
/// transcript has settled, and an average would let one large still-moving
/// transcript hide behind thousands of settled ones.
pub(crate) fn max_rel_diff(alpha_in: &[f64], alpha_out: &[f64], cutoff: f64) -> f64 {
    let mut max_d = f64::NEG_INFINITY;
    for i in 0..alpha_in.len() {
        // Skip negligible transcripts (see `alpha_check_cutoff`) and guard the
        // division against a zero denominator.
        if alpha_in[i] > cutoff && alpha_out[i] > 0.0 {
            let d = (alpha_out[i] - alpha_in[i]).abs() / alpha_out[i];
            if d > max_d {
                max_d = d;
            }
        }
    }
    // Stays -inf when no transcript qualified; callers test `is_finite()`.
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
///
/// A warm start does not change where the EM lands, only how quickly it gets
/// there: starting near the answer skips the iterations spent travelling.
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
/// through the private `run_em_counts`).
///
/// Parallelizing at only one level is deliberate: nesting a parallel M-step
/// inside parallel replicates would oversubscribe the machine and slow both down.
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
///
/// The per-nucleotide form scales the prior's strength with how much sequence a
/// transcript offers, so a flat prior does not dominate short transcripts.
pub(crate) fn prior_alphas_vec(
    opts: &EmOptions,
    eff_lens: Option<&[f64]>,
    num_txps: usize,
) -> Vec<f64> {
    match (opts.per_nucleotide_prior, eff_lens) {
        // Only usable when effective lengths were supplied and line up; the
        // catch-all arm below falls back to the flat prior otherwise.
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
    // Uniform start: spread all the observed mass evenly. Any strictly positive
    // start converges to the same place; uniform is simply the least presumptuous.
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
    // Destination buffer for each M-step; the two are swapped rather than copied.
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
    //
    // Packaging the M-step as a closure is what lets the three acceleration
    // strategies below be written once against an abstract `F`, with no knowledge
    // of EM vs VBEM or of the threading model.
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
            // The textbook loop: one M-step, check whether anything moved, repeat.
            while it < opts.max_iter {
                f(&alphas, &mut alphas_prime);
                it += 1;
                if it >= min_iter {
                    let d = max_rel_diff(&alphas, &alphas_prime, opts.alpha_check_cutoff);
                    // Swap rather than copy: the buffers just exchange roles, so
                    // no data movement is needed.
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
        EmAccel::Daarem => {
            converged = daarem::daarem_loop(&mut f, &mut alphas, num_txps, opts, min_iter, &mut it);
        }
    }
    (alphas, it, converged)
}

/// SQUAREM (SqS3) acceleration of the fixed-point iteration.
///
/// **The intuition.** Near its fixed point, EM moves along a roughly straight
/// line, each step a little shorter than the last — a geometric series. If you
/// can see the first two steps you can estimate where the whole series is
/// heading, and jump there directly. That is all SQUAREM does: two ordinary
/// steps, one extrapolation, one safety step.
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
/// Here `r` is the first step and `v` how much the second step differed from it —
/// the "acceleration" of the sequence. The clamp and the negativity fallback
/// together guarantee SQUAREM is never *worse* than taking two plain steps.
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
    // Allocated once outside the loop; these are reused every cycle.
    let mut x1 = vec![0.0f64; n];
    let mut x2 = vec![0.0f64; n];
    let mut r = vec![0.0f64; n];
    let mut v = vec![0.0f64; n];
    let mut xp = vec![0.0f64; n];

    while *it < opts.max_iter {
        // Two ordinary EM steps, which is what the extrapolation is built from.
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
        // The ratio of step size to step-size *change* estimates how far the
        // geometric series still has to run.
        let mut alpha = -(r_sq.sqrt() / v_sq.sqrt());
        if alpha > -1.0 {
            alpha = -1.0;
        }

        // Extrapolate, backing off toward the double-step if a component < 0.
        //
        // Abundances are counts and cannot be negative; an over-long jump can
        // overshoot past zero, so the step is repeatedly halved back toward the
        // plain double step until every component is valid.
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
        //
        // The extrapolated point is only a guess and need not be a legal EM
        // iterate; running one real M-step from it re-imposes non-negativity and
        // mass conservation, so the accelerated path can never drift off the
        // manifold the plain path stays on.
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
    ///
    /// Unit lengths and unit weights mean the tests below isolate the EM's
    /// apportionment logic from length correction and bias.
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

    /// With no ambiguity there is nothing to infer, so the EM must return the
    /// observed counts unchanged.
    #[test]
    fn unique_classes_recover_exact_counts() {
        // Two transcripts, only unique evidence: EM must return those counts.
        let eq = build(&[(vec![0], 30), (vec![1], 70)], 2);
        let res = optimize(&eq, 2, &EmOptions::default(), None);
        assert!((res.alphas[0] - 30.0).abs() < 1e-6);
        assert!((res.alphas[1] - 70.0).abs() < 1e-6);
    }

    /// The central behaviour: unique evidence decides how the ambiguous evidence
    /// is split.
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

    /// The EM only ever redistributes evidence, so no fragment may be created or
    /// destroyed however tangled the classes are.
    #[test]
    fn conserves_total_count() {
        let eq = build(&[(vec![0, 1, 2], 50), (vec![1, 2], 30), (vec![2], 20)], 3);
        let res = optimize(&eq, 3, &EmOptions::default(), None);
        let total: f64 = res.alphas.iter().sum();
        assert!((total - 100.0).abs() < 1e-6, "total={total}");
    }

    /// VBEM's prior perturbs the total slightly (by design); check it stays close
    /// and preserves the ordering.
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

    /// Acceleration must change the *speed*, never the *answer* — the single most
    /// important property of an opt-in accelerator.
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

    /// The same guarantee under VBEM, whose map is a different function.
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

    /// Extrapolation could in principle break mass conservation; the stabilizing
    /// M-step is what prevents it, and this pins that.
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

    /// Same fixpoint is necessary but not sufficient — the accelerator also has to
    /// actually accelerate, on a case constructed to be hard for plain EM.
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

    /// The DAAREM path carries the same four obligations as SQUAREM; this is the
    /// first (same fixpoint under plain EM).
    #[test]
    fn daarem_matches_plain_em_fixpoint() {
        let eq = ambiguous();
        let plain = optimize(&eq, 3, &EmOptions::default(), None);
        let da = optimize(
            &eq,
            3,
            &EmOptions {
                accel: EmAccel::Daarem,
                ..Default::default()
            },
            None,
        );
        for t in 0..3 {
            let rel = (da.alphas[t] - plain.alphas[t]).abs() / plain.alphas[t].max(1.0);
            assert!(
                rel < 1e-3,
                "txp {t}: daarem {} vs plain {} (rel {rel})",
                da.alphas[t],
                plain.alphas[t]
            );
        }
        let tp: f64 = plain.alphas.iter().sum();
        let td: f64 = da.alphas.iter().sum();
        assert!((tp - td).abs() < 1e-6, "totals differ: {tp} vs {td}");
        assert!(
            da.iters <= plain.iters,
            "daarem used more M-steps ({}) than plain ({})",
            da.iters,
            plain.iters
        );
    }

    /// Same fixpoint under VBEM.
    #[test]
    fn daarem_matches_plain_vbem_fixpoint() {
        let eq = ambiguous();
        let base = EmOptions {
            use_vbem: true,
            ..Default::default()
        };
        let plain = optimize(&eq, 3, &base, None);
        let da = optimize(
            &eq,
            3,
            &EmOptions {
                accel: EmAccel::Daarem,
                ..base.clone()
            },
            None,
        );
        for t in 0..3 {
            let rel = (da.alphas[t] - plain.alphas[t]).abs() / plain.alphas[t].max(1.0);
            assert!(
                rel < 1e-3,
                "vbem txp {t}: {} vs {}",
                da.alphas[t],
                plain.alphas[t]
            );
        }
    }

    /// Mass conservation through the Anderson extrapolation.
    #[test]
    fn daarem_conserves_total_count() {
        let eq = ambiguous();
        let res = optimize(
            &eq,
            3,
            &EmOptions {
                accel: EmAccel::Daarem,
                ..Default::default()
            },
            None,
        );
        let total: f64 = res.alphas.iter().sum();
        assert!((total - 1340.0).abs() < 1e-6, "total={total}");
    }

    /// And that DAAREM actually accelerates, on the same deliberately slow case.
    #[test]
    fn daarem_reduces_iters_on_slow_mixer() {
        // Same slow-mixing case as the SQUAREM test: two near-indistinguishable
        // transcripts, a huge shared class, tiny asymmetric unique evidence, tight
        // tolerance. Analytic fixpoint a0=25001, a1=75003 (see SQUAREM test).
        let ntxp = 2usize;
        let eq = build(&[(vec![0], 1), (vec![1], 3), (vec![0, 1], 100_000)], ntxp);
        let tight = EmOptions {
            rel_diff_tol: 1e-5,
            min_iter: 1,
            alpha_check_cutoff: 1e-8,
            ..Default::default()
        };
        let plain = optimize(&eq, ntxp, &tight, None);
        let da = optimize(
            &eq,
            ntxp,
            &EmOptions {
                accel: EmAccel::Daarem,
                ..tight.clone()
            },
            None,
        );
        println!(
            "daarem_reduces_iters: plain M-steps={} (converged={}) daarem M-steps={} (converged={})",
            plain.iters, plain.converged, da.iters, da.converged
        );
        assert!(da.converged, "daarem did not converge");
        assert!((da.alphas[0] - 25001.0).abs() < 1.0, "a0={}", da.alphas[0]);
        assert!((da.alphas[1] - 75003.0).abs() < 1.0, "a1={}", da.alphas[1]);
        assert!(
            da.iters * 20 < plain.iters,
            "expected DAAREM ({}) to need far fewer M-steps than plain ({})",
            da.iters,
            plain.iters
        );
    }

    /// Effective length enters through the class weights, so a longer transcript
    /// needs *more* fragments to justify the same abundance and therefore receives
    /// a smaller share of identical evidence.
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
