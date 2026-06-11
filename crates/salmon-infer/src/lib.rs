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
//! Parallelization with rayon and SQUAREM acceleration are deferred; plain
//! iteration converges to the same fixpoint.

use salmon_eqclass::CollapsedEqClasses;

mod online;
mod packed;
pub mod uncertainty;

pub use online::OnlineInference;
pub use packed::PackedEqClasses;
pub use uncertainty::{ambiguity_counts, bootstrap, gibbs_sample, GibbsOptions};

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
pub fn optimize(eq: &CollapsedEqClasses, num_txps: usize, opts: &EmOptions) -> EmResult {
    let packed = PackedEqClasses::from_collapsed(eq, num_txps);
    optimize_packed(&packed, opts, true)
}

/// Core convergence loop over a [`PackedEqClasses`]. `parallel` selects the
/// rayon M-step (for the single main run) vs. the sequential one (used by
/// bootstrap, which parallelizes across replicates instead). The per-class
/// `counts` are the packed structure's own (bootstrap passes resampled counts
/// through [`run_em_counts`]).
pub fn optimize_packed(p: &PackedEqClasses, opts: &EmOptions, parallel: bool) -> EmResult {
    let (mut alphas, iters, converged) = run_em_counts(p, &p.counts, opts, parallel, opts.min_iter);
    // truncate negligible abundances (matches salmon's cutoff)
    for a in &mut alphas {
        if *a < opts.min_alpha {
            *a = 0.0;
        }
    }
    EmResult {
        alphas,
        iters,
        converged,
    }
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
) -> (Vec<f64>, u32, bool) {
    let num_txps = p.num_txps;
    let total: u64 = counts.iter().sum();
    let init = if num_txps > 0 { total as f64 / num_txps as f64 } else { 0.0 };
    let mut alphas = vec![init; num_txps];
    let mut alphas_prime = vec![0.0f64; num_txps];
    let prior_alphas = vec![opts.vb_prior; num_txps];
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

    let mut converged = false;
    let mut it = 0u32;
    while it < opts.max_iter {
        match (opts.use_vbem, parallel) {
            (false, true) => packed::em_step_par(p, counts, &alphas, &mut alphas_prime, &mut shards),
            (false, false) => {
                packed::em_step_seq(p, counts, &alphas, &mut alphas_prime, &mut scratch)
            }
            (true, true) => packed::vbem_step_par(
                p, counts, &prior_alphas, &alphas, &mut alphas_prime, &mut exp_theta, &mut shards,
            ),
            (true, false) => packed::vbem_step_seq(
                p, counts, &prior_alphas, &alphas, &mut alphas_prime, &mut exp_theta,
                &mut scratch,
            ),
        }
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
    (alphas, it, converged)
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
        let res = optimize(&eq, 2, &EmOptions::default());
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
        let res = optimize(&eq, 2, &EmOptions::default());
        let total = res.alphas[0] + res.alphas[1];
        assert!((total - 200.0).abs() < 1e-6, "total={total}");
        assert!((res.alphas[0] - 20.0).abs() < 1e-2, "a0={}", res.alphas[0]);
        assert!((res.alphas[1] - 180.0).abs() < 1e-2, "a1={}", res.alphas[1]);
    }

    #[test]
    fn conserves_total_count() {
        let eq = build(&[(vec![0, 1, 2], 50), (vec![1, 2], 30), (vec![2], 20)], 3);
        let res = optimize(&eq, 3, &EmOptions::default());
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
        let res = optimize(&eq, 2, &opts);
        let total: f64 = res.alphas.iter().sum();
        // VBEM with a tiny prior stays very close to the EM total.
        assert!((total - 200.0).abs() < 1.0, "total={total}");
        assert!(res.alphas[1] > res.alphas[0]);
    }

    #[test]
    fn effective_length_shifts_allocation() {
        // One shared class, equal weights, but t0 is 3x longer -> more of the
        // shared mass should go to the shorter t1.
        let b = EquivalenceClassBuilder::new();
        b.add_group(TranscriptGroup::new(vec![0, 1]), vec![1.0, 1.0], 100);
        let mut eq = b.finish();
        eq.update_eff_lengths(&[300.0, 100.0]);
        let res = optimize(&eq, 2, &EmOptions::default());
        assert!(res.alphas[1] > res.alphas[0], "{:?}", res.alphas);
    }
}
