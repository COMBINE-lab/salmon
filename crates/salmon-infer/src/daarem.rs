//! DAAREM acceleration of the EM/VBEM fixed-point iteration.
//!
//! # The problem being solved
//!
//! EM is an iterative algorithm: from a guess `x` it computes a better guess
//! `F(x)`, and repeats until the guess stops changing. Mathematically it is
//! looking for a *fixed point*, a vector with `F(x) = x`. EM is reliable but
//! slow: near the answer each iteration only shaves a fixed fraction off the
//! remaining error, so thousands of iterations may be needed.
//!
//! *Acceleration* means using the history of recent iterates to guess where the
//! sequence is heading and jump ahead, rather than plodding one step at a time.
//! Each accelerated step costs one extra evaluation of `F` but can be worth
//! dozens of plain ones.
//!
//! # What DAAREM is
//!
//! Port of DAAREM (Damped Anderson Acceleration with Restarts and
//! Epsilon-Monotonicity) from Henderson & Varadhan, *J. Comput. Graph. Statist.*
//! 28(4):834-846 (2019), following the CRAN `daarem` package's
//! `daarem_base_objfn`: a proposed Anderson step is accepted unless it lowers
//! the objective by more than `MON_TOL`, otherwise the iteration falls back to a
//! plain fixed-point (EM) step.
//!
//! The objective is the one the underlying fixed-point map climbs: the
//! log-likelihood for EM, and for VBEM the evidence lower bound with the
//! assignment distribution optimized for the current Dirichlet. Neither costs a
//! separate pass — each M-step already forms every class's denominator
//! `Σ_t x_t·w_ct`, and `Σ_c n_c·ln(denominator)` is the data term; the rest is
//! `O(num_txps)` (see `run_em_counts`).
//!
//! **Why not the objective-free variant.** CRAN also ships
//! `daarem_base_noobjfn`, which judges steps by the residual norm `‖F(x) − x‖`
//! instead, and this module originally ported that. The package documents the
//! objective as "preferable", and VBEM shows why: its sparse Dirichlet prior
//! drives near-identical isoforms to winner-take-all solutions, and while a
//! losing isoform drains, the residual *grows* (the collapse accelerates). A
//! residual-monotone test therefore refuses VBEM's own path to the optimum and
//! accepts only extrapolations that jump back into the slow region. On
//! synthetic and real (GENCODE) problems that was a deterministic limit cycle —
//! one run repeated the same state every 9750 M-steps until `max_iter`; plain
//! VBEM started from the stall converged in 66. Plain EM and VBEM steps never
//! lower their objective, so an objective test cannot reject that path.
//!
//! **Local optima.** The VBEM objective under a prior below 1 has many optima
//! (one per choice of winners among near-identical isoforms), and an Anderson
//! step can cross into a different basin than plain VBEM would reach. Measured
//! over perturbed starts on three GENCODE datasets, DAAREM's optima were
//! higher on two and lower on one, and the most reproducible of plain VBEM,
//! SQUAREM and DAAREM on all three. Plain EM's objective is concave, so there
//! the question does not arise.
//!
//! Where SQUAREM extrapolates from only the last two iterates, DAAREM solves a
//! damped least-squares problem over a window of the last `nlag` residual/iterate
//! differences (a multi-secant quasi-Newton step), which converges faster on
//! high-dimensional, ill-conditioned problems — the regime of transcript-level
//! quantification. The damping (a ridge on the least-squares solve) starts strong
//! (near a plain EM step) and anneals toward the full Anderson step as accepted
//! steps accumulate; the history is periodically restarted.
//!
//! "Damping" here means blending the aggressive extrapolation back toward a
//! plain EM step. Early on the history is short and unreliable, so the step is
//! heavily damped; as accepted steps accumulate the damping relaxes.
//!
//! # Determinism
//!
//! Floating-point addition is not associative, so a sum grouped by however the
//! thread pool happened to split it would give (very slightly) different answers
//! on different thread counts, and an iterative method can amplify that into
//! visibly different output. Every O(num_txps) reduction here (the Gram matrix
//! and its right-hand side) is therefore summed over fixed-size chunks
//! ([`VEC_CHUNK`]) whose partials are added in chunk order: the grouping depends
//! on the data length only, so the result is identical sequential or parallel
//! at any thread count. Element-wise work has no grouping to worry about and is
//! simply parallel. The per-class M-step inside `F` is parallel over its own
//! fixed partition for the same reason.

use crate::max_rel_diff_mode;
use crate::EmOptions;
use rayon::prelude::*;

/// DAAREM control parameters (the CRAN `daarem` defaults).
/// Kept at the reference implementation's values so behaviour matches the
/// published method rather than a local re-tuning.
const ORDER: usize = 5; // max Anderson history depth
const A1: f64 = 1.2; // damping base
const KAPPA: f64 = 25.0; // damping schedule exponent offset
const MON_TOL: f64 = 0.01; // objective may drop by at most this on an accepted step
const CYCL_MON_TOL: f64 = 0.0; // per-restart-cycle objective slack

/// Floor on an extrapolated abundance, as a fraction of the current iterate's.
///
/// Not part of the CRAN method, which extrapolates in an unconstrained space.
/// Abundances are non-negative, and the Anderson step, being a linear
/// combination of past differences, routinely overshoots below zero: on a
/// 646k-transcript GENCODE run, a quarter of the coordinates went negative. A
/// negative iterate is not a valid input to the M-step, so its residual is
/// meaningless; the acceptance test and the convergence check then act on
/// garbage, and the loop stopped early, far from the fixpoint (Spearman 0.90
/// against a fully converged EM, with whole classes truncated away).
///
/// Rejecting any step with a negative component is no fix: in that many
/// dimensions nearly every step has one, and DAAREM degenerates into plain EM at
/// twice the cost. Projecting instead keeps the step. The floor is relative
/// because a zero floor makes a coordinate absorbing under the multiplicative
/// M-step, while `0.1 * x` lets an extrapolation shrink an abundance by at most
/// 10x per step, so a wrongly drained transcript can still recover. On the
/// GENCODE run, floors from 0.01 to 0.25 all landed within Spearman 0.997-0.999
/// of the converged fixpoint in 340-540 M-steps; 0.5 damped the steps so much
/// that the method fell back to plain-EM behaviour.
const PROJ_FLOOR: f64 = 0.1;

/// Absolute floor under [`PROJ_FLOOR`], so repeated projections cannot drive an
/// abundance to exactly zero.
///
/// The objective test lets each accepted step give back up to `MON_TOL`, which
/// a small transcript's whole contribution can fit inside; repeated 10x
/// shrinks then underflowed such transcripts to `0.0`, where the
/// multiplicative M-step can never revive them. On a cerebellum run 13
/// transcripts with real support (one with 148 reads at the optimum) ended at
/// exactly zero, 0.21 nats short of the EM optimum. The floor sits well below
/// the final `min_alpha` truncation, so it changes nothing that is reported,
/// only what can recover; it never raises a coordinate above its current value.
const PROJ_ABS_FLOOR: f64 = 1e-10;

/// Jacobi eigendecomposition of a small symmetric `n×n` matrix `a` (row-major,
/// `n <= ORDER`). Overwrites `a`; returns eigenvalues in `eval` and the
/// eigenvectors as the columns of `evec` (row-major `n×n`). Deterministic: fixed
/// cyclic sweep order and iteration cap.
///
/// **What an eigendecomposition gives us.** It rewrites the matrix in terms of
/// directions that it merely scales, which turns the least-squares solve below
/// into independent one-dimensional divisions. **Why Jacobi:** the matrix here is
/// at most 5×5, so a simple, obviously-deterministic algorithm beats calling into
/// a tuned linear-algebra library whose result could vary with build flags.
///
/// The method repeatedly picks an off-diagonal entry and applies a rotation that
/// zeroes it; other entries grow back a little, but the total off-diagonal mass
/// shrinks every sweep, so it converges.
///
/// "Row-major" means the matrix is stored as one flat array where entry `(r, c)`
/// lives at index `r * n + c`.
fn jacobi_eig(a: &mut [f64], n: usize, eval: &mut [f64], evec: &mut [f64]) {
    // evec = identity
    //
    // The eigenvector matrix accumulates every rotation applied to `a`, so it
    // starts as the transformation that does nothing.
    for r in 0..n {
        for c in 0..n {
            evec[r * n + c] = if r == c { 1.0 } else { 0.0 };
        }
    }
    // A 1×1 matrix is already diagonal.
    if n == 1 {
        eval[0] = a[0];
        return;
    }
    // Fixed iteration cap rather than "loop until converged": bounded work, and
    // one less way for the result to depend on floating-point luck.
    for _sweep in 0..50 {
        // largest off-diagonal magnitude
        //
        // Sum of squares of the upper triangle: how far from diagonal we still
        // are. Once it is negligible the decomposition is done.
        let mut off = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[p * n + q] * a[p * n + q];
            }
        }
        if off <= 1e-30 {
            break;
        }
        // One cyclic sweep: visit every off-diagonal position in fixed order.
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                // Already effectively zero; rotating would be a no-op.
                if apq.abs() <= 1e-300 {
                    continue;
                }
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                // Givens angle that zeroes a[p][q] under A' = RᵀAR with
                // R = [[c, s], [-s, c]]: tan(2θ) = 2·apq / (aqq - app).
                //
                // `atan2` is used rather than `atan` of a ratio because it stays
                // well-defined when the denominator is zero (equal diagonals).
                let phi = 0.5 * (2.0 * apq).atan2(aqq - app);
                let (s, c) = phi.sin_cos();
                // rotate rows/cols p,q
                //
                // The similarity transform RᵀAR is applied as two passes: first
                // the columns, then the rows.
                for k in 0..n {
                    let akp = a[k * n + p];
                    let akq = a[k * n + q];
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p * n + k];
                    let aqk = a[q * n + k];
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                // Accumulate the same rotation into the eigenvector matrix.
                for k in 0..n {
                    let vkp = evec[k * n + p];
                    let vkq = evec[k * n + q];
                    evec[k * n + p] = c * vkp - s * vkq;
                    evec[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    // The matrix is now (near) diagonal, so its diagonal holds the eigenvalues.
    for i in 0..n {
        eval[i] = a[i * n + i];
    }
}

/// Choose the ridge damping `lambda` for the least-squares solve, targeting the
/// DAAREM shrink schedule. Direct port of the CRAN `daarem` `DampingFind` (a
/// Newton root-find on the scalar residual-norm equation `‖s(lambda)‖ = vk`).
/// `uy_sq[i] = (Uᵀ f)_i²`, `dvec[i]` = singular values, `sk` = accepted-step count.
/// Returns `(lambda, rr)` carried into the next call.
///
/// **What "ridge damping" does.** Adding `lambda` to the denominators shrinks the
/// solution, most strongly along the directions the data constrains least. Large
/// `lambda` ⇒ nearly no extrapolation (a plain EM step); small `lambda` ⇒ the full
/// Anderson step. The schedule below picks a `lambda` that shrinks the step to a
/// target fraction of the undamped one, with that fraction rising toward 1 as
/// `sk` (the count of accepted steps, i.e. accumulated trust) grows.
///
/// This is a direct port and is best read against the paper; the local names
/// follow the reference implementation deliberately so the two can be diffed.
#[allow(clippy::too_many_arguments)]
fn damping_find(
    uy_sq: &[f64],
    dvec: &[f64],
    aa: f64,
    kappa: f64,
    sk: f64,
    ftf: f64,
    lambda_start: f64,
    r_start: f64,
) -> (f64, f64) {
    // Drop zero singular values (their directions carry no information).
    let mut d: Vec<f64> = Vec::with_capacity(dvec.len());
    let mut u: Vec<f64> = Vec::with_capacity(dvec.len());
    for i in 0..dvec.len() {
        if dvec[i] > 0.0 {
            d.push(dvec[i]);
            u.push(uy_sq[i]);
        }
    }
    // Nothing to solve: keep the caller's current damping.
    if d.is_empty() {
        return (lambda_start, r_start);
    }
    // The shrink schedule: `target` rises from near 0 toward 1 as `sk` grows,
    // annealing from "almost a plain EM step" to "the full Anderson step".
    let pow = kappa - sk;
    let target = (-0.5 * (1.0 + aa.powf(pow)).ln()).exp();
    let d_sq: Vec<f64> = d.iter().map(|x| x * x).collect();
    // Norm of the *undamped* least-squares solution, the thing being shrunk.
    let betahat_ls_norm: f64 = (0..d.len()).map(|i| u[i] / d_sq[i]).sum::<f64>().sqrt();
    let vk = target * betahat_ls_norm;
    if vk == 0.0 {
        return (lambda_start, r_start);
    }
    // Newton iteration on `‖s(lambda)‖ - vk = 0`, bracketed by [ll, uu].
    let mut lambda = lambda_start - r_start / vk;
    let denom_ll: f64 = (0..d.len()).map(|i| u[i] / (d_sq[i] * d_sq[i])).sum();
    let mut ll = (betahat_ls_norm * (betahat_ls_norm - vk)) / denom_ll;
    let mut uu = ftf / vk;
    // Acceptance band around the target: solving exactly is unnecessary, so the
    // iteration stops as soon as the shrink is within half a schedule step.
    let pow_low = pow + 0.5;
    let pow_up = pow - 0.5;
    let l_stop = (-0.5 * (1.0 + aa.powf(pow_low)).ln()).exp();
    let u_stop = (-0.5 * (1.0 + aa.powf(pow_up)).ln()).exp();

    let mut s_norm = 0.0;
    let mut phi_ratio = 0.0;
    // Capped at 10: this is an inner heuristic, not something worth converging.
    for _ in 0..10 {
        // Newton can step outside the bracket; fall back to a geometric mean,
        // which is the standard safeguard for a positive-valued parameter.
        if lambda <= ll || lambda >= uu {
            lambda = (0.0001 * uu).max((ll * uu).sqrt());
        }
        // `sn` is ‖s(lambda)‖², `der` its derivative, both accumulated per
        // singular direction.
        let mut sn = 0.0;
        let mut der = 0.0;
        for i in 0..d.len() {
            let denom = d_sq[i] + lambda;
            let dl = (d[i] / denom) * (d[i] / denom);
            sn += u[i] * dl;
            der += u[i] * (dl / denom);
        }
        s_norm = sn.sqrt();
        let phi_val = s_norm - vk;
        let phi_der = -der / s_norm;
        phi_ratio = phi_val / phi_der;
        // Inside the acceptance band: good enough.
        if s_norm <= u_stop * betahat_ls_norm && s_norm >= l_stop * betahat_ls_norm {
            break;
        }
        // Tighten the bracket, then take the (scaled) Newton step.
        uu = if phi_val >= 0.0 { uu } else { lambda };
        ll = ll.max(lambda - phi_ratio);
        lambda -= (s_norm * phi_ratio) / vk;
    }
    (lambda, s_norm * phi_ratio)
}

/// Fixed chunk for the O(num_txps) vector work below. Every reduction sums
/// per-chunk partials in chunk order, so the grouping of the floating-point sums
/// is a function of the data length alone — identical sequential or parallel,
/// at any thread count (the same scheme as `max_rel_diff_par`).
const VEC_CHUNK: usize = 8192;

/// Upper bound on the inner products one step needs: the `np(np+1)/2` Gram
/// entries plus `np` right-hand sides, at `np <= ORDER`.
const MAX_DOTS: usize = ORDER * (ORDER + 1) / 2 + ORDER;

/// `gram = Fwᵀ Fw` and `bvec = Fwᵀ fnew` over the `np` active columns of
/// `fdiff`, in one chunked pass (see [`VEC_CHUNK`]).
///
/// These were plain sequential loops, and with up to 20 dot products over every
/// transcript per step they cost several times a (parallel) M-step on a
/// transcriptome-sized problem — enough to make DAAREM slower in wall time than
/// plain EM even while taking far fewer steps.
fn gram_and_rhs(
    fdiff: &[Vec<f64>],
    fnew: &[f64],
    np: usize,
    gram: &mut [f64],
    bvec: &mut [f64],
    parallel: bool,
) {
    let n = fnew.len();
    let chunk = |c: usize| -> [f64; MAX_DOTS] {
        let lo = c * VEC_CHUNK;
        let hi = (lo + VEC_CHUNK).min(n);
        let mut out = [0.0f64; MAX_DOTS];
        let mut k = 0;
        for a in 0..np {
            for b in a..np {
                let (fa, fb) = (&fdiff[a][lo..hi], &fdiff[b][lo..hi]);
                out[k] = fa.iter().zip(fb).map(|(x, y)| x * y).sum();
                k += 1;
            }
            out[k] = fdiff[a][lo..hi]
                .iter()
                .zip(&fnew[lo..hi])
                .map(|(x, y)| x * y)
                .sum();
            k += 1;
        }
        out
    };
    let nchunks = n.div_ceil(VEC_CHUNK);
    let parts: Vec<[f64; MAX_DOTS]> = if parallel {
        (0..nchunks).into_par_iter().map(chunk).collect()
    } else {
        (0..nchunks).map(chunk).collect()
    };
    let mut k = 0;
    for a in 0..np {
        for b in a..np {
            let s: f64 = parts.iter().map(|p| p[k]).sum();
            // Symmetric: fill both halves from one dot product.
            gram[a * np + b] = s;
            gram[b * np + a] = s;
            k += 1;
        }
        bvec[a] = parts.iter().map(|p| p[k]).sum();
        k += 1;
    }
}

/// Apply `body(lo, hi, out)` to `out` in [`VEC_CHUNK`] pieces, in parallel when
/// asked. Only for element-wise work (no reductions), so the result is the same
/// either way.
fn chunked_mut(out: &mut [f64], parallel: bool, body: impl Fn(usize, &mut [f64]) + Sync + Send) {
    if parallel {
        out.par_chunks_mut(VEC_CHUNK)
            .enumerate()
            .for_each(|(c, o)| body(c * VEC_CHUNK, o));
    } else {
        out.chunks_mut(VEC_CHUNK)
            .enumerate()
            .for_each(|(c, o)| body(c * VEC_CHUNK, o));
    }
}

/// DAAREM (objective-monotonicity variant, CRAN `daarem_base_objfn`) driving
/// the fixed-point map `f`, which also returns the objective at its input. `x0` aliases the caller's abundance vector and receives the final
/// estimate. `it` accumulates the total number of `F` evaluations (M-steps), for
/// parity with the plain/SQUAREM loops. Convergence uses salmon's own
/// [`max_rel_diff`] criterion on both successive iterates and the iterate
/// versus one plain M-step from it (see the loop). Returns whether it converged.
///
/// Counting `F` evaluations rather than loop trips is what makes `--maxIter`
/// mean the same amount of work in every mode.
#[allow(clippy::too_many_arguments)]
pub(crate) fn daarem_loop(
    f: &mut impl FnMut(&[f64], &mut [f64]) -> f64,
    x0: &mut [f64],
    num_txps: usize,
    opts: &EmOptions,
    min_iter: u32,
    it: &mut u32,
    parallel: bool,
    rel_diff_partials: &mut Vec<f64>,
) -> bool {
    let n = num_txps;
    if n == 0 {
        return true;
    }
    // History depth: at most ORDER, and never more than half the problem
    // dimension (a longer window than the problem has directions is degenerate).
    let nlag = ORDER.min(n.div_ceil(2)).max(1);

    // Iterate/residual state.
    //
    // Two iterates and their residuals are needed before any extrapolation is
    // possible, so the loop is primed with two plain EM steps.
    let mut xold = x0.to_vec();
    let mut xnew = vec![0.0f64; n];
    f(&xold, &mut xnew); // xnew = F(xold)
    *it += 1;
    let mut fold: Vec<f64> = (0..n).map(|i| xnew[i] - xold[i]).collect();
    let mut fx = vec![0.0f64; n];
    // Every point the loop evaluates `F` at is a point it needs the objective
    // at, so `obj_cur` is always "the objective at `xnew`" for free.
    let mut obj_cur = f(&xnew, &mut fx); // F(xnew), L(xnew)
    *it += 1;
    let mut fnew: Vec<f64> = (0..n).map(|i| fx[i] - xnew[i]).collect();
    // Objective at the end of the previous restart cycle (`ell.star`).
    let mut ell_star = obj_cur;

    // Anderson history windows (column c is a length-n difference vector).
    //
    // These hold *differences* between successive iterates and residuals, which
    // is the information the multi-secant step extrapolates from.
    let mut fdiff: Vec<Vec<f64>> = vec![vec![0.0f64; n]; nlag];
    let mut xdiff: Vec<Vec<f64>> = vec![vec![0.0f64; n]; nlag];

    // Small (np×np) work buffers, np <= nlag.
    //
    // Allocated once outside the loop: the linear algebra is tiny but runs every
    // iteration, and per-iteration allocation would dominate it.
    let mut gram = vec![0.0f64; nlag * nlag];
    let mut bvec = vec![0.0f64; nlag];
    let mut eval = vec![0.0f64; nlag];
    let mut evec = vec![0.0f64; nlag * nlag];
    let mut dvec = vec![0.0f64; nlag];
    let mut uy = vec![0.0f64; nlag];
    let mut uy_sq = vec![0.0f64; nlag];
    let mut dd = vec![0.0f64; nlag];
    let mut gamma = vec![0.0f64; nlag];

    let mut x_prop = vec![0.0f64; n];
    let mut f_prop = vec![0.0f64; n];

    let mut count = 0usize; // active history depth since last restart
    let mut shrink_count = 0.0f64;
    // Starts enormous, i.e. maximum damping: the first steps are essentially
    // plain EM until some history has been accumulated and accepted.
    let mut lambda_ridge = 100_000.0f64;
    let mut r_penalty = 0.0f64;
    let mut converged = false;

    while *it < opts.max_iter {
        // Record this step's iterate/residual differences into the window.
        let col = count; // 0-based store slot for this step
        chunked_mut(&mut fdiff[col], parallel, |lo, o| {
            for (j, v) in o.iter_mut().enumerate() {
                *v = fnew[lo + j] - fold[lo + j];
            }
        });
        chunked_mut(&mut xdiff[col], parallel, |lo, o| {
            for (j, v) in o.iter_mut().enumerate() {
                *v = xnew[lo + j] - xold[lo + j];
            }
        });
        count += 1;
        let np = count;

        // Gram = Fwᵀ Fw  and  b = Fwᵀ fnew   over the np active columns.
        //
        // The normal equations of the least-squares problem "which combination of
        // recent residual differences best cancels the current residual?".
        // Computing the small np×np Gram matrix is far cheaper than working with
        // the n×np window directly, since n is the transcript count.
        gram_and_rhs(&fdiff, &fnew, np, &mut gram, &mut bvec, parallel);

        // Eigendecompose Gram = V diag(eval) Vᵀ  (eval = singular values squared).
        //
        // Since Gram = FᵀF, its eigenvalues are the squares of the window's
        // singular values, which is exactly what the ridge solve needs.
        jacobi_eig(
            &mut gram[..np * np],
            np,
            &mut eval[..np],
            &mut evec[..np * np],
        );
        for a in 0..np {
            // `.max(0.0)` guards against a tiny negative eigenvalue from rounding.
            dvec[a] = eval[a].max(0.0).sqrt();
        }
        // uy = Uᵀ fnew = D⁻¹ Vᵀ b ; Ftf = ‖Vᵀ b‖.
        let mut ftf_sq = 0.0;
        for a in 0..np {
            let mut vtb = 0.0;
            for r in 0..np {
                vtb += evec[r * np + a] * bvec[r];
            }
            ftf_sq += vtb * vtb;
            uy[a] = if dvec[a] > 0.0 { vtb / dvec[a] } else { 0.0 };
            uy_sq[a] = uy[a] * uy[a];
        }
        let ftf = ftf_sq.sqrt();

        // Pick how much to shrink this step (see `damping_find`).
        let (lam, rr) = damping_find(
            &uy_sq[..np],
            &dvec[..np],
            A1,
            KAPPA,
            shrink_count,
            ftf,
            lambda_ridge,
            r_penalty,
        );
        lambda_ridge = lam;
        r_penalty = rr;

        // dd = (dvec·uy)/(dvec²+λ) ; gamma = V dd.
        //
        // The ridge-regularized least-squares coefficients: `+λ` in the
        // denominator is the damping, and rotating by V returns them to the
        // original coordinates.
        for a in 0..np {
            dd[a] = (dvec[a] * uy[a]) / (dvec[a] * dvec[a] + lambda_ridge);
        }
        for r in 0..np {
            let mut g = 0.0;
            for a in 0..np {
                g += evec[r * np + a] * dd[a];
            }
            gamma[r] = g;
        }

        // x_prop = (xnew - Xw·gamma) + (fnew - Fw·gamma).
        //
        // The Anderson step: correct both the current iterate and its residual by
        // the fitted combination of recent differences, then take one implied
        // fixed-point step. This is the "jump ahead" the whole file exists for.
        // The result is projected onto the floor `PROJ_FLOOR * xnew` (but never
        // below `min(xnew, PROJ_ABS_FLOOR)`), which keeps every iterate
        // non-negative and every positive one positive (see both constants).
        let gamma_np = &gamma[..np];
        chunked_mut(&mut x_prop, parallel, |lo, o| {
            for (j, out) in o.iter_mut().enumerate() {
                let i = lo + j;
                let mut xg = 0.0;
                let mut fg = 0.0;
                for (a, &g) in gamma_np.iter().enumerate() {
                    xg += xdiff[a][i] * g;
                    fg += fdiff[a][i] * g;
                }
                let floor = (PROJ_FLOOR * xnew[i]).max(xnew[i].min(PROJ_ABS_FLOOR));
                *out = ((xnew[i] - xg) + (fnew[i] - fg)).max(floor);
            }
        });

        let obj_prop = f(&x_prop, &mut fx); // F(x_prop), L(x_prop)
        *it += 1;
        chunked_mut(&mut f_prop, parallel, |lo, o| {
            for (j, v) in o.iter_mut().enumerate() {
                *v = fx[lo + j] - x_prop[lo + j];
            }
        });

        // Track the previous iterate for the convergence check.
        std::mem::swap(&mut xold, &mut xnew); // xold = previous xnew

        // Epsilon-monotonicity on the objective: the extrapolated point is
        // accepted unless it lowers the objective by more than `MON_TOL`. Plain
        // EM/VBEM steps never lower it, so this can only refuse a step that is
        // genuinely worse — unlike the residual test it replaces, which also
        // refused VBEM's own path to the optimum (see the module docs).
        if !obj_prop.is_nan() && !obj_cur.is_nan() && obj_prop >= obj_cur - MON_TOL {
            // Accept the (damped) Anderson step.
            std::mem::swap(&mut fold, &mut fnew);
            xnew.copy_from_slice(&x_prop);
            fnew.copy_from_slice(&f_prop);
            // One more accepted step ⇒ less damping next time.
            shrink_count += 1.0;
            obj_cur = obj_prop;
        } else {
            // Reject: take a plain fixed-point (EM) step from the previous xnew.
            // The accelerated evaluation is discarded, so a rejection costs one
            // wasted `F` call — the price of the attempt.
            std::mem::swap(&mut fold, &mut fnew);
            // xnew = F(xold) = xold + fold  (fold is the residual at xold).
            chunked_mut(&mut xnew, parallel, |lo, o| {
                for (j, v) in o.iter_mut().enumerate() {
                    *v = xold[lo + j] + fold[lo + j];
                }
            });
            obj_cur = f(&xnew, &mut fx);
            *it += 1;
            chunked_mut(&mut fnew, parallel, |lo, o| {
                for (j, v) in o.iter_mut().enumerate() {
                    *v = fx[lo + j] - xnew[lo + j];
                }
            });
            // An objective that cannot be evaluated says the history led
            // somewhere meaningless; drop it, as the reference does.
            if obj_prop.is_nan() {
                count = 0;
            }
        }

        // Convergence on salmon's own `max_rel_diff` criterion, required of two
        // pairs at once:
        //   * successive iterates — the rule every mode uses; but an Anderson
        //     step can stagnate briefly, and on an EM fixture that alone ended
        //     runs 0.66 nats short of the optimum;
        //   * the iterate and one plain M-step from it, `F(xnew) = xnew + fnew`
        //     (already computed) — the question plain EM's rule asks; but where
        //     EM mixes slowly one step barely moves anything far from the
        //     fixpoint, which alone stopped a slow-mixing case ~30% off.
        // Each is fooled where the other is not, so both must agree.
        // `min_iter` prevents an early accidental "no change" from ending the run.
        if *it >= min_iter {
            let d_iter = max_rel_diff_mode(
                &xold,
                &xnew,
                opts.alpha_check_cutoff,
                parallel,
                rel_diff_partials,
            );
            // `x_prop` is free until the next proposal overwrites it.
            chunked_mut(&mut x_prop, parallel, |lo, o| {
                for (j, v) in o.iter_mut().enumerate() {
                    *v = xnew[lo + j] + fnew[lo + j];
                }
            });
            let d_step = max_rel_diff_mode(
                &xnew,
                &x_prop,
                opts.alpha_check_cutoff,
                parallel,
                rel_diff_partials,
            );
            let d = d_iter.max(d_step);
            if d.is_finite() && d < opts.rel_diff_tol {
                converged = true;
                break;
            }
        }

        // Periodic restart of the Anderson history.
        //
        // The window is finite, and stale differences describe a region the
        // iteration has left; restarting drops them wholesale, which is the "R"
        // in DAAREM.
        //
        // At each restart the cycle as a whole must not have lost objective; if
        // it did, the damping credit earned by the cycle's accepted steps is
        // taken back, so the next cycle extrapolates more cautiously.
        if count == nlag {
            count = 0;
            if obj_cur < ell_star - CYCL_MON_TOL {
                shrink_count = (shrink_count - nlag as f64).max(-2.0 * KAPPA);
            }
            ell_star = obj_cur;
        }
    }

    x0.copy_from_slice(&xnew);
    converged
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PackedEqClasses;
    use salmon_eqclass::{EquivalenceClassBuilder, TranscriptGroup};

    /// The eigensolver is the one piece here with a checkable closed-form answer,
    /// so pin it on a matrix whose spectrum is known by hand — and check the
    /// eigenvectors come out orthogonal, since the ridge solve assumes it.
    #[test]
    fn jacobi_recovers_known_spectrum() {
        // Symmetric 2x2 [[2,1],[1,2]] has eigenvalues 1 and 3.
        let mut a = vec![2.0, 1.0, 1.0, 2.0];
        let mut eval = vec![0.0; 2];
        let mut evec = vec![0.0; 4];
        jacobi_eig(&mut a, 2, &mut eval, &mut evec);
        // Order is not guaranteed, so sort before comparing.
        let mut evs = eval.clone();
        evs.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!((evs[0] - 1.0).abs() < 1e-9, "{evs:?}");
        assert!((evs[1] - 3.0).abs() < 1e-9, "{evs:?}");
        // Eigenvectors orthonormal: columns dot to ~0.
        let dot = evec[0] * evec[1] + evec[2] * evec[3];
        assert!(dot.abs() < 1e-9, "columns not orthogonal: {dot}");
    }

    /// A transcriptome-shaped problem: heavy-tailed abundances, a third of the
    /// transcripts absent, and every fragment multi-mapping to a few decoy
    /// transcripts besides its source. The absent transcripts are what matter —
    /// EM drains them toward zero slowly, and that slow drain is exactly where
    /// an Anderson extrapolation overshoots into negative abundances.
    fn skewed_fixture(num_genes: usize, num_frags: usize) -> PackedEqClasses {
        let mut seed = 0x2545_F491_4F6C_DD1Du64;
        let mut rnd = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let mut unif = move || ((rnd() >> 11) as f64) / ((1u64 << 53) as f64);
        // Genes of 1..=8 isoforms. Gene abundance is heavy-tailed (1/u^2); within
        // a gene about a third of the isoforms are absent.
        let mut genes: Vec<(u32, usize)> = Vec::with_capacity(num_genes);
        let mut iso_share: Vec<f64> = Vec::new();
        let mut cum = Vec::with_capacity(num_genes);
        let mut acc = 0.0;
        let mut next = 0u32;
        for _ in 0..num_genes {
            let n = 1 + (unif() * 8.0) as usize;
            genes.push((next, n));
            for _ in 0..n {
                iso_share.push(if unif() < 0.35 { 0.0 } else { unif() });
            }
            next += n as u32;
            let u = unif().max(1e-6);
            acc += (1.0 / (u * u)).min(1e6);
            cum.push(acc);
        }
        let num_txps = next as usize;
        let b = EquivalenceClassBuilder::new();
        for _ in 0..num_frags {
            let r = unif() * acc;
            let g = cum.partition_point(|&c| c < r).min(num_genes - 1);
            let (first, n) = genes[g];
            let shares = &iso_share[first as usize..first as usize + n];
            let tot: f64 = shares.iter().sum();
            if tot == 0.0 {
                continue;
            }
            let mut r = unif() * tot;
            let mut src = n - 1;
            for (i, &sh) in shares.iter().enumerate() {
                if r < sh {
                    src = i;
                    break;
                }
                r -= sh;
            }
            // Isoforms of a gene share most of their sequence: a fragment is
            // usually compatible with most of its siblings too.
            let mut tids: Vec<u32> = (0..n)
                .filter(|&i| i == src || unif() < 0.8)
                .map(|i| first + i as u32)
                .collect();
            if unif() < 0.1 {
                tids.push((unif() * num_txps as f64) as u32);
            }
            tids.sort_unstable();
            tids.dedup();
            let k = tids.len();
            b.add_group(TranscriptGroup::new(tids), vec![1.0; k], 1);
        }
        let mut eq = b.finish();
        let eff: Vec<f64> = (0..num_txps)
            .map(|i| 200.0 + (i % 1013) as f64 * 3.1)
            .collect();
        eq.update_eff_lengths(&eff);
        PackedEqClasses::from_collapsed(&eq, num_txps)
    }

    /// Log-likelihood up to a constant. The EM preserves total mass, so
    /// `sum_c n_c log(sum_t alpha_t w_ct)` ranks estimates of equal total.
    fn log_lik(p: &PackedEqClasses, alphas: &[f64]) -> f64 {
        (0..p.counts.len())
            .map(|c| {
                let (s, e) = (p.starts[c] as usize, p.starts[c + 1] as usize);
                let d: f64 = (s..e)
                    .map(|j| alphas[p.labels[j] as usize] * p.combined[j])
                    .sum();
                p.counts[c] as f64 * d.ln()
            })
            .sum()
    }

    /// The defect this guards: the Anderson step overshot below zero and the
    /// negative iterate was kept, so the M-step, the acceptance test and the
    /// convergence check all ran on an invalid vector. Read before truncation
    /// (`run_em_counts`), since truncation would silently clamp the evidence
    /// away. Before the projection this fixture left thousands of negative
    /// abundances.
    #[test]
    fn daarem_iterates_stay_non_negative() {
        let p = skewed_fixture(2_000, 200_000);
        let opts = EmOptions {
            accel: crate::EmAccel::Daarem,
            ..EmOptions::default()
        };
        let (alphas, _, converged) = crate::run_em_counts(
            &p,
            &p.counts,
            &opts,
            true,
            opts.min_iter,
            crate::InitAlphas::NONE,
            crate::EffLens::NONE,
        );
        assert!(converged, "did not converge");
        let bad = alphas.iter().filter(|a| a.is_nan() || **a < 0.0).count();
        assert_eq!(bad, 0, "{bad} negative or NaN abundances");
    }

    /// Acceleration must not buy speed with accuracy. DAAREM stops on the same
    /// question plain EM asks — would one more M-step move anything by more
    /// than `rel_diff_tol`? — so at a given tolerance it should land about as
    /// close to the (unique, EM being concave) optimum as plain EM, in far fewer
    /// M-steps. With the objective test but convergence judged on successive
    /// iterates alone, this fixture stopped at 71 M-steps 0.66 nats short at
    /// `1e-3`, where plain EM is 0.0015 short.
    #[test]
    fn daarem_em_stops_as_close_as_plain_em() {
        let p = skewed_fixture(2_000, 200_000);
        let run = |accel, rel_diff_tol, max_iter| {
            crate::optimize_packed(
                &p,
                &EmOptions {
                    accel,
                    rel_diff_tol,
                    max_iter,
                    ..EmOptions::default()
                },
                true,
            )
        };
        let best = log_lik(&p, &run(crate::EmAccel::Squarem, 1e-9, 200_000).alphas);
        let plain = run(crate::EmAccel::None, 1e-3, 100_000);
        let da = run(crate::EmAccel::Daarem, 1e-3, 100_000);
        assert!(da.converged);
        assert_eq!(da.dropped_mass, 0.0, "DAAREM truncated whole classes away");
        let gap_plain = best - log_lik(&p, &plain.alphas);
        let gap_da = best - log_lik(&p, &da.alphas);
        assert!(
            gap_da < 0.01 && gap_da < 5.0 * gap_plain.max(1e-3),
            "DAAREM {gap_da} nats short of the optimum, plain EM {gap_plain}"
        );
        assert!(
            da.iters * 5 < plain.iters,
            "{} vs {}",
            da.iters,
            plain.iters
        );
    }
}
