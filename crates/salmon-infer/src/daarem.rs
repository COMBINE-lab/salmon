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
//! Port of the **objective-free** DAAREM variant (Damped Anderson Acceleration
//! with Restarts and Epsilon-Monotonicity) from Henderson & Varadhan,
//! *J. Comput. Graph. Statist.* 28(4):834-846 (2019), as implemented in the CRAN
//! `daarem` package's `daarem_base_noobjfn`. Salmon's EM/VBEM does not evaluate a
//! likelihood, so the residual-norm monotonicity control (the no-objective
//! variant) is the natural fit: a proposed Anderson step is accepted only if it
//! does not increase the residual norm `‖F(x) - x‖` beyond a shrinking tolerance,
//! otherwise it falls back to a plain fixed-point (EM) step.
//!
//! The residual `F(x) - x` is how far the map still moves the point; it is zero
//! exactly at the answer, so its norm is a natural "how far from done" measure
//! even without a likelihood to evaluate.
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
//! All reductions (Gram matrix, norms, matrix-vector products) are computed
//! sequentially in index order, so the result is deterministic regardless of
//! thread count. Floating-point addition is not associative, so a parallel sum
//! would give (very slightly) different answers on different thread counts, and
//! an iterative method can amplify that into visibly different output. The
//! expensive per-class M-step stays inside the caller-supplied map `F`, which is
//! already rayon-parallel — but its parallelism is over a fixed partition, so it
//! stays reproducible.

use crate::max_rel_diff;
use crate::EmOptions;

/// DAAREM control parameters (the CRAN `daarem` no-objective defaults).
/// Kept at the reference implementation's values so behaviour matches the
/// published method rather than a local re-tuning.
const ORDER: usize = 5; // max Anderson history depth
const A1: f64 = 1.2; // damping base
const KAPPA: f64 = 25.0; // damping schedule exponent offset
const RESID_TOL: f64 = 0.95; // rho in the residual acceptance test

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

/// L2 norm of a slice, summed in index order for determinism.
///
/// Written out rather than using an iterator sum so the accumulation order is
/// unambiguously left to right; see the determinism note in the module docs.
fn norm(v: &[f64]) -> f64 {
    let mut s = 0.0;
    for &x in v {
        s += x * x;
    }
    s.sqrt()
}

/// DAAREM (no-objective, residual-monotonicity variant) driving the fixed-point
/// map `f`. `x0` aliases the caller's abundance vector and receives the final
/// estimate. `it` accumulates the total number of `F` evaluations (M-steps), for
/// parity with the plain/SQUAREM loops. Convergence uses salmon's own
/// [`max_rel_diff`] criterion (not DAAREM's residual tolerance) so all three
/// acceleration modes stop on the same rule. Returns whether it converged.
///
/// Counting `F` evaluations rather than loop trips is what makes `--maxIter`
/// mean the same amount of work in every mode.
#[allow(clippy::too_many_arguments)]
pub(crate) fn daarem_loop(
    f: &mut impl FnMut(&[f64], &mut [f64]),
    x0: &mut [f64],
    num_txps: usize,
    opts: &EmOptions,
    min_iter: u32,
    it: &mut u32,
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
    f(&xnew, &mut fx); // F(xnew)
    *it += 1;
    let mut fnew: Vec<f64> = (0..n).map(|i| fx[i] - xnew[i]).collect();
    let mut ss_resids = norm(&fnew);

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
    let mut k = 1i32;
    let mut converged = false;

    while *it < opts.max_iter {
        // Record this step's iterate/residual differences into the window.
        let col = count; // 0-based store slot for this step
        for i in 0..n {
            fdiff[col][i] = fnew[i] - fold[i];
            xdiff[col][i] = xnew[i] - xold[i];
        }
        count += 1;
        let np = count;

        // Gram = Fwᵀ Fw  and  b = Fwᵀ fnew   over the np active columns.
        //
        // The normal equations of the least-squares problem "which combination of
        // recent residual differences best cancels the current residual?".
        // Computing the small np×np Gram matrix is far cheaper than working with
        // the n×np window directly, since n is the transcript count.
        for a in 0..np {
            for b in a..np {
                let mut s = 0.0;
                for i in 0..n {
                    s += fdiff[a][i] * fdiff[b][i];
                }
                // Symmetric: fill both halves from one dot product.
                gram[a * np + b] = s;
                gram[b * np + a] = s;
            }
            let mut sb = 0.0;
            for i in 0..n {
                sb += fdiff[a][i] * fnew[i];
            }
            bvec[a] = sb;
        }

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
        for i in 0..n {
            let mut xg = 0.0;
            let mut fg = 0.0;
            for a in 0..np {
                xg += xdiff[a][i] * gamma[a];
                fg += fdiff[a][i] * gamma[a];
            }
            x_prop[i] = (xnew[i] - xg) + (fnew[i] - fg);
        }

        f(&x_prop, &mut fx); // F(x_prop)
        *it += 1;
        for i in 0..n {
            f_prop[i] = fx[i] - x_prop[i];
        }
        let ss_prop = norm(&f_prop);

        // Track the previous iterate for the convergence check.
        std::mem::swap(&mut xold, &mut xnew); // xold = previous xnew

        // Epsilon-monotonicity: the extrapolated point is accepted only if its
        // residual is no worse than the current one, allowing a slack of
        // `RESID_TOL^k` that shrinks with the iteration count. Early on, a small
        // uphill move is tolerated (it often pays off); later the test becomes
        // strict, which is what guarantees the method cannot wander.
        if ss_prop <= ss_resids * (1.0 + RESID_TOL.powi(k)) {
            // Accept the (damped) Anderson step.
            std::mem::swap(&mut fold, &mut fnew);
            xnew.copy_from_slice(&x_prop);
            fnew.copy_from_slice(&f_prop);
            // One more accepted step ⇒ less damping next time.
            shrink_count += 1.0;
            ss_resids = ss_prop;
        } else {
            // Reject: take a plain fixed-point (EM) step from the previous xnew.
            // The accelerated evaluation is discarded, so a rejection costs one
            // wasted `F` call — the price of the attempt.
            std::mem::swap(&mut fold, &mut fnew);
            // xnew = F(xold) = xold + fold  (fold is the residual at xold).
            for i in 0..n {
                xnew[i] = xold[i] + fold[i];
            }
            f(&xnew, &mut fx);
            *it += 1;
            for i in 0..n {
                fnew[i] = fx[i] - xnew[i];
            }
            ss_resids = norm(&fnew);
        }

        // Convergence on salmon's own criterion, not DAAREM's residual test, so
        // every acceleration mode stops at the same place. `min_iter` prevents an
        // early accidental "no change" from ending the run.
        if *it >= min_iter {
            let d = max_rel_diff(&xold, &xnew, opts.alpha_check_cutoff);
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
        if count == nlag {
            count = 0;
        }
        k += 1;
    }

    x0.copy_from_slice(&xnew);
    converged
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
