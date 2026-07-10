//! DAAREM acceleration of the EM/VBEM fixed-point iteration.
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
//! Where SQUAREM extrapolates from only the last two iterates, DAAREM solves a
//! damped least-squares problem over a window of the last `nlag` residual/iterate
//! differences (a multi-secant quasi-Newton step), which converges faster on
//! high-dimensional, ill-conditioned problems — the regime of transcript-level
//! quantification. The damping (a ridge on the least-squares solve) starts strong
//! (near a plain EM step) and anneals toward the full Anderson step as accepted
//! steps accumulate; the history is periodically restarted.
//!
//! All reductions (Gram matrix, norms, matrix-vector products) are computed
//! sequentially in index order, so the result is deterministic regardless of
//! thread count. The expensive per-class M-step stays inside the caller-supplied
//! map `F`, which is already rayon-parallel.

use crate::max_rel_diff;
use crate::EmOptions;

/// DAAREM control parameters (the CRAN `daarem` no-objective defaults).
const ORDER: usize = 5; // max Anderson history depth
const A1: f64 = 1.2; // damping base
const KAPPA: f64 = 25.0; // damping schedule exponent offset
const RESID_TOL: f64 = 0.95; // rho in the residual acceptance test

/// Jacobi eigendecomposition of a small symmetric `n×n` matrix `a` (row-major,
/// `n <= ORDER`). Overwrites `a`; returns eigenvalues in `eval` and the
/// eigenvectors as the columns of `evec` (row-major `n×n`). Deterministic: fixed
/// cyclic sweep order and iteration cap.
fn jacobi_eig(a: &mut [f64], n: usize, eval: &mut [f64], evec: &mut [f64]) {
    // evec = identity
    for r in 0..n {
        for c in 0..n {
            evec[r * n + c] = if r == c { 1.0 } else { 0.0 };
        }
    }
    if n == 1 {
        eval[0] = a[0];
        return;
    }
    for _sweep in 0..50 {
        // largest off-diagonal magnitude
        let mut off = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[p * n + q] * a[p * n + q];
            }
        }
        if off <= 1e-30 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = a[p * n + q];
                if apq.abs() <= 1e-300 {
                    continue;
                }
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                // Givens angle that zeroes a[p][q] under A' = RᵀAR with
                // R = [[c, s], [-s, c]]: tan(2θ) = 2·apq / (aqq - app).
                let phi = 0.5 * (2.0 * apq).atan2(aqq - app);
                let (s, c) = phi.sin_cos();
                // rotate rows/cols p,q
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
                for k in 0..n {
                    let vkp = evec[k * n + p];
                    let vkq = evec[k * n + q];
                    evec[k * n + p] = c * vkp - s * vkq;
                    evec[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    for i in 0..n {
        eval[i] = a[i * n + i];
    }
}

/// Choose the ridge damping `lambda` for the least-squares solve, targeting the
/// DAAREM shrink schedule. Direct port of the CRAN `daarem` `DampingFind` (a
/// Newton root-find on the scalar residual-norm equation `‖s(lambda)‖ = vk`).
/// `uy_sq[i] = (Uᵀ f)_i²`, `dvec[i]` = singular values, `sk` = accepted-step count.
/// Returns `(lambda, rr)` carried into the next call.
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
    if d.is_empty() {
        return (lambda_start, r_start);
    }
    let pow = kappa - sk;
    let target = (-0.5 * (1.0 + aa.powf(pow)).ln()).exp();
    let d_sq: Vec<f64> = d.iter().map(|x| x * x).collect();
    let betahat_ls_norm: f64 = (0..d.len()).map(|i| u[i] / d_sq[i]).sum::<f64>().sqrt();
    let vk = target * betahat_ls_norm;
    if vk == 0.0 {
        return (lambda_start, r_start);
    }
    let mut lambda = lambda_start - r_start / vk;
    let denom_ll: f64 = (0..d.len()).map(|i| u[i] / (d_sq[i] * d_sq[i])).sum();
    let mut ll = (betahat_ls_norm * (betahat_ls_norm - vk)) / denom_ll;
    let mut uu = ftf / vk;
    let pow_low = pow + 0.5;
    let pow_up = pow - 0.5;
    let l_stop = (-0.5 * (1.0 + aa.powf(pow_low)).ln()).exp();
    let u_stop = (-0.5 * (1.0 + aa.powf(pow_up)).ln()).exp();

    let mut s_norm = 0.0;
    let mut phi_ratio = 0.0;
    for _ in 0..10 {
        if lambda <= ll || lambda >= uu {
            lambda = (0.0001 * uu).max((ll * uu).sqrt());
        }
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
        if s_norm <= u_stop * betahat_ls_norm && s_norm >= l_stop * betahat_ls_norm {
            break;
        }
        uu = if phi_val >= 0.0 { uu } else { lambda };
        ll = ll.max(lambda - phi_ratio);
        lambda -= (s_norm * phi_ratio) / vk;
    }
    (lambda, s_norm * phi_ratio)
}

/// L2 norm of a slice, summed in index order for determinism.
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
    let nlag = ORDER.min(n.div_ceil(2)).max(1);

    // Iterate/residual state.
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
    let mut fdiff: Vec<Vec<f64>> = vec![vec![0.0f64; n]; nlag];
    let mut xdiff: Vec<Vec<f64>> = vec![vec![0.0f64; n]; nlag];

    // Small (np×np) work buffers, np <= nlag.
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
    let mut lambda_ridge = 100_000.0f64;
    let mut r_penalty = 0.0f64;
    let mut k = 1i32;
    let mut converged = false;

    while *it < opts.max_iter {
        let col = count; // 0-based store slot for this step
        for i in 0..n {
            fdiff[col][i] = fnew[i] - fold[i];
            xdiff[col][i] = xnew[i] - xold[i];
        }
        count += 1;
        let np = count;

        // Gram = Fwᵀ Fw  and  b = Fwᵀ fnew   over the np active columns.
        for a in 0..np {
            for b in a..np {
                let mut s = 0.0;
                for i in 0..n {
                    s += fdiff[a][i] * fdiff[b][i];
                }
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
        jacobi_eig(
            &mut gram[..np * np],
            np,
            &mut eval[..np],
            &mut evec[..np * np],
        );
        for a in 0..np {
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

        if ss_prop <= ss_resids * (1.0 + RESID_TOL.powi(k)) {
            // Accept the (damped) Anderson step.
            std::mem::swap(&mut fold, &mut fnew);
            xnew.copy_from_slice(&x_prop);
            fnew.copy_from_slice(&f_prop);
            shrink_count += 1.0;
            ss_resids = ss_prop;
        } else {
            // Reject: take a plain fixed-point (EM) step from the previous xnew.
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

        if *it >= min_iter {
            let d = max_rel_diff(&xold, &xnew, opts.alpha_check_cutoff);
            if d.is_finite() && d < opts.rel_diff_tol {
                converged = true;
                break;
            }
        }

        // Periodic restart of the Anderson history.
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

    #[test]
    fn jacobi_recovers_known_spectrum() {
        // Symmetric 2x2 [[2,1],[1,2]] has eigenvalues 1 and 3.
        let mut a = vec![2.0, 1.0, 1.0, 2.0];
        let mut eval = vec![0.0; 2];
        let mut evec = vec![0.0; 4];
        jacobi_eig(&mut a, 2, &mut eval, &mut evec);
        let mut evs = eval.clone();
        evs.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert!((evs[0] - 1.0).abs() < 1e-9, "{evs:?}");
        assert!((evs[1] - 3.0).abs() < 1e-9, "{evs:?}");
        // Eigenvectors orthonormal: columns dot to ~0.
        let dot = evec[0] * evec[1] + evec[2] * evec[3];
        assert!(dot.abs() < 1e-9, "columns not orthogonal: {dot}");
    }
}
