//! Cubic spline interpolation — a faithful port of the vendored `tk::spline`
//! (Tino Kluge) used by salmon's `SimplePosBias`.
//!
//! Matches salmon's default construction: **cubic** spline with **natural**
//! boundary conditions (second derivative 0 at both ends) and **quadratic**
//! extrapolation outside the knot range. The band-matrix LU solve is ported
//! verbatim (tridiagonal, one upper + one lower band) so the spline coefficients
//! match the C++ implementation.

/// A tridiagonal band matrix with one upper and one lower band (plus a saved
/// diagonal used during LU). Mirrors `tk::band_matrix` for the `n_u = n_l = 1`
/// case that cubic-spline construction needs.
struct BandMatrix {
    dim: usize,
    // m_upper[0] = main diagonal, m_upper[1] = super-diagonal
    upper: [Vec<f64>; 2],
    // m_lower[0] = saved diagonal (set during LU), m_lower[1] = sub-diagonal
    lower: [Vec<f64>; 2],
}

impl BandMatrix {
    fn new(dim: usize) -> Self {
        Self {
            dim,
            upper: [vec![0.0; dim], vec![0.0; dim]],
            lower: [vec![0.0; dim], vec![0.0; dim]],
        }
    }

    /// `A(i, j)` accessor (band index `k = j - i`; `k >= 0` -> upper, else lower).
    #[inline]
    fn get(&self, i: usize, j: usize) -> f64 {
        let k = j as isize - i as isize;
        if k >= 0 {
            self.upper[k as usize][i]
        } else {
            self.lower[(-k) as usize][i]
        }
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, v: f64) {
        let k = j as isize - i as isize;
        if k >= 0 {
            self.upper[k as usize][i] = v;
        } else {
            self.lower[(-k) as usize][i] = v;
        }
    }

    #[inline]
    fn saved_diag(&self, i: usize) -> f64 {
        self.lower[0][i]
    }
    #[inline]
    fn set_saved_diag(&mut self, i: usize, v: f64) {
        self.lower[0][i] = v;
    }

    /// LR-decomposition of the band matrix (ported from `tk::band_matrix`).
    fn lu_decompose(&mut self) {
        let dim = self.dim as isize;
        // preconditioning: normalize row i so a_ii = 1
        for i in 0..self.dim {
            let diag = self.get(i, i);
            self.set_saved_diag(i, 1.0 / diag);
            let j_min = (i as isize - 1).max(0) as usize;
            let j_max = (i + 1).min(self.dim - 1);
            let s = self.saved_diag(i);
            for j in j_min..=j_max {
                self.set(i, j, self.get(i, j) * s);
            }
            self.set(i, i, 1.0); // prevents rounding errors
        }
        // Gauss LR-decomposition
        for k in 0..self.dim {
            let i_max = ((k + 1).min(self.dim - 1)) as isize;
            let mut i = k as isize + 1;
            while i <= i_max {
                let iu = i as usize;
                let akk = self.get(k, k);
                let x = -self.get(iu, k) / akk;
                self.set(iu, k, -x); // assembly part of L
                let j_max = ((k + 1).min(self.dim - 1)) as isize;
                let mut j = k as isize + 1;
                while j <= j_max {
                    let ju = j as usize;
                    self.set(iu, ju, self.get(iu, ju) + x * self.get(k, ju));
                    j += 1;
                }
                i += 1;
            }
        }
        let _ = dim;
    }

    /// Solve `Ly = b`.
    fn l_solve(&self, b: &[f64]) -> Vec<f64> {
        let mut x = vec![0.0; self.dim];
        for i in 0..self.dim {
            let mut sum = 0.0;
            let j_start = (i as isize - 1).max(0) as usize;
            for j in j_start..i {
                sum += self.get(i, j) * x[j];
            }
            x[i] = b[i] * self.saved_diag(i) - sum;
        }
        x
    }

    /// Solve `Rx = y`.
    fn r_solve(&self, b: &[f64]) -> Vec<f64> {
        let mut x = vec![0.0; self.dim];
        for i in (0..self.dim).rev() {
            let mut sum = 0.0;
            let j_stop = (i + 1).min(self.dim - 1);
            for j in (i + 1)..=j_stop {
                if j > i {
                    sum += self.get(i, j) * x[j];
                }
            }
            x[i] = (b[i] - sum) / self.get(i, i);
        }
        x
    }

    fn lu_solve(&mut self, b: &[f64]) -> Vec<f64> {
        self.lu_decompose();
        let y = self.l_solve(b);
        self.r_solve(&y)
    }
}

/// A cubic spline `f(x) = a·(x−x_i)³ + b·(x−x_i)² + c·(x−x_i) + y_i`.
#[derive(Debug, Clone, Default)]
pub struct Spline {
    x: Vec<f64>,
    y: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    b0: f64,
    c0: f64,
}

impl Spline {
    /// Build a cubic spline with natural boundary conditions and quadratic
    /// extrapolation — salmon's default `tk::spline(xs, ys)`.
    pub fn new(xs: Vec<f64>, ys: Vec<f64>) -> Self {
        let n = xs.len();
        assert!(n >= 3, "cubic spline needs >= 3 points");
        let x = xs;
        let y = ys;

        let mut mat = BandMatrix::new(n);
        let mut rhs = vec![0.0; n];
        for i in 1..n - 1 {
            mat.set(i, i - 1, (x[i] - x[i - 1]) / 3.0);
            mat.set(i, i, 2.0 / 3.0 * (x[i + 1] - x[i - 1]));
            mat.set(i, i + 1, (x[i + 1] - x[i]) / 3.0);
            rhs[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
        }
        // natural (second_deriv = 0) boundary conditions
        mat.set(0, 0, 2.0);
        mat.set(0, 1, 0.0);
        rhs[0] = 0.0;
        mat.set(n - 1, n - 1, 2.0);
        mat.set(n - 1, n - 2, 0.0);
        rhs[n - 1] = 0.0;

        let b = mat.lu_solve(&rhs);

        let mut a = vec![0.0; n];
        let mut c = vec![0.0; n];
        for i in 0..n - 1 {
            a[i] = (b[i + 1] - b[i]) / (x[i + 1] - x[i]) / 3.0;
            c[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
                - (2.0 * b[i] + b[i + 1]) * (x[i + 1] - x[i]) / 3.0;
        }

        // quadratic-extrapolation coefficients (default)
        let b0 = 0.0;
        let c0 = c[0];
        let h = x[n - 1] - x[n - 2];
        a[n - 1] = 0.0;
        c[n - 1] = 3.0 * a[n - 2] * h * h + 2.0 * b[n - 2] * h + c[n - 2];

        Self {
            x,
            y,
            a,
            b,
            c,
            b0,
            c0,
        }
    }

    /// `m_x[idx] <= x`, with `idx = 0` even when `x < m_x[0]`.
    #[inline]
    fn closest_idx_to(&self, x: f64) -> usize {
        // lower_bound: first element >= x
        let it = self.x.partition_point(|&v| v < x);
        if it == 0 {
            0
        } else {
            it - 1
        }
    }

    /// Evaluate the spline at `x` (interpolating, or extrapolating quadratically
    /// outside `[x_0, x_{n-1}]`).
    pub fn eval(&self, x: f64) -> f64 {
        let n = self.x.len();
        let idx = self.closest_idx_to(x);
        let h = x - self.x[idx];
        if x < self.x[0] {
            (self.b0 * h + self.c0) * h + self.y[0]
        } else if x > self.x[n - 1] {
            (self.b[n - 1] * h + self.c[n - 1]) * h + self.y[n - 1]
        } else {
            ((self.a[idx] * h + self.b[idx]) * h + self.c[idx]) * h + self.y[idx]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolates_knots_exactly() {
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = vec![0.0, 1.0, 4.0, 9.0, 16.0]; // ~ x^2
        let s = Spline::new(xs.clone(), ys.clone());
        for (x, y) in xs.iter().zip(&ys) {
            assert!(
                (s.eval(*x) - y).abs() < 1e-9,
                "knot {x}: {} != {y}",
                s.eval(*x)
            );
        }
    }

    #[test]
    fn monotone_line_is_reproduced() {
        // a straight line should be reproduced exactly by a natural cubic spline
        let xs = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let ys: Vec<f64> = xs.iter().map(|x| 2.0 * x + 1.0).collect();
        let s = Spline::new(xs, ys);
        for i in 0..=10 {
            let x = i as f64 / 10.0;
            assert!((s.eval(x) - (2.0 * x + 1.0)).abs() < 1e-9, "x={x}");
        }
    }

    #[test]
    fn extrapolates_without_panic() {
        let xs = vec![0.0, 0.5, 1.0];
        let ys = vec![1.0, 2.0, 1.0];
        let s = Spline::new(xs, ys);
        let _ = s.eval(-0.3);
        let _ = s.eval(1.7);
    }
}
