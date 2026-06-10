//! Log-space arithmetic helpers.
//!
//! Salmon performs nearly all of its probability accumulation in log space to
//! avoid underflow. These helpers mirror the conventions in the C++ code where
//! `LOG_0` is negative infinity and `LOG_1` is zero.

/// log(0) = -inf
pub const LOG_0: f64 = f64::NEG_INFINITY;
/// log(1) = 0
pub const LOG_1: f64 = 0.0;
/// A very small log-space epsilon used as an effective "zero" mass in places
/// where strict -inf would propagate NaNs.
pub const LOG_EPSILON: f64 = -1.0e10;

/// Numerically stable `log(exp(x) + exp(y))`.
///
/// Handles the `-inf` identities so that `log_add(LOG_0, y) == y`.
#[inline]
pub fn log_add(x: f64, y: f64) -> f64 {
    if x == LOG_0 {
        return y;
    }
    if y == LOG_0 {
        return x;
    }
    let (hi, lo) = if x > y { (x, y) } else { (y, x) };
    // hi + log1p(exp(lo - hi)); lo - hi <= 0 so exp is in (0, 1].
    hi + (lo - hi).exp().ln_1p()
}

/// Numerically stable `log(exp(x) - exp(y))`, requires `x >= y`.
///
/// Returns `LOG_0` when `x == y`.
#[inline]
pub fn log_sub(x: f64, y: f64) -> f64 {
    debug_assert!(x >= y, "log_sub requires x >= y (x={x}, y={y})");
    if y == LOG_0 {
        return x;
    }
    if x == y {
        return LOG_0;
    }
    // x + log(1 - exp(y - x)); y - x < 0.
    x + (-((y - x).exp())).ln_1p()
}

/// Log-sum-exp over an iterator of log-space values.
pub fn log_sum<I: IntoIterator<Item = f64>>(values: I) -> f64 {
    values.into_iter().fold(LOG_0, log_add)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64) {
        if a == b {
            return; // covers the -inf == -inf case
        }
        assert!((a - b).abs() < 1e-9, "expected {b}, got {a}");
    }

    #[test]
    fn add_identities() {
        approx(log_add(LOG_0, 0.5), 0.5);
        approx(log_add(0.5, LOG_0), 0.5);
        // log(e^0 + e^0) = log(2)
        approx(log_add(LOG_1, LOG_1), 2.0_f64.ln());
    }

    #[test]
    fn add_matches_naive() {
        let (x, y) = (-3.2_f64, -1.1_f64);
        let naive = (x.exp() + y.exp()).ln();
        approx(log_add(x, y), naive);
    }

    #[test]
    fn sub_works() {
        let (x, y) = (2.0_f64, 1.0_f64);
        let naive = (x.exp() - y.exp()).ln();
        approx(log_sub(x, y), naive);
        approx(log_sub(x, LOG_0), x);
        approx(log_sub(x, x), LOG_0);
    }

    #[test]
    fn sum_works() {
        let vals = [LOG_1, LOG_1, LOG_1];
        approx(log_sum(vals), 3.0_f64.ln());
        approx(log_sum(std::iter::empty::<f64>()), LOG_0);
    }
}
