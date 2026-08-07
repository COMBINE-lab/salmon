//! Log-space arithmetic helpers.
//!
//! # Why log space
//!
//! Salmon multiplies together many probabilities: the chance a fragment has this
//! length, times the chance of this sequence context, times the chance of this
//! position, and so on, over hundreds of millions of fragments. Multiplying a
//! few dozen numbers below 1 quickly produces values smaller than the smallest
//! positive `f64` (~1e-308), at which point they collapse to exactly zero and
//! all information is lost. This is *underflow*.
//!
//! The standard fix is to store the logarithm of every probability instead.
//! Multiplication becomes addition, which cannot underflow, and the dynamic
//! range of `f64` becomes enormous (a log of -700 still represents a meaningful
//! 1e-304).
//!
//! The cost is that *adding* two probabilities is no longer trivial: given
//! `log(a)` and `log(b)` you want `log(a + b)`. Naively computing
//! `(x.exp() + y.exp()).ln()` reintroduces exactly the underflow we were
//! avoiding, so [`log_add`] uses the standard stable reformulation below.
//!
//! Salmon performs nearly all of its probability accumulation this way. These
//! helpers mirror the conventions in the C++ code where `LOG_0` is negative
//! infinity and `LOG_1` is zero.

/// log(0) = -inf
pub const LOG_0: f64 = f64::NEG_INFINITY;
/// log(1) = 0
pub const LOG_1: f64 = 0.0;
/// A very small log-space epsilon used as an effective "zero" mass in places
/// where strict -inf would propagate NaNs.
///
/// The trouble with true `-inf` is that `-inf - -inf` and `0 * -inf` are `NaN`,
/// and a single `NaN` poisons every arithmetic result it touches. A very large
/// negative *finite* number behaves like zero for every practical purpose while
/// keeping the arithmetic well-defined.
pub const LOG_EPSILON: f64 = -1.0e10;

/// Numerically stable `log(exp(x) + exp(y))`.
///
/// **How the trick works.** Factor out the larger term:
/// `log(e^hi + e^lo) = hi + log(1 + e^(lo-hi))`. Since `lo - hi <= 0`, the
/// exponential is in `(0, 1]` and cannot overflow, and the whole expression
/// stays accurate even when the two inputs differ by hundreds of orders of
/// magnitude.
///
/// Handles the `-inf` identities so that `log_add(LOG_0, y) == y`.
#[inline]
pub fn log_add(x: f64, y: f64) -> f64 {
    // Adding zero probability changes nothing — and short-circuiting also avoids
    // the `-inf - -inf = NaN` that the general path would hit.
    if x == LOG_0 {
        return y;
    }
    if y == LOG_0 {
        return x;
    }
    let (hi, lo) = if x > y { (x, y) } else { (y, x) };
    // hi + log1p(exp(lo - hi)); lo - hi <= 0 so exp is in (0, 1].
    // `ln_1p(v)` computes `log(1 + v)` accurately for tiny `v`, where computing
    // `1 + v` first would round away all the significant digits.
    hi + (lo - hi).exp().ln_1p()
}

/// Numerically stable `log(exp(x) - exp(y))`, requires `x >= y`.
///
/// Same factoring as [`log_add`], with a minus sign:
/// `log(e^x - e^y) = x + log(1 - e^(y-x))`.
///
/// The precondition is real: subtracting a larger probability from a smaller one
/// would give a negative number, which has no logarithm.
///
/// Returns `LOG_0` when `x == y`.
#[inline]
pub fn log_sub(x: f64, y: f64) -> f64 {
    // Checked in debug/test builds only; compiles away in release.
    debug_assert!(x >= y, "log_sub requires x >= y (x={x}, y={y})");
    // Subtracting zero probability changes nothing.
    if y == LOG_0 {
        return x;
    }
    // Equal probabilities cancel exactly; the general formula would give
    // log(1 - 1) = log(0), which is right but arrives via a rounding-sensitive
    // path.
    if x == y {
        return LOG_0;
    }
    // x + log(1 - exp(y - x)); y - x < 0.
    x + (-((y - x).exp())).ln_1p()
}

/// Log-sum-exp over an iterator of log-space values.
///
/// `fold` starts from `LOG_0` (probability zero, the identity for addition) and
/// folds each value in with [`log_add`], so an empty input correctly yields
/// `LOG_0` rather than an error.
pub fn log_sum<I: IntoIterator<Item = f64>>(values: I) -> f64 {
    values.into_iter().fold(LOG_0, log_add)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Compare with a tolerance, since these are floating-point results — but
    /// let exact equality through first so infinities compare cleanly.
    fn approx(a: f64, b: f64) {
        if a == b {
            return; // covers the -inf == -inf case
        }
        assert!((a - b).abs() < 1e-9, "expected {b}, got {a}");
    }

    /// Adding probability zero must be a no-op on either side, and adding two
    /// certainties must give exactly log 2.
    #[test]
    fn add_identities() {
        approx(log_add(LOG_0, 0.5), 0.5);
        approx(log_add(0.5, LOG_0), 0.5);
        // log(e^0 + e^0) = log(2)
        approx(log_add(LOG_1, LOG_1), 2.0_f64.ln());
    }

    /// Where the naive formula *is* safe (values near 1), the stable one must
    /// agree with it — this catches an algebra slip in the reformulation.
    #[test]
    fn add_matches_naive() {
        let (x, y) = (-3.2_f64, -1.1_f64);
        let naive = (x.exp() + y.exp()).ln();
        approx(log_add(x, y), naive);
    }

    /// Subtraction against the naive formula, plus both identity cases.
    #[test]
    fn sub_works() {
        let (x, y) = (2.0_f64, 1.0_f64);
        let naive = (x.exp() - y.exp()).ln();
        approx(log_sub(x, y), naive);
        approx(log_sub(x, LOG_0), x);
        approx(log_sub(x, x), LOG_0);
    }

    /// Folding must accumulate correctly and return the additive identity for an
    /// empty sequence.
    #[test]
    fn sum_works() {
        let vals = [LOG_1, LOG_1, LOG_1];
        approx(log_sum(vals), 3.0_f64.ln());
        approx(log_sum(std::iter::empty::<f64>()), LOG_0);
    }
}
