//! Positional fragment bias (`SimplePosBias`) — a port of salmon's
//! `src/model/SimplePosBias.cpp`.
//!
//! A 20-bin histogram of fragment start positions (binned uniformly along the
//! transcript) accumulated in log space, one model per *length class* and per
//! end (5'/3'). After [`finalize`](SimplePosBias::finalize) the (normalized)
//! bin masses are fit with a cubic spline over the non-uniform Cufflinks
//! position knots, and [`project_weights`](SimplePosBias::project_weights)
//! evaluates that spline at each transcript position. The observed/expected
//! per-position ratio gives the positional bias factor used in the
//! effective-length convolution.

use crate::spline::Spline;

/// salmon's default number of position bins (`SimplePosBias` `numBins`).
pub const NUM_POS_BINS: usize = 20;
/// salmon's number of transcript length classes (`numBiasBins`).
pub const NUM_LENGTH_CLASSES: usize = 5;

/// Cufflinks fractional position knots used as the spline x-coordinates
/// (`SimplePosBias::positionBins_`).
const POSITION_BINS: [f64; NUM_POS_BINS] = [
    0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.92, 0.94,
    0.96, 0.98, 1.0,
];

/// A single positional-bias distribution (20 bins).
#[derive(Debug, Clone)]
pub struct SimplePosBias {
    num_bins: usize,
    /// bin masses; log space until [`finalize`](Self::finalize), then linear.
    /// Materialized from `masses_fp` at finalize.
    masses: Vec<f64>,
    /// Fixed-point integer accumulator of LINEAR per-bin mass (`mass *
    /// BIAS_WEIGHT_SCALE`, truncated); summed order-independently across worker
    /// threads and materialized into the log-space `masses` (with the single
    /// `log(1)` pseudocount injected once) at [`finalize`]. Replaces the old
    /// non-associative log-space `log_add`. See [`crate::BIAS_WEIGHT_SCALE`].
    masses_fp: Vec<u64>,
    finalized: bool,
    spline: Option<Spline>,
}

impl Default for SimplePosBias {
    fn default() -> Self {
        Self::new(NUM_POS_BINS)
    }
}

impl SimplePosBias {
    /// A model with `num_bins` bins, masses initialized to `log(1) = 0` (salmon
    /// initializes to `LOG_1`, i.e. a pseudocount of 1 per bin).
    pub fn new(num_bins: usize) -> Self {
        Self {
            num_bins,
            masses: vec![0.0; num_bins],
            masses_fp: vec![0u64; num_bins],
            finalized: false,
            spline: None,
        }
    }

    /// A model with masses at `log(0) = -inf` — the true [`log_add`] identity, so
    /// it carries *no* pseudocount. Used as the per-thread accumulator when
    /// building an expected model in parallel: partials accumulate pure observed
    /// log-mass and merge associatively via [`combine`](Self::combine); the single
    /// `log(1)` pseudocount is then injected once by combining into a [`new`](Self::new)
    /// (or [`default`](Self::default)) model.
    pub fn new_empty(num_bins: usize) -> Self {
        Self {
            num_bins,
            masses: vec![f64::NEG_INFINITY; num_bins],
            masses_fp: vec![0u64; num_bins],
            finalized: false,
            spline: None,
        }
    }

    /// Add `mass` (LINEAR, e.g. a `[0,1]` posterior or an expected `weight *
    /// density`) to bin `bin`. Accumulated in the fixed-point integer buffer so
    /// the merge across threads is order-independent.
    #[inline]
    pub fn add_mass_bin(&mut self, bin: usize, mass: f64) {
        self.masses_fp[bin] += crate::bias_mass_to_fp(mass);
    }

    /// Add `mass` (LINEAR) to the bin for `pos` on a transcript of `length`
    /// (`bin = floor(pos / (length / numBins))`, clamped to the last bin).
    #[inline]
    pub fn add_mass(&mut self, pos: i32, length: i32, mass: f64) {
        debug_assert!(!self.finalized);
        let step = length as f64 / self.num_bins as f64;
        let mut bin = (pos as f64 / step).floor() as i32;
        if bin < 0 {
            bin = 0;
        }
        if bin as usize >= self.num_bins {
            bin = self.num_bins as i32 - 1;
        }
        self.add_mass_bin(bin as usize, mass);
    }

    /// Merge another (un-finalized) model's masses (log-space sum).
    pub fn combine(&mut self, other: &SimplePosBias) {
        debug_assert!(!self.finalized && !other.finalized);
        // Integer sum of linear masses — associative, so merging is independent
        // of the worker-thread fragment partition. The `log(1)` pseudocount is
        // injected once at `finalize`, so neither `new` nor `new_empty` carries a
        // per-partial prior (both start `masses_fp` at 0).
        for (a, b) in self.masses_fp.iter_mut().zip(&other.masses_fp) {
            *a += *b;
        }
    }

    /// Convert log masses to linear, normalize, and fit the spline over the
    /// Cufflinks position knots (duplicating the endpoints as natural knots).
    pub fn finalize(&mut self) {
        if self.finalized {
            return;
        }
        // Materialize the integer linear accumulator into the log-space `masses`,
        // injecting the single `log(1)` pseudocount once: `log(1 + Σ mass)`. This
        // matches the old serial `log_add` from a `log(1)`-seeded bin, but is
        // thread-count independent (the data sum is exact-integer).
        for (m, &fp) in self.masses.iter_mut().zip(&self.masses_fp) {
            *m = (1.0 + fp as f64 / crate::BIAS_WEIGHT_SCALE).ln();
        }
        let mut sum = 0.0;
        for m in &mut self.masses {
            *m = m.exp();
            sum += *m;
        }
        let n = self.masses.len();
        let start_knot = self.masses[0] / sum;
        let stop_knot = self.masses[n - 1] / sum;
        let spline_sum = sum + start_knot + stop_knot;

        let mut spline_mass = vec![0.0f64; n + 2];
        spline_mass[0] = start_knot;
        for i in 0..n {
            spline_mass[i + 1] = self.masses[i] / spline_sum;
            self.masses[i] /= sum;
        }
        spline_mass[n + 1] = stop_knot;

        let mut spline_bins = vec![0.0f64; n + 2];
        spline_bins[0] = 0.0;
        for i in 0..n {
            spline_bins[i + 1] = POSITION_BINS[i] - 0.01;
        }
        spline_bins[n + 1] = 1.0;

        self.spline = Some(Spline::new(spline_bins, spline_mass));
        self.finalized = true;
    }

    /// The (normalized, post-[`finalize`](Self::finalize)) bin masses — the same
    /// values salmon's `writeBinary` emits to `aux_info/*_pos.gz`.
    pub fn masses(&self) -> &[f64] {
        &self.masses
    }

    /// Project the spline weights onto `out.len()` transcript positions:
    /// `out[p] = max(0.0, spline(p / len))`. We clamp only negatives (cubic
    /// overshoot), NOT at a hard `0.001` floor: the floor turned the near-zero
    /// tails of the expected density into an artificially small divisor in the
    /// observed/expected ratio, amplifying tiny model differences. The ratio is
    /// instead additively smoothed toward 1 at the division site (see
    /// `bias::positional_factor`).
    pub fn project_weights(&self, out: &mut [f64]) {
        let spline = self
            .spline
            .as_ref()
            .expect("project_weights requires finalize()");
        let len = out.len();
        for (p, o) in out.iter_mut().enumerate() {
            let frac_p = p as f64 / len as f64;
            *o = spline.eval(frac_p).max(0.0);
        }
    }
}

/// Length-class quantile boundaries (salmon's `setTranscriptLengthClasses_`):
/// with `n > nbins`, the lengths at sorted positions `step, 2·step, …` (clamped);
/// otherwise the sorted lengths themselves.
pub fn compute_length_quantiles(lengths: &[u32], nbins: usize) -> Vec<u32> {
    let n = lengths.len();
    let mut sorted: Vec<u32> = lengths.to_vec();
    sorted.sort_unstable();
    if n > nbins {
        let step = n / nbins;
        let mut cum = 0usize;
        let mut q = Vec::with_capacity(nbins);
        for _ in 0..nbins {
            cum += step;
            let ind = cum.min(n - 1);
            q.push(sorted[ind]);
        }
        q
    } else {
        sorted
    }
}

/// The length-class index for a transcript (salmon's
/// `min(maxQuant, upper_bound(quantiles, refLen) - begin)`).
#[inline]
pub fn length_class_index(quantiles: &[u32], ref_len: u32) -> usize {
    let max_quant = quantiles.len() - 1;
    // upper_bound: first index whose quantile is > ref_len
    let ub = quantiles.partition_point(|&q| q <= ref_len);
    ub.min(max_quant)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_observations_project_near_flat() {
        // mass spread uniformly across positions -> projected weights ~ constant
        let mut m = SimplePosBias::new(NUM_POS_BINS);
        for pos in 0..1000 {
            m.add_mass(pos, 1000, 0.0); // log(1) per observation
        }
        m.finalize();
        let mut out = vec![0.0; 500];
        m.project_weights(&mut out);
        // all weights positive and finite; spread should be modest for uniform input
        let mean = out.iter().sum::<f64>() / out.len() as f64;
        assert!(out.iter().all(|&w| w >= 0.001 && w.is_finite()));
        let maxdev = out.iter().map(|w| (w - mean).abs()).fold(0.0, f64::max);
        assert!(
            maxdev / mean < 0.5,
            "uniform input not near-flat: maxdev/mean={}",
            maxdev / mean
        );
    }

    #[test]
    fn enriched_5p_has_higher_early_weight() {
        let mut m = SimplePosBias::new(NUM_POS_BINS);
        // pile mass at the 5' end
        for _ in 0..10000 {
            m.add_mass(5, 1000, 0.0);
        }
        for pos in 0..1000 {
            m.add_mass(pos, 1000, 0.0);
        }
        m.finalize();
        let mut out = vec![0.0; 1000];
        m.project_weights(&mut out);
        // early-position weight should exceed late-position weight
        assert!(
            out[10] > out[900],
            "5'-enriched: {} !> {}",
            out[10],
            out[900]
        );
    }

    #[test]
    fn length_classes_partition() {
        let lengths: Vec<u32> = (1..=1000).collect();
        let q = compute_length_quantiles(&lengths, 5);
        assert_eq!(q.len(), 5);
        // monotonic non-decreasing
        for w in q.windows(2) {
            assert!(w[1] >= w[0]);
        }
        // a tiny transcript -> class 0; a huge one -> last class
        assert_eq!(length_class_index(&q, 1), 0);
        assert_eq!(length_class_index(&q, 100000), 4);
    }
}
