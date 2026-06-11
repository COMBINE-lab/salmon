//! Fragment-length distribution.
//!
//! Direct port of salmon's `FragmentLengthDistribution`
//! (`src/model/FragmentLengthDistribution.cpp`): a log-space histogram seeded
//! with a Gaussian (or uniform) prior, updated by adding a binomial smoothing
//! kernel around each observed length. All masses and probabilities are in log
//! space. Updates are lock-free so worker threads can call [`add_val`] with a
//! shared reference, matching the C++ design.
//!
//! [`add_val`]: FragmentLengthDistribution::add_val

use salmon_core::atomic::AtomicF64;
use salmon_core::math::{log_add, LOG_0, LOG_EPSILON};
use statrs::distribution::{Binomial, ContinuousCDF, Discrete, Normal};
use std::sync::atomic::{AtomicUsize, Ordering};

/// Tracks the observed distribution of fragment lengths.
#[derive(Debug)]
pub struct FragmentLengthDistribution {
    /// logged binomial smoothing kernel
    kernel: Vec<f64>,
    /// logged observed mass per length bin
    hist: Vec<AtomicF64>,
    /// logged total observed mass (including pseudo-counts)
    tot_mass: AtomicF64,
    /// logged sum of length*mass, for fast mean computation
    sum: AtomicF64,
    /// minimum observed length (bin units)
    min: AtomicUsize,
    /// internal bin size
    bin_size: usize,

    /// cached normalized PMF, valid once [`cache`](Self::cache) is called
    cached_pmf: Vec<f64>,
    /// cached CMF
    cached_cmf: Vec<f64>,
    have_cache: bool,
}

impl FragmentLengthDistribution {
    /// Construct a distribution.
    ///
    /// * `alpha` – total pseudo-count mass (linear space).
    /// * `max_val` – maximum representable length.
    /// * `prior_mu` – Gaussian prior mean; if `<= 0`, a uniform prior is used.
    /// * `prior_sigma` – Gaussian prior standard deviation.
    /// * `kernel_n` – binomial kernel trials; must be even (after binning).
    /// * `kernel_p` – binomial kernel success probability.
    /// * `bin_size` – internal length binning (use 1 for no binning).
    pub fn new(
        alpha: f64,
        max_val: usize,
        prior_mu: f64,
        prior_sigma: f64,
        kernel_n: usize,
        kernel_p: f64,
        bin_size: usize,
    ) -> Self {
        assert!(bin_size >= 1, "bin_size must be >= 1");
        let max_val = max_val / bin_size;
        let kernel_n = kernel_n / bin_size;
        assert!(kernel_n % 2 == 0, "kernel_n must be even after binning");

        let tot = alpha.ln();
        let hist: Vec<AtomicF64>;
        let mut sum = LOG_0;
        let mut tot_mass;

        if prior_mu > 0.0 {
            let norm = Normal::new(prior_mu / bin_size as f64, prior_sigma / (bin_size * bin_size) as f64)
                .expect("valid normal prior");
            hist = (0..=max_val).map(|_| AtomicF64::new(LOG_0)).collect();
            tot_mass = LOG_0;
            for (i, slot) in hist.iter().enumerate() {
                let norm_mass = norm.cdf(i as f64 + 0.5) - norm.cdf(i as f64 - 0.5);
                let mass = if norm_mass != 0.0 {
                    tot + norm_mass.ln()
                } else {
                    LOG_EPSILON
                };
                slot.store(mass);
                sum = log_add(sum, (i as f64).ln() + mass);
                tot_mass = log_add(tot_mass, mass);
            }
        } else {
            // uniform prior
            let per = tot - (max_val as f64).ln();
            hist = (0..=max_val).map(|_| AtomicF64::new(per)).collect();
            hist[0].store(LOG_0);
            let h1 = hist.get(1).map(|a| a.load()).unwrap_or(per);
            sum = h1 + ((max_val * (max_val + 1)) as f64).ln() - 2.0_f64.ln();
            tot_mass = tot;
        }

        // binomial smoothing kernel
        let binom = Binomial::new(kernel_p, kernel_n as u64).expect("valid binomial kernel");
        let kernel: Vec<f64> = (0..=kernel_n).map(|i| binom.pmf(i as u64).ln()).collect();

        Self {
            kernel,
            hist,
            tot_mass: AtomicF64::new(tot_mass),
            sum: AtomicF64::new(sum),
            min: AtomicUsize::new(max_val),
            bin_size,
            cached_pmf: Vec::new(),
            cached_cmf: Vec::new(),
            have_cache: false,
        }
    }

    /// salmon's default fragment-length distribution: pseudo-count 1.0, max
    /// length 1000, no Gaussian prior (uniform), kernel `n=4, p=0.5`.
    pub fn default_for_paired() -> Self {
        Self::new(1.0, 1000, 0.0, 0.0, 4, 0.5, 1)
    }

    pub fn max_val(&self) -> usize {
        (self.hist.len() - 1) * self.bin_size
    }

    pub fn min_val(&self) -> usize {
        let m = self.min.load(Ordering::Relaxed);
        if m == self.hist.len() - 1 {
            1
        } else {
            m
        }
    }

    /// Add `mass` (log space) for an observed fragment of length `len`,
    /// spreading it over the smoothing kernel. Lock-free; safe to call from
    /// multiple threads. (Must not race with [`cache`](Self::cache).)
    pub fn add_val(&self, len: usize, mass: f64) {
        let mut len = len / self.bin_size;
        let max_v = self.max_val() / self.bin_size;
        if len > max_v {
            len = max_v;
        }
        self.min.fetch_min(len, Ordering::Relaxed);

        let half = self.kernel.len() / 2;
        // offset can go negative conceptually; use isize math then bound-check.
        let mut offset = len as isize - half as isize;
        for &k in &self.kernel {
            if offset > 0 && (offset as usize) < self.hist.len() {
                let o = offset as usize;
                let k_mass = mass + k;
                self.hist[o].log_add_assign(k_mass);
                self.sum.log_add_assign((o as f64).ln() + k_mass);
                self.tot_mass.log_add_assign(k_mass);
            }
            offset += 1;
        }
    }

    /// Logged probability of observing a fragment of length `len`.
    pub fn pmf(&self, len: usize) -> f64 {
        if self.have_cache {
            return *self
                .cached_pmf
                .get(len)
                .unwrap_or_else(|| self.cached_pmf.last().unwrap());
        }
        let mut l = len / self.bin_size;
        let max_v = self.max_val() / self.bin_size;
        if l > max_v {
            l = max_v;
        }
        self.hist[l].load() - self.tot_mass.load()
    }

    /// Logged cumulative mass up to and including `len`.
    pub fn cmf(&self, len: usize) -> f64 {
        if self.have_cache {
            return *self
                .cached_cmf
                .get(len)
                .unwrap_or_else(|| self.cached_cmf.last().unwrap());
        }
        let mut l = len / self.bin_size;
        let max_v = self.max_val() / self.bin_size;
        if l > max_v {
            l = max_v;
        }
        let mut cum = LOG_0;
        for i in 0..=l {
            cum = log_add(cum, self.hist[i].load());
        }
        cum - self.tot_mass.load()
    }

    /// Total observed mass (log space).
    pub fn tot_mass(&self) -> f64 {
        self.tot_mass.load()
    }

    /// Mean observed length.
    pub fn mean(&self) -> f64 {
        (self.sum.load() - self.tot_mass.load()).exp()
    }

    /// Standard deviation of the observed length distribution, computed from the
    /// cached normalized PMF (call after [`cache`](Self::cache)).
    pub fn sd(&self) -> f64 {
        let lp = self.log_pmf();
        if lp.is_empty() {
            return 0.0;
        }
        let mut mean = 0.0;
        for (l, &p) in lp.iter().enumerate() {
            mean += (l as f64) * p.exp();
        }
        let mut var = 0.0;
        for (l, &p) in lp.iter().enumerate() {
            let d = l as f64 - mean;
            var += d * d * p.exp();
        }
        var.max(0.0).sqrt()
    }

    /// Freeze the distribution and precompute normalized PMF/CMF for fast,
    /// allocation-free lookup. Call once after updates have stopped.
    pub fn cache(&mut self) {
        if self.have_cache {
            return;
        }
        let max_v = self.max_val();
        // normalized PMF over [0, max_v]
        let mut pmf = Vec::with_capacity(max_v + 1);
        let mut tot = LOG_0;
        for i in 0..=max_v {
            let p = self.pmf(i);
            pmf.push(p);
            tot = log_add(tot, p);
        }
        for p in &mut pmf {
            *p -= tot;
        }
        // CMF from the normalized PMF
        let mut cmf = Vec::with_capacity(pmf.len());
        let mut cum = LOG_0;
        for &p in &pmf {
            cum = log_add(cum, p);
            cmf.push(cum);
        }
        self.cached_pmf = pmf;
        self.cached_cmf = cmf;
        self.have_cache = true;
    }

    /// The cached, normalized log-PMF over `[0, max_val]`. Requires [`cache`](Self::cache).
    pub fn log_pmf(&self) -> &[f64] {
        debug_assert!(self.have_cache, "call cache() before log_pmf()");
        &self.cached_pmf
    }

    /// Cumulative conditional means `E[L | L ≤ i]` over `[0, max_val]`, i.e.
    /// salmon's `correctionFactorsFromMass` (`DistributionUtils.cpp`):
    /// `cm[i] = (Σ_{l≤i} l·pmf[l]) / (Σ_{l≤i} pmf[l])`.
    ///
    /// These are the per-length correction factors `computeSmoothedEffectiveLengths`
    /// subtracts from the reference length to get the base effective length. The
    /// ratio is invariant to the PMF normalization, so the cached (normalized) PMF
    /// gives the same values as salmon's `100·exp(logPMF)` mass. Requires
    /// [`cache`](Self::cache).
    pub fn conditional_means(&self) -> Vec<f64> {
        debug_assert!(self.have_cache, "call cache() before conditional_means()");
        let n = self.cached_pmf.len();
        let mut cms = vec![0.0f64; n];
        let mut vals = 0.0; // Σ l·pmf[l]
        let mut mult = 0.0; // Σ pmf[l]
        for i in 0..n {
            let p = self.cached_pmf[i].exp();
            vals += (i as f64) * p;
            mult += p;
            cms[i] = if mult > 0.0 { vals / mult } else { 0.0 };
        }
        cms
    }
}

/// salmon's base effective length (`computeSmoothedEffectiveLengths`):
/// `effLen = refLen − E[L | L ≤ refLen]`, clamped back to `refLen` if it would
/// fall below 1. `cond_means` is [`FragmentLengthDistribution::conditional_means`].
///
/// This replaces the truncated-PMF `Σ pmf(l)·(refLen−l+1)` estimate (which falls
/// back to the raw `refLen` for any transcript shorter than the FLD mean), matching
/// salmon's behaviour exactly.
pub fn smoothed_effective_length(cond_means: &[f64], ref_len: usize) -> f64 {
    if cond_means.is_empty() {
        return ref_len as f64;
    }
    let max_len = cond_means.len();
    let cf = if ref_len >= max_len {
        cond_means[max_len - 1]
    } else {
        cond_means[ref_len]
    };
    let eff = ref_len as f64 - cf;
    if eff < 1.0 {
        ref_len as f64
    } else {
        eff
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_prior_pmf_normalizes() {
        let mut fld = FragmentLengthDistribution::new(1.0, 200, 0.0, 0.0, 4, 0.5, 1);
        fld.cache();
        let total: f64 = fld.log_pmf().iter().map(|p| p.exp()).sum();
        assert!((total - 1.0).abs() < 1e-9, "pmf sums to {total}");
    }

    #[test]
    fn gaussian_prior_mean_is_near_mu() {
        let fld = FragmentLengthDistribution::new(1000.0, 1000, 250.0, 25.0, 4, 0.5, 1);
        let m = fld.mean();
        assert!((m - 250.0).abs() < 5.0, "mean {m} not near 250");
    }

    #[test]
    fn observations_shift_the_distribution() {
        let mut fld = FragmentLengthDistribution::new(1.0, 1000, 250.0, 25.0, 4, 0.5, 1);
        // pile observations around 400
        for _ in 0..100_000 {
            fld.add_val(400, 0.0); // mass = log(1) = 0
        }
        let m = fld.mean();
        assert!(m > 300.0, "mean {m} did not move toward 400");
        fld.cache();
        // length 400 should be among the most probable
        let p400 = fld.pmf(400);
        let p250 = fld.pmf(250);
        assert!(p400 > p250, "p(400)={p400} not > p(250)={p250}");
    }

    #[test]
    fn smoothed_efflen_shrinks_short_transcripts() {
        // Gaussian prior mean 250: a transcript far shorter than the mean should
        // get a heavily shrunk effective length (NOT the raw refLen the old
        // truncated-PMF estimate fell back to).
        let mut fld = FragmentLengthDistribution::new(1000.0, 1000, 250.0, 25.0, 4, 0.5, 1);
        fld.cache();
        let cm = fld.conditional_means();
        // conditional means are non-decreasing
        for w in cm.windows(2) {
            assert!(w[1] >= w[0] - 1e-9, "cond means not monotonic: {} < {}", w[1], w[0]);
        }
        let short = smoothed_effective_length(&cm, 201);
        assert!(short < 201.0 && short > 1.0, "short effLen {short} not shrunk");
        // a long transcript keeps most of its length
        let long = smoothed_effective_length(&cm, 5000);
        assert!(long > 4000.0, "long effLen {long} shrunk too much");
        // below the 1.0 barrier the raw length is returned
        let tiny = smoothed_effective_length(&cm, 2);
        assert_eq!(tiny, 2.0, "tiny transcript should fall back to refLen");
    }

    #[test]
    fn cmf_is_monotonic() {
        let mut fld = FragmentLengthDistribution::new(1000.0, 500, 200.0, 30.0, 4, 0.5, 1);
        fld.cache();
        let mut prev = f64::NEG_INFINITY;
        for l in 0..=500 {
            let c = fld.cmf(l);
            assert!(c >= prev - 1e-9, "cmf decreased at {l}: {c} < {prev}");
            prev = c;
        }
        assert!((prev - 0.0).abs() < 1e-6, "cmf endpoint {prev} != log(1)");
    }
}
