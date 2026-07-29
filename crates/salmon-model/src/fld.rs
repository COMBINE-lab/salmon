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
use salmon_core::{LibraryFormat, ReadOrientation, ReadStrandedness, ReadType};
use statrs::distribution::{Binomial, ContinuousCDF, Discrete, Normal};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, RwLock};

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

    /// Periodically-refreshed snapshot of the (un-normalized) log-PMF used during
    /// the *online* phase, indexed by raw length. Reading the live `hist`/`tot_mass`
    /// directly (two separate atomic loads on a concurrently-updated distribution)
    /// returns slightly different values for the same length across calls, which
    /// breaks the weight symmetry of exact-duplicate transcripts and is then
    /// amplified by the VBEM `α<1` prior. Mirroring C++ salmon (`cachedPMF_` +
    /// `LogCMFCache`), worker threads instead capture an immutable snapshot of this
    /// once per fragment ([`online_snapshot`](Self::online_snapshot)) so every
    /// transcript of a given length in that fragment gets an identical value;
    /// [`refresh_online`](Self::refresh_online) rebuilds it at mini-batch
    /// boundaries.
    online_pmf: RwLock<Arc<Vec<f64>>>,

    /// Periodically-refreshed snapshot of the (normalized) log-CMF, the
    /// cumulative companion to [`online_pmf`](Self::online_pmf). Used for the
    /// ambiguous (orphan / single-end) fragment-length probability and for the
    /// `pmf(flen) − cmf(txpLen)` length-conditioning of proper pairs, both of
    /// which need cumulative mass. Rebuilt alongside `online_pmf` in
    /// [`refresh_online`](Self::refresh_online); mirrors C++ salmon's
    /// `LogCMFCache`.
    online_cmf: RwLock<Arc<Vec<f64>>>,
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
        assert!(
            kernel_n.is_multiple_of(2),
            "kernel_n must be even after binning"
        );

        let tot = alpha.ln();
        let hist: Vec<AtomicF64>;
        let mut sum = LOG_0;
        let mut tot_mass;

        if prior_mu > 0.0 {
            let norm = Normal::new(
                prior_mu / bin_size as f64,
                prior_sigma / (bin_size * bin_size) as f64,
            )
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
            online_pmf: RwLock::new(Arc::new(Vec::new())),
            online_cmf: RwLock::new(Arc::new(Vec::new())),
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

    /// Rebuild the online log-PMF snapshot from the current histogram (one pass
    /// over the length bins, with a single `tot_mass` read so the snapshot is
    /// internally consistent). Call at mini-batch boundaries during the online
    /// phase; no-op once the final [`cache`](Self::cache) has been taken. Cheap
    /// relative to mapping a batch, and decouples per-fragment reads from the
    /// concurrent `add_val` writes so identical lengths read identical values.
    pub fn refresh_online(&self) {
        if self.have_cache {
            return;
        }
        let max_raw = self.max_val();
        let max_v = max_raw / self.bin_size;
        let tot = self.tot_mass.load();
        // Per-bin cumulative mass (matches `cmf()`), so the snapshot CMF at raw
        // index `raw` equals `cmf(raw)`. Built first, then both the PMF and CMF
        // snapshots are expanded over raw indices from the same `tot` read so
        // they are mutually consistent.
        let mut bin_cum = Vec::with_capacity(max_v + 1);
        let mut cum = LOG_0;
        for b in 0..=max_v {
            cum = log_add(cum, self.hist[b].load() - tot);
            bin_cum.push(cum);
        }
        let mut v = Vec::with_capacity(max_raw + 1);
        let mut c = Vec::with_capacity(max_raw + 1);
        for raw in 0..=max_raw {
            let l = (raw / self.bin_size).min(max_v);
            v.push(self.hist[l].load() - tot);
            c.push(bin_cum[l]);
        }
        *self.online_pmf.write().unwrap() = Arc::new(v);
        *self.online_cmf.write().unwrap() = Arc::new(c);
    }

    /// Cheap (one `Arc` clone) immutable handle to the current online log-PMF
    /// snapshot. Capture once per fragment and index by raw length: every
    /// transcript of a given length then reads an identical value even if another
    /// thread refreshes the shared snapshot meanwhile. Empty until the first
    /// [`refresh_online`](Self::refresh_online) (the pre-burn-in window, where this
    /// term is not folded into the eq-class weight anyway).
    pub fn online_snapshot(&self) -> Arc<Vec<f64>> {
        self.online_pmf.read().unwrap().clone()
    }

    /// Cheap (one `Arc` clone) immutable handle to the current online log-CMF
    /// snapshot, the cumulative companion to [`online_snapshot`](Self::online_snapshot).
    /// Capture once per fragment for the ambiguous (orphan / single-end)
    /// fragment-length probability and the proper-pair length-conditioning.
    /// Empty until the first [`refresh_online`](Self::refresh_online).
    pub fn online_cmf_snapshot(&self) -> Arc<Vec<f64>> {
        self.online_cmf.read().unwrap().clone()
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

    /// Reconstruct a *cached* distribution directly from a (log-space) PMF,
    /// e.g. one serialized into a RAD header during a previous run. `log_pmf` is
    /// indexed by raw length over `[0, log_pmf.len())`. The CMF, conditional
    /// means, mean and sd are all re-derived from it via [`cache`](Self::cache),
    /// so a reconstructed distribution is interchangeable with the original for
    /// every read-side use. The masses need not be pre-normalized — `cache`
    /// normalizes them — but a normalized PMF round-trips exactly.
    pub fn from_log_pmf(log_pmf: &[f64]) -> Self {
        let max_val = log_pmf.len().saturating_sub(1);
        let mut d = Self::new(1.0, max_val, 0.0, 1.0, 4, 0.5, 1);
        // Replace the prior histogram with the supplied masses and recompute the
        // aggregate statistics (so `mean`/`sd` are consistent), then cache.
        let mut tot = LOG_0;
        let mut sm = LOG_0;
        for (i, &p) in log_pmf.iter().enumerate() {
            d.hist[i].store(p);
            tot = log_add(tot, p);
            if i > 0 {
                sm = log_add(sm, (i as f64).ln() + p);
            }
        }
        d.tot_mass.store(tot);
        d.sum.store(sm);
        d.min.store(0, Ordering::Relaxed);
        d.cache();
        d
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

/// Index a length into a (log) CMF snapshot, clamping out-of-range lengths to
/// the last bin (which holds the total mass). Returns [`LOG_0`] for an empty
/// snapshot.
#[inline]
fn cmf_at(cmf: &[f64], len: i32) -> f64 {
    if cmf.is_empty() {
        return LOG_0;
    }
    let i = (len.max(0) as usize).min(cmf.len() - 1);
    cmf[i]
}

/// Logged ambiguous-fragment-length probability for an orphan / single-end
/// read, given a (log) CMF snapshot. Direct port of C++ salmon's
/// `LogCMFCache::getAmbigFragLengthProb` (`DistributionUtils.cpp`).
///
/// The mapped mate bounds the maximum possible fragment length: a forward read
/// at `pos` can extend downstream to the transcript 3' end (`txp_len − pos`); a
/// reverse read's outer (5') end sits at `pos + read_len`, bounding the upstream
/// extent toward the 5' end. The weight is the FLD mass up to that bound,
/// *conditioned* on the mass up to the full transcript length — i.e.
/// `cmf(maxFragLen) − cmf(txpLen)` — so orphan weights sit on the same
/// length-conditioned scale as proper pairs. Returns [`LOG_EPSILON`] when the
/// transcript admits no representable fragment mass, and `LOG_1` (= 0) when no
/// snapshot is available yet (pre-burn-in), leaving the weight unmodelled.
pub fn ambig_frag_log_prob(cmf: &[f64], fwd: bool, pos: i32, read_len: i32, txp_len: i32) -> f64 {
    if cmf.is_empty() {
        return 0.0; // LOG_1: no model yet
    }
    let stxp = txp_len.max(0);
    let max_frag_len = if fwd {
        stxp - pos.clamp(0, stxp)
    } else {
        (pos + read_len).clamp(0, stxp)
    };
    let ref_cm = cmf_at(cmf, stxp);
    if ref_cm <= LOG_0 {
        return LOG_EPSILON;
    }
    cmf_at(cmf, max_frag_len) - ref_cm
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

/// Order-independent accumulator for *deriving* a fragment-length distribution
/// (and library format) from uniquely-mapped proper pairs.
///
/// Unlike [`FragmentLengthDistribution`] — which accumulates in **log space** and
/// is required for the online phase (it folds in per-fragment forgetting mass and
/// is read per fragment) — this stores plain **integer counts** per length,
/// bucketed by orientation. Integer increments are commutative, so the tallies
/// are independent of thread / chunk order; the FLD is then built **once**,
/// deterministically (fixed length order), in [`finish`](Self::finish). Keeping
/// it a separate type means the online FLD's float/log-space accumulation is
/// never disturbed, and no runtime dispatch is needed.
#[derive(Debug)]
pub struct DiscreteFld {
    /// per-length counts for opposite-strand (inward/outward) proper pairs
    opp: Vec<AtomicU64>,
    /// per-length counts for same-strand proper pairs
    same: Vec<AtomicU64>,
    n_opp: AtomicU64,
    n_same: AtomicU64,
    /// observed-format tally (indexed by [`LibraryFormat::format_id`]) for
    /// order-independent `-l A` auto-detection
    fmt_counts: [AtomicU64; 12],
    max_len: usize,
}

impl DiscreteFld {
    /// Create an accumulator covering raw fragment lengths `[0, fld_max]`.
    pub fn new(fld_max: usize) -> Self {
        Self {
            opp: (0..=fld_max).map(|_| AtomicU64::new(0)).collect(),
            same: (0..=fld_max).map(|_| AtomicU64::new(0)).collect(),
            n_opp: AtomicU64::new(0),
            n_same: AtomicU64::new(0),
            fmt_counts: std::array::from_fn(|_| AtomicU64::new(0)),
            max_len: fld_max,
        }
    }

    /// Record one uniquely-mapped proper pair from its mate strands: fragment
    /// length, orientation bucket, and observed format (derived from `is_fw` /
    /// `mate_fw`, mirroring the RAD reader's `rad_frag_format` so the mapping pass
    /// and the RAD-derive path agree). Works in sketch mode, where no precomputed
    /// `LibraryFormat` is available. Thread-safe and order-independent.
    pub fn add(&self, len: usize, is_fw: bool, mate_fw: bool) {
        let l = len.min(self.max_len);
        // opposite-strand (inward/outward) vs same-strand, exactly as the RAD
        // reader classifies placements.
        let (orientation, strandedness) = if is_fw != mate_fw {
            let s = if is_fw {
                ReadStrandedness::SA
            } else {
                ReadStrandedness::AS
            };
            (ReadOrientation::Toward, s)
        } else {
            let s = if is_fw {
                ReadStrandedness::S
            } else {
                ReadStrandedness::A
            };
            (ReadOrientation::Same, s)
        };
        if orientation == ReadOrientation::Same {
            self.same[l].fetch_add(1, Ordering::Relaxed);
            self.n_same.fetch_add(1, Ordering::Relaxed);
        } else {
            self.opp[l].fetch_add(1, Ordering::Relaxed);
            self.n_opp.fetch_add(1, Ordering::Relaxed);
        }
        let fmt = LibraryFormat::new(ReadType::PairedEnd, orientation, strandedness);
        self.fmt_counts[fmt.format_id() as usize].fetch_add(1, Ordering::Relaxed);
    }

    /// Total unique proper pairs seen across both orientation buckets.
    pub fn count(&self) -> u64 {
        self.n_opp.load(Ordering::Relaxed) + self.n_same.load(Ordering::Relaxed)
    }

    /// Build the FLD from the majority-orientation bucket (deterministically, in
    /// fixed length order) and infer the library format from the format tally.
    /// `fld_mean`/`fld_sd` seed the prior. Returns the cached FLD and the detected
    /// format (`None` when no proper pairs were seen).
    pub fn finish(
        &self,
        fld_mean: f64,
        fld_sd: f64,
    ) -> (FragmentLengthDistribution, Option<LibraryFormat>) {
        let n_opp = self.n_opp.load(Ordering::Relaxed);
        let n_same = self.n_same.load(Ordering::Relaxed);
        let chosen = if n_opp >= n_same {
            &self.opp
        } else {
            &self.same
        };
        let mut fld =
            FragmentLengthDistribution::new(1.0, self.max_len, fld_mean, fld_sd, 4, 0.5, 1);
        // `add_val(len, ln count)` equals `count` unit `add_val(len, 0)` calls in
        // log space, so this reproduces the per-fragment FLD — but deterministically.
        for (len, c) in chosen.iter().enumerate() {
            let c = c.load(Ordering::Relaxed);
            if c > 0 {
                fld.add_val(len, (c as f64).ln());
            }
        }
        fld.cache();
        let tally: Vec<u64> = self
            .fmt_counts
            .iter()
            .map(|a| a.load(Ordering::Relaxed))
            .collect();
        let detected = if tally.iter().sum::<u64>() > 0 {
            Some(crate::infer_format_from_counts(&tally, ReadType::PairedEnd))
        } else {
            None
        };
        (fld, detected)
    }
}

/// Where a run's fragment-length distribution came from, reported as
/// `frag_length_source` in `aux_info/meta_info.json`.
///
/// `--fldMean`/`--fldSD` are priors, so which of these applies decides whether
/// they influenced the result at all: they fully determine [`Self::Prior`], seed
/// [`Self::Reads`]/[`Self::Alignments`]/[`Self::RadDerived`] with weight 1
/// against the observations, and are not consulted at all for the two baked
/// variants.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FragLengthSource {
    /// Trained during the read-mapping pass (reads mode).
    Reads,
    /// Derived from the input alignments (alignment mode, `-a`).
    Alignments,
    /// Read from the RAD header, where it was observed by the paired-end run
    /// that wrote the file.
    RadBaked,
    /// Read from the RAD header, but written by a *single-end* run. No fragment
    /// lengths existed to observe, so the baked distribution is that run's
    /// `--fldMean`/`--fldSD` prior rather than an empirical distribution.
    RadBakedPrior,
    /// Derived at read time from this RAD's uniquely-mapped proper pairs
    /// (a piscem RAD, or `--fldPolicy derive`).
    RadDerived,
    /// The `--fldMean`/`--fldSD` prior alone, with no observations folded in.
    Prior,
}

impl FragLengthSource {
    /// The `meta_info.json` spelling.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Reads => "reads",
            Self::Alignments => "alignments",
            Self::RadBaked => "rad_baked",
            Self::RadBakedPrior => "rad_baked_prior",
            Self::RadDerived => "rad_derived",
            Self::Prior => "prior",
        }
    }

    /// Whether `--fldMean`/`--fldSD`/`--fldMax` were consulted at all. False
    /// only for the baked variants, which take the distribution verbatim from
    /// the RAD header.
    pub fn uses_fld_prior_args(self) -> bool {
        !matches!(self, Self::RadBaked | Self::RadBakedPrior)
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
            assert!(
                w[1] >= w[0] - 1e-9,
                "cond means not monotonic: {} < {}",
                w[1],
                w[0]
            );
        }
        let short = smoothed_effective_length(&cm, 201);
        assert!(
            short < 201.0 && short > 1.0,
            "short effLen {short} not shrunk"
        );
        // a long transcript keeps most of its length
        let long = smoothed_effective_length(&cm, 5000);
        assert!(long > 4000.0, "long effLen {long} shrunk too much");
        // below the 1.0 barrier the raw length is returned
        let tiny = smoothed_effective_length(&cm, 2);
        assert_eq!(tiny, 2.0, "tiny transcript should fall back to refLen");
    }

    #[test]
    fn ambig_frag_prob_bounds_and_orientation() {
        let mut fld = FragmentLengthDistribution::new(1000.0, 1000, 250.0, 25.0, 4, 0.5, 1);
        fld.cache();
        // Use the cached (frozen) CMF as a stand-in for the online snapshot.
        let cmf = fld.cached_cmf.clone();
        let txp_len = 2000i32;
        // A forward read with ample downstream space (mate fits at the typical
        // insert) should be near log(1) ≈ 0, since cmf(maxFrag) ≈ cmf(txpLen).
        let ample = ambig_frag_log_prob(&cmf, true, 100, 75, txp_len);
        assert!(ample > -0.01, "ample-space orphan logProb {ample} not ~0");
        // A forward read crammed against the 3' end (little downstream space)
        // implies an implausibly short fragment -> much smaller probability.
        let crammed = ambig_frag_log_prob(&cmf, true, txp_len - 50, 75, txp_len);
        assert!(
            crammed < ample - 1.0,
            "crammed orphan {crammed} not << ample {ample}"
        );
        // Reverse-strand orientation uses pos + read_len for the upstream bound:
        // a reverse read whose outer end is near the 5' start is likewise crammed.
        let rc_crammed = ambig_frag_log_prob(&cmf, false, 0, 50, txp_len);
        assert!(
            rc_crammed < ample - 1.0,
            "rc crammed orphan {rc_crammed} not << ample {ample}"
        );
        // Empty snapshot -> unmodelled (LOG_1 = 0).
        assert_eq!(ambig_frag_log_prob(&[], true, 100, 75, txp_len), 0.0);
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
