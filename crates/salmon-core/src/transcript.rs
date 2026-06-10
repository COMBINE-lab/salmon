//! Transcript record: identity, lengths, accumulated mass/counts, and the
//! effective-length computation.
//!
//! Mirrors the modeling-relevant parts of salmon's `Transcript`
//! (`include/.../model/Transcript.hpp`). Mutable count/mass state is stored in
//! atomics so worker threads can update it concurrently during a quantification
//! round, matching the C++ design. Per-transcript sequence storage and bias
//! state are added in a later phase (bias models).

use crate::atomic::AtomicF64;
use crate::math::{log_add, LOG_0};
use std::sync::atomic::{AtomicU64, Ordering};

/// A reference transcript and its mutable quantification state.
#[derive(Debug)]
pub struct Transcript {
    /// dense transcript id (index into the experiment's transcript vector)
    pub id: u32,
    /// reference name as it appears in the index / output
    pub name: String,
    /// length used for modeling (may be poly-A clipped)
    pub ref_length: u32,
    /// full reference length (reported in `quant.sf` `Length` column)
    pub complete_length: u32,
    /// true if this is a decoy sequence
    pub is_decoy: bool,

    /// cached log effective length
    log_eff_length: AtomicF64,
    /// fragments mapping uniquely to this transcript
    unique_count: AtomicU64,
    /// total fragments touching this transcript (unique + shared)
    total_count: AtomicU64,
    /// accumulated (fractional) mass assigned to this transcript, log space
    mass: AtomicF64,
}

impl Transcript {
    pub fn new(id: u32, name: impl Into<String>, complete_length: u32) -> Self {
        Self {
            id,
            name: name.into(),
            ref_length: complete_length,
            complete_length,
            is_decoy: false,
            log_eff_length: AtomicF64::new(LOG_0),
            unique_count: AtomicU64::new(0),
            total_count: AtomicU64::new(0),
            mass: AtomicF64::new(LOG_0),
        }
    }

    pub fn with_ref_length(mut self, ref_length: u32) -> Self {
        self.ref_length = ref_length;
        self
    }

    pub fn as_decoy(mut self) -> Self {
        self.is_decoy = true;
        self
    }

    // --- count / mass accessors ------------------------------------------------

    pub fn unique_count(&self) -> u64 {
        self.unique_count.load(Ordering::Relaxed)
    }
    pub fn total_count(&self) -> u64 {
        self.total_count.load(Ordering::Relaxed)
    }
    pub fn add_unique(&self, n: u64) {
        self.unique_count.fetch_add(n, Ordering::Relaxed);
    }
    pub fn add_total(&self, n: u64) {
        self.total_count.fetch_add(n, Ordering::Relaxed);
    }

    /// Current accumulated mass (log space).
    pub fn mass(&self) -> f64 {
        self.mass.load()
    }
    pub fn set_mass(&self, log_mass: f64) {
        self.mass.store(log_mass);
    }
    /// Atomically add `log_inc` (a log-space value) to the accumulated mass.
    pub fn add_mass(&self, log_inc: f64) {
        self.mass.log_add_assign(log_inc);
    }

    // --- effective length ------------------------------------------------------

    pub fn log_effective_length(&self) -> f64 {
        self.log_eff_length.load()
    }
    pub fn set_log_effective_length(&self, v: f64) {
        self.log_eff_length.store(v);
    }
    pub fn effective_length(&self) -> f64 {
        self.log_eff_length.load().exp()
    }

    /// Compute the log effective length given a log-space fragment-length PMF
    /// covering lengths `[min_val, min_val + log_pmf.len())`.
    ///
    /// Mirrors salmon's `Transcript::computeLogEffectiveLength`:
    /// `effLen = sum_{l=min}^{min(refLen,max)} pmf(l) * (refLen - l + 1)` in
    /// log space. If the transcript is shorter than the minimum fragment
    /// length, it falls back to `log(refLen)`.
    pub fn compute_log_effective_length(&self, log_pmf: &[f64], min_val: usize) -> f64 {
        let ref_len = self.ref_length as usize;
        let max_val = min_val + log_pmf.len().saturating_sub(1);
        let upper = ref_len.min(max_val);

        let mut eff = LOG_0;
        if ref_len >= min_val {
            let mut l = min_val;
            while l <= upper {
                let pmf = log_pmf[l - min_val];
                // log(refLen - l + 1)
                let span = ((ref_len - l + 1) as f64).ln();
                eff = log_add(eff, pmf + span);
                l += 1;
            }
        }

        // Guard: a too-short transcript or empty support falls back to log(refLen).
        if !(eff.is_finite()) || eff < 1.0_f64.ln() {
            (ref_len.max(1) as f64).ln()
        } else {
            eff
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::LOG_1;

    #[test]
    fn mass_accumulates_in_log_space() {
        let t = Transcript::new(0, "t0", 1000);
        assert_eq!(t.mass(), LOG_0);
        t.add_mass(LOG_1); // add 1.0
        t.add_mass(LOG_1); // add 1.0 -> log(2)
        assert!((t.mass() - 2.0_f64.ln()).abs() < 1e-9);
    }

    #[test]
    fn counts_accumulate() {
        let t = Transcript::new(1, "t1", 500);
        t.add_unique(3);
        t.add_total(5);
        assert_eq!(t.unique_count(), 3);
        assert_eq!(t.total_count(), 5);
    }

    #[test]
    fn eff_length_point_mass_is_intuitive() {
        // A point-mass FLD at length 100 on a 1000nt transcript gives
        // effLen = 1000 - 100 + 1 = 901.
        let t = Transcript::new(0, "t", 1000);
        let min_val = 100usize;
        let log_pmf = vec![LOG_1]; // pmf(100) = 1
        let eff = t.compute_log_effective_length(&log_pmf, min_val).exp();
        assert!((eff - 901.0).abs() < 1e-6, "got {eff}");
    }

    #[test]
    fn eff_length_short_transcript_falls_back() {
        // transcript shorter than the minimum fragment length -> log(refLen)
        let t = Transcript::new(0, "t", 50);
        let log_pmf = vec![LOG_1];
        let eff = t.compute_log_effective_length(&log_pmf, 100).exp();
        assert!((eff - 50.0).abs() < 1e-6, "got {eff}");
    }
}
