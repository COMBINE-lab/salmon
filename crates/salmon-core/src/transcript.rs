//! Transcript record: identity, lengths, accumulated mass/counts, and the
//! effective-length computation.
//!
//! Mirrors the modeling-relevant parts of salmon's `Transcript`
//! (`include/.../model/Transcript.hpp`). Mutable count/mass state is stored in
//! atomics so worker threads can update it concurrently during a quantification
//! round, matching the C++ design. Per-transcript sequence storage and bias
//! state are added in a later phase (bias models).
//!
//! # Effective length, in one paragraph
//!
//! A transcript's *effective* length is the number of distinct positions a
//! fragment could have started at, given how long fragments actually are. If
//! fragments are 200 nt and a transcript is 1000 nt, only the first 801 positions
//! can start a fragment that fits, so its effective length is 801, not 1000. This
//! matters because a longer transcript yields proportionally more fragments at
//! the same molar abundance; dividing by effective length is what converts
//! "fragments seen" into "molecules present". A transcript shorter than a typical
//! fragment is the awkward case, handled explicitly below.

use crate::atomic::AtomicF64;
use crate::math::{log_add, LOG_0};
use std::sync::atomic::{AtomicU64, Ordering};

/// A reference transcript and its mutable quantification state.
///
/// Note every mutable field is atomic and every mutating method takes `&self`
/// rather than `&mut self`: many mapping threads share one `Transcript` and
/// update it simultaneously, so there is no exclusive access to hand out.
#[derive(Debug)]
pub struct Transcript {
    /// dense transcript id (index into the experiment's transcript vector)
    pub id: u32,
    /// reference name as it appears in the index / output
    pub name: String,
    /// length used for modeling (may be poly-A clipped)
    ///
    /// Poly-A tails are added during library prep and are not part of the
    /// sequence reads come from, so modeling uses the clipped length while
    /// reporting keeps the original.
    pub ref_length: u32,
    /// full reference length (reported in `quant.sf` `Length` column)
    pub complete_length: u32,
    /// true if this is a decoy sequence (see [`crate::RefProvider::is_decoy`])
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
    /// Create a transcript with no accumulated state.
    ///
    /// `impl Into<String>` means the caller may pass either a `String` or a
    /// `&str`; the conversion happens here.
    pub fn new(id: u32, name: impl Into<String>, complete_length: u32) -> Self {
        Self {
            id,
            name: name.into(),
            // Default to unclipped; `with_ref_length` overrides when poly-A
            // clipping applies.
            ref_length: complete_length,
            complete_length,
            is_decoy: false,
            // LOG_0 is -inf, i.e. zero mass in log space — nothing observed yet.
            log_eff_length: AtomicF64::new(LOG_0),
            unique_count: AtomicU64::new(0),
            total_count: AtomicU64::new(0),
            mass: AtomicF64::new(LOG_0),
        }
    }

    /// Builder-style override of the modeling length. Takes and returns `self`
    /// so calls can be chained: `Transcript::new(..).with_ref_length(..)`.
    pub fn with_ref_length(mut self, ref_length: u32) -> Self {
        self.ref_length = ref_length;
        self
    }

    /// Builder-style marker for decoy sequences.
    pub fn as_decoy(mut self) -> Self {
        self.is_decoy = true;
        self
    }

    // --- count / mass accessors ------------------------------------------------
    //
    // All of these use `Ordering::Relaxed`: the counters are independent tallies
    // read only after the threads have joined, so no cross-thread ordering needs
    // to be enforced and the cheapest atomic operation is the right one.

    /// Fragments that mapped to this transcript and nowhere else.
    pub fn unique_count(&self) -> u64 {
        self.unique_count.load(Ordering::Relaxed)
    }
    /// Fragments compatible with this transcript, whether or not also with others.
    pub fn total_count(&self) -> u64 {
        self.total_count.load(Ordering::Relaxed)
    }
    /// Add to the unique tally from any thread.
    pub fn add_unique(&self, n: u64) {
        self.unique_count.fetch_add(n, Ordering::Relaxed);
    }
    /// Add to the total tally from any thread.
    pub fn add_total(&self, n: u64) {
        self.total_count.fetch_add(n, Ordering::Relaxed);
    }

    /// Current accumulated mass (log space).
    ///
    /// "Mass" is fractional read assignment: an ambiguous fragment splits itself
    /// across the transcripts it is compatible with, so this is generally not a
    /// whole number.
    pub fn mass(&self) -> f64 {
        self.mass.load()
    }
    /// Overwrite the accumulated mass (used when a round resets it).
    pub fn set_mass(&self, log_mass: f64) {
        self.mass.store(log_mass);
    }
    /// Atomically add `log_inc` (a log-space value) to the accumulated mass.
    ///
    /// This is a *probability* addition performed in log space, not an addition
    /// of logs; see [`crate::math::log_add`].
    pub fn add_mass(&self, log_inc: f64) {
        self.mass.log_add_assign(log_inc);
    }

    // --- effective length ------------------------------------------------------

    /// Cached effective length, in log space (the form the models use).
    pub fn log_effective_length(&self) -> f64 {
        self.log_eff_length.load()
    }
    /// Store a freshly computed effective length.
    pub fn set_log_effective_length(&self, v: f64) {
        self.log_eff_length.store(v);
    }
    /// Cached effective length in ordinary linear space (for reporting).
    pub fn effective_length(&self) -> f64 {
        self.log_eff_length.load().exp()
    }

    /// Compute the log effective length given a log-space fragment-length PMF
    /// covering lengths `[min_val, min_val + log_pmf.len())`.
    ///
    /// **The formula.** Mirrors salmon's `Transcript::computeLogEffectiveLength`:
    /// `effLen = sum_{l=min}^{min(refLen,max)} pmf(l) * (refLen - l + 1)` in
    /// log space. Read it as: for each possible fragment length `l`, weight the
    /// number of start positions that fit (`refLen - l + 1`) by how probable that
    /// length is (`pmf(l)`), and sum. Lengths longer than the transcript
    /// contribute nothing and are excluded by the `upper` bound.
    ///
    /// "PMF" is a probability mass function: `pmf(l)` is the probability that a
    /// fragment has exactly length `l`.
    ///
    /// If the transcript is shorter than the minimum fragment length, it falls
    /// back to `log(refLen)`.
    pub fn compute_log_effective_length(&self, log_pmf: &[f64], min_val: usize) -> f64 {
        let ref_len = self.ref_length as usize;
        // Largest length the PMF describes. `saturating_sub` avoids underflow on
        // an empty PMF (0 - 1 would wrap around on an unsigned integer).
        let max_val = min_val + log_pmf.len().saturating_sub(1);
        // No fragment longer than the transcript can fit in it.
        let upper = ref_len.min(max_val);

        let mut eff = LOG_0;
        // Skip the loop entirely for a transcript shorter than any fragment;
        // there is no valid `l` and the fallback below applies.
        if ref_len >= min_val {
            let mut l = min_val;
            while l <= upper {
                let pmf = log_pmf[l - min_val];
                // log(refLen - l + 1)
                let span = ((ref_len - l + 1) as f64).ln();
                // Adding logs = multiplying probabilities; `log_add` then sums
                // the products in log space.
                eff = log_add(eff, pmf + span);
                l += 1;
            }
        }

        // Guard: a too-short transcript or empty support falls back to log(refLen).
        //
        // `eff < log(1)` means an effective length below one position, which
        // would make the downstream division explode; `max(1)` likewise keeps
        // `log` defined for a zero-length reference.
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

    /// Two unit masses must accumulate to a mass of 2, expressed as log 2 — the
    /// check that log-space accumulation really is probability addition.
    #[test]
    fn mass_accumulates_in_log_space() {
        let t = Transcript::new(0, "t0", 1000);
        assert_eq!(t.mass(), LOG_0);
        t.add_mass(LOG_1); // add 1.0
        t.add_mass(LOG_1); // add 1.0 -> log(2)
        assert!((t.mass() - 2.0_f64.ln()).abs() < 1e-9);
    }

    /// The unique and total tallies are independent.
    #[test]
    fn counts_accumulate() {
        let t = Transcript::new(1, "t1", 500);
        t.add_unique(3);
        t.add_total(5);
        assert_eq!(t.unique_count(), 3);
        assert_eq!(t.total_count(), 5);
    }

    /// With all fragments at exactly one length, the effective length collapses
    /// to the closed form `refLen - fragLen + 1`, which is easy to verify by
    /// hand and pins the off-by-one.
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

    /// A transcript no fragment can fit in must fall back to its raw length
    /// rather than producing zero or a negative effective length.
    #[test]
    fn eff_length_short_transcript_falls_back() {
        // transcript shorter than the minimum fragment length -> log(refLen)
        let t = Transcript::new(0, "t", 50);
        let log_pmf = vec![LOG_1];
        let eff = t.compute_log_effective_length(&log_pmf, 100).exp();
        assert!((eff - 50.0).abs() < 1e-6, "got {eff}");
    }
}
