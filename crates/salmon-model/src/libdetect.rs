//! Automatic library-type detection.
//!
//! Port of salmon's `LibraryTypeDetector`
//! (`include/.../model/LibraryTypeDetector.hpp`). During the first reads of a
//! run, the observed [`LibraryFormat`] of each confidently mapped fragment is
//! tallied; once enough samples are seen, the most likely orientation and
//! strandedness are inferred from the count ratios using salmon's 30%/70%
//! thresholds.

use salmon_core::{LibraryFormat, ReadOrientation, ReadStrandedness, ReadType};
use std::sync::atomic::{AtomicBool, AtomicI64, AtomicU64, AtomicU8, Ordering};

/// Default number of samples to collect before guessing (matches salmon).
pub const DEFAULT_SAMPLES_NEEDED: i64 = 50_000;

/// Sentinel for `resolved` meaning "not yet locked in" (no valid format has this
/// id; `MAX_FORMAT_ID` is 11).
const UNSET_FORMAT: u8 = 0xFF;

/// Accumulates observed library formats and infers the most likely type.
#[derive(Debug)]
pub struct LibraryTypeDetector {
    /// still sampling (gates [`add_sample`]); cleared once the format locks in
    active: AtomicBool,
    read_type: ReadType,
    samples_needed: AtomicI64,
    counts: Vec<AtomicU64>,
    /// the locked-in format id once detection completes, else [`UNSET_FORMAT`]
    resolved: AtomicU8,
}

impl LibraryTypeDetector {
    pub fn new(read_type: ReadType) -> Self {
        let counts = (0..=LibraryFormat::MAX_FORMAT_ID)
            .map(|_| AtomicU64::new(0))
            .collect();
        Self {
            active: AtomicBool::new(true),
            read_type,
            samples_needed: AtomicI64::new(DEFAULT_SAMPLES_NEEDED),
            counts,
            resolved: AtomicU8::new(UNSET_FORMAT),
        }
    }

    pub fn is_active(&self) -> bool {
        self.active.load(Ordering::Relaxed)
    }

    /// True once enough samples have been collected to guess.
    pub fn can_guess(&self) -> bool {
        self.samples_needed.load(Ordering::Relaxed) <= 0
    }

    /// Record one confidently mapped fragment's observed format. Only formats
    /// matching the detector's read type are counted, and only until the sample
    /// budget is exhausted. Thread-safe.
    pub fn add_sample(&self, f: LibraryFormat) {
        if f.read_type == self.read_type && self.samples_needed.load(Ordering::Relaxed) >= 0 {
            self.counts[f.format_id() as usize].fetch_add(1, Ordering::Relaxed);
            self.samples_needed.fetch_sub(1, Ordering::Relaxed);
        }
    }

    /// Mid-run resolution (salmon's prefix-detect-then-apply): once enough
    /// samples have been collected ([`can_guess`]), infer and **lock in** the
    /// library format (one writer wins the CAS), stop sampling, and return it;
    /// idempotent thereafter. Returns `None` while still sampling. The caller
    /// applies the returned format as a strand-compatibility filter for the rest
    /// of the run.
    pub fn resolved_format(&self) -> Option<LibraryFormat> {
        let r = self.resolved.load(Ordering::Acquire);
        if r != UNSET_FORMAT {
            return Some(LibraryFormat::from_format_id(r));
        }
        if !self.can_guess() {
            return None;
        }
        let f = self.infer_format();
        match self.resolved.compare_exchange(
            UNSET_FORMAT,
            f.format_id(),
            Ordering::AcqRel,
            Ordering::Acquire,
        ) {
            Ok(_) => {
                self.active.store(false, Ordering::Release);
                Some(f)
            }
            // Another thread locked in first; use its result.
            Err(existing) => Some(LibraryFormat::from_format_id(existing)),
        }
    }

    /// The final library format to report at end of run: the locked-in format if
    /// resolution happened mid-run, else inferred from whatever samples were
    /// collected (recorded so repeat calls agree). Always returns a format.
    pub fn final_format(&self) -> LibraryFormat {
        let r = self.resolved.load(Ordering::Acquire);
        if r != UNSET_FORMAT {
            return LibraryFormat::from_format_id(r);
        }
        let f = self.infer_format();
        let _ = self.resolved.compare_exchange(
            UNSET_FORMAT,
            f.format_id(),
            Ordering::AcqRel,
            Ordering::Acquire,
        );
        LibraryFormat::from_format_id(self.resolved.load(Ordering::Acquire))
    }

    /// Pure inference of the most likely library format from the accumulated
    /// counts (no state change). Falls back to inward/unstranded when there are
    /// no usable samples.
    fn infer_format(&self) -> LibraryFormat {
        let count = |id: u8| self.counts[id as usize].load(Ordering::Relaxed);

        match self.read_type {
            ReadType::SingleEnd => {
                let mut nf = 0u64;
                let mut nr = 0u64;
                for id in 0..=LibraryFormat::MAX_FORMAT_ID {
                    let f = LibraryFormat::from_format_id(id);
                    let c = count(id);
                    nf += if f.strandedness == ReadStrandedness::S {
                        c
                    } else {
                        0
                    };
                    nr += if f.strandedness == ReadStrandedness::A {
                        c
                    } else {
                        0
                    };
                }
                let strandedness = if nf + nr == 0 {
                    ReadStrandedness::U
                } else {
                    // Single-end uses the matching (S/A) encoding, like a
                    // paired "same"-orientation library.
                    strandedness_from_fw_ratio(nf as f64 / (nf + nr) as f64, true)
                };
                LibraryFormat::new(ReadType::SingleEnd, ReadOrientation::None, strandedness)
            }
            ReadType::PairedEnd => {
                let (mut nsf, mut nsr) = (0u64, 0u64);
                let (mut ninward, mut noutward, mut nsame) = (0u64, 0u64, 0u64);
                for id in 0..=LibraryFormat::MAX_FORMAT_ID {
                    let f = LibraryFormat::from_format_id(id);
                    let c = count(id);
                    nsf += matches!(f.strandedness, ReadStrandedness::S | ReadStrandedness::SA)
                        .then_some(c)
                        .unwrap_or(0);
                    nsr += matches!(f.strandedness, ReadStrandedness::A | ReadStrandedness::AS)
                        .then_some(c)
                        .unwrap_or(0);
                    match f.orientation {
                        ReadOrientation::Toward => ninward += c,
                        ReadOrientation::Away => noutward += c,
                        ReadOrientation::Same => nsame += c,
                        ReadOrientation::None => {}
                    }
                }

                let num_orient = ninward + noutward + nsame;
                if num_orient > 0 && (nsf + nsr) > 0 {
                    let ratio_in = ninward as f64 / num_orient as f64;
                    let ratio_out = noutward as f64 / num_orient as f64;
                    let ratio_same = nsame as f64 / num_orient as f64;

                    let (orientation, same) = if ratio_in >= ratio_out && ratio_in >= ratio_same {
                        (ReadOrientation::Toward, false)
                    } else if ratio_out >= ratio_in && ratio_out >= ratio_same {
                        (ReadOrientation::Away, false)
                    } else {
                        (ReadOrientation::Same, true)
                    };

                    let ratio_fw = nsf as f64 / (nsf + nsr) as f64;
                    let strandedness = strandedness_from_fw_ratio(ratio_fw, same);
                    LibraryFormat::new(ReadType::PairedEnd, orientation, strandedness)
                } else {
                    LibraryFormat::new(
                        ReadType::PairedEnd,
                        ReadOrientation::Toward,
                        ReadStrandedness::U,
                    )
                }
            }
        }
    }
}

/// Map a forward-strand fraction to a strandedness using salmon's 30%/70%
/// thresholds. `same` selects between the matching (S/A) and opposite (SA/AS)
/// stranded encodings for paired-end "same"-orientation libraries; for
/// single-end pass `false`.
fn strandedness_from_fw_ratio(ratio_fw: f64, same: bool) -> ReadStrandedness {
    if ratio_fw < 0.3 {
        if same {
            ReadStrandedness::A
        } else {
            ReadStrandedness::AS
        }
    } else if ratio_fw < 0.7 {
        ReadStrandedness::U
    } else if same {
        ReadStrandedness::S
    } else {
        ReadStrandedness::SA
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_end_detects_sense() {
        let d = LibraryTypeDetector::new(ReadType::SingleEnd);
        let sf = LibraryFormat::parse("SF").unwrap();
        let sr = LibraryFormat::parse("SR").unwrap();
        for _ in 0..90 {
            d.add_sample(sf);
        }
        for _ in 0..10 {
            d.add_sample(sr);
        }
        assert_eq!(d.infer_format().canonical(), "SF");
        // final_format records and returns the same result idempotently
        assert_eq!(d.final_format().canonical(), "SF");
        assert_eq!(d.final_format().canonical(), "SF");
    }

    #[test]
    fn single_end_detects_unstranded() {
        let d = LibraryTypeDetector::new(ReadType::SingleEnd);
        let sf = LibraryFormat::parse("SF").unwrap();
        let sr = LibraryFormat::parse("SR").unwrap();
        for _ in 0..50 {
            d.add_sample(sf);
            d.add_sample(sr);
        }
        assert_eq!(d.infer_format().canonical(), "U");
    }

    #[test]
    fn paired_end_detects_isr() {
        let d = LibraryTypeDetector::new(ReadType::PairedEnd);
        let isr = LibraryFormat::parse("ISR").unwrap();
        for _ in 0..100 {
            d.add_sample(isr);
        }
        // ISR: inward + antisense -> toward + AS
        assert_eq!(d.infer_format().canonical(), "ISR");
    }

    #[test]
    fn paired_end_detects_iu() {
        let d = LibraryTypeDetector::new(ReadType::PairedEnd);
        let isf = LibraryFormat::parse("ISF").unwrap();
        let isr = LibraryFormat::parse("ISR").unwrap();
        for _ in 0..50 {
            d.add_sample(isf);
            d.add_sample(isr);
        }
        // balanced strandedness -> unstranded, inward -> IU
        assert_eq!(d.infer_format().canonical(), "IU");
    }

    #[test]
    fn resolved_format_locks_in_after_prefix() {
        let d = LibraryTypeDetector::new(ReadType::PairedEnd);
        let isr = LibraryFormat::parse("ISR").unwrap();
        // Before the sample budget is consumed: no resolution yet, so the caller
        // applies no filter, and the detector keeps sampling.
        assert!(d.resolved_format().is_none());
        assert!(d.is_active());
        // Feed the full prefix budget.
        for _ in 0..DEFAULT_SAMPLES_NEEDED {
            d.add_sample(isr);
        }
        assert!(d.can_guess());
        // Now it locks in to the inferred type, stops sampling, and is idempotent.
        assert_eq!(d.resolved_format().unwrap().canonical(), "ISR");
        assert!(!d.is_active());
        assert_eq!(d.resolved_format().unwrap().canonical(), "ISR");
        assert_eq!(d.final_format().canonical(), "ISR");
    }

    #[test]
    fn final_format_without_lockin_infers_from_partial() {
        // Fewer than the budget: never locks in mid-run, but end-of-run reporting
        // still returns a best-guess format from the partial samples.
        let d = LibraryTypeDetector::new(ReadType::PairedEnd);
        let isf = LibraryFormat::parse("ISF").unwrap();
        for _ in 0..100 {
            d.add_sample(isf);
        }
        assert!(d.resolved_format().is_none()); // not enough to lock in mid-run
        assert_eq!(d.final_format().canonical(), "ISF");
    }

    #[test]
    fn sample_budget_is_respected() {
        let d = LibraryTypeDetector::new(ReadType::SingleEnd);
        assert!(!d.can_guess());
        let sf = LibraryFormat::parse("SF").unwrap();
        // exhaust the budget
        let mut n = DEFAULT_SAMPLES_NEEDED + 5;
        while n > 0 {
            d.add_sample(sf);
            n -= 1;
        }
        assert!(d.can_guess());
    }
}
