//! Automatic library-type detection.
//!
//! Port of salmon's `LibraryTypeDetector`
//! (`include/.../model/LibraryTypeDetector.hpp`). During the first reads of a
//! run, the observed [`LibraryFormat`] of each confidently mapped fragment is
//! tallied; once enough samples are seen, the most likely orientation and
//! strandedness are inferred from the count ratios using salmon's 30%/70%
//! thresholds.

use salmon_core::{LibraryFormat, ReadOrientation, ReadStrandedness, ReadType};
use std::sync::atomic::{AtomicBool, AtomicI64, AtomicU64, Ordering};

/// Default number of samples to collect before guessing (matches salmon).
pub const DEFAULT_SAMPLES_NEEDED: i64 = 50_000;

/// Accumulates observed library formats and infers the most likely type.
#[derive(Debug)]
pub struct LibraryTypeDetector {
    active: AtomicBool,
    read_type: ReadType,
    samples_needed: AtomicI64,
    counts: Vec<AtomicU64>,
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

    /// Infer the most likely library format from the accumulated counts. Marks
    /// the detector inactive (so this is effectively done once). Returns `None`
    /// if the detector was already inactive.
    pub fn most_likely_type(&self) -> Option<LibraryFormat> {
        if !self.active.swap(false, Ordering::AcqRel) {
            return None;
        }

        let count = |id: u8| self.counts[id as usize].load(Ordering::Relaxed);

        let fmt = match self.read_type {
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
        };

        Some(fmt)
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
        let f = d.most_likely_type().unwrap();
        assert_eq!(f.canonical(), "SF");
        // detector becomes inactive after guessing
        assert!(d.most_likely_type().is_none());
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
        assert_eq!(d.most_likely_type().unwrap().canonical(), "U");
    }

    #[test]
    fn paired_end_detects_isr() {
        let d = LibraryTypeDetector::new(ReadType::PairedEnd);
        let isr = LibraryFormat::parse("ISR").unwrap();
        for _ in 0..100 {
            d.add_sample(isr);
        }
        // ISR: inward + antisense -> toward + AS
        assert_eq!(d.most_likely_type().unwrap().canonical(), "ISR");
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
        assert_eq!(d.most_likely_type().unwrap().canonical(), "IU");
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
