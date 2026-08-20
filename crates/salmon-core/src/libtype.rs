//! Library format types and parsing.
//!
//! # What a "library type" is
//!
//! An RNA-seq protocol makes specific promises about how the sequenced reads
//! relate to the original RNA molecule, and salmon uses those promises to reject
//! placements that could not have arisen. There are three independent facts:
//!
//! 1. **Read type** — one read per fragment (single-end) or two (paired-end).
//! 2. **Orientation** — for pairs, how the two mates point relative to each
//!    other: toward each other (`I`, the Illumina default), away (`O`), or the
//!    same way (`M`).
//! 3. **Strandedness** — whether the protocol preserved which strand the RNA came
//!    from. An *unstranded* protocol (`U`) loses that information; a *stranded*
//!    one records it, so read 1 is expected on the sense strand (`SF`) or the
//!    antisense strand (`SR`).
//!
//! Written together these give the familiar salmon strings: `IU` (inward,
//! unstranded), `ISR` (inward, read 1 antisense), `SF` (single-end, sense), and
//! so on. Getting this wrong is one of the most common causes of a bad
//! quantification, which is why salmon can also detect it automatically (`-l A`)
//! and warns when an explicit setting disagrees with what it observes.
//!
//! Mirrors salmon's `LibraryFormat` (defined in the vendored RapMap/pufferfish
//! headers) and the `parseLibraryFormatStringNew` logic in
//! `src/util/LibraryTypeUtils.cpp`. A library format is the triple of read
//! type, relative orientation, and strandedness; only the 12 combinations
//! salmon accepts are representable through [`LibraryFormat::parse`].

use crate::error::SalmonError;
use crate::mate::MateStatus;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Single- vs paired-end.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ReadType {
    SingleEnd,
    PairedEnd,
}

/// Relative orientation of the two mates of a fragment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ReadOrientation {
    /// mates point toward each other (`I`, "inward") — the Illumina default
    Toward,
    /// mates point away from each other (`O`, "outward")
    Away,
    /// mates map to the same strand (`M`, "matching")
    Same,
    /// not applicable (single-end)
    None,
}

/// Strand from which fragments are expected to originate.
///
/// "Sense" means the same strand as the transcript; "antisense" the opposite.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ReadStrandedness {
    /// read 1 sense, read 2 antisense
    SA,
    /// read 1 antisense, read 2 sense
    AS,
    /// sense (single-end / matching)
    S,
    /// antisense (single-end / matching)
    A,
    /// unstranded — the protocol did not preserve strand information, so both
    /// orientations are acceptable
    U,
}

/// The full library format: read type + orientation + strandedness.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct LibraryFormat {
    pub read_type: ReadType,
    pub orientation: ReadOrientation,
    pub strandedness: ReadStrandedness,
}

impl LibraryFormat {
    /// Assemble a format from its three parts.
    ///
    /// `const fn` means this can be evaluated at compile time, which is what
    /// lets the defaults below be constants.
    pub const fn new(
        read_type: ReadType,
        orientation: ReadOrientation,
        strandedness: ReadStrandedness,
    ) -> Self {
        Self {
            read_type,
            orientation,
            strandedness,
        }
    }

    /// Default paired-end format used by salmon (`IU`): inward-facing mates,
    /// no strand information assumed. The safest default, since it rejects the
    /// least.
    pub const fn paired_default() -> Self {
        Self::new(
            ReadType::PairedEnd,
            ReadOrientation::Toward,
            ReadStrandedness::U,
        )
    }

    /// Default single-end format used by salmon (`U`).
    pub const fn single_default() -> Self {
        Self::new(
            ReadType::SingleEnd,
            ReadOrientation::None,
            ReadStrandedness::U,
        )
    }

    /// Parse a library-type string (case-insensitive), e.g. `"IU"`, `"ISR"`,
    /// `"SF"`. Mirrors `parseLibraryFormatStringNew`.
    ///
    /// Only these 12 strings are accepted. The exhaustive table is deliberate:
    /// the alternative — parsing the letters independently — would accept
    /// nonsense like `OSS` and silently build a format no protocol produces.
    pub fn parse(s: &str) -> Result<Self, SalmonError> {
        // `use` inside a function pulls the variant names into scope just for
        // this block, so the table below reads as `Toward` rather than
        // `ReadOrientation::Toward` twelve times over.
        use ReadOrientation::*;
        use ReadStrandedness::*;
        use ReadType::*;
        let up = s.to_ascii_uppercase();
        let fmt = match up.as_str() {
            "IU" => Self::new(PairedEnd, Toward, U),
            "ISF" => Self::new(PairedEnd, Toward, SA),
            "ISR" => Self::new(PairedEnd, Toward, AS),
            "OU" => Self::new(PairedEnd, Away, U),
            "OSF" => Self::new(PairedEnd, Away, SA),
            "OSR" => Self::new(PairedEnd, Away, AS),
            "MU" => Self::new(PairedEnd, Same, U),
            "MSF" => Self::new(PairedEnd, Same, S),
            "MSR" => Self::new(PairedEnd, Same, A),
            "U" => Self::new(SingleEnd, None, U),
            "SF" => Self::new(SingleEnd, None, S),
            "SR" => Self::new(SingleEnd, None, A),
            _ => return Err(SalmonError::UnknownLibraryFormat(s.to_string())),
        };
        Ok(fmt)
    }

    /// Is this the `A`/`a` auto-detect sentinel? (handled by the caller, not a
    /// real [`LibraryFormat`]).
    ///
    /// `A` is not a format at all — it means "work it out from the data" — so it
    /// cannot be represented by this type and is tested for separately.
    pub fn is_auto(s: &str) -> bool {
        s.eq_ignore_ascii_case("a")
    }

    /// The canonical string for this format, e.g. `"ISR"`.
    ///
    /// The exact inverse of [`Self::parse`] for the 12 real formats. Matching on
    /// the triple at once is what makes the mapping total.
    pub fn canonical(&self) -> &'static str {
        use ReadOrientation::*;
        use ReadStrandedness::*;
        use ReadType::*;
        match (self.read_type, self.orientation, self.strandedness) {
            (PairedEnd, Toward, U) => "IU",
            (PairedEnd, Toward, SA) => "ISF",
            (PairedEnd, Toward, AS) => "ISR",
            (PairedEnd, Away, U) => "OU",
            (PairedEnd, Away, SA) => "OSF",
            (PairedEnd, Away, AS) => "OSR",
            (PairedEnd, Same, U) => "MU",
            (PairedEnd, Same, S) => "MSF",
            (PairedEnd, Same, A) => "MSR",
            // Single-end has no meaningful orientation, so `_` ignores it.
            (SingleEnd, _, U) => "U",
            (SingleEnd, _, S) => "SF",
            (SingleEnd, _, A) => "SR",
            // Any non-canonical combination collapses to unstranded of its type.
            // Reachable only if a format is assembled field-by-field rather than
            // parsed; degrading to the permissive default beats panicking.
            (PairedEnd, _, _) => "IU",
            (SingleEnd, _, _) => "U",
        }
    }

    /// Stable integer id for the 12 canonical formats. Useful as a dense index
    /// (e.g. for `lib_format_counts.json` accounting).
    ///
    /// "Dense" means the ids are exactly `0..=11` with no gaps, so they can index
    /// straight into a fixed-size array instead of needing a hash map. The ids
    /// are also written into RAD files, so they must never be renumbered.
    pub fn format_id(&self) -> u8 {
        match self.canonical() {
            "IU" => 0,
            "ISF" => 1,
            "ISR" => 2,
            "OU" => 3,
            "OSF" => 4,
            "OSR" => 5,
            "MU" => 6,
            "MSF" => 7,
            "MSR" => 8,
            "U" => 9,
            "SF" => 10,
            "SR" => 11,
            // `unreachable!` documents an invariant: `canonical()` returns one of
            // the strings above by construction, so reaching here would mean the
            // two functions had drifted apart.
            _ => unreachable!("canonical() only returns the 12 known formats"),
        }
    }

    /// Whether this format describes paired-end reads.
    pub fn is_paired(&self) -> bool {
        matches!(self.read_type, ReadType::PairedEnd)
    }

    /// Largest valid [`format_id`](Self::format_id) (there are 12 formats, ids 0..=11).
    pub const MAX_FORMAT_ID: u8 = 11;

    /// Inverse of [`format_id`](Self::format_id). Mirrors `formatFromID`.
    ///
    /// Panics on an unknown id: these come from salmon's own files, so an
    /// out-of-range value means the input is corrupt rather than merely unusual.
    pub fn from_format_id(id: u8) -> Self {
        let s = match id {
            0 => "IU",
            1 => "ISF",
            2 => "ISR",
            3 => "OU",
            4 => "OSF",
            5 => "OSR",
            6 => "MU",
            7 => "MSF",
            8 => "MSR",
            9 => "U",
            10 => "SF",
            11 => "SR",
            _ => panic!("invalid library format id: {id}"),
        };
        Self::parse(s).expect("canonical format string parses")
    }
}

// `Display` is the standard trait for user-facing text, so a format can be
// dropped straight into a log line or a `format!` string.
impl fmt::Display for LibraryFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.canonical())
    }
}

/// Whether a properly-paired mapping with the `observed` format is compatible
/// with the `expected` library format. Ports salmon's `compatibleHit(expected,
/// observed)`: orientations must match, and the expected strandedness must be
/// unstranded or equal to the observed.
///
/// The asymmetry is the point: an *unstranded expectation* accepts any observed
/// strandedness (the protocol never promised one), while a *stranded
/// expectation* accepts only its own. Orientation, by contrast, is a physical
/// property of the protocol and must always agree.
pub fn compatible_paired(expected: LibraryFormat, observed: LibraryFormat) -> bool {
    if expected.orientation != observed.orientation {
        return false;
    }
    expected.strandedness == ReadStrandedness::U || expected.strandedness == observed.strandedness
}

/// Whether a single-end or orphan mapping (mapping forward = `is_forward`, with
/// the given mate status) is compatible with the `expected` format. Ports
/// salmon's `compatibleHit(expected, start, isForward, ms)`.
///
/// With only one end observed there is no pair orientation to compare, so
/// compatibility reduces to: does this read's strand match what the protocol
/// promises *for that mate*? Note read 1 and read 2 have mirrored expectations,
/// which is why the two orphan cases below swap `SA` and `AS`.
pub fn compatible_single(expected: LibraryFormat, is_forward: bool, status: MateStatus) -> bool {
    use ReadStrandedness::*;
    let es = expected.strandedness;
    let same = expected.orientation == ReadOrientation::Same;
    match status {
        // Single-end: sense expects forward, antisense expects reverse,
        // unstranded accepts either.
        MateStatus::SingleEnd => {
            if is_forward {
                es == U || es == S
            } else {
                es == U || es == A
            }
        }
        // Left (read 1) orphan.
        MateStatus::PairedEndLeft => {
            if same {
                // `M` protocols use the single-ended S/A vocabulary, since both
                // mates share a strand.
                es == U || (es == S && is_forward) || (es == A && !is_forward)
            } else if is_forward {
                // Forward read 1 is what `SA` ("read 1 sense") promises.
                es == U || es == SA
            } else {
                es == U || es == AS
            }
        }
        // Right (read 2) orphan: the mirror image of the left case.
        MateStatus::PairedEndRight => {
            if same {
                es == U || (es == S && is_forward) || (es == A && !is_forward)
            } else if is_forward {
                // A forward read 2 means read 1 was antisense, i.e. `AS`.
                es == U || es == AS
            } else {
                es == U || es == SA
            }
        }
        // a properly-paired mapping should use `compatible_paired`
        MateStatus::PairedEndPaired => true,
    }
}

/// The observed library format of a properly-paired mapping, from the two mates'
/// reference strands. Opposite strands (`is_fw != mate_fw`) are inward
/// (`Toward`); same strands are `Same`. The strandedness follows read-1's strand
/// (`SA`/`AS` for inward, `S`/`A` for same). Shared by fragment-length-derivation
/// and naive-eq orientation tagging so they classify orientation identically.
///
/// This is the "what did we actually see" side of the comparison, versus the
/// "what did the user promise" side held in the expected format.
///
/// Note that `Away` is never produced: outward-facing pairs are indistinguishable
/// from inward ones by strand alone (the difference is which mate is leftmost),
/// so this function reports the inward reading and callers that care about `O`
/// protocols compare positions themselves.
pub fn observed_paired_format(is_fw: bool, mate_fw: bool) -> LibraryFormat {
    let (orientation, strandedness) = if is_fw != mate_fw {
        (
            ReadOrientation::Toward,
            if is_fw {
                ReadStrandedness::SA
            } else {
                ReadStrandedness::AS
            },
        )
    } else {
        (
            ReadOrientation::Same,
            if is_fw {
                ReadStrandedness::S
            } else {
                ReadStrandedness::A
            },
        )
    };
    LibraryFormat::new(ReadType::PairedEnd, orientation, strandedness)
}

/// Unified compatibility check: dispatches to [`compatible_paired`] for a
/// properly-paired mapping (using its `observed` format) and to
/// [`compatible_single`] otherwise.
///
/// `observed == None` on a proper pair means the caller did not record an
/// orientation; rather than guess, the mapping is accepted.
pub fn is_compatible(
    expected: LibraryFormat,
    observed: Option<LibraryFormat>,
    is_forward: bool,
    status: MateStatus,
) -> bool {
    if status == MateStatus::PairedEndPaired {
        match observed {
            Some(o) => compatible_paired(expected, o),
            None => true,
        }
    } else {
        compatible_single(expected, is_forward, status)
    }
}

/// [`is_compatible`], but for *accounting only*: a proper pair with no recorded
/// `observed` format is judged on its two mates' strands instead of being
/// accepted.
///
/// # Why this is a separate function
///
/// Sketch mode does not fill in a paired mapping's `observed` format — the two
/// mates' orientations are known, but nothing assembles them into a
/// [`LibraryFormat`]. [`is_compatible`] therefore accepts every sketch pair,
/// which is deliberate: making the *quantification* filter judge them would
/// change which fragments are assigned and what a sketch run reports as
/// abundance, not merely how it is described.
///
/// Counting is a different question. "How many fragments landed on the wrong
/// strand" is answerable from `is_forward`/`mate_forward` alone, and answering
/// it changes nothing about the run — so `lib_format_counts.json` can report a
/// real number in sketch mode while the filter keeps its behaviour. Using this
/// function anywhere that decides what is quantified would silently make that
/// change; use [`is_compatible`] there.
///
/// The derivation is [`observed_paired_format`], the same one the
/// deterministic pass already uses to infer a sketch run's library type.
pub fn is_compatible_for_counting(
    expected: LibraryFormat,
    observed: Option<LibraryFormat>,
    is_forward: bool,
    mate_forward: bool,
    status: MateStatus,
) -> bool {
    if status == MateStatus::PairedEndPaired {
        let o = observed.unwrap_or_else(|| observed_paired_format(is_forward, mate_forward));
        compatible_paired(expected, o)
    } else {
        compatible_single(expected, is_forward, status)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The counting variant must judge a pair that carries no observed format —
    /// which is every paired mapping in sketch mode — instead of accepting it.
    #[test]
    fn counting_judges_a_pair_with_no_recorded_format() {
        let isf = LibraryFormat::parse("ISF").unwrap();
        let isr = LibraryFormat::parse("ISR").unwrap();
        // Inward, read 1 forward: an `ISF` fragment.
        let (is_fw, mate_fw) = (true, false);
        let st = MateStatus::PairedEndPaired;

        // The quantification filter accepts it under either expectation, which
        // is the behaviour sketch mode relies on and must keep.
        assert!(is_compatible(isf, None, is_fw, st));
        assert!(is_compatible(isr, None, is_fw, st));

        // The counting variant derives the orientation and tells them apart.
        assert!(is_compatible_for_counting(isf, None, is_fw, mate_fw, st));
        assert!(!is_compatible_for_counting(isr, None, is_fw, mate_fw, st));
        // The mirrored fragment flips both verdicts.
        assert!(!is_compatible_for_counting(isf, None, false, true, st));
        assert!(is_compatible_for_counting(isr, None, false, true, st));
    }

    /// With a recorded format, or on anything that is not a proper pair, the two
    /// must agree — the counting variant only ever fills in a missing format.
    #[test]
    fn counting_matches_the_filter_wherever_the_format_is_known() {
        let formats: Vec<LibraryFormat> = ["IU", "ISF", "ISR", "OSF", "MSR", "U", "SF", "SR"]
            .iter()
            .map(|s| LibraryFormat::parse(s).unwrap())
            .collect();
        for &exp in &formats {
            for &obs in &formats {
                for is_fw in [true, false] {
                    for mate_fw in [true, false] {
                        for st in [
                            MateStatus::PairedEndPaired,
                            MateStatus::PairedEndLeft,
                            MateStatus::PairedEndRight,
                            MateStatus::SingleEnd,
                        ] {
                            assert_eq!(
                                is_compatible(exp, Some(obs), is_fw, st),
                                is_compatible_for_counting(exp, Some(obs), is_fw, mate_fw, st),
                                "diverged for expected={} observed={} is_fw={is_fw} status={st:?}",
                                exp.canonical(),
                                obs.canonical()
                            );
                            if st != MateStatus::PairedEndPaired {
                                assert_eq!(
                                    is_compatible(exp, None, is_fw, st),
                                    is_compatible_for_counting(exp, None, is_fw, mate_fw, st),
                                    "diverged on a non-pair for expected={}",
                                    exp.canonical()
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    /// `parse` and `canonical` must be exact inverses for all 12 formats,
    /// otherwise a format could change identity on a round trip through text.
    #[test]
    fn roundtrip_all_formats() {
        for s in [
            "IU", "ISF", "ISR", "OU", "OSF", "OSR", "MU", "MSF", "MSR", "U", "SF", "SR",
        ] {
            let fmt = LibraryFormat::parse(s).unwrap();
            assert_eq!(fmt.canonical(), s, "canonical mismatch for {s}");
        }
    }

    /// Users type `-l isr` as often as `-l ISR`.
    #[test]
    fn parse_is_case_insensitive() {
        assert_eq!(
            LibraryFormat::parse("isr").unwrap(),
            LibraryFormat::parse("ISR").unwrap()
        );
    }

    /// A typo must be an error, not a silently-accepted default.
    #[test]
    fn unknown_format_errors() {
        assert!(LibraryFormat::parse("ZZ").is_err());
    }

    /// `A` is the auto-detect request and must not be confused with a format.
    #[test]
    fn auto_detection_sentinel() {
        assert!(LibraryFormat::is_auto("A"));
        assert!(LibraryFormat::is_auto("a"));
        assert!(!LibraryFormat::is_auto("IU"));
    }

    /// The ids index arrays directly and are written into RAD files, so they
    /// must stay unique and gap-free.
    #[test]
    fn format_ids_are_unique_and_dense() {
        let formats = [
            "IU", "ISF", "ISR", "OU", "OSF", "OSR", "MU", "MSF", "MSR", "U", "SF", "SR",
        ];
        let mut ids: Vec<u8> = formats
            .iter()
            .map(|s| LibraryFormat::parse(s).unwrap().format_id())
            .collect();
        ids.sort_unstable();
        assert_eq!(ids, (0..12).collect::<Vec<_>>());
    }

    /// Pin the defaults, since they decide what an unspecified run assumes.
    #[test]
    fn defaults() {
        assert_eq!(LibraryFormat::paired_default().canonical(), "IU");
        assert_eq!(LibraryFormat::single_default().canonical(), "U");
    }

    /// The three rules of paired compatibility: unstranded accepts anything of
    /// its orientation, stranded accepts only its own strandedness, and a
    /// different orientation is always fatal.
    #[test]
    fn paired_compatibility() {
        let isf = LibraryFormat::parse("ISF").unwrap();
        let isr = LibraryFormat::parse("ISR").unwrap();
        let iu = LibraryFormat::parse("IU").unwrap();
        let osf = LibraryFormat::parse("OSF").unwrap();
        // unstranded expected accepts any same-orientation observation
        assert!(compatible_paired(iu, isf));
        assert!(compatible_paired(iu, isr));
        // stranded expected requires matching strandedness
        assert!(compatible_paired(isf, isf));
        assert!(!compatible_paired(isf, isr));
        // different orientation is always incompatible
        assert!(!compatible_paired(isf, osf));
    }

    /// Single-end compatibility is purely a strand test.
    #[test]
    fn single_compatibility() {
        let sf = LibraryFormat::parse("SF").unwrap(); // single sense
        let sr = LibraryFormat::parse("SR").unwrap(); // single antisense
        let u = LibraryFormat::parse("U").unwrap();
        // SF expects forward reads
        assert!(compatible_single(sf, true, MateStatus::SingleEnd));
        assert!(!compatible_single(sf, false, MateStatus::SingleEnd));
        // SR expects reverse reads
        assert!(compatible_single(sr, false, MateStatus::SingleEnd));
        assert!(!compatible_single(sr, true, MateStatus::SingleEnd));
        // U accepts both
        assert!(compatible_single(u, true, MateStatus::SingleEnd));
        assert!(compatible_single(u, false, MateStatus::SingleEnd));
    }
}
