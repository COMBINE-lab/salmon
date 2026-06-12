//! Library format types and parsing.
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
    /// unstranded
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

    /// Default paired-end format used by salmon (`IU`).
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
    pub fn parse(s: &str) -> Result<Self, SalmonError> {
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
    pub fn is_auto(s: &str) -> bool {
        s.eq_ignore_ascii_case("a")
    }

    /// The canonical string for this format, e.g. `"ISR"`.
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
            (SingleEnd, _, U) => "U",
            (SingleEnd, _, S) => "SF",
            (SingleEnd, _, A) => "SR",
            // Any non-canonical combination collapses to unstranded of its type.
            (PairedEnd, _, _) => "IU",
            (SingleEnd, _, _) => "U",
        }
    }

    /// Stable integer id for the 12 canonical formats. Useful as a dense index
    /// (e.g. for `lib_format_counts.json` accounting).
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
            _ => unreachable!("canonical() only returns the 12 known formats"),
        }
    }

    pub fn is_paired(&self) -> bool {
        matches!(self.read_type, ReadType::PairedEnd)
    }

    /// Largest valid [`format_id`](Self::format_id) (there are 12 formats, ids 0..=11).
    pub const MAX_FORMAT_ID: u8 = 11;

    /// Inverse of [`format_id`](Self::format_id). Mirrors `formatFromID`.
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

impl fmt::Display for LibraryFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.canonical())
    }
}

/// Whether a properly-paired mapping with the `observed` format is compatible
/// with the `expected` library format. Ports salmon's `compatibleHit(expected,
/// observed)`: orientations must match, and the expected strandedness must be
/// unstranded or equal to the observed.
pub fn compatible_paired(expected: LibraryFormat, observed: LibraryFormat) -> bool {
    if expected.orientation != observed.orientation {
        return false;
    }
    expected.strandedness == ReadStrandedness::U || expected.strandedness == observed.strandedness
}

/// Whether a single-end or orphan mapping (mapping forward = `is_forward`, with
/// the given mate status) is compatible with the `expected` format. Ports
/// salmon's `compatibleHit(expected, start, isForward, ms)`.
pub fn compatible_single(expected: LibraryFormat, is_forward: bool, status: MateStatus) -> bool {
    use ReadStrandedness::*;
    let es = expected.strandedness;
    let same = expected.orientation == ReadOrientation::Same;
    match status {
        MateStatus::SingleEnd => {
            if is_forward {
                es == U || es == S
            } else {
                es == U || es == A
            }
        }
        MateStatus::PairedEndLeft => {
            if same {
                es == U || (es == S && is_forward) || (es == A && !is_forward)
            } else if is_forward {
                es == U || es == SA
            } else {
                es == U || es == AS
            }
        }
        MateStatus::PairedEndRight => {
            if same {
                es == U || (es == S && is_forward) || (es == A && !is_forward)
            } else if is_forward {
                es == U || es == AS
            } else {
                es == U || es == SA
            }
        }
        // a properly-paired mapping should use `compatible_paired`
        MateStatus::PairedEndPaired => true,
    }
}

/// Unified compatibility check: dispatches to [`compatible_paired`] for a
/// properly-paired mapping (using its `observed` format) and to
/// [`compatible_single`] otherwise.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_all_formats() {
        for s in [
            "IU", "ISF", "ISR", "OU", "OSF", "OSR", "MU", "MSF", "MSR", "U", "SF", "SR",
        ] {
            let fmt = LibraryFormat::parse(s).unwrap();
            assert_eq!(fmt.canonical(), s, "canonical mismatch for {s}");
        }
    }

    #[test]
    fn parse_is_case_insensitive() {
        assert_eq!(
            LibraryFormat::parse("isr").unwrap(),
            LibraryFormat::parse("ISR").unwrap()
        );
    }

    #[test]
    fn unknown_format_errors() {
        assert!(LibraryFormat::parse("ZZ").is_err());
    }

    #[test]
    fn auto_detection_sentinel() {
        assert!(LibraryFormat::is_auto("A"));
        assert!(LibraryFormat::is_auto("a"));
        assert!(!LibraryFormat::is_auto("IU"));
    }

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

    #[test]
    fn defaults() {
        assert_eq!(LibraryFormat::paired_default().canonical(), "IU");
        assert_eq!(LibraryFormat::single_default().canonical(), "U");
    }

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
