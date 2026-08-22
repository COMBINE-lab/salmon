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

/// Number of distinct library formats ([`LibraryFormat::format_id`] is dense in
/// `0..NUM_LIB_FORMATS`). Sizes the observed-format histogram every quant path
/// tallies for `lib_format_counts.json`.
pub const NUM_LIB_FORMATS: usize = LibraryFormat::MAX_FORMAT_ID as usize + 1;

/// Per-fragment observed-format counts: `counts[f.format_id()]` is the number of
/// mapped fragments with at least one placement observed as format `f`. A
/// fragment with placements of several formats counts once under each (ports
/// salmon's per-fragment `libTypeCountsPerFrag[i] > 0` accumulation), so the
/// histogram can sum to more than the mapped total.
pub type LibFormatCountsArray = [u64; NUM_LIB_FORMATS];

/// The observed library format of a single-end mapping or an orphaned mate, for
/// the observed-format tally: a single-end observation on its strand. The pair
/// counterpart is [`observed_paired_format`].
pub fn observed_single_format(is_forward: bool) -> LibraryFormat {
    LibraryFormat::new(
        ReadType::SingleEnd,
        ReadOrientation::None,
        if is_forward {
            ReadStrandedness::S
        } else {
            ReadStrandedness::A
        },
    )
}

/// `lib_format_counts.json` fields derived from the observed-format histogram.
/// See [`summarize_lib_format_counts`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LibFormatSummary {
    /// Fragments whose observed format agrees with the expected one — exactly
    /// the expected format for a stranded expectation, either strandedness of
    /// the expected orientation for an unstranded one.
    pub num_concordant_consistent: u64,
    /// Fragments observed only in formats that disagree (wrong orientation,
    /// wrong strand under a stranded expectation, or orphaned mates).
    pub num_inconsistent_or_orphan: u64,
    /// Of the two strandedness variants of the expected orientation, the
    /// fraction observed as the *sense* variant. Near 0.5 for a genuinely
    /// unstranded library; near 0 or 1 for a stranded one.
    pub strand_mapping_bias: f64,
}

/// Derive the histogram-summary fields of `lib_format_counts.json` from the
/// observed-format counts and the (resolved) expected format. Ports the
/// derivation in salmon's `ReadExperiment::summarizeLibraryTypeCounts`, so the
/// fields keep their C++ meaning: *concordant and consistent* is agreement of
/// observed format with expected — not merely "both mates placed" — and
/// `strand_mapping_bias` is measured, which is what makes the classic
/// "potential strand bias in an unstranded protocol" check possible.
pub fn summarize_lib_format_counts(
    expected: LibraryFormat,
    counts: &LibFormatCountsArray,
) -> LibFormatSummary {
    // The two strandedness variants of the expected orientation: sense-first
    // and antisense-first (e.g. ISF/ISR for an inward expectation, SF/SR for
    // single-end or matching).
    let (s1, s2) = match expected.orientation {
        ReadOrientation::Same | ReadOrientation::None => (ReadStrandedness::S, ReadStrandedness::A),
        ReadOrientation::Away | ReadOrientation::Toward => {
            (ReadStrandedness::SA, ReadStrandedness::AS)
        }
    };
    let fmt1 = LibraryFormat::new(expected.read_type, expected.orientation, s1);
    let fmt2 = LibraryFormat::new(expected.read_type, expected.orientation, s2);
    let n1 = counts[fmt1.format_id() as usize];
    let n2 = counts[fmt2.format_id() as usize];

    let (agree, disagree) = if expected.strandedness == ReadStrandedness::U {
        // Unstranded: both strandedness variants of the orientation agree.
        let agree = n1 + n2;
        let disagree: u64 = counts
            .iter()
            .enumerate()
            .filter(|&(i, _)| i != fmt1.format_id() as usize && i != fmt2.format_id() as usize)
            .map(|(_, &c)| c)
            .sum();
        (agree, disagree)
    } else {
        // Stranded: only the expected format itself agrees.
        let eid = expected.format_id() as usize;
        let agree = counts[eid];
        let disagree: u64 = counts
            .iter()
            .enumerate()
            .filter(|&(i, _)| i != eid)
            .map(|(_, &c)| c)
            .sum();
        (agree, disagree)
    };
    let strand_mapping_bias = if agree > 0 && n1 + n2 > 0 {
        n1 as f64 / (n1 + n2) as f64
    } else {
        0.0
    };
    LibFormatSummary {
        num_concordant_consistent: agree,
        num_inconsistent_or_orphan: disagree,
        strand_mapping_bias,
    }
}

impl LibFormatSummary {
    /// The classic end-of-run library-format warnings (ported from
    /// `summarizeLibraryTypeCounts`) as structured [`crate::Diagnostic`]s, so
    /// they reach `meta_info.json`'s machine-readable `diagnostics` array and
    /// not just whoever reads the log (#1140): no concordant-consistent
    /// mappings at all; a strand bias beyond 1% under an unstranded
    /// expectation; and more than 5% of judged fragments incompatible with the
    /// expected format.
    pub fn diagnostics(
        &self,
        expected: LibraryFormat,
        compatible_fragment_ratio: f64,
    ) -> Vec<crate::Diagnostic> {
        let mut out = Vec::new();
        if expected.strandedness == ReadStrandedness::U {
            if self.num_concordant_consistent == 0 {
                out.push(crate::Diagnostic::new(
                    "no_concordant_mappings",
                    "warning",
                    "found no concordant and consistent mappings; if this is a \
                     paired-end library, are you sure the reads are properly \
                     paired? See lib_format_counts.json for details"
                        .to_string(),
                ));
            } else if (self.strand_mapping_bias - 0.5).abs() > 0.01 {
                out.push(crate::Diagnostic::new(
                    "strand_bias_unstranded",
                    "warning",
                    format!(
                        "detected a *potential* strand bias > 1% in an unstranded \
                         protocol (strand_mapping_bias = {:.4}); see \
                         lib_format_counts.json for details",
                        self.strand_mapping_bias
                    ),
                ));
            }
        }
        if 1.0 - compatible_fragment_ratio > 0.05 {
            out.push(crate::Diagnostic::new(
                "high_incompatible_fraction",
                "warning",
                format!(
                    "greater than 5% of the fragments disagreed with the expected \
                     library type ({}); see lib_format_counts.json for details",
                    expected.canonical()
                ),
            ));
        }
        out
    }

    /// Log [`Self::diagnostics`] as warnings — one condition set, spoken and
    /// recorded from the same source so the two cannot drift.
    pub fn log_warnings(&self, expected: LibraryFormat, compatible_fragment_ratio: f64) {
        for d in self.diagnostics(expected, compatible_fragment_ratio) {
            tracing::warn!("{}", d.message);
        }
    }
}

/// The full contents of `lib_format_counts.json`, shared by every writer (reads
/// mode, alignment mode, RAD requant) so the file cannot drift between paths.
/// Field names and meanings match C++ salmon's `summarizeLibraryTypeCounts`,
/// with `num_incompatible_fragments` as the one (documented) extension (#1130);
/// the trailing 12 keys are the raw observed-format histogram C++ appends.
#[derive(Debug, Serialize)]
pub struct LibFormatCountsFile {
    pub read_files: Vec<String>,
    pub expected_format: String,
    pub compatible_fragment_ratio: f64,
    pub num_compatible_fragments: u64,
    pub num_incompatible_fragments: u64,
    pub num_assigned_fragments: u64,
    pub num_frags_with_concordant_consistent_mappings: u64,
    pub num_frags_with_inconsistent_or_orphan_mappings: u64,
    pub strand_mapping_bias: f64,
    #[serde(rename = "IU")]
    pub iu: u64,
    #[serde(rename = "ISF")]
    pub isf: u64,
    #[serde(rename = "ISR")]
    pub isr: u64,
    #[serde(rename = "OU")]
    pub ou: u64,
    #[serde(rename = "OSF")]
    pub osf: u64,
    #[serde(rename = "OSR")]
    pub osr: u64,
    #[serde(rename = "MU")]
    pub mu: u64,
    #[serde(rename = "MSF")]
    pub msf: u64,
    #[serde(rename = "MSR")]
    pub msr: u64,
    #[serde(rename = "U")]
    pub u: u64,
    #[serde(rename = "SF")]
    pub sf: u64,
    #[serde(rename = "SR")]
    pub sr: u64,
}

impl LibFormatCountsFile {
    /// Assemble the file from a run's end-of-pass aggregates.
    ///
    /// `expected_format` is the *resolved* library-type string (the detected
    /// type under `-l A` when detection resolved, else the user's). The ratio is
    /// over the fragments the strand filter actually judged
    /// (`num_compatible + num_incompatible`), with the historical
    /// `num_mapped` / ratio-1.0 fallback when nothing ever was — under `-l A`
    /// the fragments consumed before detection resolves are in neither tally,
    /// and an unstranded type judges every fragment compatible. The
    /// histogram-derived fields need a concrete expected format; when the
    /// string does not name one (an unresolved `-l A` on an empty run) they are
    /// zero.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        read_files: Vec<String>,
        expected_format: String,
        num_mapped: u64,
        num_compatible: u64,
        num_incompatible: u64,
        counts: &LibFormatCountsArray,
    ) -> Self {
        let judged = num_compatible + num_incompatible;
        let ratio = if judged > 0 {
            num_compatible as f64 / judged as f64
        } else {
            1.0
        };
        let summary = LibraryFormat::parse(&expected_format)
            .ok()
            .map(|exp| summarize_lib_format_counts(exp, counts))
            .unwrap_or(LibFormatSummary {
                num_concordant_consistent: 0,
                num_inconsistent_or_orphan: 0,
                strand_mapping_bias: 0.0,
            });
        let c = |s: &str| counts[LibraryFormat::parse(s).unwrap().format_id() as usize];
        Self {
            read_files,
            expected_format,
            compatible_fragment_ratio: ratio,
            num_compatible_fragments: if judged > 0 {
                num_compatible
            } else {
                num_mapped
            },
            num_incompatible_fragments: num_incompatible,
            num_assigned_fragments: num_mapped,
            num_frags_with_concordant_consistent_mappings: summary.num_concordant_consistent,
            num_frags_with_inconsistent_or_orphan_mappings: summary.num_inconsistent_or_orphan,
            strand_mapping_bias: summary.strand_mapping_bias,
            iu: c("IU"),
            isf: c("ISF"),
            isr: c("ISR"),
            ou: c("OU"),
            osf: c("OSF"),
            osr: c("OSR"),
            mu: c("MU"),
            msf: c("MSF"),
            msr: c("MSR"),
            u: c("U"),
            sf: c("SF"),
            sr: c("SR"),
        }
    }

    /// The classic end-of-run warnings for this file's contents, as structured
    /// diagnostics (see [`LibFormatSummary::diagnostics`]). Empty when the
    /// expected format never resolved, or when nothing was assigned *and*
    /// nothing was judged incompatible — a truly empty run has nothing to warn
    /// about, but a run whose every fragment was dropped as wrong-strand
    /// (assigned 0, incompatible many) is exactly the run these warnings exist
    /// for.
    pub fn diagnostics(&self) -> Vec<crate::Diagnostic> {
        if self.num_assigned_fragments == 0 && self.num_incompatible_fragments == 0 {
            return Vec::new();
        }
        let Ok(exp) = LibraryFormat::parse(&self.expected_format) else {
            return Vec::new();
        };
        let summary = LibFormatSummary {
            num_concordant_consistent: self.num_frags_with_concordant_consistent_mappings,
            num_inconsistent_or_orphan: self.num_frags_with_inconsistent_or_orphan_mappings,
            strand_mapping_bias: self.strand_mapping_bias,
        };
        summary.diagnostics(exp, self.compatible_fragment_ratio)
    }

    /// Log [`Self::diagnostics`] as warnings; the writers also merge the same
    /// diagnostics into `meta_info.json`, so the log and the report agree.
    pub fn log_warnings(&self) {
        for d in self.diagnostics() {
            tracing::warn!("{}", d.message);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

    /// The histogram summary must keep the C++ field meanings: agreement is of
    /// observed format vs expected, and the bias is the sense fraction within
    /// the expected orientation's two variants.
    #[test]
    fn lib_format_summary_matches_cpp_derivation() {
        let mut counts: LibFormatCountsArray = [0; NUM_LIB_FORMATS];
        let isf = LibraryFormat::parse("ISF").unwrap();
        let isr = LibraryFormat::parse("ISR").unwrap();
        let sf = LibraryFormat::parse("SF").unwrap();
        counts[isf.format_id() as usize] = 30;
        counts[isr.format_id() as usize] = 20;
        counts[sf.format_id() as usize] = 5; // orphans, observed single-end

        // Unstranded inward expectation: both inward variants agree; the
        // orphan observations disagree; bias is the sense share.
        let iu = LibraryFormat::parse("IU").unwrap();
        let s = summarize_lib_format_counts(iu, &counts);
        assert_eq!(s.num_concordant_consistent, 50);
        assert_eq!(s.num_inconsistent_or_orphan, 5);
        assert!((s.strand_mapping_bias - 0.6).abs() < 1e-12);

        // Stranded expectation: only the exact format agrees; everything else
        // (wrong strand and orphans) disagrees. Bias is unchanged — it is a
        // property of the orientation's two variants, not of the expectation.
        let s = summarize_lib_format_counts(isf, &counts);
        assert_eq!(s.num_concordant_consistent, 30);
        assert_eq!(s.num_inconsistent_or_orphan, 25);
        assert!((s.strand_mapping_bias - 0.6).abs() < 1e-12);

        // Single-end expectation uses the S/A variants.
        let mut se: LibFormatCountsArray = [0; NUM_LIB_FORMATS];
        se[sf.format_id() as usize] = 8;
        se[LibraryFormat::parse("SR").unwrap().format_id() as usize] = 2;
        let u = LibraryFormat::parse("U").unwrap();
        let s = summarize_lib_format_counts(u, &se);
        assert_eq!(s.num_concordant_consistent, 10);
        assert_eq!(s.num_inconsistent_or_orphan, 0);
        assert!((s.strand_mapping_bias - 0.8).abs() < 1e-12);

        // Nothing observed at all: everything zero, bias 0 by convention.
        let empty: LibFormatCountsArray = [0; NUM_LIB_FORMATS];
        let s = summarize_lib_format_counts(iu, &empty);
        assert_eq!(s.num_concordant_consistent, 0);
        assert_eq!(s.num_inconsistent_or_orphan, 0);
        assert_eq!(s.strand_mapping_bias, 0.0);
    }

    /// Orphan/single-end observations are a strand on a single end.
    #[test]
    fn observed_single_format_is_se_strand() {
        assert_eq!(observed_single_format(true).canonical(), "SF");
        assert_eq!(observed_single_format(false).canonical(), "SR");
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
