//! RAD-format schema, record types, and writer for salmon (libradicl-backed).
//!
//! salmon can optionally emit its per-fragment mappings as a RAD file, in one of
//! two profiles:
//!   * [`RadProfile::Sketch`] — equivalent to piscem `map-bulk` (`bulk_with_pos`):
//!     per-hit `(tid, orientation, pos, frag_len)`, no alignment scores, decoys
//!     excluded.
//!   * [`RadProfile::SelectiveAlignment`] — the same, plus a per-hit alignment
//!     score so the placements can be re-weighted on requant exactly like
//!     internal selective alignments.
//!
//! Determinism is achieved structurally (a fixed FLD before order-independent
//! equivalence-class assembly), so no per-read name hash is stored; a salmon RAD
//! is identified instead by its `rad_type` file tag.
//!
//! The per-hit `(compressed_ori_ref, pos, frag_len)` triple is byte-identical to
//! piscem's `bulk_with_pos` layout (`tid` in the low 30 bits, `is_fw` in bit 31,
//! `mate_fw` in bit 30), which is what libradicl's reader decodes. (Note:
//! libradicl 0.13's `PiscemBulkReadRecord::write` uses a different, inconsistent
//! orientation encoding; salmon's writer matches the *reader* / piscem-rs.)

use salmon_core::mate::MateStatus;

pub mod record;
pub mod schema;
pub mod writer;

/// Re-exported so consumers can select a RAD chunk compression codec without a
/// direct libradicl dependency.
pub use libradicl::ChunkCodec;
pub use record::{SalmonBulkContext, SalmonBulkRecord};
pub use writer::{FragmentChunkBuf, RadOutputWriter};

/// The profile of a RAD file presented for *input*, as detected from its
/// prelude tag signature (see [`detect_input_profile`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RadInputProfile {
    /// piscem `map-bulk` (`bulk_with_pos`): `frag_map_type` read tag + the
    /// `(compressed_ori_ref, pos, frag_len)` aln triple, no scores.
    PiscemBulk,
    /// salmon-written RAD (carries the `rad_type` file tag); the inner
    /// [`RadProfile`] distinguishes sketch vs. selective-alignment.
    Salmon(RadProfile),
}

/// Detect the input profile of a RAD file from its prelude's tag signature.
///
/// libradicl has no single "kind" field, so detection is by which tags are
/// present: salmon RAD is identified by the `rad_type` file tag (and
/// selective-alignment by the `alignment_score` aln tag); a piscem bulk file is
/// identified by `frag_map_type` + the `compressed_ori_ref`/`pos`/`frag_len`
/// alignment triple.
pub fn detect_input_profile(
    prelude: &libradicl::header::RadPrelude,
) -> anyhow::Result<RadInputProfile> {
    let ft = &prelude.file_tags;
    let rt = &prelude.read_tags;
    let at = &prelude.aln_tags;
    // salmon RAD is identified by its `rad_type` file tag; piscem uses a
    // differently-named `known_rad_type` tag, so this cleanly distinguishes them.
    if ft.has_tag("rad_type") {
        let p = if at.has_tag("alignment_score") {
            RadProfile::SelectiveAlignment
        } else {
            RadProfile::Sketch
        };
        return Ok(RadInputProfile::Salmon(p));
    }
    if rt.has_tag("frag_map_type")
        && at.has_tag("compressed_ori_ref")
        && at.has_tag("pos")
        && at.has_tag("frag_len")
    {
        return Ok(RadInputProfile::PiscemBulk);
    }
    anyhow::bail!(
        "unrecognized RAD profile: expected a salmon RAD (with a `rad_type` file tag) \
         or a piscem bulk RAD (with `frag_map_type` + `compressed_ori_ref`/`pos`/`frag_len`)"
    )
}

/// Which flavor of salmon RAD this is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RadProfile {
    /// piscem-`map-bulk`-compatible mappings, no scores.
    Sketch,
    /// mappings carrying a per-hit alignment score.
    SelectiveAlignment,
}

impl RadProfile {
    /// The `rad_type` file-tag string salmon writes for this profile.
    pub fn rad_type_str(self) -> &'static str {
        match self {
            RadProfile::Sketch => "salmon_bulk_sketch",
            RadProfile::SelectiveAlignment => "salmon_bulk_sa",
        }
    }

    /// True if records in this profile carry a per-hit `alignment_score` tag.
    pub fn has_scores(self) -> bool {
        matches!(self, RadProfile::SelectiveAlignment)
    }
}

/// Orientation/reference bit layout of `compressed_ori_ref`, matching piscem
/// (`piscem-rs::io::rad::write_bulk_record`) and libradicl's bulk *reader*.
pub const ORI_FW: u32 = 0x8000_0000;
/// Mate-forward bit.
pub const MATE_FW: u32 = 0x4000_0000;
/// Transcript-id mask (low 30 bits).
pub const TID_MASK: u32 = 0x3FFF_FFFF;

/// Sentinel `frag_len` for non-proper-pair hits (orphan / single-end), matching
/// piscem.
pub const FRAG_LEN_UNPAIRED: u16 = u16::MAX;

/// File-tag name: a `u8` bitfield recording which reserved values were actually
/// baked at finalize ([`BAKED_FLD`] / [`BAKED_ABUND`]); `0` ⇒ none (placeholders).
pub const BAKED_FLAGS_TAG: &str = "baked_flags";
/// File-tag name: the fragment-length distribution as a log-PMF over raw lengths
/// `[0, fld_len)`, baked at write time so a reader can quantify in a single pass
/// with exact FLD parity instead of deriving it.
pub const FRAG_LENGTH_DIST_TAG: &str = "frag_length_dist";
/// File-tag name: initial per-reference abundance estimates (one `f64` each),
/// baked when the write run quantified; a prior for future bias-aware requant.
pub const INITIAL_ABUNDANCES_TAG: &str = "initial_abundances";
/// File-tag name: the resolved library format as a `u8` `format_id`
/// ([`salmon_core::LibraryFormat::format_id`]), baked so a reader can apply
/// concordance filtering under `-l A` without re-inferring the type (the write
/// run already detected it). Absent ⇒ the reader auto-detects from the RAD.
pub const LIBRARY_FORMAT_TAG: &str = "library_format";
/// [`BAKED_FLAGS_TAG`] bit: a fragment-length distribution is present.
pub const BAKED_FLD: u8 = 0x01;
/// [`BAKED_FLAGS_TAG`] bit: initial abundance estimates are present.
pub const BAKED_ABUND: u8 = 0x02;
/// [`BAKED_FLAGS_TAG`] bit: a resolved library format is present.
pub const BAKED_LIBFMT: u8 = 0x04;
/// [`BAKED_FLAGS_TAG`] bit: a non-default [`SCORE_KIND_TAG`] is present (the
/// per-hit score is a quantized log-weight, not a raw aligner `AS`).
pub const BAKED_SCORE_KIND: u8 = 0x08;

/// File-tag name: how a reader should interpret the per-hit `alignment_score`
/// (present only in the selective-alignment profile). A `u8`:
/// [`SCORE_KIND_AS`] (`0`, default/absent) ⇒ the raw aligner score, weighted as
/// `exp(-scoreExp·(best−score))`; [`SCORE_KIND_LOGWEIGHT`] (`1`) ⇒ a quantized
/// alignment-model log-likelihood ratio `round(Σ(fg−bg)·`[`SCORE_LOG_SCALE`]`)`,
/// used directly (÷ scale) as the eq-class log-weight basis. Absent ⇒ `0`.
pub const SCORE_KIND_TAG: &str = "score_kind";
/// [`SCORE_KIND_TAG`] value: raw aligner score (`AS`), soft-weighted by score
/// difference. The default when the tag is absent (all pre-existing salmon RAD).
pub const SCORE_KIND_AS: u8 = 0;
/// [`SCORE_KIND_TAG`] value: a quantized alignment-model log-likelihood ratio
/// (deterministic error model); the stored `i32` is `Σ(fg−bg)·`[`SCORE_LOG_SCALE`].
pub const SCORE_KIND_LOGWEIGHT: u8 = 1;
/// Fixed-point scale for a [`SCORE_KIND_LOGWEIGHT`] score: the per-hit log-LR is
/// stored as `round(logLR · SCORE_LOG_SCALE)` and read back as `score / scale`.
/// 1000 keeps ~3 decimal digits of the log-weight, well within an `i32`.
pub const SCORE_LOG_SCALE: f64 = 1000.0;

/// piscem `frag_map_type` (a.k.a. `MappingType`) code for an unmapped fragment.
pub const FRAG_MAP_TYPE_UNMAPPED: u8 = 0;

/// A single salmon mapping placement, in the form written to / read from RAD.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RadHit {
    /// transcript id (must fit in 30 bits).
    pub tid: u32,
    /// read-1 / fragment-forward orientation on the reference.
    pub is_fw: bool,
    /// mate (read-2) forward orientation; meaningful only for proper pairs.
    pub mate_fw: bool,
    /// leftmost reference position of the fragment.
    pub pos: u32,
    /// fragment length, or [`FRAG_LEN_UNPAIRED`] for orphan / single-end.
    pub frag_len: u16,
    /// alignment score (only meaningful / written for [`RadProfile::SelectiveAlignment`]).
    pub score: i32,
}

impl RadHit {
    /// Pack `(tid, is_fw, mate_fw)` into the `compressed_ori_ref` u32, matching
    /// piscem's bulk layout and libradicl's bulk reader.
    #[inline]
    pub fn compressed_ori_ref(&self) -> u32 {
        (self.tid & TID_MASK)
            | if self.is_fw { ORI_FW } else { 0 }
            | if self.mate_fw { MATE_FW } else { 0 }
    }

    /// Decode a `compressed_ori_ref` u32 into `(tid, is_fw, mate_fw)`.
    #[inline]
    pub fn decode_ori_ref(v: u32) -> (u32, bool, bool) {
        (v & TID_MASK, (v & ORI_FW) != 0, (v & MATE_FW) != 0)
    }

    /// Build a hit for one placement following piscem's bulk `frag_len`
    /// convention, shared by every RAD producer (reads, alignment-BAM, genome
    /// projection): a proper pair stores its real `fragment_len`; an orphan /
    /// single-end stores the mapped mate's `read_len` in the slot (clamped below
    /// [`FRAG_LEN_UNPAIRED`]) so the reader can recover the bounded-CMF orphan
    /// fragment-length weight rather than a flat penalty. `pos` is clamped to
    /// `>= 0`. `mate_fw` is meaningful only for proper pairs (callers pass
    /// `false` otherwise).
    pub fn for_placement(
        tid: u32,
        is_fw: bool,
        mate_fw: bool,
        pos: i32,
        paired: bool,
        fragment_len: i32,
        read_len: i32,
        score: i32,
    ) -> Self {
        RadHit {
            tid,
            is_fw,
            mate_fw: paired && mate_fw,
            pos: pos.max(0) as u32,
            frag_len: if paired {
                fragment_len.clamp(0, u16::MAX as i32) as u16
            } else {
                read_len.clamp(0, (u16::MAX - 1) as i32) as u16
            },
            score,
        }
    }
}

/// piscem `MappingType` codes for the read-level `frag_map_type` tag.
///
/// salmon's [`MateStatus`] maps onto the subset piscem uses for bulk data.
pub mod frag_map_type {
    use super::MateStatus;

    /// both mates mapped as a proper pair.
    pub const MAPPED_PAIR: u8 = 4;
    /// only the left/first mate mapped (orphan).
    pub const LEFT_ORPHAN: u8 = 2;
    /// only the right/second mate mapped (orphan).
    pub const RIGHT_ORPHAN: u8 = 3;
    /// single-end mapping.
    pub const SINGLE: u8 = 1;
    /// no mapping.
    pub const UNMAPPED: u8 = 0;

    /// Per-hit `frag_map_type` for a single placement's [`MateStatus`].
    #[inline]
    pub fn from_mate_status(s: MateStatus) -> u8 {
        match s {
            MateStatus::PairedEndPaired => MAPPED_PAIR,
            MateStatus::PairedEndLeft => LEFT_ORPHAN,
            MateStatus::PairedEndRight => RIGHT_ORPHAN,
            MateStatus::SingleEnd => SINGLE,
        }
    }

    /// "Completeness" rank used to pick a single fragment-level `frag_map_type`
    /// when a fragment's placements disagree (a proper pair is more complete
    /// than an orphan, which is more complete than a single-end / unmapped).
    #[inline]
    fn rank(code: u8) -> u8 {
        match code {
            MAPPED_PAIR => 3,
            LEFT_ORPHAN | RIGHT_ORPHAN => 2,
            SINGLE => 1,
            _ => 0,
        }
    }

    /// Fragment-level `frag_map_type` = the most complete status across hits.
    pub fn fragment_level(statuses: impl IntoIterator<Item = MateStatus>) -> u8 {
        statuses
            .into_iter()
            .map(from_mate_status)
            .max_by_key(|&c| rank(c))
            .unwrap_or(UNMAPPED)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ori_ref_roundtrip() {
        for &(tid, fw, mfw) in &[
            (0u32, false, false),
            (123, true, false),
            (TID_MASK, true, true),
        ] {
            let h = RadHit {
                tid,
                is_fw: fw,
                mate_fw: mfw,
                pos: 7,
                frag_len: 11,
                score: -3,
            };
            let v = h.compressed_ori_ref();
            assert_eq!(RadHit::decode_ori_ref(v), (tid, fw, mfw));
        }
    }

    #[test]
    fn fragment_level_picks_most_complete() {
        use frag_map_type::*;
        assert_eq!(
            fragment_level([MateStatus::SingleEnd, MateStatus::PairedEndPaired]),
            MAPPED_PAIR
        );
        assert_eq!(
            fragment_level([MateStatus::PairedEndLeft, MateStatus::SingleEnd]),
            LEFT_ORPHAN
        );
        assert_eq!(fragment_level([MateStatus::SingleEnd]), SINGLE);
        assert_eq!(fragment_level([]), UNMAPPED);
    }

    #[test]
    fn for_placement_frag_len_convention() {
        // Proper pair: real fragment length, mate strand kept, pos clamped >= 0.
        let p = RadHit::for_placement(5, true, true, -2, true, 247, 100, -3);
        assert_eq!(
            p,
            RadHit {
                tid: 5,
                is_fw: true,
                mate_fw: true,
                pos: 0,
                frag_len: 247,
                score: -3
            }
        );
        // Orphan / single-end: the slot carries the mate read length (clamped
        // below the unpaired sentinel), and mate_fw is forced false.
        let o = RadHit::for_placement(5, false, true, 10, false, 247, 100, 7);
        assert_eq!(
            o,
            RadHit {
                tid: 5,
                is_fw: false,
                mate_fw: false,
                pos: 10,
                frag_len: 100,
                score: 7
            }
        );
        // Clamps: fragment length saturates at u16::MAX; orphan read length at MAX-1.
        assert_eq!(
            RadHit::for_placement(0, true, false, 0, true, 1 << 20, 0, 0).frag_len,
            u16::MAX
        );
        assert_eq!(
            RadHit::for_placement(0, true, false, 0, false, 0, 1 << 20, 0).frag_len,
            u16::MAX - 1
        );
    }
}
