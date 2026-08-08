//! RAD-format schema, record types, and writer for salmon (libradicl-backed).
//!
//! # What RAD is
//!
//! RAD ("Reduced Alignment Data") is a compact binary file format from the
//! COMBINE-lab tool family. Instead of a full SAM/BAM alignment record per read
//! — read name, sequence, qualities, CIGAR string, all of it — a RAD file stores
//! only what a quantifier actually needs: for each fragment, the list of places
//! it could have come from. That is roughly an order of magnitude smaller and
//! much faster to re-read.
//!
//! # Why salmon writes one
//!
//! Mapping reads is the expensive half of a salmon run; the statistical half
//! (EM) is comparatively cheap. Writing a RAD file lets you map once and then
//! re-quantify many times (different bias settings, different priors) without
//! touching the FASTQ again. It is also how `--deterministic` works: map to a
//! RAD, then quantify from it with a fixed fragment-length distribution so the
//! answer is byte-identical across runs and thread counts.
//!
//! # The two profiles
//!
//! salmon can emit its per-fragment mappings in one of two profiles:
//!   * [`RadProfile::Sketch`] — equivalent to piscem `map-bulk`
//!     (`bulk_with_pos`): per-hit `(tid, orientation, pos, frag_len)`, no
//!     alignment scores, decoys excluded.
//!   * [`RadProfile::SelectiveAlignment`] — the same, plus a per-hit alignment
//!     score so the placements can be re-weighted on requant exactly like
//!     internal selective alignments.
//!
//! Determinism is achieved structurally (a fixed FLD before order-independent
//! equivalence-class assembly), so no per-read name hash is stored; a salmon RAD
//! is identified instead by its `rad_type` file tag.
//!
//! # Binary-compatibility note
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
///
/// salmon can quantify from RAD files it wrote itself *and* from files piscem
/// wrote, so on the reading side it must first work out which it is looking at.
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
/// **How a RAD file describes itself.** Its header ("prelude") carries three
/// lists of named tags: file-level, per-read, and per-alignment. There is no
/// single "kind" field, so detection is by *which tags are present*.
///
/// salmon RAD is identified by the `rad_type` file tag (and selective-alignment
/// by the `alignment_score` aln tag); a piscem bulk file is identified by
/// `frag_map_type` + the `compressed_ori_ref`/`pos`/`frag_len` alignment triple.
pub fn detect_input_profile(
    prelude: &libradicl::header::RadPrelude,
) -> anyhow::Result<RadInputProfile> {
    let ft = &prelude.file_tags;
    let rt = &prelude.read_tags;
    let at = &prelude.aln_tags;
    // salmon RAD is identified by its `rad_type` file tag; piscem uses a
    // differently-named `known_rad_type` tag, so this cleanly distinguishes them.
    if ft.has_tag("rad_type") {
        // Within salmon's own files, the presence of per-hit scores is what
        // separates the two profiles.
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
    // Neither signature matched: fail loudly rather than guess, since guessing
    // wrong would decode the binary payload with the wrong layout.
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
    /// The `rad_type` file-tag string salmon writes for this profile. This
    /// string is what [`detect_input_profile`] later reads back.
    pub fn rad_type_str(self) -> &'static str {
        match self {
            RadProfile::Sketch => "salmon_bulk_sketch",
            RadProfile::SelectiveAlignment => "salmon_bulk_sa",
        }
    }

    /// True if records in this profile carry a per-hit `alignment_score` tag.
    pub fn has_scores(self) -> bool {
        // `matches!` is a compact "does this value have this shape" test.
        matches!(self, RadProfile::SelectiveAlignment)
    }
}

// ---------------------------------------------------------------------------
// Bit layout of `compressed_ori_ref`
//
// Three pieces of information are packed into one 32-bit integer to keep the
// per-hit record small (a big run writes billions of hits, so every byte
// counts):
//
//   bit 31        : is_fw    — fragment/read-1 aligned to the forward strand
//   bit 30        : mate_fw  — read-2 aligned to the forward strand
//   bits 29..0    : tid      — transcript id, up to ~1.07 billion references
//
// A "mask" is a constant with 1-bits exactly where a field lives, so `v & MASK`
// extracts the field and `v | BIT` sets a flag.
// ---------------------------------------------------------------------------

/// Orientation/reference bit layout of `compressed_ori_ref`, matching piscem
/// (`piscem-rs::io::rad::write_bulk_record`) and libradicl's bulk *reader*.
pub const ORI_FW: u32 = 0x8000_0000;
/// Mate-forward bit.
pub const MATE_FW: u32 = 0x4000_0000;
/// Transcript-id mask (low 30 bits).
pub const TID_MASK: u32 = 0x3FFF_FFFF;

/// Sentinel `frag_len` for non-proper-pair hits (orphan / single-end), matching
/// piscem. A sentinel is an otherwise-impossible value used to mean "not
/// applicable" without spending an extra field.
pub const FRAG_LEN_UNPAIRED: u16 = u16::MAX;

// ---------------------------------------------------------------------------
// File tags salmon bakes into its own RAD files.
//
// "Baking" means storing a value the writing run already computed so a later
// requant does not have to recompute (or worse, recompute differently). This is
// what makes a RAD requant reproduce the original run exactly.
// ---------------------------------------------------------------------------

/// File-tag name: a `u8` bitfield recording which reserved values were actually
/// baked at finalize ([`BAKED_FLD`] / [`BAKED_ABUND`]); `0` ⇒ none (placeholders).
///
/// The tag slots are written up front at a fixed size so they can be patched in
/// place at the end of the run; this bitfield says which slots hold real data.
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
// ---------------------------------------------------------------------------
// Mapping-pass provenance.
//
// The counters and index identity below are observations of the *mapping* pass,
// not of the fragments it emitted, so a reader cannot recover them from the
// records: a RAD holds only the fragments that mapped. Without them a requant
// can say how many fragments it quantified but not how many were seen, so it
// reports a 100% mapping rate by construction.
//
// Two groups, because they become known at different times. The index identity
// is known before the first record is written, so it is written directly. The
// counters are only final once the pass ends, so their slots are reserved and
// backpatched at finalize, exactly like the FLD (see [`BAKED_MAP_COUNTERS`]).
// ---------------------------------------------------------------------------

/// File-tag name: total fragments the mapping pass *observed*, mapped or not.
///
/// The one quantity a RAD reader can never derive, and the one every mapping-rate
/// figure depends on.
pub const NUM_PROCESSED_TAG: &str = "num_processed";
/// File-tag name: fragments whose mappings were dovetailed.
pub const NUM_DOVETAIL_TAG: &str = "num_dovetail_fragments";
/// File-tag name: fragments dropped by the mapping-score filter.
pub const NUM_FILTERED_VM_TAG: &str = "num_fragments_filtered_vm";
/// File-tag name: alignments scoring below threshold among mapped fragments.
pub const NUM_BELOW_THRESH_VM_TAG: &str = "num_alignments_below_threshold_vm";
/// File-tag name: fragments whose best mapping was to a decoy.
///
/// Decoy mappings are filtered before records are written, so this count is
/// unrecoverable from the file itself.
pub const NUM_DECOY_FRAGMENTS_TAG: &str = "num_decoy_fragments";
/// File-tag names identifying the index the mappings were made against.
///
/// Recorded from the mapping pass rather than read from whatever index the
/// requant happens to be given: these describe the index that actually produced
/// the records, which is what the provenance in `meta_info.json` claims.
pub const INDEX_SEQ_HASH_TAG: &str = "index_seq_hash";
/// See [`INDEX_SEQ_HASH_TAG`].
pub const INDEX_NAME_HASH_TAG: &str = "index_name_hash";
/// See [`INDEX_SEQ_HASH_TAG`].
pub const INDEX_SEQ_HASH512_TAG: &str = "index_seq_hash512";
/// See [`INDEX_SEQ_HASH_TAG`].
pub const INDEX_NAME_HASH512_TAG: &str = "index_name_hash512";
/// See [`INDEX_SEQ_HASH_TAG`].
pub const INDEX_DECOY_SEQ_HASH_TAG: &str = "index_decoy_seq_hash";
/// See [`INDEX_SEQ_HASH_TAG`].
pub const INDEX_DECOY_NAME_HASH_TAG: &str = "index_decoy_name_hash";
/// File-tag name: whether the index was built with `--keepDuplicates`.
///
/// Written only when the index actually recorded it; absent means unknown,
/// which a reader must not collapse to `false`.
pub const KEEP_DUPLICATES_TAG: &str = "keep_duplicates";
/// File-tag name: the `@PG` header records of the BAM these fragments came from,
/// verbatim, one per line.
///
/// Stored as the raw SAM lines rather than a parsed structure so nothing the
/// aligner recorded is dropped in transit; a reader parses them back out for
/// `meta_info.json`. Present only for a RAD derived from an alignment file.
pub const SOURCE_PROGRAMS_TAG: &str = "source_programs";

/// File-tag name: how the fragments in this RAD were produced.
///
/// Recorded explicitly rather than inferred from the record profile, because the
/// two do not agree: a BAM-derived RAD is written in the selective-alignment
/// profile but its fragments came from an aligner, not from salmon's mapper.
pub const MAPPING_TYPE_TAG: &str = "mapping_type";

/// How the fragments in a RAD were produced; `meta_info.json`'s `mapping_type`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MappingType {
    /// salmon's selective-alignment mapper
    Mapping,
    /// salmon's sketch (pseudoalignment) mapper
    Pseudo,
    /// an external aligner, via a BAM
    Alignment,
}

impl MappingType {
    /// The string salmon reports in `meta_info.json`.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Mapping => "mapping",
            Self::Pseudo => "pseudo",
            Self::Alignment => "alignment",
        }
    }

    /// Parse the tag written by [`Self::as_str`]; `None` for anything else.
    pub fn from_str_tag(s: &str) -> Option<Self> {
        match s {
            "mapping" => Some(Self::Mapping),
            "pseudo" => Some(Self::Pseudo),
            "alignment" => Some(Self::Alignment),
            _ => None,
        }
    }
}

/// What a writer records about the pass that produced a RAD.
///
/// Grouped rather than passed as loose arguments because the parts are not
/// independent: a BAM-derived RAD has a `mapping_type` but no index behind it,
/// and `None` there must reach the reader as "unknown" rather than as empty
/// strings that look like real hashes.
#[derive(Clone, Debug)]
pub struct WriterProvenance {
    /// how the fragments were produced
    pub mapping_type: MappingType,
    /// identity of the index the mappings were made against, when there is one
    pub index: Option<IndexProvenance>,
    /// verbatim `@PG` lines from the source BAM, when the fragments came from
    /// one; empty otherwise. Propagates what the aligner said about itself, so a
    /// requant reports how the alignments were produced and not merely that they
    /// were alignments.
    pub source_programs: Vec<String>,
}

/// Identity of the index a mapping pass ran against, recorded in the RAD.
///
/// Lets a requant report the same index provenance the original run did. It has
/// to be recorded rather than read back from the index at requant time, both
/// because the `--rad` path loads no index at all and because the user may hand
/// it a different one than the mappings were made against.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IndexProvenance {
    /// hash of the reference sequences
    pub seq_hash: String,
    /// hash of the reference names
    pub name_hash: String,
    /// 512-bit form of [`Self::seq_hash`]
    pub seq_hash512: String,
    /// 512-bit form of [`Self::name_hash`]
    pub name_hash512: String,
    /// hash of the decoy sequences (empty when the index has no decoys)
    pub decoy_seq_hash: String,
    /// hash of the decoy names (empty when the index has no decoys)
    pub decoy_name_hash: String,
    /// whether the index was built with `--keepDuplicates`; `None` when the
    /// index predates recording it, which must not be reported as `false`
    pub keep_duplicates: Option<bool>,
}

/// Counters describing what a mapping pass *saw*, as opposed to what it wrote.
///
/// Every field here is unrecoverable from the records: a RAD contains only the
/// fragments that mapped, so without these a reader cannot report a mapping rate
/// (see [`NUM_PROCESSED_TAG`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MapCounters {
    /// total fragments observed, mapped or not
    pub num_processed: u64,
    /// fragments whose mappings were dovetailed
    pub num_dovetail: u64,
    /// fragments dropped by the mapping-score filter
    pub num_filtered_vm: u64,
    /// alignments below threshold among mapped fragments
    pub num_below_threshold_vm: u64,
    /// fragments whose best mapping was to a decoy
    pub num_decoy_fragments: u64,
}

/// [`BAKED_FLAGS_TAG`] bit: a fragment-length distribution is present.
pub const BAKED_FLD: u8 = 0x01;
/// [`BAKED_FLAGS_TAG`] bit: initial abundance estimates are present.
pub const BAKED_ABUND: u8 = 0x02;
/// [`BAKED_FLAGS_TAG`] bit: a resolved library format is present.
pub const BAKED_LIBFMT: u8 = 0x04;
/// [`BAKED_FLAGS_TAG`] bit: a non-default [`SCORE_KIND_TAG`] is present (the
/// per-hit score is a quantized log-weight, not a raw aligner `AS`).
pub const BAKED_SCORE_KIND: u8 = 0x08;
/// [`BAKED_FLAGS_TAG`] bit: the writer reached finalize, so the file is complete.
///
/// This is **confirmatory only**. Its presence proves the file was fully
/// written; its absence proves nothing, because a conforming producer may stream
/// a RAD file without ever backpatching (the same reason `num_chunks == 0` is a
/// valid "unknown chunk count" signal). Readers may use it as a fast path, but
/// must never require it, or they would reject legitimate files from other
/// tools.
pub const WRITE_COMPLETE: u8 = 0x10;
/// [`BAKED_FLAGS_TAG`] bit: the mapping-pass counters hold real values.
///
/// Distinguishes "the pass observed zero of these" from "nobody filled the
/// slot", which a zero alone cannot: the slots are reserved as zeros.
pub const BAKED_MAP_COUNTERS: u8 = 0x20;
/// [`BAKED_FLAGS_TAG`] bit: the index-identity tags are present.
pub const BAKED_INDEX_PROV: u8 = 0x40;

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
///
/// Storing a scaled integer rather than a float keeps the record fixed-width and
/// makes the value bit-reproducible across platforms.
pub const SCORE_LOG_SCALE: f64 = 1000.0;

/// piscem `frag_map_type` (a.k.a. `MappingType`) code for an unmapped fragment.
pub const FRAG_MAP_TYPE_UNMAPPED: u8 = 0;

/// A single salmon mapping placement, in the form written to / read from RAD.
///
/// One fragment usually has several of these: "this fragment could have come
/// from transcript 12 at position 340, or transcript 88 at position 12, or …".
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
    ///
    /// `#[inline]` asks the compiler to paste the body into the caller: this is
    /// called once per written hit, so even a function call would be visible.
    #[inline]
    pub fn compressed_ori_ref(&self) -> u32 {
        // Mask the id down to 30 bits, then OR in each flag bit if set.
        (self.tid & TID_MASK)
            | if self.is_fw { ORI_FW } else { 0 }
            | if self.mate_fw { MATE_FW } else { 0 }
    }

    /// Decode a `compressed_ori_ref` u32 into `(tid, is_fw, mate_fw)`, the exact
    /// inverse of [`Self::compressed_ori_ref`].
    #[inline]
    pub fn decode_ori_ref(v: u32) -> (u32, bool, bool) {
        (v & TID_MASK, (v & ORI_FW) != 0, (v & MATE_FW) != 0)
    }

    /// Build a hit for one placement following piscem's bulk `frag_len`
    /// convention, shared by every RAD producer (reads, alignment-BAM, genome
    /// projection).
    ///
    /// **The `frag_len` convention.** A proper pair stores its real
    /// `fragment_len`. An orphan or single-end read has no fragment length —
    /// only one end was observed — so the slot instead stores the mapped mate's
    /// `read_len` (clamped below the [`FRAG_LEN_UNPAIRED`] sentinel). That lets
    /// the reader recover the bounded-CMF orphan fragment-length weight, i.e.
    /// "the true fragment was at least this long", rather than falling back to a
    /// flat penalty.
    ///
    /// `pos` is clamped to `>= 0` (a soft-clipped alignment can compute a
    /// negative start). `mate_fw` is meaningful only for proper pairs, and is
    /// forced to `false` otherwise so the encoded bit is deterministic.
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
                // Saturate rather than wrap: an absurd fragment length must not
                // silently become a small one.
                fragment_len.clamp(0, u16::MAX as i32) as u16
            } else {
                // MAX - 1 so a real read length can never be mistaken for the
                // "unpaired" sentinel.
                read_len.clamp(0, (u16::MAX - 1) as i32) as u16
            },
            score,
        }
    }
}

/// piscem `MappingType` codes for the read-level `frag_map_type` tag.
///
/// salmon's [`MateStatus`] maps onto the subset piscem uses for bulk data. The
/// numeric codes are fixed by the format, not chosen here.
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
    ///
    /// Note this is *not* the same ordering as the raw codes above, which is
    /// exactly why the rank function exists: `MAPPED_PAIR` is 4 but
    /// `RIGHT_ORPHAN` is 3, so comparing codes directly would work here by luck,
    /// while `SINGLE` (1) vs `LEFT_ORPHAN` (2) would not generalize.
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
    ///
    /// The tag is per read, but hits may disagree (the same fragment can be a
    /// proper pair on one transcript and an orphan on another), so one summary
    /// value has to be chosen. `unwrap_or(UNMAPPED)` covers the empty case.
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

    /// The packing is a wire format shared with piscem, so encode/decode must be
    /// an exact round trip, including at the extreme transcript id.
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

    /// When a fragment's hits disagree, the summary tag must report the most
    /// complete mapping seen, not the first or the last.
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

    /// Pin the paired/orphan `frag_len` convention and both saturation limits,
    /// since a reader decodes these bytes with no way to detect a wrong rule.
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
