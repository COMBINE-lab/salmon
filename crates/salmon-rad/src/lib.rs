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
//! Both profiles additionally carry a read-level tag — a deterministic, seeded
//! `xxh3` hash of the read name — so a future opt-in determinism feature can
//! induce a fixed input order by sorting on it. Emitting the tag has no
//! behavioral effect on its own.
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

pub use record::{SalmonBulkContext, SalmonBulkRecord};
pub use writer::{FragmentChunkBuf, RadOutputWriter};

/// The profile of a RAD file presented for *input*, as detected from its
/// prelude tag signature (see [`detect_input_profile`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RadInputProfile {
    /// piscem `map-bulk` (`bulk_with_pos`): `frag_map_type` read tag + the
    /// `(compressed_ori_ref, pos, frag_len)` aln triple, no scores.
    PiscemBulk,
    /// salmon-written RAD (carries the `frag_name_hash` read tag); the inner
    /// [`RadProfile`] distinguishes sketch vs. selective-alignment.
    Salmon(RadProfile),
}

/// Detect the input profile of a RAD file from its prelude's tag signature.
///
/// libradicl has no single "kind" field, so detection is by which tags are
/// present: salmon RAD is identified by the `frag_name_hash` read tag (and
/// selective-alignment by the `alignment_score` aln tag); a piscem bulk file is
/// identified by `frag_map_type` + the `compressed_ori_ref`/`pos`/`frag_len`
/// alignment triple.
pub fn detect_input_profile(
    prelude: &libradicl::header::RadPrelude,
) -> anyhow::Result<RadInputProfile> {
    let rt = &prelude.read_tags;
    let at = &prelude.aln_tags;
    if rt.has_tag("frag_name_hash") {
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
        "unrecognized RAD profile: expected a salmon RAD (with a `frag_name_hash` read tag) \
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

/// Seed for the read-name hash tag. Fixed so the hash is reproducible across
/// runs/machines; recorded in `meta_info.json` as `rad_hash_seed`. ("SalmonRD")
pub const SALMON_RAD_HASH_SEED: u64 = 0x5361_6c6d_6f6e_5244;

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

/// The seeded `xxh3`-128 hash of a read name, used as the `frag_name_hash`
/// read-level tag. `name` should be the read id trimmed at the first ASCII
/// whitespace (the same slice salmon uses elsewhere); raw bytes are hashed, no
/// UTF-8 decoding.
#[inline]
pub fn name_hash(name: &[u8]) -> u128 {
    xxhash_rust::xxh3::xxh3_128_with_seed(name, SALMON_RAD_HASH_SEED)
}

/// Trim a read id at the first ASCII whitespace, matching salmon's read-name
/// handling (so the hash is stable regardless of FASTQ comment fields).
#[inline]
pub fn trim_read_name(id: &[u8]) -> &[u8] {
    match id.iter().position(|b| b.is_ascii_whitespace()) {
        Some(i) => &id[..i],
        None => id,
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
    fn name_hash_is_seeded_and_trimmed() {
        let full = b"READ_42 1:N:0:ACGT";
        let trimmed = trim_read_name(full);
        assert_eq!(trimmed, b"READ_42");
        assert_eq!(
            name_hash(trimmed),
            xxhash_rust::xxh3::xxh3_128_with_seed(b"READ_42", SALMON_RAD_HASH_SEED)
        );
        // stable + distinct from the unseeded hash
        assert_ne!(
            name_hash(b"READ_42"),
            xxhash_rust::xxh3::xxh3_128(b"READ_42")
        );
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
}
