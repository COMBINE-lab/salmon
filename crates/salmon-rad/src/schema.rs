//! Building the [`RadPrelude`] (header + tag sections) and matching file-tag
//! [`TagMap`] for a salmon RAD file.
//!
//! # What a "prelude" is
//!
//! Every RAD file starts with a self-describing header: how many references
//! there are and what they are called, plus three *tag sections* that declare
//! the fields each later record will contain, in order. A reader parses the tag
//! sections first and then knows exactly how to decode the binary payload; there
//! is no per-record field naming, which is what keeps the format compact.
//!
//! The three sections are:
//!   * **file tags** — one value per file (written once, right after the header);
//!   * **read tags** — one value per fragment;
//!   * **alignment tags** — one value per placement of a fragment.
//!
//! This module only *declares* the read/alignment tags; the values are written
//! per record by [`crate::record`].

use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{
    RadIntId, RadType, TagDesc, TagMap, TagSection, TagSectionLabel, TagValue,
};
use libradicl::{ChunkCodec, CHUNK_CODEC_TAG};

use crate::RadProfile;

/// Shorthand for declaring an integer-typed tag: a name plus a width
/// (`RadIntId::U32` and friends). Used purely to keep the declarations below
/// readable.
fn int_desc(name: &str, id: RadIntId) -> TagDesc {
    TagDesc {
        name: name.to_string(),
        typeid: RadType::Int(id),
    }
}

/// Build the prelude and file-tag values for a salmon RAD file.
///
/// `ref_names` / `ref_lengths` cover the **full** reference table (including any
/// decoys); decoys simply never appear in records (salmon filters them at write
/// time). `ref_lengths` is preserved so a reader can compute effective lengths
/// without the original index.
///
/// `fld_len` reserves a fixed-size `frag_length_dist` slot (a log-PMF over raw
/// lengths `[0, fld_len)`) and `initial_abundances` reserves one f64 per
/// reference. Both are written as zero placeholders and backpatched at finalize
/// (see [`crate::RadOutputWriter`]); a `baked_flags` byte records which were
/// actually filled, so a reader can do a single-pass quant when the FLD (and,
/// for bias, abundances) are present and fall back to deriving them otherwise.
///
/// **Why placeholders + backpatching.** The fragment-length distribution and the
/// abundance estimates are only known once the whole run has finished, but the
/// header must be written before any records. So the slots are reserved at the
/// correct size up front and overwritten in place at the end — which only works
/// because their size does not depend on the data.
pub fn build_prelude(
    profile: RadProfile,
    is_paired: bool,
    ref_names: &[&str],
    ref_lengths: &[u32],
    fld_len: usize,
    codec: ChunkCodec,
    provenance: &crate::WriterProvenance,
) -> (RadPrelude, TagMap) {
    let num_refs = ref_names.len();
    // File tags carry concrete values (written once, after the header). The
    // reserved placeholders are backpatched once the FLD / abundances /
    // library format are known at the end of the pass.
    let mut file_tag_values: Vec<(&str, TagValue)> = vec![
        // Identifies the file as salmon's and says which profile it is.
        (
            "rad_type",
            TagValue::String(profile.rad_type_str().to_string()),
        ),
        ("ref_lengths", TagValue::ArrayU32(ref_lengths.to_vec())),
        // Zeroed bitfield: "nothing baked yet".
        (crate::BAKED_FLAGS_TAG, TagValue::U8(0)),
        (
            crate::FRAG_LENGTH_DIST_TAG,
            TagValue::ArrayF64(vec![0.0; fld_len]),
        ),
        (
            crate::INITIAL_ABUNDANCES_TAG,
            TagValue::ArrayF64(vec![0.0; num_refs]),
        ),
        (crate::LIBRARY_FORMAT_TAG, TagValue::U8(0)),
        // Mapping-pass counters. Reserved as zeros and backpatched at finalize;
        // the `BAKED_MAP_COUNTERS` flag says whether they were actually filled,
        // since a real count of zero is indistinguishable from a placeholder.
        (crate::NUM_PROCESSED_TAG, TagValue::U64(0)),
        (crate::NUM_DOVETAIL_TAG, TagValue::U64(0)),
        (crate::NUM_FILTERED_VM_TAG, TagValue::U64(0)),
        (crate::NUM_BELOW_THRESH_VM_TAG, TagValue::U64(0)),
        (crate::NUM_DECOY_FRAGMENTS_TAG, TagValue::U64(0)),
    ];
    // Index identity, known before the first record, so written outright rather
    // than reserved. Omitted entirely when the producer has no index (a BAM-derived
    // RAD), because an absent tag says "unknown" and an empty string does not.
    file_tag_values.push((
        crate::MAPPING_TYPE_TAG,
        TagValue::String(provenance.mapping_type.as_str().to_string()),
    ));
    if let Some(prov) = &provenance.index {
        file_tag_values.extend([
            (
                crate::INDEX_SEQ_HASH_TAG,
                TagValue::String(prov.seq_hash.clone()),
            ),
            (
                crate::INDEX_NAME_HASH_TAG,
                TagValue::String(prov.name_hash.clone()),
            ),
            (
                crate::INDEX_SEQ_HASH512_TAG,
                TagValue::String(prov.seq_hash512.clone()),
            ),
            (
                crate::INDEX_NAME_HASH512_TAG,
                TagValue::String(prov.name_hash512.clone()),
            ),
            (
                crate::INDEX_DECOY_SEQ_HASH_TAG,
                TagValue::String(prov.decoy_seq_hash.clone()),
            ),
            (
                crate::INDEX_DECOY_NAME_HASH_TAG,
                TagValue::String(prov.decoy_name_hash.clone()),
            ),
            (
                crate::KEEP_DUPLICATES_TAG,
                TagValue::U8(prov.keep_duplicates as u8),
            ),
        ]);
    }
    // Score interpretation, only meaningful for the scored (SA) profile. Reserved
    // as the default (raw AS) and backpatched at finalize if the write run scored
    // by the deterministic error model instead (see [`crate::SCORE_KIND_TAG`]).
    if profile.has_scores() {
        file_tag_values.push((crate::SCORE_KIND_TAG, TagValue::U8(crate::SCORE_KIND_AS)));
    }
    // Advertise chunk compression only when enabled; omitting the tag for
    // uncompressed output keeps the file byte-compatible with readers that
    // predate the feature (absent tag ⇒ ChunkCodec::None).
    if codec != ChunkCodec::None {
        file_tag_values.push((CHUNK_CODEC_TAG, TagValue::U8(codec.as_u8())));
    }
    // Returns both the serializable section and a name-keyed map, so the writer
    // can later look a tag up by name in order to backpatch it.
    let (file_tags, file_tag_map) =
        TagSection::from_tag_values(TagSectionLabel::FileTags, &file_tag_values);

    // Read-level tags (descriptors only; values are per-record in each chunk).
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    read_tags.add_tag_desc(int_desc("frag_map_type", RadIntId::U8));

    // Alignment-level tags. The first three are byte-identical to piscem
    // `bulk_with_pos`; SA adds a per-hit score. Declaration order *is* the
    // on-disk field order, so these lines must match the record writer exactly.
    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    aln_tags.add_tag_desc(int_desc("compressed_ori_ref", RadIntId::U32));
    aln_tags.add_tag_desc(int_desc("pos", RadIntId::U32));
    aln_tags.add_tag_desc(int_desc("frag_len", RadIntId::U16));
    if profile.has_scores() {
        aln_tags.add_tag_desc(int_desc("alignment_score", RadIntId::I32));
    }

    let hdr = RadHeader {
        is_paired: is_paired as u8,
        ref_count: ref_names.len() as u64,
        ref_names: ref_names.iter().map(|s| s.to_string()).collect(),
        num_chunks: 0, // backpatched by RadFileWriter::finalize
    };

    (
        RadPrelude {
            hdr,
            file_tags,
            read_tags,
            aln_tags,
        },
        file_tag_map,
    )
}
