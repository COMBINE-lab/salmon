//! Building the [`RadPrelude`] (header + tag sections) and matching file-tag
//! [`TagMap`] for a salmon RAD file.

use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{
    RadIntId, RadType, TagDesc, TagMap, TagSection, TagSectionLabel, TagValue,
};
use libradicl::{ChunkCodec, CHUNK_CODEC_TAG};

use crate::RadProfile;

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
pub fn build_prelude(
    profile: RadProfile,
    is_paired: bool,
    ref_names: &[&str],
    ref_lengths: &[u32],
    fld_len: usize,
    codec: ChunkCodec,
) -> (RadPrelude, TagMap) {
    let num_refs = ref_names.len();
    // File tags carry concrete values (written once, after the header). The
    // reserved placeholders are backpatched once the FLD / abundances /
    // library format are known at the end of the pass.
    let mut file_tag_values: Vec<(&str, TagValue)> = vec![
        (
            "rad_type",
            TagValue::String(profile.rad_type_str().to_string()),
        ),
        ("ref_lengths", TagValue::ArrayU32(ref_lengths.to_vec())),
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
    ];
    // Advertise chunk compression only when enabled; omitting the tag for
    // uncompressed output keeps the file byte-compatible with readers that
    // predate the feature (absent tag ⇒ ChunkCodec::None).
    if codec != ChunkCodec::None {
        file_tag_values.push((CHUNK_CODEC_TAG, TagValue::U8(codec.as_u8())));
    }
    let (file_tags, file_tag_map) =
        TagSection::from_tag_values(TagSectionLabel::FileTags, &file_tag_values);

    // Read-level tags (descriptors only; values are per-record in each chunk).
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    read_tags.add_tag_desc(int_desc("frag_map_type", RadIntId::U8));

    // Alignment-level tags. The first three are byte-identical to piscem
    // `bulk_with_pos`; SA adds a per-hit score.
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
