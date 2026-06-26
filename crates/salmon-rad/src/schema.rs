//! Building the [`RadPrelude`] (header + tag sections) and matching file-tag
//! [`TagMap`] for a salmon RAD file.

use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{
    RadIntId, RadType, TagDesc, TagMap, TagSection, TagSectionLabel, TagValue,
};

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
pub fn build_prelude(
    profile: RadProfile,
    is_paired: bool,
    ref_names: &[&str],
    ref_lengths: &[u32],
) -> (RadPrelude, TagMap) {
    // File tags carry concrete values (written once, after the header).
    let (file_tags, file_tag_map) = TagSection::from_tag_values(
        TagSectionLabel::FileTags,
        &[
            (
                "rad_type",
                TagValue::String(profile.rad_type_str().to_string()),
            ),
            ("ref_lengths", TagValue::ArrayU32(ref_lengths.to_vec())),
        ],
    );

    // Read-level tags (descriptors only; values are per-record in each chunk).
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    read_tags.add_tag_desc(int_desc("frag_map_type", RadIntId::U8));
    read_tags.add_tag_desc(int_desc("frag_name_hash", RadIntId::U128));

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
