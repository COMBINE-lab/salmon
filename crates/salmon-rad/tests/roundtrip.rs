//! Round-trip and concurrency tests for salmon RAD output.

use std::io::{Cursor, Read};

use libradicl::chunk::Chunk;
use libradicl::header::RadPrelude;
use libradicl::rad_types::TagValue;

use salmon_rad::{
    detect_input_profile, frag_map_type, ChunkCodec, FragmentChunkBuf, RadHit, RadInputProfile,
    RadOutputWriter, RadProfile, SalmonBulkContext, SalmonBulkRecord, FRAG_LEN_UNPAIRED,
};

/// Writer provenance standing in for a mapping pass, so the round-trip covers
/// the provenance tags rather than only the shapes that omit them.
fn test_prov() -> salmon_rad::WriterProvenance {
    salmon_rad::WriterProvenance {
        mapping_type: salmon_rad::MappingType::Mapping,
        index: Some(salmon_rad::IndexProvenance {
            seq_hash: "seqhash".into(),
            name_hash: "namehash".into(),
            seq_hash512: "seqhash512".into(),
            name_hash512: "namehash512".into(),
            decoy_seq_hash: String::new(),
            decoy_name_hash: String::new(),
            first_decoy_index: None,
            num_decoys: None,
            keep_duplicates: Some(false),
        }),
        source_programs: vec![],
    }
}

fn sample_records(profile: RadProfile) -> Vec<SalmonBulkRecord> {
    let score = |s: i32| if profile.has_scores() { s } else { 0 };
    vec![
        SalmonBulkRecord::new(
            frag_map_type::MAPPED_PAIR,
            vec![
                RadHit {
                    tid: 0,
                    is_fw: true,
                    mate_fw: false,
                    pos: 10,
                    frag_len: 250,
                    score: score(-3),
                },
                RadHit {
                    tid: 2,
                    is_fw: false,
                    mate_fw: true,
                    pos: 77,
                    frag_len: 312,
                    score: score(-5),
                },
            ],
        ),
        SalmonBulkRecord::new(
            frag_map_type::LEFT_ORPHAN,
            vec![RadHit {
                tid: 1,
                is_fw: true,
                mate_fw: false,
                pos: 5,
                frag_len: FRAG_LEN_UNPAIRED,
                score: score(-1),
            }],
        ),
    ]
}

fn write_file(path: &std::path::Path, profile: RadProfile, recs: &[SalmonBulkRecord]) {
    write_file_codec(path, profile, recs, ChunkCodec::None);
}

fn write_file_codec(
    path: &std::path::Path,
    profile: RadProfile,
    recs: &[SalmonBulkRecord],
    codec: ChunkCodec,
) {
    let names = ["t0", "t1", "t2"];
    let lens = [100u32, 200, 300];
    let w = RadOutputWriter::create(
        path,
        &names,
        &lens,
        true,
        profile,
        1024,
        codec,
        &test_prov(),
    )
    .unwrap();
    let ctx: SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity_codec(1024, codec);
    for r in recs {
        cb.write(r, &ctx).unwrap();
    }
    w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
    w.finalize().unwrap();
}

fn read_all(
    path: &std::path::Path,
) -> (
    RadPrelude,
    TagValue,
    SalmonBulkContext,
    Vec<SalmonBulkRecord>,
) {
    let mut data = Vec::new();
    std::fs::File::open(path)
        .unwrap()
        .read_to_end(&mut data)
        .unwrap();
    let mut rc = Cursor::new(data);
    let prelude = RadPrelude::from_bytes(&mut rc).unwrap();
    let ftm = prelude.file_tags.parse_tags_from_bytes(&mut rc).unwrap();
    let ref_lengths = ftm.get("ref_lengths").unwrap().clone();
    let ctx: SalmonBulkContext = prelude.get_record_context().unwrap();
    let mut recs = Vec::new();
    for _ in 0..prelude.hdr.num_chunks {
        let chunk = Chunk::<SalmonBulkRecord>::from_bytes(&mut rc, &ctx);
        recs.extend(chunk.reads);
    }
    (prelude, ref_lengths, ctx, recs)
}

/// Read one file tag by name from a written RAD file.
fn read_file_tag(path: &std::path::Path, name: &str) -> TagValue {
    read_file_tag_opt(path, name).unwrap_or_else(|| panic!("missing file tag {name}"))
}

/// Like `read_file_tag`, but `None` when the tag is absent — the distinction the
/// unknown-vs-false cases turn on.
fn read_file_tag_opt(path: &std::path::Path, name: &str) -> Option<TagValue> {
    let mut data = Vec::new();
    std::fs::File::open(path)
        .unwrap()
        .read_to_end(&mut data)
        .unwrap();
    let mut rc = Cursor::new(data);
    let prelude = RadPrelude::from_bytes(&mut rc).unwrap();
    let ftm = prelude.file_tags.parse_tags_from_bytes(&mut rc).unwrap();
    ftm.get(name).cloned()
}

fn check_roundtrip(profile: RadProfile, expected_detect: RadInputProfile) {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let recs = sample_records(profile);
    write_file(tmp.path(), profile, &recs);

    let (prelude, ref_lengths, ctx, got) = read_all(tmp.path());

    // header
    assert_eq!(prelude.hdr.is_paired, 1);
    assert_eq!(prelude.hdr.num_chunks, 1, "num_chunks must be backpatched");
    assert_eq!(prelude.hdr.ref_names, vec!["t0", "t1", "t2"]);
    // file tag: ref_lengths preserved
    assert_eq!(ref_lengths, TagValue::ArrayU32(vec![100, 200, 300]));
    // profile detection + record context
    assert_eq!(detect_input_profile(&prelude).unwrap(), expected_detect);
    assert_eq!(ctx.profile, profile);
    // records: exact equality (scores are 0 in sketch by construction)
    assert_eq!(got, recs);
}

#[test]
/// What is written must read back identically: the writer and reader are separate
/// implementations of one binary layout.
fn sketch_roundtrip() {
    check_roundtrip(
        RadProfile::Sketch,
        RadInputProfile::Salmon(RadProfile::Sketch),
    );
}

/// The `chunk_codec` file tag is absent for uncompressed output (so pre-feature
/// readers are unaffected) and carries the codec id (1=lz4, 2=zstd) otherwise.
#[test]
/// A compressed file must announce its codec, or a reader cannot decode it.
fn chunk_codec_tag_advertised() {
    use libradicl::CHUNK_CODEC_TAG;
    let recs = sample_records(RadProfile::Sketch);
    for (codec, expect) in [
        (ChunkCodec::None, None),
        (ChunkCodec::Lz4, Some(1u8)),
        (ChunkCodec::Zstd, Some(2u8)),
    ] {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        write_file_codec(tmp.path(), RadProfile::Sketch, &recs, codec);
        // Parse only the prelude + file-tag values (the chunks may be compressed).
        let mut data = Vec::new();
        std::fs::File::open(tmp.path())
            .unwrap()
            .read_to_end(&mut data)
            .unwrap();
        let mut rc = Cursor::new(data);
        let prelude = RadPrelude::from_bytes(&mut rc).unwrap();
        let ftm = prelude.file_tags.parse_tags_from_bytes(&mut rc).unwrap();
        match (ftm.get(CHUNK_CODEC_TAG), expect) {
            (None, None) => {}
            (Some(TagValue::U8(v)), Some(e)) => assert_eq!(*v, e, "codec {codec:?}"),
            (got, _) => panic!("codec {codec:?}: unexpected chunk_codec tag {got:?}"),
        }
    }
}

#[test]
/// The selective-alignment profile adds a per-hit score field, which must survive
/// the round trip.
fn sa_roundtrip_carries_scores() {
    check_roundtrip(
        RadProfile::SelectiveAlignment,
        RadInputProfile::Salmon(RadProfile::SelectiveAlignment),
    );
    // explicitly confirm a non-zero score survives
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let recs = sample_records(RadProfile::SelectiveAlignment);
    write_file(tmp.path(), RadProfile::SelectiveAlignment, &recs);
    let (_, _, _, got) = read_all(tmp.path());
    assert_eq!(got[0].hits[0].score, -3);
    assert_eq!(got[0].hits[1].score, -5);
}

#[test]
/// Many threads append concurrently, so no chunk may be lost or double-counted.
fn concurrent_writer_counts_chunks() {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let names = ["t0", "t1", "t2"];
    let lens = [100u32, 200, 300];
    let w = RadOutputWriter::create(
        tmp.path(),
        &names,
        &lens,
        true,
        RadProfile::Sketch,
        1024,
        ChunkCodec::None,
        &test_prov(),
    )
    .unwrap();
    let ctx: SalmonBulkContext = *w.context();

    std::thread::scope(|s| {
        for t in 0..4u32 {
            let w = &w;
            s.spawn(move || {
                let mut cb = FragmentChunkBuf::with_capacity(256);
                let rec = SalmonBulkRecord::new(
                    frag_map_type::SINGLE,
                    vec![RadHit {
                        tid: t,
                        is_fw: true,
                        mate_fw: false,
                        pos: t,
                        frag_len: FRAG_LEN_UNPAIRED,
                        score: 0,
                    }],
                );
                cb.write(&rec, &ctx).unwrap();
                w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
            });
        }
    });
    w.finalize().unwrap();

    let (prelude, _, _, got) = read_all(tmp.path());
    assert_eq!(prelude.hdr.num_chunks, 4, "one chunk per thread");
    assert_eq!(got.len(), 4, "all four records present");
}

/// A RAD file must only take the name the user asked for once it is complete.
///
/// Regression test for salmon#1105: a run that died mid-write (a full disk) left
/// a truncated file at the requested path, and a later `--rad` run consumed it
/// as if it were whole.
#[test]
fn abandoned_write_never_claims_the_final_path() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out.rad");
    let names = ["t0", "t1", "t2"];
    let lens = [100u32, 200, 300];

    let w = RadOutputWriter::create(
        &path,
        &names,
        &lens,
        true,
        RadProfile::SelectiveAlignment,
        1024,
        ChunkCodec::None,
        &test_prov(),
    )
    .unwrap();
    let ctx: SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity_codec(1024, ChunkCodec::None);
    for r in &sample_records(RadProfile::SelectiveAlignment) {
        cb.write(r, &ctx).unwrap();
    }
    w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
    // Simulates the run dying before finalize: chunks are on disk, but the
    // header was never backpatched.
    drop(w);

    assert!(
        !path.exists(),
        "an unfinalized write must not appear at the requested path"
    );
    let partials: Vec<_> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n.starts_with("out.rad.partial-"))
        .collect();
    assert_eq!(
        partials.len(),
        1,
        "the incomplete bytes should remain under a `.partial-*` name, found {partials:?}"
    );
}

/// `finalize` renames the partial into place and marks the file complete.
#[test]
fn finalize_renames_and_marks_complete() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out.rad");
    write_file_codec(
        &path,
        RadProfile::SelectiveAlignment,
        &sample_records(RadProfile::SelectiveAlignment),
        ChunkCodec::None,
    );

    assert!(
        path.exists(),
        "finalize should rename the partial into place"
    );
    let leftovers: Vec<_> = std::fs::read_dir(dir.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|n| n.contains(".partial-"))
        .collect();
    assert!(leftovers.is_empty(), "partial left behind: {leftovers:?}");

    let flags = match read_file_tag(&path, salmon_rad::BAKED_FLAGS_TAG) {
        TagValue::U8(v) => v,
        other => panic!("baked_flags should be a u8, got {other:?}"),
    };
    assert_ne!(
        flags & salmon_rad::WRITE_COMPLETE,
        0,
        "a finalized file should carry the confirmatory completeness bit"
    );
}

/// A requant can only reproduce the mapping run's `meta_info.json` if the RAD
/// records what the mapping pass *saw*, which its records cannot show: a RAD
/// holds only the fragments that mapped.
#[test]
fn mapping_pass_provenance_round_trips() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out.rad");
    let names = ["t0", "t1", "t2"];
    let lens = [100u32, 200, 300];

    let mut w = RadOutputWriter::create(
        &path,
        &names,
        &lens,
        true,
        RadProfile::SelectiveAlignment,
        1024,
        ChunkCodec::None,
        &test_prov(),
    )
    .unwrap();
    let ctx: SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity_codec(1024, ChunkCodec::None);
    for r in &sample_records(RadProfile::SelectiveAlignment) {
        cb.write(r, &ctx).unwrap();
    }
    w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
    // Only two of these fragments were written, out of ten the pass observed.
    w.set_map_counters(salmon_rad::MapCounters {
        num_processed: 10,
        num_dovetail: 3,
        num_filtered_vm: 4,
        num_below_threshold_vm: 5,
        num_decoy_fragments: 6,
    });
    w.finalize().unwrap();

    let u64_tag = |name: &str| match read_file_tag(&path, name) {
        TagValue::U64(v) => v,
        other => panic!("{name} should be a u64, got {other:?}"),
    };
    assert_eq!(u64_tag(salmon_rad::NUM_PROCESSED_TAG), 10);
    assert_eq!(u64_tag(salmon_rad::NUM_DOVETAIL_TAG), 3);
    assert_eq!(u64_tag(salmon_rad::NUM_FILTERED_VM_TAG), 4);
    assert_eq!(u64_tag(salmon_rad::NUM_BELOW_THRESH_VM_TAG), 5);
    assert_eq!(u64_tag(salmon_rad::NUM_DECOY_FRAGMENTS_TAG), 6);

    let flags = match read_file_tag(&path, salmon_rad::BAKED_FLAGS_TAG) {
        TagValue::U8(v) => v,
        other => panic!("baked_flags should be a u8, got {other:?}"),
    };
    assert_ne!(flags & salmon_rad::BAKED_MAP_COUNTERS, 0);
    assert_ne!(flags & salmon_rad::BAKED_INDEX_PROV, 0);

    match read_file_tag(&path, salmon_rad::INDEX_SEQ_HASH_TAG) {
        TagValue::String(v) => assert_eq!(v, "seqhash"),
        other => panic!("index_seq_hash should be a string, got {other:?}"),
    }
    match read_file_tag(&path, salmon_rad::MAPPING_TYPE_TAG) {
        TagValue::String(v) => assert_eq!(v, "mapping"),
        other => panic!("mapping_type should be a string, got {other:?}"),
    }
}

/// A writer that records no counters must leave the flag clear, so a reader can
/// tell "the pass observed none of these" from "nobody filled the slot" — the
/// slots are reserved as zeros, so the values alone cannot say.
#[test]
fn unbaked_counters_are_distinguishable_from_zero() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out.rad");
    write_file_codec(
        &path,
        RadProfile::SelectiveAlignment,
        &sample_records(RadProfile::SelectiveAlignment),
        ChunkCodec::None,
    );
    let flags = match read_file_tag(&path, salmon_rad::BAKED_FLAGS_TAG) {
        TagValue::U8(v) => v,
        other => panic!("baked_flags should be a u8, got {other:?}"),
    };
    assert_eq!(
        flags & salmon_rad::BAKED_MAP_COUNTERS,
        0,
        "counters were never set, so the flag must stay clear"
    );
    assert!(matches!(
        read_file_tag(&path, salmon_rad::NUM_PROCESSED_TAG),
        TagValue::U64(0)
    ));
}

/// `keep_duplicates` must survive as three distinct states, not two: an index
/// that predates recording it is *unknown*, and reporting that as `false` would
/// state as fact something no one observed.
#[test]
fn unknown_keep_duplicates_is_not_written_as_false() {
    for (keep, expect_tag) in [
        (Some(true), Some(1u8)),
        (Some(false), Some(0)),
        (None, None),
    ] {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("out.rad");
        let mut prov = test_prov();
        prov.index.as_mut().unwrap().keep_duplicates = keep;

        let w = RadOutputWriter::create(
            &path,
            &["t0", "t1", "t2"],
            &[100u32, 200, 300],
            true,
            RadProfile::SelectiveAlignment,
            1024,
            ChunkCodec::None,
            &prov,
        )
        .unwrap();
        w.finalize().unwrap();

        let got = match read_file_tag_opt(&path, salmon_rad::KEEP_DUPLICATES_TAG) {
            Some(TagValue::U8(v)) => Some(v),
            None => None,
            other => panic!("keep_duplicates should be a u8 or absent, got {other:?}"),
        };
        assert_eq!(got, expect_tag, "for keep_duplicates = {keep:?}");
    }
}

/// A BAM-derived RAD carries the aligner's `@PG` chain verbatim, so a requant
/// can report how the alignments were produced rather than only that they were
/// alignments. Tab-separated fields must survive the join/split round trip.
#[test]
fn source_programs_round_trip() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("out.rad");
    let lines = vec![
        "@PG\tID:bowtie2\tPN:bowtie2\tVN:2.4.5\tCL:bowtie2-align-s --wrapper basic-0 -x idx"
            .to_string(),
        "@PG\tID:samtools\tPN:samtools\tPP:bowtie2\tVN:1.17\tCL:samtools view -b".to_string(),
    ];
    let prov = salmon_rad::WriterProvenance {
        mapping_type: salmon_rad::MappingType::Alignment,
        index: None,
        source_programs: lines.clone(),
    };
    let w = RadOutputWriter::create(
        &path,
        &["t0", "t1", "t2"],
        &[100u32, 200, 300],
        true,
        RadProfile::SelectiveAlignment,
        1024,
        ChunkCodec::None,
        &prov,
    )
    .unwrap();
    w.finalize().unwrap();

    match read_file_tag(&path, salmon_rad::SOURCE_PROGRAMS_TAG) {
        TagValue::String(v) => {
            let got: Vec<&str> = v.split('\n').collect();
            assert_eq!(got, lines, "the @PG chain must survive verbatim");
        }
        other => panic!("source_programs should be a string, got {other:?}"),
    }
    match read_file_tag(&path, salmon_rad::MAPPING_TYPE_TAG) {
        TagValue::String(v) => assert_eq!(v, "alignment"),
        other => panic!("mapping_type should be a string, got {other:?}"),
    }
}
