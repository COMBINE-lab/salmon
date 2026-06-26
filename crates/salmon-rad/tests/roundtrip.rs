//! Round-trip and concurrency tests for salmon RAD output.

use std::io::{Cursor, Read};

use libradicl::chunk::Chunk;
use libradicl::header::RadPrelude;
use libradicl::rad_types::TagValue;

use salmon_rad::{
    detect_input_profile, frag_map_type, name_hash, FragmentChunkBuf, RadHit, RadInputProfile,
    RadOutputWriter, RadProfile, SalmonBulkContext, SalmonBulkRecord, FRAG_LEN_UNPAIRED,
};

fn sample_records(profile: RadProfile) -> Vec<SalmonBulkRecord> {
    let score = |s: i32| if profile.has_scores() { s } else { 0 };
    vec![
        SalmonBulkRecord::new(
            frag_map_type::MAPPED_PAIR,
            name_hash(b"read0"),
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
            name_hash(b"read1"),
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
    let names = ["t0", "t1", "t2"];
    let lens = [100u32, 200, 300];
    let w = RadOutputWriter::create(path, &names, &lens, true, profile).unwrap();
    let ctx: SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity(1024);
    for r in recs {
        cb.write(r, &ctx).unwrap();
    }
    w.append_chunk_bytes(&cb.take_bytes()).unwrap();
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
fn sketch_roundtrip() {
    check_roundtrip(
        RadProfile::Sketch,
        RadInputProfile::Salmon(RadProfile::Sketch),
    );
}

#[test]
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
fn concurrent_writer_counts_chunks() {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let names = ["t0", "t1", "t2"];
    let lens = [100u32, 200, 300];
    let w = RadOutputWriter::create(tmp.path(), &names, &lens, true, RadProfile::Sketch).unwrap();
    let ctx: SalmonBulkContext = *w.context();

    std::thread::scope(|s| {
        for t in 0..4u32 {
            let w = &w;
            s.spawn(move || {
                let mut cb = FragmentChunkBuf::with_capacity(256);
                let rec = SalmonBulkRecord::new(
                    frag_map_type::SINGLE,
                    name_hash(format!("read{t}").as_bytes()),
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
                w.append_chunk_bytes(&cb.take_bytes()).unwrap();
            });
        }
    });
    w.finalize().unwrap();

    let (prelude, _, _, got) = read_all(tmp.path());
    assert_eq!(prelude.hdr.num_chunks, 4, "one chunk per thread");
    assert_eq!(got.len(), 4, "all four records present");
}
