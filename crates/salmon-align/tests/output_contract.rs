//! The requant writer's output contract (#1140): decoy references stay out of
//! `quant.sf` and the files positionally aligned with it, and `meta_info.json`
//! reports what the run actually did (`samp_type`, `serialized_eq_classes`,
//! `num_decoy_fragments`) instead of hardcoded placeholders.
//!
//! A reads-mode RAD's header must carry every reference — records index into
//! it — so the *writer* is where decoys are excluded, exactly as reads mode
//! excludes them. Before #1140's audit, every `quantify_rad` path (including
//! the deterministic default) emitted decoy rows the docs promise are absent.

use std::path::Path;

use salmon_align::{quantify_rad, AlignQuantOptions};
use salmon_rad::{frag_map_type, FragmentChunkBuf, RadHit, RadOutputWriter, RadProfile};

/// A RAD whose reference table is `[t0, t1, dA]` with `dA` marked as the decoy
/// block via the baked index provenance, plus baked mapping counters. Proper
/// pairs land only on the real transcripts, as the mapper guarantees.
fn write_decoy_rad(path: &Path) {
    let prov = salmon_rad::WriterProvenance {
        mapping_type: salmon_rad::MappingType::Mapping,
        index: Some(salmon_rad::IndexProvenance {
            seq_hash: "sh".into(),
            name_hash: "nh".into(),
            seq_hash512: String::new(),
            name_hash512: String::new(),
            decoy_seq_hash: "dsh".into(),
            decoy_name_hash: "dnh".into(),
            keep_duplicates: Some(false),
            first_decoy_index: Some(2),
            num_decoys: Some(1),
        }),
        source_programs: vec![],
    };
    let mut w = RadOutputWriter::create(
        path,
        &["t0", "t1", "dA"],
        &[4000u32, 4000, 4000],
        true,
        RadProfile::Sketch,
        1024,
        salmon_rad::ChunkCodec::None,
        &prov,
    )
    .unwrap();
    let ctx: salmon_rad::SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity_codec(4096, salmon_rad::ChunkCodec::None);
    for (tid, n) in [(0u32, 30), (1u32, 20)] {
        for _ in 0..n {
            let hit = RadHit {
                tid,
                is_fw: true,
                mate_fw: false,
                pos: 0,
                frag_len: 200,
                score: 0,
            };
            let rec = salmon_rad::SalmonBulkRecord::new(frag_map_type::MAPPED_PAIR, vec![hit]);
            cb.write(&rec, &ctx).unwrap();
        }
    }
    w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
    w.set_map_counters(salmon_rad::MapCounters {
        num_processed: 60,
        num_dovetail: 0,
        num_filtered_vm: 0,
        num_below_threshold_vm: 0,
        num_decoy_fragments: 7,
    });
    w.finalize().unwrap();
}

#[test]
fn requant_outputs_honor_the_contract() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("decoy.rad");
    write_decoy_rad(&rad);

    let out = tmp.path().join("out");
    let mut o = AlignQuantOptions::new(rad.clone(), out.clone());
    o.lib_type = "IU".to_string();
    o.fld_mean = 200.0;
    o.fld_sd = 20.0;
    o.num_bootstraps = 2;
    o.dump_eq = true;
    let res = quantify_rad(&o, &rad).expect("quantify_rad");
    assert_eq!(res.num_mapped, 50);

    // quant.sf: the decoy block is excluded, exactly as reads mode excludes it.
    let sf = std::fs::read_to_string(out.join("quant.sf")).unwrap();
    let names: Vec<&str> = sf
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert_eq!(names, ["t0", "t1"], "decoy dA must not be a quant.sf row");

    // ambig_info.tsv stays row-aligned with quant.sf.
    let ambig = std::fs::read_to_string(out.join("aux_info").join("ambig_info.tsv")).unwrap();
    assert_eq!(
        ambig.lines().count(),
        1 + 2,
        "header + one row per quant row"
    );

    // Bootstrap files carry the same row set, positionally aligned.
    {
        use std::io::Read as _;
        let f = std::fs::File::open(out.join("aux_info").join("bootstrap").join("names.tsv.gz"))
            .unwrap();
        let mut txt = String::new();
        flate2::read::GzDecoder::new(f)
            .read_to_string(&mut txt)
            .unwrap();
        assert_eq!(txt.trim_end(), "t0\tt1");

        let f = std::fs::File::open(out.join("aux_info").join("bootstrap").join("bootstraps.gz"))
            .unwrap();
        let mut bytes = Vec::new();
        flate2::read::GzDecoder::new(f)
            .read_to_end(&mut bytes)
            .unwrap();
        assert_eq!(
            bytes.len(),
            2 * 2 * 8,
            "2 samples x 2 quant rows x 8 bytes; decoy columns would make it 2 x 3 x 8"
        );
    }

    // meta_info.json reports what the run did, not placeholders.
    let meta: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(meta["num_valid_targets"], 2);
    assert_eq!(
        meta["samp_type"], "bootstrap",
        "bootstraps ran and were written; tximport keys off this"
    );
    assert_eq!(
        meta["serialized_eq_classes"], true,
        "the eq-class dump was requested and written"
    );
    assert_eq!(
        meta["num_decoy_fragments"], 7,
        "from the RAD-baked mapping counters, like its siblings"
    );
    assert_eq!(meta["num_processed"], 60);
    assert!(out.join("aux_info").join("eq_classes.txt.gz").is_file());
}

/// A RAD without the decoy boundary (piscem, BAM-derived, or pre-boundary
/// salmon) behaves exactly as before: every header reference is a row.
#[test]
fn rad_without_decoy_boundary_emits_every_reference() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("plain.rad");
    let prov = salmon_rad::WriterProvenance {
        mapping_type: salmon_rad::MappingType::Mapping,
        index: None,
        source_programs: vec![],
    };
    let w = RadOutputWriter::create(
        &rad,
        &["t0", "t1"],
        &[4000u32, 4000],
        true,
        RadProfile::Sketch,
        1024,
        salmon_rad::ChunkCodec::None,
        &prov,
    )
    .unwrap();
    let ctx: salmon_rad::SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity_codec(4096, salmon_rad::ChunkCodec::None);
    for _ in 0..10 {
        let hit = RadHit {
            tid: 0,
            is_fw: true,
            mate_fw: false,
            pos: 0,
            frag_len: 200,
            score: 0,
        };
        let rec = salmon_rad::SalmonBulkRecord::new(frag_map_type::MAPPED_PAIR, vec![hit]);
        cb.write(&rec, &ctx).unwrap();
    }
    w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
    w.finalize().unwrap();

    let out = tmp.path().join("out");
    let mut o = AlignQuantOptions::new(rad.clone(), out.clone());
    o.lib_type = "IU".to_string();
    o.fld_mean = 200.0;
    o.fld_sd = 20.0;
    quantify_rad(&o, &rad).expect("quantify_rad");
    let sf = std::fs::read_to_string(out.join("quant.sf")).unwrap();
    assert_eq!(sf.lines().count(), 1 + 2);
}
