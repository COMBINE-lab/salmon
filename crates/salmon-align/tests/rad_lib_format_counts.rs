//! RAD-input quantification must report measured strand-compatibility counts in
//! `lib_format_counts.json` (#1130), not the historical hardcoded `1.0`.
//!
//! This is what makes the full `--deterministic` flow correct: its phase 2 is a
//! RAD requant into the final output directory, whose `write_outputs` rewrites
//! the mapping pass's `lib_format_counts.json`. If the requant did not measure,
//! every deterministic run would end with placeholder values no matter what the
//! mapping pass tallied — which is exactly the bug this test would catch.

use std::path::Path;

use salmon_align::{quantify_rad, AlignQuantOptions};
use salmon_rad::{frag_map_type, FragmentChunkBuf, RadHit, RadOutputWriter, RadProfile};

/// Write a salmon RAD of uniquely-mapped proper pairs on one transcript:
/// `n_isf` fragments observed inward-FR (read1 forward — `ISF`-compatible) and
/// `n_isr` observed inward-RF (read1 reverse — `ISR`-compatible).
fn write_mixed_orientation_rad(path: &Path, n_isf: usize, n_isr: usize) {
    let prov = salmon_rad::WriterProvenance {
        mapping_type: salmon_rad::MappingType::Mapping,
        index: None,
        source_programs: vec![],
    };
    let w = RadOutputWriter::create(
        path,
        &["t0"],
        &[4000u32],
        true,
        RadProfile::Sketch,
        1024,
        salmon_rad::ChunkCodec::None,
        &prov,
    )
    .unwrap();
    let ctx: salmon_rad::SalmonBulkContext = *w.context();
    let mut cb = FragmentChunkBuf::with_capacity_codec(4096, salmon_rad::ChunkCodec::None);
    let mut emit = |is_fw: bool, n: usize| {
        for _ in 0..n {
            let hit = RadHit {
                tid: 0,
                is_fw,
                mate_fw: !is_fw,
                pos: 0,
                frag_len: 200,
                score: 0,
            };
            let rec = salmon_rad::SalmonBulkRecord::new(frag_map_type::MAPPED_PAIR, vec![hit]);
            cb.write(&rec, &ctx).unwrap();
        }
    };
    emit(true, n_isf);
    emit(false, n_isr);
    w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
    w.finalize().unwrap();
}

fn lib_counts(dir: &Path) -> serde_json::Value {
    let s = std::fs::read_to_string(dir.join("lib_format_counts.json")).unwrap();
    serde_json::from_str(&s).unwrap()
}

fn u64_field(v: &serde_json::Value, k: &str) -> u64 {
    v[k].as_u64().unwrap()
}

#[test]
fn rad_requant_lib_format_counts_report_incompatible_fragments() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("mixed.rad");
    let (n_isf, n_isr) = (30usize, 20usize);
    let total = (n_isf + n_isr) as u64;
    write_mixed_orientation_rad(&rad, n_isf, n_isr);

    let run = |lib_type: &str, out: &Path| -> salmon_align::AlignQuantResult {
        let mut o = AlignQuantOptions::new(rad.clone(), out.to_path_buf());
        o.lib_type = lib_type.to_string();
        o.fld_mean = 200.0;
        o.fld_sd = 20.0;
        quantify_rad(&o, &rad).expect("quantify_rad")
    };

    // Stranded ISF: the ISR-observed fragments are judged and found wrong-strand;
    // with the default `--incompatPrior 0` they are dropped, so they are counted
    // incompatible yet excluded from the assigned total.
    let out_isf = tmp.path().join("isf");
    let res = run("ISF", &out_isf);
    let v = lib_counts(&out_isf);
    assert_eq!(u64_field(&v, "num_compatible_fragments"), n_isf as u64);
    assert_eq!(u64_field(&v, "num_incompatible_fragments"), n_isr as u64);
    let ratio = v["compatible_fragment_ratio"].as_f64().unwrap();
    let expect = n_isf as f64 / total as f64;
    assert!((ratio - expect).abs() < 1e-12, "ratio {ratio} != {expect}");
    assert_eq!(u64_field(&v, "num_assigned_fragments"), res.num_mapped);
    assert_eq!(res.num_mapped, n_isf as u64);

    // The opposite stranded type flips the split exactly.
    let out_isr = tmp.path().join("isr");
    run("ISR", &out_isr);
    let v = lib_counts(&out_isr);
    assert_eq!(u64_field(&v, "num_compatible_fragments"), n_isr as u64);
    assert_eq!(u64_field(&v, "num_incompatible_fragments"), n_isf as u64);

    // Unstranded: judged like any other type, and everything passes.
    let out_iu = tmp.path().join("iu");
    run("IU", &out_iu);
    let v = lib_counts(&out_iu);
    assert_eq!(u64_field(&v, "num_compatible_fragments"), total);
    assert_eq!(u64_field(&v, "num_incompatible_fragments"), 0);
    assert_eq!(v["compatible_fragment_ratio"].as_f64().unwrap(), 1.0);
}
