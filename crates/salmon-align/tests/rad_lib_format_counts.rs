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
/// `baked_fmt` bakes a library format id into the header, as a mapping run
/// resolves and bakes its own; `None` leaves the header without one.
fn write_mixed_orientation_rad_baked(
    path: &Path,
    n_isf: usize,
    n_isr: usize,
    baked_fmt: Option<salmon_core::LibraryFormat>,
) {
    let prov = salmon_rad::WriterProvenance {
        mapping_type: salmon_rad::MappingType::Mapping,
        index: None,
        source_programs: vec![],
    };
    let mut w = RadOutputWriter::create(
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
    if let Some(f) = baked_fmt {
        w.set_library_format(f.format_id());
    }
    w.append_chunk_bytes(&cb.take_bytes().unwrap()).unwrap();
    w.finalize().unwrap();
}

fn write_mixed_orientation_rad(path: &Path, n_isf: usize, n_isr: usize) {
    write_mixed_orientation_rad_baked(path, n_isf, n_isr, None)
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
    // The observed-format histogram is measured from the raw placements
    // (pre-filter), and the C++-semantics derived fields follow from it: under
    // a stranded ISF expectation only the ISF observations are concordant, and
    // the strand bias is the sense share of the inward orientation.
    assert_eq!(u64_field(&v, "ISF"), n_isf as u64);
    assert_eq!(u64_field(&v, "ISR"), n_isr as u64);
    assert_eq!(u64_field(&v, "IU"), 0);
    assert_eq!(u64_field(&v, "SF"), 0);
    assert_eq!(
        u64_field(&v, "num_frags_with_concordant_consistent_mappings"),
        n_isf as u64
    );
    assert_eq!(
        u64_field(&v, "num_frags_with_inconsistent_or_orphan_mappings"),
        n_isr as u64
    );
    let bias = v["strand_mapping_bias"].as_f64().unwrap();
    let expect_bias = n_isf as f64 / total as f64;
    assert!(
        (bias - expect_bias).abs() < 1e-12,
        "bias {bias} != {expect_bias}"
    );
    assert_eq!(v["expected_format"].as_str().unwrap(), "ISF");

    // The opposite stranded type flips the split exactly.
    let out_isr = tmp.path().join("isr");
    run("ISR", &out_isr);
    let v = lib_counts(&out_isr);
    assert_eq!(u64_field(&v, "num_compatible_fragments"), n_isr as u64);
    assert_eq!(u64_field(&v, "num_incompatible_fragments"), n_isf as u64);

    // Unstranded: judged like any other type, and everything passes — but the
    // histogram-derived fields still expose the strand split, which is exactly
    // the "potential strand bias in an unstranded protocol" signal.
    let out_iu = tmp.path().join("iu");
    run("IU", &out_iu);
    let v = lib_counts(&out_iu);
    assert_eq!(u64_field(&v, "num_compatible_fragments"), total);
    assert_eq!(u64_field(&v, "num_incompatible_fragments"), 0);
    assert_eq!(v["compatible_fragment_ratio"].as_f64().unwrap(), 1.0);
    assert_eq!(
        u64_field(&v, "num_frags_with_concordant_consistent_mappings"),
        total
    );
    assert_eq!(
        u64_field(&v, "num_frags_with_inconsistent_or_orphan_mappings"),
        0
    );
    let bias = v["strand_mapping_bias"].as_f64().unwrap();
    let expect_bias = n_isf as f64 / total as f64;
    assert!(
        (bias - expect_bias).abs() < 1e-12,
        "bias {bias} != {expect_bias}"
    );
}

/// `-l A` on a RAD with no baked library format must detect the type from the
/// records during the FLD-derive pass and then *filter* the quant pass with it
/// — whole-run detect-then-filter, the requant analog of the reads-mode
/// detector.
#[test]
fn rad_requant_auto_detects_and_filters() {
    let tmp = tempfile::tempdir().unwrap();

    // Heavily skewed: 48 ISF-observed vs 2 ISR-observed is a 96% forward
    // ratio, past salmon's 70% strandedness threshold, so detection must name
    // `ISF` — and the quant pass must then drop the 2 wrong-strand fragments.
    let rad = tmp.path().join("skewed.rad");
    write_mixed_orientation_rad(&rad, 48, 2);
    let out = tmp.path().join("auto_skewed");
    let mut o = AlignQuantOptions::new(rad.clone(), out.clone());
    o.lib_type = "A".to_string();
    o.fld_mean = 200.0;
    o.fld_sd = 20.0;
    let res = quantify_rad(&o, &rad).expect("quantify_rad");
    assert_eq!(res.detected_library_type.as_deref(), Some("ISF"));
    assert_eq!(res.num_mapped, 48);
    let v = lib_counts(&out);
    assert_eq!(v["expected_format"].as_str().unwrap(), "ISF");
    assert_eq!(u64_field(&v, "num_compatible_fragments"), 48);
    assert_eq!(u64_field(&v, "num_incompatible_fragments"), 2);
    assert_eq!(u64_field(&v, "num_assigned_fragments"), 48);
    assert_eq!(u64_field(&v, "ISF"), 48);
    assert_eq!(u64_field(&v, "ISR"), 2);

    // A balanced mix (60% forward, inside the 30–70% unstranded band) detects
    // `IU`: nothing is filtered, and the strand split is still reported.
    let rad = tmp.path().join("balanced.rad");
    write_mixed_orientation_rad(&rad, 30, 20);
    let out = tmp.path().join("auto_balanced");
    let mut o = AlignQuantOptions::new(rad.clone(), out.clone());
    o.lib_type = "A".to_string();
    o.fld_mean = 200.0;
    o.fld_sd = 20.0;
    let res = quantify_rad(&o, &rad).expect("quantify_rad");
    assert_eq!(res.detected_library_type.as_deref(), Some("IU"));
    assert_eq!(res.num_mapped, 50);
    let v = lib_counts(&out);
    assert_eq!(v["expected_format"].as_str().unwrap(), "IU");
    assert_eq!(u64_field(&v, "num_compatible_fragments"), 50);
    assert_eq!(u64_field(&v, "num_incompatible_fragments"), 0);
    let bias = v["strand_mapping_bias"].as_f64().unwrap();
    assert!((bias - 0.6).abs() < 1e-12, "bias {bias} != 0.6");
}

/// A library format baked into the RAD header is authoritative under `-l A`:
/// the producing run already resolved it, and the requant must filter with it
/// and report it — even when re-deriving from the records would say otherwise
/// (this mix re-derives as `IU`). This is the deterministic two-phase flow's
/// phase 2, where phase 1 bakes its `-l A` resolution.
#[test]
fn rad_requant_baked_library_format_is_authoritative() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("baked.rad");
    let isr = salmon_core::LibraryFormat::parse("ISR").unwrap();
    write_mixed_orientation_rad_baked(&rad, 30, 20, Some(isr));
    let out = tmp.path().join("auto_baked");
    let mut o = AlignQuantOptions::new(rad.clone(), out.clone());
    o.lib_type = "A".to_string();
    o.fld_mean = 200.0;
    o.fld_sd = 20.0;
    let res = quantify_rad(&o, &rad).expect("quantify_rad");
    assert_eq!(res.detected_library_type.as_deref(), Some("ISR"));
    // Only the 20 ISR-observed pairs survive the baked type's filter.
    assert_eq!(res.num_mapped, 20);
    let v = lib_counts(&out);
    assert_eq!(v["expected_format"].as_str().unwrap(), "ISR");
    assert_eq!(u64_field(&v, "num_compatible_fragments"), 20);
    assert_eq!(u64_field(&v, "num_incompatible_fragments"), 30);
    assert_eq!(u64_field(&v, "num_assigned_fragments"), 20);
    assert_eq!(u64_field(&v, "ISF"), 30);
    assert_eq!(u64_field(&v, "ISR"), 20);
}

/// A two-phase driver hands its phase-1 timing to the requant, and
/// `meta_info.json` must report the whole run plus which inference path ran
/// (#1140): `total_time_seconds` previously covered only the second pass
/// (0.4s reported against minutes of wall on real data), and nothing recorded
/// that the deterministic path produced the output.
#[test]
fn requant_reports_run_spanning_time_and_inference_path() {
    let tmp = tempfile::tempdir().unwrap();
    let rad = tmp.path().join("timed.rad");
    write_mixed_orientation_rad(&rad, 30, 20);
    let out = tmp.path().join("timed_out");
    let mut o = AlignQuantOptions::new(rad.clone(), out.clone());
    o.lib_type = "IU".to_string();
    o.fld_mean = 200.0;
    o.fld_sd = 20.0;
    o.prior_seconds = 100.0;
    o.external_start_time = Some("Thu Jan  1 00:00:00 2026".to_string());
    let res = quantify_rad(&o, &rad).expect("quantify_rad");
    assert!(
        res.total_seconds >= 100.0,
        "total_seconds {} must include the driver's prior phase",
        res.total_seconds
    );
    let meta: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap(),
    )
    .unwrap();
    assert!(meta["total_time_seconds"].as_f64().unwrap() >= 100.0);
    assert_eq!(
        meta["start_time"].as_str().unwrap(),
        "Thu Jan  1 00:00:00 2026"
    );
    assert_eq!(meta["inference_path"].as_str().unwrap(), "deterministic");
}
