//! Regression test for the alignment-mode analog of issue #1025: on a stranded
//! library, `num_mapped` must count only fragments with a surviving strand-
//! compatible placement, so the mass-conservation invariant `num_mapped == Σ
//! counts` holds (and `percent_mapped` is not inflated). Before the fix,
//! alignment mode counted every aligned fragment as mapped, regardless of the
//! strand-compatibility filter applied during assignment.

use salmon_align::{quantify_alignments, AlignQuantOptions};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Write a tiny transcriptome SAM with mixed-orientation FR fragments to `dir`.
/// `n_fwd` fragments have read1 on the forward strand (inward, ISF-compatible);
/// `n_rev` have read1 on the reverse strand (ISR-compatible, ISF-incompatible).
/// Records of a fragment are consecutive (grouped by read name), as an aligner
/// emits them. Returns the SAM path.
fn write_sam(dir: &Path, n_fwd: usize, n_rev: usize) -> PathBuf {
    let path = dir.join("aln.sam");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:unsorted").unwrap();
    writeln!(f, "@SQ\tSN:tx0\tLN:5000").unwrap();
    let emit = |f: &mut std::fs::File, name: &str, fwd: bool, p: usize| {
        let (l, r) = (p + 1, p + 201); // 1-based leftmost of the two mates
        let tlen = 300i32;
        if fwd {
            // read1 forward upstream, read2 reverse downstream -> inward, ISF
            writeln!(f, "{name}\t99\ttx0\t{l}\t60\t100M\t=\t{r}\t{tlen}\t*\t*").unwrap();
            writeln!(f, "{name}\t147\ttx0\t{r}\t60\t100M\t=\t{l}\t-{tlen}\t*\t*").unwrap();
        } else {
            // read1 reverse downstream, read2 forward upstream -> inward, ISR
            writeln!(f, "{name}\t83\ttx0\t{r}\t60\t100M\t=\t{l}\t-{tlen}\t*\t*").unwrap();
            writeln!(f, "{name}\t163\ttx0\t{l}\t60\t100M\t=\t{r}\t{tlen}\t*\t*").unwrap();
        }
    };
    let mut pos = 1usize;
    for i in 0..n_fwd {
        emit(&mut f, &format!("fwd{i}"), true, pos);
        pos += 400;
    }
    for i in 0..n_rev {
        emit(&mut f, &format!("rev{i}"), false, pos);
        pos += 400;
    }
    path
}

/// Write a tiny transcriptome SAM of *genuine single-end* records (BAM `0x1`
/// unset, so neither `0x40` nor `0x80`) to `dir`. `n_fwd` reads align to the
/// forward strand (flag `0x0`); `n_rev` to the reverse strand (flag `0x10`) —
/// the flag-16-heavy shape reported in issue #1057. Returns the SAM path.
fn write_single_end_sam(dir: &Path, n_fwd: usize, n_rev: usize) -> PathBuf {
    let path = dir.join("se.sam");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:unsorted").unwrap();
    writeln!(f, "@SQ\tSN:tx0\tLN:5000").unwrap();
    // 0x0 = mapped forward single-end; 0x10 = mapped reverse single-end.
    let emit = |f: &mut std::fs::File, name: &str, flag: u16, p: usize| {
        writeln!(f, "{name}\t{flag}\ttx0\t{}\t60\t100M\t*\t0\t0\t*\t*", p + 1).unwrap();
    };
    let mut pos = 1usize;
    for i in 0..n_fwd {
        emit(&mut f, &format!("fwd{i}"), 0, pos);
        pos += 200;
    }
    for i in 0..n_rev {
        emit(&mut f, &format!("rev{i}"), 0x10, pos);
        pos += 200;
    }
    path
}

fn run(sam: &Path, out: &Path, lib_type: &str) -> salmon_align::AlignQuantResult {
    let mut opts = AlignQuantOptions::new(sam.to_path_buf(), out.to_path_buf());
    opts.lib_type = lib_type.to_string();
    opts.transcripts = None; // lengths from @SQ; no error model / bias
    opts.no_error_model = true;
    quantify_alignments(&opts).expect("quantify_alignments")
}

#[test]
/// The regression itself: `num_mapped` must count fragments that survived the
/// strand filter, not every fragment that aligned, or the reported mapping rate
/// is inflated and mass conservation breaks.
fn stranded_align_num_mapped_matches_quantified_mass() {
    let tmp = tempfile::tempdir().unwrap();
    let (n_fwd, n_rev) = (5usize, 3usize);
    let sam = write_sam(tmp.path(), n_fwd, n_rev);
    let total = (n_fwd + n_rev) as u64;

    let mass_conserved = |res: &salmon_align::AlignQuantResult, tag: &str| {
        let sum: f64 = res.counts.iter().sum();
        assert!(
            (res.num_mapped as f64 - sum).abs() <= 1e-6,
            "{tag}: mass not conserved: num_mapped={} but Σ counts={sum:.6}",
            res.num_mapped
        );
    };

    // Unstranded: every fragment is compatible -> all mapped, mass conserved.
    let iu = run(&sam, &tmp.path().join("iu"), "IU");
    assert_eq!(iu.num_processed, total, "IU: num_processed");
    assert_eq!(iu.num_mapped, total, "IU: every fragment should map");
    mass_conserved(&iu, "IU");

    // Stranded ISF: the reverse-strand fragments are strand-incompatible and are
    // dropped during assignment. They must NOT be counted as mapped, and mass
    // must stay conserved. (Before the fix, num_mapped was `total` for both.)
    let isf = run(&sam, &tmp.path().join("isf"), "ISF");
    assert_eq!(
        isf.num_processed, total,
        "ISF: num_processed counts all aligned"
    );
    assert!(
        isf.num_mapped < total,
        "ISF: strand-incompatible fragments counted as mapped: num_mapped={} of {total}",
        isf.num_mapped
    );
    assert_eq!(
        isf.num_mapped, n_fwd as u64,
        "ISF: only the forward (ISF-compatible) fragments should map"
    );
    mass_conserved(&isf, "ISF");
}

/// Regression for issue #1057: genuine single-end records (BAM `0x1` unset) must
/// be classified `SingleEnd`, so the single-end strand filters accept them by
/// their own alignment orientation. Before the fix they were mislabeled as
/// paired-end right orphans and dropped under `SF`/`SR` (0 fragments mapped),
/// leaving `U` as the only library type that worked.
#[test]
/// Single-end records must be filtered on their own orientation, which is the
/// only strand information they carry.
fn single_end_strand_filters_respect_alignment_orientation() {
    let tmp = tempfile::tempdir().unwrap();
    let (n_fwd, n_rev) = (12usize, 48usize); // flag-16-heavy, as reported
    let sam = write_single_end_sam(tmp.path(), n_fwd, n_rev);
    let total = (n_fwd + n_rev) as u64;

    let mass_conserved = |res: &salmon_align::AlignQuantResult, tag: &str| {
        let sum: f64 = res.counts.iter().sum();
        assert!(
            (res.num_mapped as f64 - sum).abs() <= 1e-6,
            "{tag}: mass not conserved: num_mapped={} but Σ counts={sum:.6}",
            res.num_mapped
        );
    };

    // U (unstranded single-end): accepts both orientations -> all map.
    let u = run(&sam, &tmp.path().join("u"), "U");
    assert_eq!(u.num_processed, total, "U: num_processed");
    assert_eq!(u.num_mapped, total, "U: every single-end read should map");
    mass_conserved(&u, "U");

    // SF (forward single-end): accepts only the forward reads.
    let sf = run(&sam, &tmp.path().join("sf"), "SF");
    assert_eq!(
        sf.num_mapped, n_fwd as u64,
        "SF: only forward-strand single-end reads should map (was 0 before #1057 fix)"
    );
    mass_conserved(&sf, "SF");

    // SR (reverse single-end): accepts only the reverse reads (the reporter's case).
    let sr = run(&sam, &tmp.path().join("sr"), "SR");
    assert_eq!(
        sr.num_mapped, n_rev as u64,
        "SR: only reverse-strand single-end reads should map (was 0 before #1057 fix)"
    );
    mass_conserved(&sr, "SR");
}

/// `lib_format_counts.json` must report measured strand-compatibility counts in
/// alignment mode too, not the historical hardcoded `1.0` (#1130). The BAM
/// streaming pass judges every reported alignment against the expected format
/// already; this pins that the judgment reaches the output file: exact
/// `num_compatible_fragments` / `num_incompatible_fragments`, the ratio over
/// judged fragments, and `num_assigned_fragments` still excluding the dropped
/// wrong-strand fragments (default `--incompatPrior 0`).
#[test]
fn stranded_align_lib_format_counts_report_incompatible_fragments() {
    let tmp = tempfile::tempdir().unwrap();
    let (n_fwd, n_rev) = (12usize, 7usize);
    let sam = write_sam(tmp.path(), n_fwd, n_rev);
    let total = (n_fwd + n_rev) as u64;

    let lib_counts = |dir: &Path| -> serde_json::Value {
        let s = std::fs::read_to_string(dir.join("lib_format_counts.json")).unwrap();
        serde_json::from_str(&s).unwrap()
    };
    let u64_field = |v: &serde_json::Value, k: &str| v[k].as_u64().unwrap();

    // Stranded ISF: the reverse fragments aligned but only on the wrong strand.
    let out_isf = tmp.path().join("isf_counts");
    let isf_res = run(&sam, &out_isf, "ISF");
    let isf = lib_counts(&out_isf);
    assert_eq!(u64_field(&isf, "num_compatible_fragments"), n_fwd as u64);
    assert_eq!(u64_field(&isf, "num_incompatible_fragments"), n_rev as u64);
    let ratio = isf["compatible_fragment_ratio"].as_f64().unwrap();
    let expect = n_fwd as f64 / total as f64;
    assert!((ratio - expect).abs() < 1e-12, "ratio {ratio} != {expect}");
    // Dropped, so not assigned: the two counts answer different questions.
    assert_eq!(
        u64_field(&isf, "num_assigned_fragments"),
        isf_res.num_mapped
    );
    assert_eq!(isf_res.num_mapped, n_fwd as u64);
    // The observed-format histogram (pre-filter) and its C++-semantics derived
    // fields: under a stranded ISF expectation only ISF observations are
    // concordant, and the strand bias is the sense share of the inward
    // orientation.
    assert_eq!(u64_field(&isf, "ISF"), n_fwd as u64);
    assert_eq!(u64_field(&isf, "ISR"), n_rev as u64);
    assert_eq!(
        u64_field(&isf, "num_frags_with_concordant_consistent_mappings"),
        n_fwd as u64
    );
    assert_eq!(
        u64_field(&isf, "num_frags_with_inconsistent_or_orphan_mappings"),
        n_rev as u64
    );
    let bias = isf["strand_mapping_bias"].as_f64().unwrap();
    let expect_bias = n_fwd as f64 / total as f64;
    assert!(
        (bias - expect_bias).abs() < 1e-12,
        "bias {bias} != {expect_bias}"
    );
    assert_eq!(isf["expected_format"].as_str().unwrap(), "ISF");

    // Unstranded: judged like any other type, and everything passes — but the
    // histogram still reports the strand split.
    let out_iu = tmp.path().join("iu_counts");
    run(&sam, &out_iu, "IU");
    let iu = lib_counts(&out_iu);
    assert_eq!(u64_field(&iu, "num_compatible_fragments"), total);
    assert_eq!(u64_field(&iu, "num_incompatible_fragments"), 0);
    assert_eq!(iu["compatible_fragment_ratio"].as_f64().unwrap(), 1.0);
    assert_eq!(
        u64_field(&iu, "num_frags_with_concordant_consistent_mappings"),
        total
    );
    assert_eq!(
        u64_field(&iu, "num_frags_with_inconsistent_or_orphan_mappings"),
        0
    );
    let bias = iu["strand_mapping_bias"].as_f64().unwrap();
    assert!(
        (bias - expect_bias).abs() < 1e-12,
        "bias {bias} != {expect_bias}"
    );
}

/// `-l A` in online alignment mode must infer the library type from the
/// alignments (the reads-mode prefix detector, fed from the BAM/SAM stream) and
/// report it — where before it silently skipped detection and filtering
/// entirely, leaving `expected_format: "A"` and a fictitious ratio.
///
/// Both fixtures are far below the 50k-sample lock-in budget, so no live
/// filtering happens (same as reads mode on a small run); the end-of-run best
/// guess is still reported and names the file's expected format.
#[test]
fn align_auto_library_type_detects_from_alignments() {
    let tmp = tempfile::tempdir().unwrap();
    let lib_counts = |dir: &Path| -> serde_json::Value {
        let s = std::fs::read_to_string(dir.join("lib_format_counts.json")).unwrap();
        serde_json::from_str(&s).unwrap()
    };

    // Paired, 12 inward-FR vs 7 inward-RF: a 63% forward ratio is inside
    // salmon's 30–70% unstranded band, so detection names `IU`.
    let (n_fwd, n_rev) = (12usize, 7usize);
    let sam = write_sam(tmp.path(), n_fwd, n_rev);
    let out = tmp.path().join("auto_pe");
    let res = run(&sam, &out, "A");
    assert_eq!(res.detected_library_type.as_deref(), Some("IU"));
    assert_eq!(res.num_mapped, (n_fwd + n_rev) as u64);
    let v = lib_counts(&out);
    assert_eq!(v["expected_format"].as_str().unwrap(), "IU");
    assert_eq!(v["ISF"].as_u64().unwrap(), n_fwd as u64);
    assert_eq!(v["ISR"].as_u64().unwrap(), n_rev as u64);

    // Single-end, 48 forward vs 2 reverse: 96% forward is past the 70%
    // threshold, so detection names `SF` (the read type comes from peeking the
    // first record's flags, C++'s `peekBAMIsPaired`). Below the lock-in budget
    // nothing was judged, so the counts keep the historical fallback form.
    let se = write_single_end_sam(tmp.path(), 48, 2);
    let out = tmp.path().join("auto_se");
    let res = run(&se, &out, "A");
    assert_eq!(res.detected_library_type.as_deref(), Some("SF"));
    assert_eq!(res.num_mapped, 50);
    let v = lib_counts(&out);
    assert_eq!(v["expected_format"].as_str().unwrap(), "SF");
    assert_eq!(v["SF"].as_u64().unwrap(), 48);
    assert_eq!(v["SR"].as_u64().unwrap(), 2);
    assert_eq!(v["num_incompatible_fragments"].as_u64().unwrap(), 0);
    assert_eq!(v["num_assigned_fragments"].as_u64().unwrap(), 50);
    let bias = v["strand_mapping_bias"].as_f64().unwrap();
    assert!((bias - 0.96).abs() < 1e-12, "bias {bias} != 0.96");
}
