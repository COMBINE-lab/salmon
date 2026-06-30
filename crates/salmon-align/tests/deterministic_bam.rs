//! Deterministic alignment mode: a transcriptome BAM/SAM projected to a salmon
//! RAD via `write_alignment_rad`, then quantified with `quantify_rad`. Verifies
//! the RAD is produced + fully baked (so requant is single-pass), mass is
//! conserved, and the whole pipeline is reproducible run-to-run.

use salmon_align::{quantify_rad, write_alignment_rad, AlignQuantOptions};
use salmon_rad::ChunkCodec;
use std::io::Write;
use std::path::{Path, PathBuf};

/// Tiny two-transcript transcriptome SAM, name-grouped (queryname), with inward
/// (ISF) paired fragments: unique to txA, unique to txB, and ambiguous (a
/// primary placement on txA + a secondary on txB) so an eq-class spans both.
fn write_sam(dir: &Path) -> PathBuf {
    let path = dir.join("aln.sam");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:queryname").unwrap();
    writeln!(f, "@SQ\tSN:txA\tLN:5000").unwrap();
    writeln!(f, "@SQ\tSN:txB\tLN:5000").unwrap();
    let pair = |f: &mut std::fs::File, name: &str, tx: &str, p: usize, second: bool| {
        let (l, r) = (p + 1, p + 201);
        let sf = if second { 256 } else { 0 };
        writeln!(f, "{name}\t{}\t{tx}\t{l}\t60\t100M\t=\t{r}\t300\t*\t*", 99 | sf).unwrap();
        writeln!(f, "{name}\t{}\t{tx}\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*", 147 | sf).unwrap();
    };
    let mut pos = 1usize;
    for i in 0..40 {
        pair(&mut f, &format!("a{i}"), "txA", pos, false);
        pos += 400;
    }
    for i in 0..20 {
        pair(&mut f, &format!("b{i}"), "txB", pos, false);
        pos += 400;
    }
    for i in 0..30 {
        pair(&mut f, &format!("amb{i}"), "txA", pos, false);
        pair(&mut f, &format!("amb{i}"), "txB", pos, true);
        pos += 400;
    }
    path
}

fn opts_for(sam: &Path, out: &Path) -> AlignQuantOptions {
    let mut o = AlignQuantOptions::new(sam.to_path_buf(), out.to_path_buf());
    o.lib_type = "IU".to_string();
    o.no_error_model = true; // lengths from @SQ; deterministic mode is score-based
    o
}

fn quant_counts(out: &Path) -> Vec<(String, f64)> {
    let txt = std::fs::read_to_string(out.join("quant.sf")).unwrap();
    txt.lines()
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[0].to_string(), c[4].parse::<f64>().unwrap())
        })
        .collect()
}

#[test]
fn deterministic_bam_writes_rad_and_quantifies() {
    let dir = tempfile::tempdir().unwrap();
    let sam = write_sam(dir.path());

    // Phase 1: BAM -> RAD.
    let out = dir.path().join("out");
    std::fs::create_dir_all(&out).unwrap();
    let opts = opts_for(&sam, &out);
    let rad = dir.path().join("mappings.rad");
    let summary = write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
    assert!(rad.exists(), "RAD must be written");
    assert_eq!(summary.num_mapped, 90, "40 + 20 + 30 fragments");

    // Phase 2: quantify from the RAD.
    let res = quantify_rad(&opts, &rad).unwrap();
    let counts = quant_counts(&out);
    assert_eq!(counts.len(), 2, "two transcripts");
    let total: f64 = counts.iter().map(|(_, c)| c).sum();
    assert!(
        (total - 90.0).abs() < 1e-6,
        "mass conserved: Σ NumReads = {total}, expected 90"
    );
    assert_eq!(res.num_mapped, 90);
}

#[test]
fn deterministic_bam_is_reproducible() {
    let dir = tempfile::tempdir().unwrap();
    let sam = write_sam(dir.path());

    let run = |tag: &str| -> Vec<(String, f64)> {
        let out = dir.path().join(tag);
        std::fs::create_dir_all(&out).unwrap();
        let opts = opts_for(&sam, &out);
        let rad = dir.path().join(format!("{tag}.rad"));
        write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
        quantify_rad(&opts, &rad).unwrap();
        quant_counts(&out)
    };

    // Two independent write+quant runs must agree exactly (the producer + the
    // baked-FLD requant are both order-independent).
    assert_eq!(run("r1"), run("r2"), "deterministic alignment quant must reproduce");
}
