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

fn run(sam: &Path, out: &Path, lib_type: &str) -> salmon_align::AlignQuantResult {
    let mut opts = AlignQuantOptions::new(sam.to_path_buf(), out.to_path_buf());
    opts.lib_type = lib_type.to_string();
    opts.transcripts = None; // lengths from @SQ; no error model / bias
    opts.no_error_model = true;
    quantify_alignments(&opts).expect("quantify_alignments")
}

#[test]
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
