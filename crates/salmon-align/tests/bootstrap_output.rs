//! Alignment mode supports posterior uncertainty (bootstrap / Gibbs) through the
//! same packed equivalence-class samplers as reads mode. This test confirms the
//! `aux_info/bootstrap/{names.tsv.gz, bootstraps.gz}` files are emitted, align
//! positionally and by name with `quant.sf`, and that each replicate conserves
//! mass (sums to the number of quantified fragments).

use flate2::read::GzDecoder;
use salmon_align::{quantify_alignments, AlignQuantOptions};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

/// Tiny two-transcript transcriptome SAM with inward (ISF) paired fragments,
/// some uniquely mapping and some ambiguous (aligning to both transcripts) so
/// the bootstrap actually has reassignment to resample over.
fn write_sam(dir: &Path) -> PathBuf {
    let path = dir.join("aln.sam");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:unsorted").unwrap();
    writeln!(f, "@SQ\tSN:txA\tLN:5000").unwrap();
    writeln!(f, "@SQ\tSN:txB\tLN:5000").unwrap();
    let pair = |f: &mut std::fs::File, name: &str, tx: &str, p: usize| {
        let (l, r) = (p + 1, p + 201);
        writeln!(f, "{name}\t99\t{tx}\t{l}\t60\t100M\t=\t{r}\t300\t*\t*").unwrap();
        writeln!(f, "{name}\t147\t{tx}\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*").unwrap();
    };
    let mut pos = 1usize;
    // 40 unique to txA, 20 unique to txB, 30 ambiguous (both txA and txB).
    for i in 0..40 {
        pair(&mut f, &format!("a{i}"), "txA", pos);
        pos += 400;
    }
    for i in 0..20 {
        pair(&mut f, &format!("b{i}"), "txB", pos);
        pos += 400;
    }
    for i in 0..30 {
        // ambiguous: a primary placement on txA and a secondary on txB
        let (l, r) = (pos + 1, pos + 201);
        let name = format!("amb{i}");
        writeln!(f, "{name}\t99\ttxA\t{l}\t60\t100M\t=\t{r}\t300\t*\t*").unwrap();
        writeln!(f, "{name}\t147\ttxA\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*").unwrap();
        writeln!(f, "{name}\t353\ttxB\t{l}\t60\t100M\t=\t{r}\t300\t*\t*").unwrap();
        writeln!(f, "{name}\t401\ttxB\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*").unwrap();
        pos += 400;
    }
    path
}

fn run(sam: &Path, out: &Path, configure: impl FnOnce(&mut AlignQuantOptions)) {
    let mut opts = AlignQuantOptions::new(sam.to_path_buf(), out.to_path_buf());
    opts.lib_type = "IU".to_string();
    opts.transcripts = None; // lengths from @SQ; no error model / bias
    opts.no_error_model = true;
    configure(&mut opts);
    quantify_alignments(&opts).expect("quantify_alignments");
}

fn read_quant(out: &Path) -> (Vec<String>, Vec<f64>) {
    let txt = std::fs::read_to_string(out.join("quant.sf")).unwrap();
    let mut names = Vec::new();
    let mut counts = Vec::new();
    for line in txt.lines().skip(1) {
        let cols: Vec<&str> = line.split('\t').collect();
        names.push(cols[0].to_string());
        counts.push(cols[4].parse::<f64>().unwrap());
    }
    (names, counts)
}

fn read_replicates(out: &Path, n_txp: usize) -> (Vec<String>, Vec<Vec<f64>>) {
    let bdir = out.join("aux_info").join("bootstrap");
    let mut s = String::new();
    GzDecoder::new(std::fs::File::open(bdir.join("names.tsv.gz")).unwrap())
        .read_to_string(&mut s)
        .unwrap();
    let names: Vec<String> = s.trim().split('\t').map(|x| x.to_string()).collect();

    let mut bytes = Vec::new();
    GzDecoder::new(std::fs::File::open(bdir.join("bootstraps.gz")).unwrap())
        .read_to_end(&mut bytes)
        .unwrap();
    assert_eq!(bytes.len() % (n_txp * 8), 0, "ragged bootstrap stream");
    let reps: Vec<Vec<f64>> = bytes
        .chunks_exact(n_txp * 8)
        .map(|rep| {
            let (words, _) = rep.as_chunks::<8>();
            words.iter().copied().map(f64::from_le_bytes).collect()
        })
        .collect();
    (names, reps)
}

fn check(out: &Path, requested: usize) {
    let (qnames, qcounts) = read_quant(out);
    let (bnames, reps) = read_replicates(out, qnames.len());
    // names.tsv.gz aligns positionally and by name with quant.sf
    assert_eq!(bnames, qnames, "bootstrap names must match quant.sf rows");
    assert_eq!(reps.len(), requested, "replicate count");
    let total: f64 = qcounts.iter().sum();
    for (i, rep) in reps.iter().enumerate() {
        assert_eq!(rep.len(), qnames.len(), "replicate {i} length");
        let s: f64 = rep.iter().sum();
        assert!(
            (s - total).abs() < 1e-3,
            "replicate {i} mass {s} != quantified total {total}"
        );
    }
}

#[test]
/// Bootstrap output must line up with `quant.sf` positionally *and* by name —
/// tximport reads the two together — and every replicate must conserve mass.
fn align_bootstrap_output_aligns_with_quant_and_conserves_mass() {
    let tmp = tempfile::tempdir().unwrap();
    let sam = write_sam(tmp.path());
    let out = tmp.path().join("boot");
    run(&sam, &out, |o| o.num_bootstraps = 16);
    check(&out, 16);
}

#[test]
/// The same obligations for Gibbs samples, which take a different code path to
/// the same file layout.
fn align_gibbs_output_aligns_with_quant_and_conserves_mass() {
    let tmp = tempfile::tempdir().unwrap();
    let sam = write_sam(tmp.path());
    let out = tmp.path().join("gibbs");
    run(&sam, &out, |o| {
        o.num_gibbs_samples = 16;
        o.thinning_factor = 8;
    });
    check(&out, 16);
}
