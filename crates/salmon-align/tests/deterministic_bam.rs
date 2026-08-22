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
        writeln!(
            f,
            "{name}\t{}\t{tx}\t{l}\t60\t100M\t=\t{r}\t300\t*\t*",
            99 | sf
        )
        .unwrap();
        writeln!(
            f,
            "{name}\t{}\t{tx}\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*",
            147 | sf
        )
        .unwrap();
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
/// The two-pass contract: the first pass must produce a RAD with the models
/// fully baked, so the second pass is single-shot and conserves mass.
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

// ---- error-model path (`--deterministic` with `-t`) -----------------------

/// Two independent 50 kb pseudo-random transcripts. A read drawn from txA matches
/// txA and mismatches an (independent) txB at ~3/4 of positions, so once the
/// model has learned that matches dominate, the correct placement scores higher.
fn tx_seqs() -> (String, String) {
    const BASES: [u8; 4] = *b"ACGT";
    let gen = |mut s: u64| -> String {
        (0..50_000)
            .map(|_| {
                // SplitMix64-ish LCG; top bits are the most mixed.
                s = s
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                BASES[(s >> 33) as usize & 3] as char
            })
            .collect()
    };
    (gen(0x1234_5678), gen(0x9E37_79B9_7F4A_7C15))
}

/// Sibling of [`write_sam`] carrying real 100 bp `SEQ` (so the error model can
/// walk the CIGAR) + a matching transcriptome FASTA. Ambiguous fragments' reads
/// are drawn from txA, with the secondary placement on the mismatching txB.
fn write_seq_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let (sa, sb) = tx_seqs();
    let fasta = dir.join("txome.fa");
    {
        let mut f = std::fs::File::create(&fasta).unwrap();
        writeln!(f, ">txA\n{sa}\n>txB\n{sb}").unwrap();
    }
    let path = dir.join("aln_seq.sam");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:queryname").unwrap();
    writeln!(f, "@SQ\tSN:txA\tLN:50000").unwrap();
    writeln!(f, "@SQ\tSN:txB\tLN:50000").unwrap();
    // SEQ is always in forward-reference orientation, so a perfect match uses the
    // transcript substring for both mates. `src` is the sequence the read is
    // drawn from; `tx` is where this placement maps (differs for the amb decoy).
    let pair = |f: &mut std::fs::File, name: &str, tx: &str, src: &str, p: usize, second: bool| {
        let (l, r) = (p + 1, p + 201);
        let sf = if second { 256 } else { 0 };
        let s1 = &src[p..p + 100];
        let s2 = &src[p + 200..p + 300];
        writeln!(
            f,
            "{name}\t{}\t{tx}\t{l}\t60\t100M\t=\t{r}\t300\t{s1}\t*",
            99 | sf
        )
        .unwrap();
        writeln!(
            f,
            "{name}\t{}\t{tx}\t{r}\t60\t100M\t=\t{l}\t-300\t{s2}\t*",
            147 | sf
        )
        .unwrap();
    };
    let mut pos = 1usize;
    for i in 0..40 {
        pair(&mut f, &format!("a{i}"), "txA", &sa, pos, false);
        pos += 400;
    }
    for i in 0..20 {
        pair(&mut f, &format!("b{i}"), "txB", &sb, pos, false);
        pos += 400;
    }
    for i in 0..30 {
        // Read is from txA; primary maps txA (match), secondary maps txB (mismatch).
        pair(&mut f, &format!("amb{i}"), "txA", &sa, pos, false);
        pair(&mut f, &format!("amb{i}"), "txB", &sa, pos, true);
        pos += 400;
    }
    (path, fasta)
}

#[test]
/// The error model must prefer the alignment that really matches, and must do so
/// identically on every run.
fn deterministic_bam_error_model_prefers_correct_and_reproduces() {
    let dir = tempfile::tempdir().unwrap();
    let (sam, fasta) = write_seq_fixture(dir.path());

    let run = |tag: &str, error_model: bool| -> Vec<(String, f64)> {
        let out = dir.path().join(tag);
        std::fs::create_dir_all(&out).unwrap();
        let mut opts = AlignQuantOptions::new(sam.clone(), out.clone());
        opts.lib_type = "IU".to_string();
        // The error model is opt-in (`--errorModel`) and needs a transcriptome.
        opts.deterministic_error_model = error_model;
        if error_model {
            opts.transcripts = Some(fasta.clone());
        }
        let rad = dir.path().join(format!("{tag}.rad"));
        write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
        quantify_rad(&opts, &rad).unwrap();
        let mut c = quant_counts(&out);
        c.sort_by(|a, b| a.0.cmp(&b.0));
        c
    };

    // The error model must be reproducible run-to-run (order-independent counting
    // + fixed-model scoring).
    assert_eq!(
        run("em1", true),
        run("em2", true),
        "error-model quant must reproduce"
    );

    // Ambiguous fragments' reads are from txA, so the error model should pull them
    // toward txA versus the AS-less default (which weighs both placements equally).
    let em = run("em", true);
    let asm = run("as", false);
    let get = |v: &[(String, f64)], name: &str| v.iter().find(|(n, _)| n == name).unwrap().1;
    assert!(
        get(&em, "txA") > get(&asm, "txA") + 1.0,
        "error model should assign more amb mass to txA: em txA={}, as txA={}",
        get(&em, "txA"),
        get(&asm, "txA")
    );
    // Mass is still conserved.
    let total: f64 = em.iter().map(|(_, c)| c).sum();
    assert!(
        (total - 90.0).abs() < 1e-6,
        "Σ NumReads = {total}, expected 90"
    );
}

#[test]
/// The whole point of the mode: identical input must give identical output.
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
    assert_eq!(
        run("r1"),
        run("r2"),
        "deterministic alignment quant must reproduce"
    );
}

// ---- AS-tag diagnostics ---------------------------------------------------

/// Deterministic alignment mode without `--errorModel` scores each placement by
/// the BAM's `AS`. These tests cover what the run says when that tag cannot
/// carry the weight, because the failure is otherwise invisible: a missing tag
/// reads as score 0, every placement ties, and multireads end up apportioned by
/// effective length alone. On a stripped-`AS` BAM that measured three times the
/// MARD of the same file with its tags intact (COMBINE-lab/salmon#1140).
///
/// `policy` decides the tag written on each record: `None` writes none.
fn write_sam_with_as(dir: &Path, name: &str, policy: impl Fn(usize) -> Option<i32>) -> PathBuf {
    let path = dir.join(name);
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:queryname").unwrap();
    writeln!(f, "@SQ\tSN:txA\tLN:5000").unwrap();
    writeln!(f, "@SQ\tSN:txB\tLN:5000").unwrap();
    let mut pos = 1usize;
    for i in 0..40 {
        let (l, r) = (pos + 1, pos + 201);
        let tag = |idx: usize| match policy(idx) {
            Some(score) => format!("\tAS:i:{score}"),
            None => String::new(),
        };
        writeln!(
            f,
            "a{i}\t99\ttxA\t{l}\t60\t100M\t=\t{r}\t300\t*\t*{}",
            tag(2 * i)
        )
        .unwrap();
        writeln!(
            f,
            "a{i}\t147\ttxA\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*{}",
            tag(2 * i + 1)
        )
        .unwrap();
        pos += 400;
    }
    path
}

fn warning_for(dir: &Path, name: &str, policy: impl Fn(usize) -> Option<i32>) -> Option<String> {
    let sam = write_sam_with_as(dir, name, policy);
    let out = dir.join(format!("{name}.out"));
    std::fs::create_dir_all(&out).unwrap();
    let opts = opts_for(&sam, &out);
    let rad = dir.join(format!("{name}.rad"));
    write_alignment_rad(&opts, &rad, ChunkCodec::None)
        .unwrap()
        .score_tag_warning
}

#[test]
/// A BAM with no `AS` at all cannot be scored, and the run has to say so rather
/// than reporting numbers derived from an all-zero score column.
fn an_as_less_bam_is_reported() {
    let dir = tempfile::tempdir().unwrap();
    let warning = warning_for(dir.path(), "no_as.sam", |_| None).expect("a warning");
    assert!(
        warning.contains("no alignment in this BAM carries an AS tag"),
        "{warning}"
    );
    assert!(
        warning.contains("--errorModel"),
        "the warning has to name the way out: {warning}"
    );
}

#[test]
/// One constant `AS` conveys exactly as much as none, and has to be caught at
/// the record level: a proper pair sums to twice an orphan's score, so judging
/// constancy on the per-placement sum would call a constant BAM varied.
fn a_constant_as_is_reported_too() {
    let dir = tempfile::tempdir().unwrap();
    let warning = warning_for(dir.path(), "const_as.sam", |_| Some(100)).expect("a warning");
    assert!(
        warning.contains("every alignment in this BAM reports the same AS (100)"),
        "{warning}"
    );
}

#[test]
/// A partially tagged BAM is worse than either: the untagged records score 0
/// against real scores, which is not a comparison between them.
fn a_partially_tagged_bam_is_reported() {
    let dir = tempfile::tempdir().unwrap();
    let warning = warning_for(dir.path(), "partial_as.sam", |i| {
        (i % 2 == 0).then_some(200 - i as i32)
    })
    .expect("a warning");
    assert!(warning.contains("carry no AS tag"), "{warning}");
}

#[test]
/// Varying tags are exactly what the mode expects, so nothing is said.
fn usable_as_tags_draw_no_complaint() {
    let dir = tempfile::tempdir().unwrap();
    let warning = warning_for(dir.path(), "good_as.sam", |i| Some(150 + (i as i32 % 7)));
    assert!(warning.is_none(), "unexpected complaint: {warning:?}");
}

#[test]
/// With `--errorModel` the scores come from the alignments' own bases, so the
/// state of the `AS` tags is irrelevant and saying anything about them would be
/// noise.
fn the_error_model_path_says_nothing_about_as() {
    let dir = tempfile::tempdir().unwrap();
    let (sam, fasta) = write_seq_fixture(dir.path());
    let out = dir.path().join("em.out");
    std::fs::create_dir_all(&out).unwrap();
    let mut opts = AlignQuantOptions::new(sam, out);
    opts.lib_type = "IU".to_string();
    opts.deterministic_error_model = true;
    opts.transcripts = Some(fasta);
    let rad = dir.path().join("em.rad");
    let summary = write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
    assert!(summary.score_tag_warning.is_none());
}
