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

/// The baked library-format tag records what the pass *observed* (detected
/// wins), matching the reads-mode deterministic writer — so an explicit `-l`
/// that contradicts the alignments fires `library_type_mismatch` in phase 2.
/// The writer used to bake the user's assertion (expected wins), which both
/// silenced that diagnostic on this path and echoed the assertion back as
/// `detected_library_type` (#1140, audit F18).
#[test]
fn a_wrong_explicit_lib_type_is_diagnosed() {
    let dir = tempfile::tempdir().unwrap();
    let sam = write_sam(dir.path()); // fixture pairs are ISF

    let out = dir.path().join("wrong");
    let rad = dir.path().join("wrong.rad");
    let mut opts = opts_for(&sam, &out);
    opts.lib_type = "ISR".to_string();
    write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
    quantify_rad(&opts, &rad).unwrap();
    let meta = std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap();
    assert!(
        meta.contains("library_type_mismatch"),
        "an ISR assertion over ISF alignments must be diagnosed: {meta}"
    );
    assert!(
        meta.contains("observed format 'ISF'"),
        "the diagnostic must name what was actually observed: {meta}"
    );

    // Control: the matching -l draws no complaint.
    let out_ok = dir.path().join("right");
    let rad_ok = dir.path().join("right.rad");
    let mut opts_ok = opts_for(&sam, &out_ok);
    opts_ok.lib_type = "ISF".to_string();
    write_alignment_rad(&opts_ok, &rad_ok, ChunkCodec::None).unwrap();
    quantify_rad(&opts_ok, &rad_ok).unwrap();
    let meta_ok = std::fs::read_to_string(out_ok.join("aux_info").join("meta_info.json")).unwrap();
    assert!(
        !meta_ok.contains("library_type_mismatch"),
        "a correct explicit -l must not be flagged: {meta_ok}"
    );
}

/// `extra_diagnostics` is documented as appended on every path, but only
/// `quantify_rad` honoured it — the online alignment path dropped the
/// driver's diagnostics silently (#1140, audit F20). The append now lives in
/// the shared output writer; this pins the alignment path.
#[test]
fn driver_diagnostics_reach_meta_on_the_alignment_path() {
    let dir = tempfile::tempdir().unwrap();
    let sam = write_sam(dir.path());
    let out = dir.path().join("q");
    let mut opts = opts_for(&sam, &out);
    opts.extra_diagnostics.push(salmon_core::Diagnostic::new(
        "test_probe_diagnostic",
        "warning",
        "carried across from the driver".to_string(),
    ));
    salmon_align::quantify_alignments(&opts).unwrap();
    let meta = std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap();
    assert!(
        meta.contains("test_probe_diagnostic"),
        "driver-supplied diagnostics must reach meta_info.json on the online -a path: {meta}"
    );
}

// ---- flags that were inert on this path (#1140, audit batches C and D) -----

/// A SAM with proper pairs on txA plus half-mapped fragments on txB, so
/// `--discardOrphans` has something to discard and something to keep.
fn write_sam_with_orphans(dir: &Path) -> PathBuf {
    let path = dir.join("orphans.sam");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:queryname").unwrap();
    writeln!(f, "@SQ\tSN:txA\tLN:5000").unwrap();
    writeln!(f, "@SQ\tSN:txB\tLN:5000").unwrap();
    let mut pos = 1usize;
    for i in 0..40 {
        let (l, r) = (pos + 1, pos + 201);
        writeln!(f, "p{i}\t99\ttxA\t{l}\t60\t100M\t=\t{r}\t300\t*\t*").unwrap();
        writeln!(f, "p{i}\t147\ttxA\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*").unwrap();
        pos += 400;
    }
    // Half-mapped: only read 1 is placed (0x1 paired, 0x8 mate unmapped).
    let mut opos = 1usize;
    for i in 0..25 {
        let l = opos + 1;
        writeln!(f, "o{i}\t73\ttxB\t{l}\t60\t100M\t*\t0\t0\t*\t*").unwrap();
        opos += 400;
    }
    path
}

#[test]
/// `--discardOrphans` was read only by the online path, so on the 2.6.0
/// default (`-a` deterministic) it was accepted and silently did nothing
/// (#1140, audit C10). Dropping the half-mapped fragments before they reach
/// the RAD makes both alignment paths answer the flag the same way.
fn discard_orphans_is_honoured_on_the_deterministic_path() {
    let dir = tempfile::tempdir().unwrap();
    let sam = write_sam_with_orphans(dir.path());

    let run = |tag: &str, discard: bool| -> (u64, Vec<(String, f64)>) {
        let out = dir.path().join(tag);
        std::fs::create_dir_all(&out).unwrap();
        let mut opts = opts_for(&sam, &out);
        opts.discard_orphans = discard;
        let rad = dir.path().join(format!("{tag}.rad"));
        let summary = write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
        quantify_rad(&opts, &rad).unwrap();
        let mut counts = quant_counts(&out);
        counts.sort_by(|a, b| a.0.cmp(&b.0));
        (summary.num_mapped, counts)
    };

    let (kept_mapped, kept) = run("keep", false);
    let (dropped_mapped, dropped) = run("drop", true);

    assert_eq!(kept_mapped, 65, "40 proper pairs + 25 orphans");
    assert_eq!(
        dropped_mapped, 40,
        "--discardOrphans must leave only the proper pairs"
    );

    let txb_kept = kept.iter().find(|(n, _)| n == "txB").unwrap().1;
    let txb_dropped = dropped.iter().find(|(n, _)| n == "txB").unwrap().1;
    assert!(
        txb_kept > 1.0,
        "the orphan-only transcript must carry mass when orphans are kept: {txb_kept}"
    );
    assert!(
        txb_dropped < 1e-6,
        "and none when they are discarded: {txb_dropped}"
    );
}

/// Proper pairs on txC (which train the FLD) plus orphan-only fragments
/// contested between a long transcript and the 3' end of a short one.
///
/// The contest has to be built this way. The RAD carries one map type per
/// *fragment*, so a fragment cannot be a proper pair on one transcript and an
/// orphan on another; and an orphan's modelled weight only differs between
/// transcripts when the length it bounds falls inside the FLD's support — at
/// position 3850 of a 4000-base transcript only ~150 bases remain, which
/// cannot host a ~300-base fragment, while the long transcript is unbounded.
/// The two lengths are close enough (effective lengths within 2x) that the
/// orphan term — which caps an orphan's penalty at `ln 2`, since its mate may
/// lie off the transcript — still decides the contest.
fn write_sam_contested_orphans(dir: &Path) -> PathBuf {
    let path = dir.join("contested.sam");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:queryname").unwrap();
    writeln!(f, "@SQ\tSN:txA\tLN:5000").unwrap();
    writeln!(f, "@SQ\tSN:txB\tLN:4000").unwrap();
    writeln!(f, "@SQ\tSN:txC\tLN:5000").unwrap();
    let mut pos = 1usize;
    for i in 0..60 {
        let (l, r) = (pos + 1, pos + 201);
        writeln!(f, "p{i}\t99\ttxC\t{l}\t60\t100M\t=\t{r}\t300\t*\t*").unwrap();
        writeln!(f, "p{i}\t147\ttxC\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*").unwrap();
        pos += 400;
    }
    for i in 0..60 {
        // flag 73 / 329: paired, mate unmapped, first segment (329 adds
        // secondary) — a genuine orphan placement on each transcript.
        let a = (i * 13) % 500 + 1;
        let b = 3850 + (i % 40);
        writeln!(f, "o{i}\t73\ttxA\t{a}\t60\t100M\t*\t0\t0\t*\t*").unwrap();
        writeln!(f, "o{i}\t329\ttxB\t{b}\t60\t100M\t*\t0\t0\t*\t*").unwrap();
    }
    path
}

#[test]
/// `--noSingleFragProb` stopped at the mapping pass, so it was inert on every
/// RAD path — including the 2.6.0 default (#1140, audit D14).
///
/// Modelled, an orphan's bounded-CMF weight knows that txB's remaining ~150
/// bases cannot hold this library's ~300-base fragments, so the contested mass
/// goes to txA (the penalty is capped at `ln 2`, which still outweighs the
/// two transcripts' effective-length difference). With the flag, both placements take the same flat penalty and
/// the split falls back to effective length, which favours the shorter txB.
/// The flag flipping the answer this completely is the point: it was silently
/// doing nothing.
fn no_single_frag_prob_changes_orphan_weighting_on_the_rad_path() {
    let dir = tempfile::tempdir().unwrap();
    let sam = write_sam_contested_orphans(dir.path());
    let rad = dir.path().join("shared.rad");
    let base = opts_for(&sam, &dir.path().join("unused"));
    write_alignment_rad(&base, &rad, ChunkCodec::None).unwrap();

    let quant = |tag: &str, model: bool| -> Vec<(String, f64)> {
        let out = dir.path().join(tag);
        std::fs::create_dir_all(&out).unwrap();
        let mut opts = opts_for(&sam, &out);
        opts.model_single_frag_prob = model;
        quantify_rad(&opts, &rad).unwrap();
        let mut c = quant_counts(&out);
        c.sort_by(|a, b| a.0.cmp(&b.0));
        c
    };
    let get = |v: &[(String, f64)], name: &str| v.iter().find(|(n, _)| n == name).unwrap().1;

    let modelled = quant("modelled", true);
    let flat = quant("flat", false);

    assert!(
        get(&modelled, "txA") > 59.0 && get(&modelled, "txB") < 1.0,
        "modelled: the bounded-CMF weight must reject the 3'-end placement: {modelled:?}"
    );
    assert!(
        get(&flat, "txB") > 59.0 && get(&flat, "txA") < 1.0,
        "flat: equal orphan weights leave effective length to decide: {flat:?}"
    );
}

/// Outward pairs (read 1 forward downstream of read 2) under `-l OU`. A RAD hit
/// keeps only the fragment's leftmost position, so the reader cannot see they
/// point outward and used to call them inward, which `-l OU` then dropped whole.
#[test]
fn an_outward_library_maps_on_the_rad_path() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("ou.sam");
    let mut f = std::fs::File::create(&sam).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:queryname").unwrap();
    writeln!(f, "@SQ\tSN:txA\tLN:5000").unwrap();
    for i in 0..200 {
        let l = 1 + (i * 7) % 3000;
        let r = l + 200;
        writeln!(
            f,
            "r{i}\t97\ttxA\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:-10"
        )
        .unwrap();
        writeln!(
            f,
            "r{i}\t145\ttxA\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:-10"
        )
        .unwrap();
    }
    drop(f);
    let out = dir.path().join("q");
    std::fs::create_dir_all(&out).unwrap();
    let rad = out.join("aln.rad");
    let mut opts = opts_for(&sam, &out);
    opts.lib_type = "OU".to_string();
    write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
    assert_eq!(quantify_rad(&opts, &rad).unwrap().num_mapped, 200);
    // `--hardFilter` judges strand the same way, or no placement sets its bar.
    opts.hard_filter = true;
    assert_eq!(quantify_rad(&opts, &rad).unwrap().num_mapped, 200);
}

/// `--hardFilter` keeps a fragment's best placements. A wrong-strand placement
/// that `-l` drops anyway must not set that bar, or the fragment is lost whole.
#[test]
fn hard_filter_ignores_the_score_of_a_wrong_strand_placement() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("hf.sam");
    let mut f = std::fs::File::create(&sam).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:queryname").unwrap();
    writeln!(f, "@SQ\tSN:txA\tLN:5000").unwrap();
    writeln!(f, "@SQ\tSN:txB\tLN:5000").unwrap();
    for i in 0..1000 {
        let l = 1 + (i * 3) % 3000;
        let r = l + 200;
        // ISF on txA, worse score; ISR on txB, better score.
        writeln!(
            f,
            "r{i}\t99\ttxA\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:-10"
        )
        .unwrap();
        writeln!(
            f,
            "r{i}\t147\ttxA\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:-10"
        )
        .unwrap();
        writeln!(
            f,
            "r{i}\t339\ttxB\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:0"
        )
        .unwrap();
        writeln!(
            f,
            "r{i}\t419\ttxB\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:0"
        )
        .unwrap();
    }
    drop(f);
    let out = dir.path().join("q");
    std::fs::create_dir_all(&out).unwrap();
    let rad = out.join("aln.rad");
    let mut opts = opts_for(&sam, &out);
    opts.lib_type = "ISF".to_string();
    opts.hard_filter = true;
    write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
    let res = quantify_rad(&opts, &rad).unwrap();
    let counts = quant_counts(&out);
    let get = |n: &str| counts.iter().find(|(x, _)| x == n).unwrap().1;
    assert_eq!(res.num_mapped, 1000);
    assert!(
        (get("txA") - 1000.0).abs() < 1e-6 && get("txB") < 1e-6,
        "{counts:?}"
    );
}

/// Quantify `body` (SAM records over txA and txB) under `lib`, returning the
/// counts of txA and txB and `num_mapped`.
fn quant_records(dir: &Path, lib: &str, body: &str) -> (f64, f64, u64) {
    let sam = dir.join("mix.sam");
    std::fs::write(
        &sam,
        format!("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:txA\tLN:5000\n@SQ\tSN:txB\tLN:5000\n{body}"),
    )
    .unwrap();
    let out = dir.join("q");
    std::fs::create_dir_all(&out).unwrap();
    let rad = out.join("aln.rad");
    let mut opts = opts_for(&sam, &out);
    opts.lib_type = lib.to_string();
    write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
    let res = quantify_rad(&opts, &rad).unwrap();
    let counts = quant_counts(&out);
    let get = |n: &str| counts.iter().find(|(x, _)| x == n).unwrap().1;
    (get("txA"), get("txB"), res.num_mapped)
}

/// R1 alone on txA and R2 alone on txB: each placement must be judged on its
/// own mate's strand, whatever the other one is. Both agree with ISR here, so
/// neither may be filtered and the two transcripts are symmetric.
#[test]
fn left_and_right_orphans_of_one_fragment_keep_their_own_side() {
    let dir = tempfile::tempdir().unwrap();
    let mut body = String::new();
    for i in 0..1000 {
        // Mid-transcript, so neither side's weight is cut by a transcript end.
        let p = 2000 + (i * 3) % 1000;
        body += &format!("r{i}\t89\ttxA\t{p}\t60\t100M\t*\t0\t0\t*\t*\tAS:i:0\n");
        body += &format!("r{i}\t137\ttxB\t{p}\t60\t100M\t*\t0\t0\t*\t*\tAS:i:0\n");
    }
    let (a, b, mapped) = quant_records(dir.path(), "ISR", &body);
    assert_eq!(mapped, 1000);
    assert!((a - 500.0).abs() < 1.0, "txA={a} txB={b}");
    assert!((b - 500.0).abs() < 1.0, "txA={a} txB={b}");
}

/// A proper pair on txA beats an orphan on txB: the orphan's read length must
/// not be read back as a fragment length. It is the pairs' 300 here, so that
/// misreading would give txB as much weight as txA.
#[test]
fn an_orphan_beside_a_proper_pair_takes_no_mass() {
    let dir = tempfile::tempdir().unwrap();
    let mut body = String::new();
    for i in 0..500 {
        let l = 1 + (i * 3) % 3000;
        let r = l + 200;
        body += &format!("r{i}\t99\ttxA\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:0\n");
        body += &format!("r{i}\t147\ttxA\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:0\n");
        body += &format!("r{i}\t73\ttxB\t{l}\t60\t300M\t*\t0\t0\t*\t*\tAS:i:0\n");
    }
    let (a, b, mapped) = quant_records(dir.path(), "IU", &body);
    assert_eq!(mapped, 500);
    assert!((a - 500.0).abs() < 1.0 && b < 1.0, "txA={a} txB={b}");
}

/// Single-end SAM: `n_a` forward reads on txA, `n_b` reverse reads on txB, so
/// the library is stranded (SF) and txB's reads are all antisense.
fn write_single_end_sam(dir: &Path, n_a: usize, n_b: usize) -> PathBuf {
    let path = dir.join("se.sam");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:queryname").unwrap();
    writeln!(f, "@SQ\tSN:txA\tLN:5000").unwrap();
    writeln!(f, "@SQ\tSN:txB\tLN:5000").unwrap();
    for i in 0..n_a + n_b {
        let (tx, flag) = if i < n_a { ("txA", 0) } else { ("txB", 16) };
        let p = 1 + (i * 7) % 4000;
        writeln!(
            f,
            "r{i}\t{flag}\t{tx}\t{p}\t60\t100M\t*\t0\t0\t*\t*\tAS:i:0"
        )
        .unwrap();
    }
    path
}

fn meta_info(out: &Path) -> serde_json::Value {
    let text = std::fs::read_to_string(out.join("aux_info/meta_info.json")).unwrap();
    serde_json::from_str(&text).unwrap()
}

/// `-l A` on a stranded single-end BAM must find the strand from the reads, as
/// `-l SF` is told it, and filter antisense reads the same way.
#[test]
fn auto_detects_strand_of_a_single_end_bam() {
    let dir = tempfile::tempdir().unwrap();
    let sam = write_single_end_sam(dir.path(), 1500, 500);
    let quant = |lib: &str| {
        let out = dir.path().join(lib);
        std::fs::create_dir_all(&out).unwrap();
        let rad = out.join("aln.rad");
        let mut opts = opts_for(&sam, &out);
        opts.lib_type = lib.to_string();
        write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
        let res = quantify_rad(&opts, &rad).unwrap();
        (res.num_mapped, quant_counts(&out), meta_info(&out))
    };

    let (auto_mapped, auto_counts, auto_meta) = quant("A");
    let (sf_mapped, sf_counts, sf_meta) = quant("SF");

    assert_eq!(auto_meta["detected_library_type"], "SF");
    assert_eq!(auto_mapped, sf_mapped);
    assert_eq!(auto_counts, sf_counts);
    assert_eq!(
        auto_meta["frag_length_source"],
        sf_meta["frag_length_source"]
    );
}

/// A single-end read is not an orphan, so `--discardOrphans` under `-l A` must
/// not empty the run.
#[test]
fn discard_orphans_keeps_a_single_end_bam_under_auto() {
    let dir = tempfile::tempdir().unwrap();
    let sam = write_single_end_sam(dir.path(), 300, 0);
    let out = dir.path().join("q");
    std::fs::create_dir_all(&out).unwrap();
    let rad = out.join("aln.rad");
    let mut opts = opts_for(&sam, &out);
    opts.lib_type = "A".to_string();
    opts.discard_orphans = true;
    write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
    let res = quantify_rad(&opts, &rad).unwrap();
    assert_eq!(res.num_mapped, 300);
}

/// The genome path takes the same answer from the projection: reads running
/// with the transcript are sense, on either genomic strand.
#[test]
fn auto_detects_strand_of_a_single_end_genome_bam() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = dir.path().join("a.gtf");
    let mut f = std::fs::File::create(&gtf).unwrap();
    for (tx, start, end, strand) in [("txA", 1001, 3000, '+'), ("txB", 5001, 7000, '-')] {
        let attrs = format!("gene_id \"g{tx}\"; transcript_id \"{tx}\";");
        for kind in ["transcript", "exon"] {
            writeln!(
                f,
                "chr1\tt\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"
            )
            .unwrap();
        }
    }
    let sam = dir.path().join("g.sam");
    let mut f = std::fs::File::create(&sam).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:queryname").unwrap();
    writeln!(f, "@SQ\tSN:chr1\tLN:10000").unwrap();
    for i in 0..200 {
        let (flag, p) = if i % 2 == 0 { (0, 1001) } else { (16, 5001) };
        writeln!(
            f,
            "r{i}\t{flag}\tchr1\t{}\t60\t100M\t*\t0\t0\t*\t*\tAS:i:0",
            p + i
        )
        .unwrap();
    }

    let out = dir.path().join("q");
    std::fs::create_dir_all(&out).unwrap();
    let rad = out.join("aln.rad");
    let gp = salmon_align::GenomeProjectOptions {
        bam: sam.clone(),
        annotation: gtf,
        genome_fasta: None,
        lib_type: "A".to_string(),
        junc_miss_discount: None,
        bias: false,
        fld_mean: 250.0,
        fld_sd: 25.0,
        fld_max: 1000,
        bias_seed_em_iters: 0,
        em: Default::default(),
        rad_codec: ChunkCodec::None,
    };
    salmon_align::project_genome_bam_to_rad(&gp, &rad).unwrap();
    let mut opts = opts_for(&sam, &out);
    opts.lib_type = "A".to_string();
    let res = quantify_rad(&opts, &rad).unwrap();
    assert_eq!(res.num_mapped, 200);
    assert_eq!(meta_info(&out)["detected_library_type"], "SF");
}

/// A secondary record often has no bases (SEQ `*`). Scored as nothing, it got
/// log-likelihood 0, which beats any real alignment with a mismatch, so reads moved
/// to whichever transcript was reported second. The bases of the primary record
/// are the same read, so they score it: here txA really is the better placement.
#[test]
fn a_secondary_without_bases_is_scored_from_its_primary() {
    let dir = tempfile::tempdir().unwrap();
    let random = |mut x: u64, n: usize| -> Vec<u8> {
        (0..n)
            .map(|_| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize & 3]
            })
            .collect()
    };
    // Unrelated transcripts: a read from txA does not fit txB at all.
    let tx_a = random(5, 3000);
    let tx_b = random(19, 3000);
    let flip = |b: u8| if b == b'A' { b'C' } else { b'A' };
    let fasta = dir.path().join("t.fa");
    std::fs::write(
        &fasta,
        format!(
            ">txA\n{}\n>txB\n{}\n",
            String::from_utf8_lossy(&tx_a),
            String::from_utf8_lossy(&tx_b)
        ),
    )
    .unwrap();

    let sam = dir.path().join("a.sam");
    let mut f = std::fs::File::create(&sam).unwrap();
    writeln!(f, "@HD\tVN:1.6\tSO:unsorted").unwrap();
    writeln!(f, "@SQ\tSN:txA\tLN:3000").unwrap();
    writeln!(f, "@SQ\tSN:txB\tLN:3000").unwrap();
    for i in 0..1000 {
        let p = (i * 3) % 2800;
        // Two sequencing errors, so the read does not fit even txA perfectly.
        let mut read = tx_a[p..p + 100].to_vec();
        read[30] = flip(read[30]);
        read[70] = flip(read[70]);
        let read = String::from_utf8_lossy(&read);
        let (q, one) = ("I".repeat(100), p + 1);
        writeln!(f, "r{i}\t0\ttxA\t{one}\t60\t100M\t*\t0\t0\t{read}\t{q}").unwrap();
        writeln!(f, "r{i}\t256\ttxB\t{one}\t0\t100M\t*\t0\t0\t*\t*").unwrap();
    }
    drop(f);

    let out = dir.path().join("q");
    let mut opts = AlignQuantOptions::new(sam, out.clone());
    opts.lib_type = "U".to_string();
    opts.deterministic_error_model = true;
    opts.transcripts = Some(fasta);
    let rad = dir.path().join("a.rad");
    write_alignment_rad(&opts, &rad, ChunkCodec::None).unwrap();
    quantify_rad(&opts, &rad).unwrap();
    let counts = quant_counts(&out);
    assert!(
        counts[0].1 > counts[1].1,
        "txA fits the reads better, so it must get more of them: {counts:?}"
    );
}
