//! `--skipQuant --deterministic` maps and stops.
//!
//! It used to quantify anyway: phase 1 sets `skip_quant` for its own purposes
//! (map, write the RAD, do not run the online EM), so the user's request was
//! overwritten there and never reached phase 2. The run produced a full
//! `quant.sf` from a flag asking for no quantification, and then deleted the RAD
//! it had just been asked to produce (COMBINE-lab/salmon#1140).
//!
//! Driven through the binary because the bug was in how the two passes were
//! wired rather than in either of them.

use std::path::{Path, PathBuf};
use std::process::Command;

const SALMON: &str = env!("CARGO_BIN_EXE_salmon");

fn sequence(mut seed: u64, len: usize) -> String {
    (0..len)
        .map(|_| {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            "ACGT".as_bytes()[(seed >> 33) as usize & 3] as char
        })
        .collect()
}

fn fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let a = sequence(7, 3000);
    let b = sequence(23, 3000);
    let fasta = dir.join("txp.fa");
    std::fs::write(&fasta, format!(">txA\n{a}\n>txB\n{b}\n")).unwrap();

    let (mut r1, mut r2) = (String::new(), String::new());
    let qual = "I".repeat(100);
    for i in 0..200 {
        let source = if i % 2 == 0 { &a } else { &b };
        let start = (i * 11) % (source.len() - 300);
        let mate1 = &source[start..start + 100];
        let mate2: String = source[start + 200..start + 300]
            .chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                _ => 'A',
            })
            .collect();
        r1.push_str(&format!("@f{i}\n{mate1}\n+\n{qual}\n"));
        r2.push_str(&format!("@f{i}\n{mate2}\n+\n{qual}\n"));
    }
    let (p1, p2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    std::fs::write(&p1, r1).unwrap();
    std::fs::write(&p2, r2).unwrap();
    (fasta, p1, p2)
}

#[test]
fn skip_quant_under_deterministic_maps_and_stops() {
    let dir = tempfile::tempdir().unwrap();
    let (fasta, r1, r2) = fixture(dir.path());
    let index = dir.path().join("idx");
    assert!(Command::new(SALMON)
        .args(["index", "-t"])
        .arg(&fasta)
        .arg("-i")
        .arg(&index)
        .args(["-p", "2"])
        .status()
        .unwrap()
        .success());

    let run = |tag: &str, skip: bool, scratch: Option<&Path>| -> PathBuf {
        let out = dir.path().join(tag);
        let mut cmd = Command::new(SALMON);
        cmd.args(["quant", "-i"])
            .arg(&index)
            .args(["-l", "A", "-1"])
            .arg(&r1)
            .arg("-2")
            .arg(&r2)
            .args(["-p", "2", "-o"])
            .arg(&out)
            .arg("--deterministic");
        if skip {
            cmd.arg("--skipQuant");
        }
        if let Some(d) = scratch {
            cmd.arg("--radScratchDir").arg(d);
        }
        assert!(cmd.status().unwrap().success(), "{tag} run failed");
        out
    };

    let full = run("full", false, None);
    assert!(
        full.join("quant.sf").is_file(),
        "a full run still quantifies"
    );
    assert!(
        !full.join("intermediate_mappings.rad").exists(),
        "and still cleans up its intermediate"
    );

    let mapped = run("mapped", true, None);
    assert!(
        !mapped.join("quant.sf").exists(),
        "--skipQuant must not produce abundances"
    );
    let rad = mapped.join("intermediate_mappings.rad");
    assert!(
        rad.is_file(),
        "the RAD is the deliverable of a mapping-only run, not an intermediate to delete"
    );
    assert!(
        std::fs::metadata(&rad).unwrap().len() > 0,
        "and it has to have something in it"
    );

    // Scratch placement is for intermediates. Once --skipQuant promotes the RAD
    // to the run's output, sending it to a scratch volume under a pid-suffixed
    // name would be the wrong place and the wrong name for something the user
    // asked to keep.
    let scratch = dir.path().join("scratch");
    let scratched = run("scratched", true, Some(&scratch));
    assert!(
        scratched.join("intermediate_mappings.rad").is_file(),
        "a --skipQuant deliverable belongs in the output directory"
    );
    let in_scratch: Vec<_> = std::fs::read_dir(&scratch)
        .map(|rd| rd.flatten().map(|e| e.file_name()).collect())
        .unwrap_or_default();
    assert!(
        in_scratch.is_empty(),
        "nothing should be left in --radScratchDir: {in_scratch:?}"
    );
}

/// A tiny name-grouped transcriptome SAM: inward pairs, unique placements, no
/// `AS` (which alignment mode tolerates and warns about; irrelevant here).
fn write_sam(dir: &Path) -> PathBuf {
    let path = dir.join("aln.sam");
    let mut text = String::from("@HD\tVN:1.6\tSO:queryname\n");
    text.push_str("@SQ\tSN:txA\tLN:3000\n@SQ\tSN:txB\tLN:3000\n");
    for i in 0..60 {
        let tx = if i % 2 == 0 { "txA" } else { "txB" };
        let (l, r) = (1 + i * 20, 201 + i * 20);
        text.push_str(&format!(
            "p{i}\t99\t{tx}\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\n"
        ));
        text.push_str(&format!(
            "p{i}\t147\t{tx}\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\n"
        ));
    }
    std::fs::write(&path, text).unwrap();
    path
}

#[test]
/// The `-a` spelling has the same deliverable, so it needs the same placement
/// rule: `-a --deterministic --skipQuant` writes its RAD to the output
/// directory even when `--radScratchDir` points somewhere else.
fn skip_quant_under_deterministic_alignment_keeps_the_rad_in_the_output_dir() {
    let dir = tempfile::tempdir().unwrap();
    let sam = write_sam(dir.path());
    let out = dir.path().join("out");
    let scratch = dir.path().join("scratch");

    assert!(Command::new(SALMON)
        .args(["quant", "-a"])
        .arg(&sam)
        .args(["-l", "IU", "-p", "2", "-o"])
        .arg(&out)
        .args(["--deterministic", "--skipQuant", "--radScratchDir"])
        .arg(&scratch)
        .status()
        .unwrap()
        .success());

    assert!(
        !out.join("quant.sf").exists(),
        "--skipQuant must not produce abundances"
    );
    assert!(
        out.join("intermediate_alignments.rad").is_file(),
        "the deliverable belongs in the output directory under its stable name"
    );
    let in_scratch: Vec<_> = std::fs::read_dir(&scratch)
        .map(|rd| rd.flatten().map(|e| e.file_name()).collect())
        .unwrap_or_default();
    assert!(
        in_scratch.is_empty(),
        "nothing should be left in --radScratchDir: {in_scratch:?}"
    );
}

#[test]
/// Genome projection under `--skipQuant` has the same deliverable, so a phase 2
/// that fails afterwards must not delete it. Its error path cleaned the RAD up
/// as an intermediate, while its success path (and the `-a` driver) kept it.
fn skip_quant_genome_projection_keeps_the_rad_when_phase_2_fails() {
    let dir = tempfile::tempdir().unwrap();
    let gtf = dir.path().join("a.gtf");
    let mut gtf_text = String::new();
    for (tx, start, end) in [("txA", 1001, 3000), ("txB", 5001, 7000)] {
        let attrs = format!("gene_id \"g{tx}\"; transcript_id \"{tx}\";");
        for kind in ["transcript", "exon"] {
            gtf_text.push_str(&format!(
                "chr1\tt\t{kind}\t{start}\t{end}\t.\t+\t.\t{attrs}\n"
            ));
        }
    }
    std::fs::write(&gtf, gtf_text).unwrap();
    let sam = dir.path().join("g.sam");
    let mut sam_text = String::from("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:10000\n");
    for i in 0..100 {
        let p = if i % 2 == 0 { 1001 } else { 5001 } + i;
        sam_text.push_str(&format!(
            "r{i}\t0\tchr1\t{p}\t60\t100M\t*\t0\t0\t*\t*\tAS:i:0\n"
        ));
    }
    std::fs::write(&sam, sam_text).unwrap();
    let out = dir.path().join("out");
    std::fs::create_dir_all(&out).unwrap();
    // A file where phase 2 wants its `aux_info` directory makes phase 2 fail
    // after phase 1 has written the RAD.
    std::fs::write(out.join("aux_info"), "").unwrap();

    let ok = Command::new(SALMON)
        .args(["quant", "-a"])
        .arg(&sam)
        .arg("--annotation")
        .arg(&gtf)
        .args(["-l", "SF", "-p", "2", "--skipQuant", "-o"])
        .arg(&out)
        .status()
        .unwrap()
        .success();
    assert!(!ok, "phase 2 was set up to fail");
    assert!(
        out.join("genome_projected.rad").is_file(),
        "the --skipQuant deliverable must survive a failing phase 2"
    );
}

/// With an explicit `-l ISR`, the full run drops fragments whose every placement
/// is wrong-strand, so they are not mapped. `--skipQuant` reports the mapping
/// pass's own totals, which used to count them as mapped and so reported a
/// higher mapping rate than the full run for the same input
/// (COMBINE-lab/salmon#1314).
#[test]
fn skip_quant_does_not_count_wrong_strand_fragments_as_mapped() {
    let dir = tempfile::tempdir().unwrap();
    let (fasta, f1, f2) = fixture(dir.path());
    // The fixture is ISF; appending it with the mates swapped adds as many ISR
    // fragments, so half the input is wrong-strand for `-l ISR`.
    let (a, b) = (std::fs::read(&f1).unwrap(), std::fs::read(&f2).unwrap());
    let (r1, r2) = (dir.path().join("m1.fq"), dir.path().join("m2.fq"));
    std::fs::write(&r1, [a.as_slice(), &b].concat()).unwrap();
    std::fs::write(&r2, [b.as_slice(), &a].concat()).unwrap();
    let index = dir.path().join("idx");
    assert!(Command::new(SALMON)
        .args(["index", "-t"])
        .arg(&fasta)
        .arg("-i")
        .arg(&index)
        .args(["-p", "2"])
        .status()
        .unwrap()
        .success());

    let num_mapped = |tag: &str, skip: bool| -> u64 {
        let out = dir.path().join(tag);
        let mut cmd = Command::new(SALMON);
        cmd.args(["quant", "-i"])
            .arg(&index)
            .args(["-l", "ISR", "-1"])
            .arg(&r1)
            .arg("-2")
            .arg(&r2)
            .args(["-p", "2", "-o"])
            .arg(&out);
        if skip {
            cmd.arg("--skipQuant");
        }
        assert!(cmd.status().unwrap().success(), "{tag} run failed");
        let meta: serde_json::Value =
            serde_json::from_slice(&std::fs::read(out.join("aux_info/meta_info.json")).unwrap())
                .unwrap();
        meta["num_mapped"].as_u64().unwrap()
    };

    let full = num_mapped("full", false);
    assert_eq!(full, 200, "only the ISR half maps");
    assert_eq!(num_mapped("skip", true), full, "--skipQuant must agree");
}

/// `-g` under `--skipQuant` writes no gene table because there are no
/// abundances. The BAM, RAD, projection and online drivers say so; the default
/// reads driver returned before reaching that message and dropped `-g`
/// silently (COMBINE-lab/salmon#1314).
#[test]
fn skip_quant_in_default_reads_mode_says_the_gene_map_is_skipped() {
    let dir = tempfile::tempdir().unwrap();
    let (fasta, r1, r2) = fixture(dir.path());
    let index = dir.path().join("idx");
    assert!(Command::new(SALMON)
        .args(["index", "-t"])
        .arg(&fasta)
        .arg("-i")
        .arg(&index)
        .args(["-p", "2"])
        .status()
        .unwrap()
        .success());
    let gm = dir.path().join("t2g.tsv");
    std::fs::write(&gm, "txA\tgeneA\ntxB\tgeneB\n").unwrap();
    let out = dir.path().join("out");

    let o = Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "A", "-1"])
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .args(["-p", "2", "--skipQuant", "-g"])
        .arg(&gm)
        .arg("-o")
        .arg(&out)
        .output()
        .unwrap();
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(o.status.success(), "{err}");
    assert!(!out.join("quant.genes.sf").exists());
    assert!(
        err.contains("gene-level aggregation (-g) skipped"),
        "-g must not be dropped silently:\n{err}"
    );
}
