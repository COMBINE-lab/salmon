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
