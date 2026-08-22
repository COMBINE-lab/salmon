//! `--dumpEq` has to describe the classes the run actually used, in every mode.
//!
//! Under `--deterministic` it did not: the mapping pass wrote the file before any
//! class existed, declaring zero of them, and the quantification pass never
//! replaced it (COMBINE-lab/salmon#1140). A well-formed file saying "0 classes"
//! is worse than a missing one, since a reader has no way to tell it is wrong,
//! which is why this test drives the binary end to end rather than any single
//! function: the bug lived in how the two passes were wired, not in the writer.

use std::io::Read as _;
use std::path::{Path, PathBuf};
use std::process::Command;

const SALMON: &str = env!("CARGO_BIN_EXE_salmon");

/// Deterministic pseudo-random bases, so the fixture is reproducible without a
/// data file.
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

/// Three transcripts, two of which share a long stretch so some fragments are
/// ambiguous and the classes are not all singletons.
fn fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let a = sequence(11, 4000);
    let b = sequence(29, 4000);
    let shared = format!("{}{}", &a[..1500], sequence(41, 2500));

    let fasta = dir.join("txp.fa");
    std::fs::write(&fasta, format!(">txA\n{a}\n>txB\n{b}\n>txC\n{shared}\n")).unwrap();

    let (mut r1, mut r2) = (String::new(), String::new());
    let qual = "I".repeat(100);
    for (i, source) in [&a, &b, &shared].iter().cycle().take(300).enumerate() {
        let start = (i * 37) % (source.len() - 300);
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

struct Dump {
    names: Vec<String>,
    declared: usize,
    classes: Vec<(Vec<u32>, u64)>,
}

fn read_dump(out: &Path) -> Dump {
    let f = std::fs::File::open(out.join("aux_info").join("eq_classes.txt.gz"))
        .expect("the dump has to exist");
    let mut text = String::new();
    flate2::read::GzDecoder::new(f)
        .read_to_string(&mut text)
        .unwrap();
    let mut lines = text.lines();
    let num_txps: usize = lines.next().unwrap().parse().unwrap();
    let declared: usize = lines.next().unwrap().parse().unwrap();
    let names: Vec<String> = (0..num_txps)
        .map(|_| lines.next().unwrap().to_string())
        .collect();
    let classes = lines
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let size: usize = f[0].parse().unwrap();
            let txps = f[1..=size].iter().map(|t| t.parse().unwrap()).collect();
            (txps, f[size + 1].parse().unwrap())
        })
        .collect();
    Dump {
        names,
        declared,
        classes,
    }
}

fn num_mapped(out: &Path) -> u64 {
    let text = std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap();
    let key = "\"num_mapped\":";
    let start = text.find(key).unwrap() + key.len();
    text[start..]
        .trim_start()
        .split(|c: char| !c.is_ascii_digit())
        .next()
        .unwrap()
        .parse()
        .unwrap()
}

fn quant(index: &Path, r1: &Path, r2: &Path, out: &Path, deterministic: bool) {
    let mut cmd = Command::new(SALMON);
    cmd.args(["quant", "-i"])
        .arg(index)
        .args(["-l", "A", "-1"])
        .arg(r1)
        .arg("-2")
        .arg(r2)
        .args(["-p", "2", "-o"])
        .arg(out)
        .arg("--dumpEq");
    if deterministic {
        cmd.arg("--deterministic");
    }
    let status = cmd.status().unwrap();
    assert!(status.success(), "quantification failed: {status}");
}

#[test]
fn dump_eq_describes_the_classes_in_both_modes() {
    let dir = tempfile::tempdir().unwrap();
    let (fasta, r1, r2) = fixture(dir.path());
    let index = dir.path().join("idx");
    let status = Command::new(SALMON)
        .args(["index", "-t"])
        .arg(&fasta)
        .arg("-i")
        .arg(&index)
        .args(["-p", "2"])
        .status()
        .unwrap();
    assert!(status.success(), "indexing failed");

    let online = dir.path().join("online");
    let det = dir.path().join("det");
    quant(&index, &r1, &r2, &online, false);
    quant(&index, &r1, &r2, &det, true);

    for (label, out) in [("one-pass", &online), ("deterministic", &det)] {
        let dump = read_dump(out);
        // The regression: the deterministic dump used to declare zero classes,
        // written by a pass that had none to declare.
        assert!(dump.declared > 0, "{label} dump declares no classes at all");
        assert_eq!(
            dump.declared,
            dump.classes.len(),
            "{label} dump declares {} classes and lists {}",
            dump.declared,
            dump.classes.len()
        );
        let total: u64 = dump.classes.iter().map(|(_, c)| c).sum();
        assert_eq!(
            total,
            num_mapped(out),
            "{label} dump has to account for every mapped fragment"
        );
        assert_eq!(dump.names, ["txA", "txB", "txC"], "{label} name block");
    }
}

#[test]
/// `--dumpEq` with `--skipQuant --deterministic` has nothing to write, and has
/// to say so.
///
/// The classes are built by the quantification pass, which `--skipQuant` skips,
/// so no file is possible. One-pass `--skipQuant --dumpEq` does produce a dump,
/// so silence here would look exactly like the flag being dropped, which is the
/// failure mode this whole batch exists to remove.
fn dump_eq_with_skip_quant_says_it_has_nothing_to_dump() {
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

    let out = dir.path().join("skipped");
    let result = Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "A", "-1"])
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .args(["-p", "2", "-o"])
        .arg(&out)
        .args(["--deterministic", "--skipQuant", "--dumpEq"])
        .output()
        .unwrap();
    assert!(result.status.success());

    assert!(
        !out.join("aux_info").join("eq_classes.txt.gz").exists(),
        "there are no classes to dump, so there must be no file pretending otherwise"
    );
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(
        stderr.contains("--dumpEq/--dumpEqWeights has nothing to write"),
        "the run has to say why the file is absent: {stderr}"
    );
    assert!(
        stderr.contains("--rad"),
        "and point at the way to get the classes from the RAD it kept: {stderr}"
    );
    assert!(
        out.join("intermediate_mappings.rad").is_file(),
        "the RAD the warning points at has to be there"
    );
}
