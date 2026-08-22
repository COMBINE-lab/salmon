//! 2.6.0 flips the default quantification path to deterministic (#1140).
//!
//! These tests pin the flip's contract at the binary level: a bare run takes
//! the deterministic two-phase path (auditable via `meta_info.json`'s
//! `inference_path`), `--online` selects the deprecated one-pass path with a
//! deprecation warning, `--deterministic` stays accepted as a no-op (it is in
//! scripts and CI configs already), the contradictory spelling errors out, and
//! online-only knobs warn instead of being silently inert under the default.

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
    let a = sequence(3, 3000);
    let b = sequence(17, 3000);
    let fasta = dir.join("txp.fa");
    std::fs::write(&fasta, format!(">txA\n{a}\n>txB\n{b}\n")).unwrap();

    let (mut r1, mut r2) = (String::new(), String::new());
    let qual = "I".repeat(100);
    for i in 0..200 {
        let source = if i % 2 == 0 { &a } else { &b };
        let start = (i * 13) % (source.len() - 300);
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

fn inference_path(out: &Path) -> String {
    let text = std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap();
    let key = "\"inference_path\":";
    let start = text.find(key).unwrap() + key.len();
    text[start..]
        .split('"')
        .nth(1)
        .expect("inference_path value")
        .to_string()
}

#[test]
fn deterministic_is_the_default_and_online_is_the_deprecated_escape_hatch() {
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

    let run = |tag: &str, extra: &[&str]| -> (PathBuf, std::process::Output) {
        let out = dir.path().join(tag);
        let mut cmd = Command::new(SALMON);
        cmd.args(["quant", "-i"])
            .arg(&index)
            .args(["-l", "IU", "-1"])
            .arg(&r1)
            .arg("-2")
            .arg(&r2)
            .args(["-p", "2", "-o"])
            .arg(&out)
            .args(extra);
        (out.clone(), cmd.output().unwrap())
    };

    // A bare run is deterministic: the auditable field says so, and the
    // intermediate RAD was cleaned up on success.
    let (out, res) = run("bare", &[]);
    assert!(res.status.success());
    assert_eq!(inference_path(&out), "deterministic");
    assert!(
        !out.join("intermediate_mappings.rad").exists(),
        "the intermediate is cleaned up on success"
    );
    assert!(out.join("quant.sf").is_file());

    // `--deterministic` stays accepted as a no-op: it is in scripts already,
    // and erroring on it would break exactly the users who adopted early.
    let (out, res) = run("explicit", &["--deterministic"]);
    assert!(res.status.success(), "--deterministic must stay accepted");
    assert_eq!(inference_path(&out), "deterministic");

    // `--online` selects the one-pass path and says it is deprecated.
    let (out, res) = run("online", &["--online"]);
    assert!(res.status.success());
    assert_eq!(inference_path(&out), "online");
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(
        err.contains("deprecated"),
        "--online must announce its deprecation: {err}"
    );

    // The contradictory spelling is an error, not a guess.
    let (_, res) = run("conflict", &["--online", "--deterministic"]);
    assert!(
        !res.status.success(),
        "--online --deterministic must not pick a side silently"
    );

    // An online-only knob without `--online` warns instead of being silently
    // inert — the failure mode this whole release cycle existed to eradicate.
    let (out, res) = run("inert", &["--forgettingFactor", "0.8"]);
    assert!(res.status.success());
    assert_eq!(inference_path(&out), "deterministic");
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(
        err.contains("--forgettingFactor"),
        "an inert online knob must be named in a warning: {err}"
    );

    // ...and the same knob under `--online` is used without complaint about
    // inertness (the only warning is the deprecation one).
    let (_, res) = run("online_knob", &["--online", "--forgettingFactor", "0.8"]);
    assert!(res.status.success());
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(
        !err.contains("ignored"),
        "an online knob under --online is not inert: {err}"
    );
}
