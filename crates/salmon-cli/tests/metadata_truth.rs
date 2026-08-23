//! Groups 4 and 5 of the pre-2.6.0 readiness sweep: `meta_info.json` fields and
//! log text that asserted something the run did not do, and errors that arrived
//! too late or named the wrong flag. All machine-read or user-facing, all on
//! default paths.

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

/// Decoy-aware index plus ISF-by-construction read pairs.
fn fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let (a, b) = (sequence(5, 3000), sequence(19, 3000));
    let decoy = sequence(101, 4000);
    let fasta = dir.join("gentrome.fa");
    std::fs::write(&fasta, format!(">txA\n{a}\n>txB\n{b}\n>decoy1\n{decoy}\n")).unwrap();
    let decoys = dir.join("decoys.txt");
    std::fs::write(&decoys, "decoy1\n").unwrap();
    let index = dir.join("idx");
    assert!(Command::new(SALMON)
        .args(["index", "-t"])
        .arg(&fasta)
        .arg("-d")
        .arg(&decoys)
        .arg("-i")
        .arg(&index)
        .args(["-p", "2"])
        .status()
        .unwrap()
        .success());

    let rc = |x: &str| -> String {
        x.chars()
            .rev()
            .map(|c| match c {
                'A' => 'T',
                'C' => 'G',
                'G' => 'C',
                _ => 'A',
            })
            .collect()
    };
    let (mut r1, mut r2) = (String::new(), String::new());
    let q = "I".repeat(100);
    for i in 0..200 {
        let s = if i % 2 == 0 { &a } else { &b };
        let st = (i * 11) % (s.len() - 350);
        r1.push_str(&format!("@f{i}\n{}\n+\n{q}\n", &s[st..st + 100]));
        r2.push_str(&format!("@f{i}\n{}\n+\n{q}\n", rc(&s[st + 200..st + 300])));
    }
    let (p1, p2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    std::fs::write(&p1, r1).unwrap();
    std::fs::write(&p2, r2).unwrap();
    (index, p1, p2)
}

fn meta(out: &Path) -> serde_json::Value {
    let s = std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap();
    serde_json::from_str(&s).unwrap()
}

/// `library_types` reported what the DETECTOR observed rather than what the run
/// applied, so an explicit `-l` contradicted `cmd_info.json` and
/// `lib_format_counts.json` in the same directory. Both facts are wanted, so
/// they now live in separate fields: applied in `library_types`, observed in
/// `detected_library_type`.
#[test]
fn library_types_reports_the_applied_format_and_detected_the_observed_one() {
    let dir = tempfile::tempdir().unwrap();
    let (index, r1, r2) = fixture(dir.path());
    let run = |lib: &str, out: &Path| {
        assert!(Command::new(SALMON)
            .args(["quant", "-i"])
            .arg(&index)
            .args(["-l", lib, "-1"])
            .arg(&r1)
            .arg("-2")
            .arg(&r2)
            .args(["-p", "2", "-o"])
            .arg(out)
            .status()
            .unwrap()
            .success());
    };

    // The reads are ISF. Declaring ISR must be reported as ISR — the run really
    // did filter that way (and mapped nothing as a result).
    let wrong = dir.path().join("isr");
    run("ISR", &wrong);
    let m = meta(&wrong);
    assert_eq!(m["library_types"][0], "ISR", "applied format");
    assert_eq!(m["detected_library_type"], "ISF", "observed format");
    assert_eq!(m["num_mapped"], 0, "ISR really was applied");

    // Under `-l A` detection *is* the declaration, so both agree.
    let auto = dir.path().join("auto");
    run("A", &auto);
    let m = meta(&auto);
    assert_eq!(m["library_types"][0], "ISF");
    assert_eq!(m["detected_library_type"], "ISF");
}

/// `num_decoy_targets` was hardcoded 0 while `num_decoy_fragments` beside it was
/// real — self-contradictory, with `meta_info_complete: true`.
#[test]
fn decoy_targets_are_counted() {
    let dir = tempfile::tempdir().unwrap();
    let (index, r1, r2) = fixture(dir.path());
    let out = dir.path().join("q");
    assert!(Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "IU", "-1"])
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .args(["-p", "2", "-o"])
        .arg(&out)
        .status()
        .unwrap()
        .success());
    let m = meta(&out);
    assert_eq!(m["num_valid_targets"], 2, "txA + txB");
    assert_eq!(
        m["num_decoy_targets"], 1,
        "decoy1 must be counted, not zeroed"
    );
}

/// A plain `-a` run reported `frag_length_source: "rad_baked"`, naming an
/// internal intermediate the user never asked for. Genuine `--rad` input keeps
/// that spelling, because there the RAD really is the source.
#[test]
fn frag_length_source_names_the_users_input() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("a.sam");
    let mut s = String::from("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:txA\tLN:3000\n");
    for i in 0..40 {
        let (l, r) = (i * 20 + 1, i * 20 + 201);
        s.push_str(&format!(
            "p{i}\t99\ttxA\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:200\n"
        ));
        s.push_str(&format!(
            "p{i}\t147\ttxA\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:200\n"
        ));
    }
    std::fs::write(&sam, s).unwrap();

    let aln = dir.path().join("aln");
    let rad = dir.path().join("m.rad");
    assert!(Command::new(SALMON)
        .args(["quant", "-a"])
        .arg(&sam)
        .args(["-l", "IU", "--writeRad"])
        .arg(&rad)
        .arg("-o")
        .arg(&aln)
        .status()
        .unwrap()
        .success());
    assert_eq!(
        meta(&aln)["frag_length_source"],
        "alignments",
        "alignment input must not report the internal RAD handoff"
    );

    let from_rad = dir.path().join("rad");
    assert!(Command::new(SALMON)
        .args(["quant", "--rad"])
        .arg(&rad)
        .args(["-l", "IU", "-o"])
        .arg(&from_rad)
        .status()
        .unwrap()
        .success());
    assert_eq!(
        meta(&from_rad)["frag_length_source"],
        "rad_baked",
        "genuine RAD input keeps rad_baked — the RAD really is the source"
    );
}

/// `logs/salmon_quant.log` labelled a default reads run "alignment mode",
/// dropped the documented mapping-type/rate lines, and called unmapped
/// fragments strand-incompatible on a run with none.
#[test]
fn the_run_log_describes_the_mode_the_user_chose() {
    let dir = tempfile::tempdir().unwrap();
    let (index, r1, r2) = fixture(dir.path());
    let out = dir.path().join("q");
    assert!(Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "IU", "-1"])
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .args(["-p", "2", "-o"])
        .arg(&out)
        .status()
        .unwrap()
        .success());
    let log = std::fs::read_to_string(out.join("logs").join("salmon_quant.log")).unwrap();
    assert!(log.contains("reads mode"), "{log}");
    assert!(!log.contains("alignment mode"), "{log}");
    assert!(
        log.contains("mapping type:"),
        "documented line missing:\n{log}"
    );
    assert!(
        log.contains("mapping rate:"),
        "documented line missing:\n{log}"
    );
    assert!(
        !log.contains("strand-compatible"),
        "reads mode has no strand-compatible/aligned distinction:\n{log}"
    );
}

/// `-a` with bias but no `-t` used to read the entire BAM before failing, with a
/// message naming `--rad` — a flag the invocation did not contain.
#[test]
fn alignment_bias_without_targets_fails_immediately_and_names_the_right_flags() {
    let dir = tempfile::tempdir().unwrap();
    let sam = dir.path().join("a.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:txA\tLN:3000\n\
         r0\t99\ttxA\t1\t60\t100M\t=\t201\t300\t*\t*\tAS:i:200\n\
         r0\t147\ttxA\t201\t60\t100M\t=\t1\t-300\t*\t*\tAS:i:200\n",
    )
    .unwrap();
    let o = Command::new(SALMON)
        .args(["quant", "-a"])
        .arg(&sam)
        .args(["-l", "IU", "--gcBias", "-o"])
        .arg(dir.path().join("q"))
        .output()
        .unwrap();
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("--gcBias") && err.contains("-t"), "{err}");
    assert!(
        !err.contains("--rad"),
        "the message must not name a flag the user did not pass:\n{err}"
    );
    assert!(
        !dir.path().join("q").join("quant.sf").exists(),
        "it must fail before doing the work"
    );
}

/// A FASTQ parse error named neither the file nor anything else useful, though
/// `-1/-2/-r` accept many files.
#[test]
fn a_fastq_parse_error_names_the_file() {
    let dir = tempfile::tempdir().unwrap();
    let (index, _r1, _r2) = fixture(dir.path());
    // Valid records, then one with a bad separator mid-file, then more — so the
    // failure is a genuine parse error rather than the completeness preflight.
    let bad = dir.path().join("mid.fq");
    let mut s = String::new();
    let seq = sequence(3, 60);
    for i in 0..40 {
        s.push_str(&format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(60)));
    }
    s.push_str(&format!("@bad\n{seq}\nN\n{}\n", "I".repeat(60)));
    for i in 40..60 {
        s.push_str(&format!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(60)));
    }
    std::fs::write(&bad, s).unwrap();

    let o = Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "U", "-r"])
        .arg(&bad)
        .args(["-p", "2", "-o"])
        .arg(dir.path().join("q"))
        .output()
        .unwrap();
    assert!(!o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("mid.fq"),
        "the error must name the offending file:\n{err}"
    );
}

/// The `--decoyThreshold` inertness warning sat before mode dispatch, so it
/// fired in `-a`/`--rad` — asserting sketch was running on a run that never
/// mapped, and where `--sketch` itself was already reported inert.
#[test]
fn the_sketch_decoy_threshold_warning_only_fires_in_reads_mode() {
    let dir = tempfile::tempdir().unwrap();
    let (index, r1, r2) = fixture(dir.path());
    let sam = dir.path().join("a.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:txA\tLN:3000\n\
         r0\t99\ttxA\t1\t60\t100M\t=\t201\t300\t*\t*\tAS:i:200\n\
         r0\t147\ttxA\t201\t60\t100M\t=\t1\t-300\t*\t*\tAS:i:200\n",
    )
    .unwrap();
    let out = dir.path().join("q");

    let o = Command::new(SALMON)
        .args(["quant", "-a"])
        .arg(&sam)
        .args(["-l", "IU", "--sketch", "--decoyThreshold", "0.5", "-o"])
        .arg(&out)
        .output()
        .unwrap();
    assert!(
        !String::from_utf8_lossy(&o.stderr).contains("no effect in --sketch mode"),
        "must not claim sketch is running in alignment mode"
    );

    let o = Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "IU", "-1"])
        .arg(&r1)
        .arg("-2")
        .arg(&r2)
        .args(["--sketch", "--decoyThreshold", "0.5", "-p", "2", "-o"])
        .arg(&out)
        .output()
        .unwrap();
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("no effect in --sketch mode"),
        "it must still fire where sketch really runs"
    );
}
