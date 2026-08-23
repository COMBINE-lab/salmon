//! Regressions for the four defects the pre-2.6.0 release-readiness sweep found
//! on `develop` (#1140 follow-up). Three of them were introduced by the audit
//! fixes themselves, which is why each is pinned end to end through the binary
//! rather than at the unit level: every one of them type-checked, passed the
//! existing suite, and still produced wrong or missing output.

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

fn rc(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect()
}

/// Index over three transcripts, plus inward FR read pairs drawn from them.
/// `mate1_len` shortens R1: below `k` it cannot seed, so R2 becomes the anchor
/// and R1 has to be recovered — the case B3 got backwards.
fn fixture(dir: &Path, mate1_len: usize) -> (PathBuf, PathBuf, PathBuf) {
    let txs = [sequence(5, 3000), sequence(19, 3000), sequence(101, 2000)];
    let fasta = dir.join("txome.fa");
    std::fs::write(
        &fasta,
        format!(">txA\n{}\n>txB\n{}\n>txC\n{}\n", txs[0], txs[1], txs[2]),
    )
    .unwrap();
    let index = dir.join("idx");
    assert!(Command::new(SALMON)
        .args(["index", "-t"])
        .arg(&fasta)
        .arg("-i")
        .arg(&index)
        .args(["-p", "2"])
        .status()
        .unwrap()
        .success());

    let (mut r1, mut r2) = (String::new(), String::new());
    for i in 0..300 {
        let s = &txs[i % 3];
        let st = (i * 7) % (s.len() - 350);
        // ISF by construction: mate 1 forward, mate 2 reverse-complement, inward.
        r1.push_str(&format!(
            "@x{i}\n{}\n+\n{}\n",
            &s[st..st + mate1_len],
            "I".repeat(mate1_len)
        ));
        r2.push_str(&format!(
            "@x{i}\n{}\n+\n{}\n",
            rc(&s[st + 200..st + 300]),
            "I".repeat(100)
        ));
    }
    let (p1, p2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    std::fs::write(&p1, r1).unwrap();
    std::fs::write(&p2, r2).unwrap();
    (index, p1, p2)
}

fn meta_number(out: &Path, key: &str) -> f64 {
    let m = std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap();
    m.lines()
        .find(|l| l.contains(&format!("\"{key}\"")))
        .and_then(|l| l.rsplit(':').next())
        .map(|v| v.trim().trim_end_matches(',').parse().unwrap())
        .unwrap_or_else(|| panic!("{key} missing from meta_info.json:\n{m}"))
}

/// B1: `-g` produced no `quant.genes.sf` at all in the DEFAULT reads mode — the
/// single most common salmon invocation — while `-a`, `--rad` and `--online`
/// all wrote one. The gene-level guard tested `map_opts.skip_quant`, which
/// phase 1 forces true regardless of what the user asked, so the write was dead
/// code and the run explained itself with a `--skipQuant` message on a command
/// line that never contained `--skipQuant`.
#[test]
fn gene_map_is_honoured_in_the_default_reads_mode() {
    let dir = tempfile::tempdir().unwrap();
    let (index, p1, p2) = fixture(dir.path(), 100);
    let gmap = dir.path().join("t2g.tsv");
    std::fs::write(&gmap, "txA\tgeneA\ntxB\tgeneB\ntxC\tgeneC\n").unwrap();

    let out = dir.path().join("q");
    let run = Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "IU", "-1"])
        .arg(&p1)
        .arg("-2")
        .arg(&p2)
        .args(["-p", "2", "-g"])
        .arg(&gmap)
        .arg("-o")
        .arg(&out)
        .output()
        .unwrap();
    assert!(run.status.success());

    let genes = out.join("quant.genes.sf");
    assert!(
        genes.exists(),
        "-g must write quant.genes.sf on the default path; stderr:\n{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let body = std::fs::read_to_string(&genes).unwrap();
    assert_eq!(body.lines().count(), 4, "header + 3 genes:\n{body}");

    // The run must not claim --skipQuant on a command line that lacks it.
    let err = String::from_utf8_lossy(&run.stderr);
    assert!(
        !err.contains("--skipQuant"),
        "no --skipQuant explanation belongs here:\n{err}"
    );
}

/// B3: `--recoverOrphans` took the fragment's orientation from the *anchor*
/// rather than from mate 1, so whenever mate 2 was the anchor the recorded
/// strand was inverted. A stranded library then judged every recovered pair
/// backwards — mapping 0% under the correct `-l` and 100% under the wrong one.
/// Latent until the audit made the recorded format actually drive the filter.
#[test]
fn recovered_pairs_keep_their_strand_when_mate_two_anchors() {
    let dir = tempfile::tempdir().unwrap();
    // Mate 1 is 25 bp: below k=31, so it cannot seed and mate 2 must anchor.
    let (index, p1, p2) = fixture(dir.path(), 25);

    let run = |lib: &str, out: &Path| {
        assert!(Command::new(SALMON)
            .args(["quant", "-i"])
            .arg(&index)
            .args(["-l", lib, "-1"])
            .arg(&p1)
            .arg("-2")
            .arg(&p2)
            .args(["--recoverOrphans", "-p", "2", "-o"])
            .arg(out)
            .status()
            .unwrap()
            .success());
    };

    let isf = dir.path().join("isf");
    let isr = dir.path().join("isr");
    run("ISF", &isf);
    run("ISR", &isr);

    // The fixture is ISF by construction, so ISF must keep the fragments and
    // ISR must reject them. Inverted, this assertion reads exactly backwards.
    assert!(
        meta_number(&isf, "percent_mapped") > 90.0,
        "ISF is the true library type and must map"
    );
    assert!(
        meta_number(&isr, "percent_mapped") < 10.0,
        "ISR is the wrong strand and must be rejected"
    );

    let counts = std::fs::read_to_string(isf.join("lib_format_counts.json")).unwrap();
    let isf_obs = counts
        .lines()
        .find(|l| l.contains("\"ISF\""))
        .expect("ISF key in lib_format_counts.json");
    assert!(
        !isf_obs.contains(": 0"),
        "the observed format histogram must record ISF, not its opposite:\n{counts}"
    );
}

/// B4: an unsatisfiable `-p` aborted the process with a raw rayon panic (exit
/// 134) whose message named neither `-p` nor the cause. An unset
/// `$SLURM_CPUS_PER_TASK` or a typo'd digit is an ordinary user error.
#[test]
fn an_unsatisfiable_thread_count_is_clamped_not_a_panic() {
    let dir = tempfile::tempdir().unwrap();
    let (index, p1, p2) = fixture(dir.path(), 100);
    let out = dir.path().join("q");
    let run = Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "IU", "-1"])
        .arg(&p1)
        .arg("-2")
        .arg(&p2)
        .args(["-p", "999999", "-o"])
        .arg(&out)
        .output()
        .unwrap();

    assert!(
        run.status.success(),
        "an over-large -p must not abort the run; exit {:?}\n{}",
        run.status.code(),
        String::from_utf8_lossy(&run.stderr)
    );
    let err = String::from_utf8_lossy(&run.stderr);
    assert!(
        err.contains("-p 999999") && err.contains("hardware thread"),
        "the warning must name the flag and the real limit:\n{err}"
    );
    assert!(
        !err.contains("panicked"),
        "no panic may reach the user:\n{err}"
    );
    assert!(
        out.join("quant.sf").exists(),
        "the run must still produce output"
    );
}
