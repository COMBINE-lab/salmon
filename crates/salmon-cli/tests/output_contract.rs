//! End-to-end output-contract checks for the deterministic flow (#1140):
//! decoy exclusion survives the full index → RAD → requant round trip,
//! phase 2 preserves phase 1's invocation record, and single-end provenance
//! tells the truth about where the FLD came from.

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

/// A gentrome-style fixture: two transcripts plus one decoy sequence (listed
/// in decoys.txt, appended last as `salmon index -d` expects), and inward
/// paired reads drawn from the transcripts only.
fn fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf, PathBuf) {
    let a = sequence(5, 3000);
    let b = sequence(19, 3000);
    let decoy = sequence(101, 4000);
    let fasta = dir.join("gentrome.fa");
    std::fs::write(&fasta, format!(">txA\n{a}\n>txB\n{b}\n>decoy1\n{decoy}\n")).unwrap();
    let decoys = dir.join("decoys.txt");
    std::fs::write(&decoys, "decoy1\n").unwrap();

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
    (fasta, decoys, p1, p2)
}

fn quant_names(out: &Path) -> Vec<String> {
    std::fs::read_to_string(out.join("quant.sf"))
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect()
}

#[test]
fn deterministic_flow_honors_the_output_contract_with_decoys() {
    let dir = tempfile::tempdir().unwrap();
    let (fasta, decoys, r1, r2) = fixture(dir.path());
    let index = dir.path().join("idx");
    assert!(Command::new(SALMON)
        .args(["index", "-t"])
        .arg(&fasta)
        .args(["-d"])
        .arg(&decoys)
        .arg("-i")
        .arg(&index)
        .args(["-p", "2"])
        .status()
        .unwrap()
        .success());

    let run = |tag: &str, extra: &[&str]| -> PathBuf {
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
        assert!(cmd.status().unwrap().success(), "{tag} run failed");
        out
    };

    // The decoy round trip: phase 1 bakes the boundary into the RAD, phase 2
    // excludes the block — so the deterministic flow's quant.sf carries the
    // same row set as the one-pass path's, decoy-free.
    let det = run("det", &["--deterministic"]);
    let online = run("online", &[]);
    let det_names = quant_names(&det);
    assert!(
        !det_names.iter().any(|n| n == "decoy1"),
        "decoy row leaked into deterministic quant.sf: {det_names:?}"
    );
    assert_eq!(
        det_names,
        quant_names(&online),
        "deterministic and one-pass runs must emit the same quant.sf row set"
    );

    // Phase 2 keeps phase 1's invocation record: the index and the read files
    // stay in cmd_info.json instead of being clobbered by a RAD-centric one.
    let cmd_info = std::fs::read_to_string(det.join("cmd_info.json")).unwrap();
    assert!(
        cmd_info.contains("\"index\"") && cmd_info.contains("r1.fq"),
        "cmd_info.json must record the real invocation: {cmd_info}"
    );
    // ...and lib_format_counts.json names the reads rather than [].
    let lfc = std::fs::read_to_string(det.join("lib_format_counts.json")).unwrap();
    assert!(
        lfc.contains("r1.fq"),
        "lib_format_counts.json read_files must name the inputs: {lfc}"
    );
}

#[test]
fn single_end_frag_length_source_is_the_prior() {
    let dir = tempfile::tempdir().unwrap();
    let (fasta, _decoys, r1, _r2) = fixture(dir.path());
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

    // Single-end (one-pass, where the field was wrong): no fragment lengths
    // exist to observe, so the FLD is the --fldMean/--fldSD prior and
    // meta_info.json must say so.
    let out = dir.path().join("se");
    assert!(Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "U", "-r"])
        .arg(&r1)
        .args(["-p", "2", "-o"])
        .arg(&out)
        .status()
        .unwrap()
        .success());
    let meta = std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap();
    assert!(
        meta.contains("\"frag_length_source\": \"prior\""),
        "a single-end run's FLD is the prior, and provenance must say so: {meta}"
    );
}

fn revcomp(s: &str) -> String {
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

/// An orphan's SAM record must carry its own mate's strand. Every orphan keeps
/// its strand in `is_fw`, but a right orphan's record read `r2_fw`, which only a
/// proper pair sets, so each was written as reverse (FLAG 0x10) with its
/// sequence reverse-complemented, whatever its orientation (#1272).
#[test]
fn an_orphan_record_carries_its_own_strand() {
    let dir = tempfile::tempdir().unwrap();
    let (fasta, _, _, _) = fixture(dir.path());
    let index = dir.path().join("idx");
    assert!(Command::new(SALMON)
        .args(["index", "-t"])
        .arg(&fasta)
        .arg("-i")
        .arg(&index)
        .args(["-p", "2"])
        .output()
        .unwrap()
        .status
        .success());

    // One mate from txA, on the strand named, and one that maps nowhere.
    let a = sequence(5, 3000);
    let junk = sequence(777, 100);
    let fragments = [
        ("r1_fw", a[500..600].to_string(), junk.clone()),
        ("r1_rc", revcomp(&a[900..1000]), junk.clone()),
        ("r2_fw", junk.clone(), a[1300..1400].to_string()),
        ("r2_rc", junk.clone(), revcomp(&a[1700..1800])),
    ];
    let qual = "I".repeat(100);
    let (mut r1, mut r2) = (String::new(), String::new());
    for (name, m1, m2) in &fragments {
        r1.push_str(&format!("@{name}\n{m1}\n+\n{qual}\n"));
        r2.push_str(&format!("@{name}\n{m2}\n+\n{qual}\n"));
    }
    let (p1, p2) = (dir.path().join("o1.fq"), dir.path().join("o2.fq"));
    std::fs::write(&p1, r1).unwrap();
    std::fs::write(&p2, r2).unwrap();

    for mode in [None, Some("--sketch")] {
        let out = dir.path().join(format!("q{}", mode.is_some()));
        let sam = dir.path().join(format!("o{}.sam", mode.is_some()));
        let o = Command::new(SALMON)
            .args(["quant", "-l", "IU", "-p", "2", "-i"])
            .arg(&index)
            .arg("-1")
            .arg(&p1)
            .arg("-2")
            .arg(&p2)
            .arg("-o")
            .arg(&out)
            .arg(format!("--writeMappings={}", sam.display()))
            .args(mode)
            .output()
            .unwrap();
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        let sam = std::fs::read_to_string(&sam).unwrap();
        for (name, m1, m2) in &fragments {
            let rec: Vec<&str> = sam
                .lines()
                .find(|l| l.split('\t').next() == Some(*name))
                .unwrap_or_else(|| panic!("{mode:?}: no record for {name}:\n{sam}"))
                .split('\t')
                .collect();
            let flag: u16 = rec[1].parse().unwrap();
            let (read2, given) = if name.starts_with("r2") {
                (true, m2)
            } else {
                (false, m1)
            };
            let reverse = name.ends_with("rc");
            assert_eq!(flag & 0x80 != 0, read2, "{mode:?} {name}: FLAG {flag}");
            assert_eq!(flag & 0x10 != 0, reverse, "{mode:?} {name}: FLAG {flag}");
            // SEQ is on the reference's forward strand.
            let seq = if reverse {
                revcomp(given)
            } else {
                given.to_string()
            };
            assert_eq!(rec[9], seq, "{mode:?} {name}: SEQ");
        }
    }
}
