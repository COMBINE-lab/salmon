//! Group 2 of the pre-2.6.0 readiness sweep: the mode x flag applicability
//! check existed for sketch-vs-selective-alignment and for the online-only
//! knobs, but had never been written for reads-vs-`-a`/`--rad`. ~30 mapping
//! flags were accepted there in total silence.
//!
//! Two of them turned out to be missing *features* rather than missing
//! warnings, and those are pinned by result rather than by message.

use std::path::{Path, PathBuf};
use std::process::Command;

const SALMON: &str = env!("CARGO_BIN_EXE_salmon");

fn os(s: &str) -> &std::ffi::OsStr {
    std::ffi::OsStr::new(s)
}

/// A SAM whose fragments are contested between two transcripts of equal length,
/// txA scoring one point better, plus unique anchors for both so neither is
/// driven to zero. This is the only shape where soft-weighting and hard
/// filtering give visibly different answers.
fn contested_sam(dir: &Path) -> PathBuf {
    let path = dir.join("contested.sam");
    let mut s =
        String::from("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:txA\tLN:3000\n@SQ\tSN:txB\tLN:3000\n");
    let mut pos = 1usize;
    for i in 0..200 {
        let (l, r) = (pos % 900 + 1, pos % 900 + 201);
        s.push_str(&format!(
            "m{i}\t99\ttxA\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:200\n"
        ));
        s.push_str(&format!(
            "m{i}\t147\ttxA\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:200\n"
        ));
        s.push_str(&format!(
            "m{i}\t355\ttxB\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:199\n"
        ));
        s.push_str(&format!(
            "m{i}\t403\ttxB\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:199\n"
        ));
        pos += 7;
    }
    for tx in ["txA", "txB"] {
        for i in 0..60 {
            let (l, r) = (pos % 900 + 1, pos % 900 + 201);
            s.push_str(&format!(
                "u{tx}{i}\t99\t{tx}\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:200\n"
            ));
            s.push_str(&format!(
                "u{tx}{i}\t147\t{tx}\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:200\n"
            ));
            pos += 11;
        }
    }
    std::fs::write(&path, s).unwrap();
    path
}

/// A SAM with proper pairs plus half-mapped fragments, for `--discardOrphans`.
fn orphan_sam(dir: &Path) -> PathBuf {
    let path = dir.join("orphans.sam");
    let mut s =
        String::from("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:txA\tLN:5000\n@SQ\tSN:txB\tLN:5000\n");
    let mut pos = 1usize;
    for i in 0..40 {
        let (l, r) = (pos + 1, pos + 201);
        s.push_str(&format!(
            "p{i}\t99\ttxA\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:200\n"
        ));
        s.push_str(&format!(
            "p{i}\t147\ttxA\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:200\n"
        ));
        pos += 400;
    }
    for i in 0..25 {
        s.push_str(&format!(
            "o{i}\t73\ttxB\t{}\t60\t100M\t*\t0\t0\t*\t*\tAS:i:150\n",
            i * 61 + 1
        ));
    }
    std::fs::write(&path, s).unwrap();
    path
}

fn quant(args: &[&std::ffi::OsStr]) -> std::process::Output {
    let mut v = vec![os("quant")];
    v.extend_from_slice(args);
    Command::new(SALMON).args(v).output().unwrap()
}

fn counts(out: &Path) -> Vec<(String, f64)> {
    std::fs::read_to_string(out.join("quant.sf"))
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[0].to_string(), c[4].parse().unwrap())
        })
        .collect()
}

/// Reads-mode mapping knobs cannot apply to input that arrives already aligned
/// or already mapped, and must say so rather than being accepted in silence.
#[test]
fn reads_mapping_flags_warn_in_alignment_and_rad_modes() {
    let dir = tempfile::tempdir().unwrap();
    let sam = contested_sam(dir.path());
    let out = dir.path().join("q");
    let rad = dir.path().join("m.rad");

    let o = quant(&[
        os("-a"),
        sam.as_os_str(),
        os("-l"),
        os("IU"),
        os("--sketch"),
        os("--ma"),
        os("3"),
        os("--recoverOrphans"),
        os("--writeRad"),
        rad.as_os_str(),
        os("-o"),
        out.as_os_str(),
    ]);
    assert!(o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("no effect in alignment mode"), "{err}");
    for f in ["--sketch", "--ma", "--recoverOrphans"] {
        assert!(err.contains(f), "the warning must name {f}:\n{err}");
    }

    let o = quant(&[
        os("--rad"),
        rad.as_os_str(),
        os("-l"),
        os("IU"),
        os("--minScoreFraction"),
        os("0.8"),
        os("-o"),
        out.as_os_str(),
    ]);
    assert!(o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("no effect in RAD-input mode") && err.contains("--minScoreFraction"),
        "{err}"
    );
    // Singular reads as singular: "has", not "have".
    assert!(
        err.contains("--minScoreFraction has no effect"),
        "one flag must take a singular verb:\n{err}"
    );

    // Control: a run that passes none of them says nothing about applicability.
    let o = quant(&[
        os("--rad"),
        rad.as_os_str(),
        os("-l"),
        os("IU"),
        os("-o"),
        out.as_os_str(),
    ]);
    assert!(!String::from_utf8_lossy(&o.stderr).contains("no effect in"));
}

/// `--hardFilter` is the other half of the `--scoreExp` policy, and `--scoreExp`
/// was already honoured on these paths — so accepting `--hardFilter` and
/// ignoring it left the policy half-applied. Pinned by result: soft weighting
/// shares contested mass, hard filtering gives it to the best-scoring target.
#[test]
fn hard_filter_takes_effect_on_alignment_input() {
    let dir = tempfile::tempdir().unwrap();
    let sam = contested_sam(dir.path());
    let (soft, hard) = (dir.path().join("soft"), dir.path().join("hard"));

    assert!(quant(&[
        os("-a"),
        sam.as_os_str(),
        os("-l"),
        os("IU"),
        os("-o"),
        soft.as_os_str()
    ])
    .status
    .success());
    assert!(quant(&[
        os("-a"),
        sam.as_os_str(),
        os("-l"),
        os("IU"),
        os("--hardFilter"),
        os("-o"),
        hard.as_os_str()
    ])
    .status
    .success());

    let get = |v: &[(String, f64)], n: &str| v.iter().find(|(x, _)| x == n).unwrap().1;
    let (s, h) = (counts(&soft), counts(&hard));
    assert!(
        get(&h, "txA") > get(&s, "txA"),
        "hard filtering must give the contested mass to the better-scoring transcript: \
         soft {:?} vs hard {:?}",
        s,
        h
    );
    // Every contested fragment goes to txA: 200 contested + 60 unique.
    assert!(
        (get(&h, "txA") - 260.0).abs() < 1e-6 && (get(&h, "txB") - 60.0).abs() < 1e-6,
        "hard: {h:?}"
    );
}

/// A RAD records orphan placements, so the flag that governs them means
/// something in `--rad` — without it, requantifying a RAD contradicted the
/// `-a` run that wrote it.
#[test]
fn discard_orphans_takes_effect_on_rad_input() {
    let dir = tempfile::tempdir().unwrap();
    let sam = orphan_sam(dir.path());
    let rad = dir.path().join("o.rad");
    assert!(quant(&[
        os("-a"),
        sam.as_os_str(),
        os("-l"),
        os("IU"),
        os("--writeRad"),
        rad.as_os_str(),
        os("-o"),
        dir.path().join("seed").as_os_str(),
    ])
    .status
    .success());

    let (keep, drop) = (dir.path().join("keep"), dir.path().join("drop"));
    assert!(quant(&[
        os("--rad"),
        rad.as_os_str(),
        os("-l"),
        os("IU"),
        os("-o"),
        keep.as_os_str()
    ])
    .status
    .success());
    assert!(quant(&[
        os("--rad"),
        rad.as_os_str(),
        os("-l"),
        os("IU"),
        os("--discardOrphans"),
        os("-o"),
        drop.as_os_str()
    ])
    .status
    .success());

    let get = |v: &[(String, f64)], n: &str| v.iter().find(|(x, _)| x == n).unwrap().1;
    let (k, d) = (counts(&keep), counts(&drop));
    assert!(
        get(&k, "txB") > get(&d, "txB") + 1.0,
        "the orphan-only transcript must lose its mass when orphans are discarded: \
         keep {k:?} vs discard {d:?}"
    );
}

/// Alignment-input knobs in reads mode, and pair-only knobs under `-r`, are
/// equally inert and equally silent before this.
#[test]
fn wrong_mode_flags_warn_in_reads_and_single_end() {
    let dir = tempfile::tempdir().unwrap();
    // Minimal reads fixture.
    let mut seed = 7u64;
    let mut base = || {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        "ACGT".as_bytes()[(seed >> 33) as usize & 3] as char
    };
    let tx: String = (0..3000).map(|_| base()).collect();
    let fasta = dir.path().join("t.fa");
    std::fs::write(&fasta, format!(">txA\n{tx}\n")).unwrap();
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
    let se = dir.path().join("se.fq");
    let mut s = String::new();
    for i in 0..100 {
        let st = (i * 13) % 2800;
        s.push_str(&format!(
            "@r{i}\n{}\n+\n{}\n",
            &tx[st..st + 100],
            "I".repeat(100)
        ));
    }
    std::fs::write(&se, s).unwrap();
    let out = dir.path().join("q");

    let o = quant(&[
        os("-i"),
        index.as_os_str(),
        os("-l"),
        os("U"),
        os("-r"),
        se.as_os_str(),
        os("--errorModel"),
        os("--numErrorBins"),
        os("6"),
        os("--recoverOrphans"),
        os("--allowDovetail"),
        os("-p"),
        os("2"),
        os("-o"),
        out.as_os_str(),
    ]);
    assert!(o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("no effect in reads mode") && err.contains("--errorModel"),
        "alignment-only flags must warn in reads mode:\n{err}"
    );
    assert!(
        err.contains("no effect in single-end mode") && err.contains("--recoverOrphans"),
        "pair-only flags must warn under -r:\n{err}"
    );
}
