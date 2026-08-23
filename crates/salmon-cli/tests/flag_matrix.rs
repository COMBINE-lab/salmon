//! Mode×flag applicability, end to end (#1140, audit batch E): a flag that
//! does nothing in the selected mode must say so, and the flags the audit
//! found silently dropped must now act. One index build feeds every check.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

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

/// Two transcripts and inward paired reads drawn from them.
fn fixture(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let a = sequence(5, 3000);
    let b = sequence(19, 3000);
    let fasta = dir.join("txome.fa");
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

fn salmon(args: &[&std::ffi::OsStr]) -> Output {
    Command::new(SALMON).args(args).output().unwrap()
}

fn os(s: &str) -> &std::ffi::OsStr {
    std::ffi::OsStr::new(s)
}

#[test]
fn mode_flag_matrix_end_to_end() {
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
    let quant_reads = |extra: &[&std::ffi::OsStr], out: &Path| -> Output {
        let mut v: Vec<&std::ffi::OsStr> = vec![
            os("quant"),
            os("-i"),
            index.as_os_str(),
            os("-l"),
            os("IU"),
            os("-1"),
            r1.as_os_str(),
            os("-2"),
            r2.as_os_str(),
            os("-p"),
            os("2"),
            os("-o"),
            out.as_os_str(),
        ];
        v.extend_from_slice(extra);
        salmon(&v)
    };

    // --dumpBiasModels must produce bias_models.txt on the DEFAULT
    // (deterministic) path — it silently vanished with the 2.6.0 flip because
    // only the one-pass writer had a dump site. Also produce the RAD the later
    // --rad checks quantify.
    let out = dir.path().join("det");
    let rad = dir.path().join("mappings.rad");
    let o = quant_reads(
        &[
            os("--gcBias"),
            os("--dumpBiasModels"),
            os("--writeRad"),
            rad.as_os_str(),
        ],
        &out,
    );
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(
        out.join("bias_models.txt").exists(),
        "--dumpBiasModels must write bias_models.txt under the deterministic default"
    );
    assert!(rad.exists());

    // Sketch mode computes no alignment scores: the selective-alignment knobs
    // must warn by name instead of being silently inert.
    let o = quant_reads(
        &[
            os("--sketch"),
            os("--ma"),
            os("3"),
            os("--minScoreFraction"),
            os("0.8"),
            os("--fullLengthAlignment"),
        ],
        &dir.path().join("sk_warn"),
    );
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let err = String::from_utf8_lossy(&o.stderr).into_owned();
    for needle in ["--ma", "--minScoreFraction", "--fullLengthAlignment"] {
        assert!(
            err.contains(needle),
            "sketch run must name inert flag {needle}: {err}"
        );
    }
    assert!(err.contains("no effect under --sketch"), "{err}");

    // Control: a plain sketch run draws no such complaint, and a plain
    // selective-alignment run does not warn about sketch knobs it honours.
    let o = quant_reads(&[os("--sketch")], &dir.path().join("sk_plain"));
    assert!(
        !String::from_utf8_lossy(&o.stderr).contains("no effect under --sketch"),
        "a plain --sketch run must not warn"
    );
    let o = quant_reads(&[os("--ma"), os("3")], &dir.path().join("sa_ma"));
    assert!(
        !String::from_utf8_lossy(&o.stderr).contains("no effect under --sketch"),
        "--ma is honoured in selective-alignment mode and must not warn"
    );

    // --sketchStrictOrphans without --sketch is the reverse orphan: warn.
    let o = quant_reads(&[os("--sketchStrictOrphans")], &dir.path().join("sso"));
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("--sketchStrictOrphans has no effect"),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );

    // RAD-input mode: --online has no online path to select, and it must not
    // suppress the inert-knob warning for the online knobs it used to mask;
    // the RAD-shaping flags have nothing to act on either.
    let o = salmon(&[
        os("quant"),
        os("--rad"),
        rad.as_os_str(),
        os("-l"),
        os("IU"),
        os("--online"),
        os("--forgettingFactor"),
        os("0.9"),
        os("--keepRad"),
        os("-o"),
        dir.path().join("rad_online").as_os_str(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let err = String::from_utf8_lossy(&o.stderr).into_owned();
    assert!(
        err.contains("--online has no effect in RAD-input mode"),
        "{err}"
    );
    assert!(
        err.contains("--forgettingFactor"),
        "--online must not suppress the inert-knob warning in --rad mode: {err}"
    );
    assert!(
        !err.contains("deprecated one-pass quantification path"),
        "--rad --online must not claim the one-pass path will run: {err}"
    );
    assert!(err.contains("--keepRad"), "{err}");

    // --genome / --juncMissDiscount are genome-projection knobs; without
    // --annotation clap rejects them up front instead of silently dropping.
    let o = salmon(&[
        os("quant"),
        os("-i"),
        index.as_os_str(),
        os("-l"),
        os("IU"),
        os("-1"),
        r1.as_os_str(),
        os("-2"),
        r2.as_os_str(),
        os("--genome"),
        fasta.as_os_str(),
        os("-o"),
        dir.path().join("g").as_os_str(),
    ]);
    assert!(
        !o.status.success(),
        "--genome without --annotation must be rejected"
    );

    // -g under --skipQuant: no abundances exist, so no quant.genes.sf full of
    // zeros — the run says why instead.
    let gm = dir.path().join("t2g.tsv");
    std::fs::write(&gm, "txA\tgeneA\ntxB\tgeneB\n").unwrap();
    let out_sq = dir.path().join("rad_sq");
    let o = salmon(&[
        os("quant"),
        os("--rad"),
        rad.as_os_str(),
        os("-l"),
        os("IU"),
        os("--skipQuant"),
        os("-g"),
        gm.as_os_str(),
        os("-o"),
        out_sq.as_os_str(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(
        !out_sq.join("quant.genes.sf").exists(),
        "--skipQuant must not write a zero-filled quant.genes.sf"
    );
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("gene-level aggregation (-g) skipped"),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
}
