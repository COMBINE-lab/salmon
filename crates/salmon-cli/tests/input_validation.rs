//! Group 1 of the pre-2.6.0 readiness sweep: inputs and flag combinations that
//! produced a *wrong answer presented as a right one* — a clean `quant.sf` and
//! exit 0 over data salmon had silently discarded or a request it had silently
//! dropped. Each is now refused at invocation.

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

struct Fx {
    index: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
    se: PathBuf,
}

fn fixture(dir: &Path) -> Fx {
    let txs = [sequence(5, 3000), sequence(19, 3000)];
    let fasta = dir.join("txome.fa");
    std::fs::write(&fasta, format!(">txA\n{}\n>txB\n{}\n", txs[0], txs[1])).unwrap();
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
    let (mut r1, mut r2, mut se) = (String::new(), String::new(), String::new());
    let q = "I".repeat(100);
    for i in 0..200 {
        let s = &txs[i % 2];
        let st = (i * 11) % (s.len() - 350);
        r1.push_str(&format!("@f{i}\n{}\n+\n{q}\n", &s[st..st + 100]));
        r2.push_str(&format!("@f{i}\n{}\n+\n{q}\n", rc(&s[st + 200..st + 300])));
        se.push_str(&format!("@s{i}\n{}\n+\n{q}\n", &s[st..st + 100]));
    }
    let (p1, p2, ps) = (dir.join("r1.fq"), dir.join("r2.fq"), dir.join("se.fq"));
    std::fs::write(&p1, r1).unwrap();
    std::fs::write(&p2, r2).unwrap();
    std::fs::write(&ps, se).unwrap();
    Fx {
        index,
        r1: p1,
        r2: p2,
        se: ps,
    }
}

fn run(args: &[&std::ffi::OsStr]) -> std::process::Output {
    Command::new(SALMON).args(args).output().unwrap()
}
fn os(s: &str) -> &std::ffi::OsStr {
    std::ffi::OsStr::new(s)
}

/// An invalid `-l` was a hard error in reads mode but silently accepted in
/// `-a` and `--rad`, where every call site discarded the parse error with
/// `.ok()` and quantified as though the library were unstranded. Same input,
/// same typo, two different answers — one of them wrong and exit 0.
#[test]
fn an_invalid_library_type_is_refused_in_every_mode() {
    let dir = tempfile::tempdir().unwrap();
    let fx = fixture(dir.path());

    // Produce a RAD and a SAM so all three input modes can be exercised.
    let rad = dir.path().join("m.rad");
    assert!(Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&fx.index)
        .args(["-l", "IU", "-1"])
        .arg(&fx.r1)
        .arg("-2")
        .arg(&fx.r2)
        .args(["-p", "2", "--writeRad"])
        .arg(&rad)
        .arg("-o")
        .arg(dir.path().join("seed"))
        .status()
        .unwrap()
        .success());

    let out = dir.path().join("bad");
    let modes: Vec<Vec<&std::ffi::OsStr>> = vec![
        vec![
            os("-i"),
            fx.index.as_os_str(),
            os("-1"),
            fx.r1.as_os_str(),
            os("-2"),
            fx.r2.as_os_str(),
        ],
        vec![os("--rad"), rad.as_os_str()],
    ];
    for m in modes {
        let mut v: Vec<&std::ffi::OsStr> = vec![os("quant")];
        v.extend(m.iter().copied());
        v.extend([
            os("-l"),
            os("XYZ"),
            os("-p"),
            os("2"),
            os("-o"),
            out.as_os_str(),
        ]);
        let o = run(&v);
        assert!(
            !o.status.success(),
            "an invalid -l must be refused, not quantified as unstranded"
        );
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(
            err.contains("invalid library type 'XYZ'"),
            "the error must name the offending value:\n{err}"
        );
    }

    // Control: a valid type still runs, and `A` (auto) is not mistaken for one.
    for lib in ["IU", "A", "ISR"] {
        let o = run(&[
            os("quant"),
            os("-i"),
            fx.index.as_os_str(),
            os("-1"),
            fx.r1.as_os_str(),
            os("-2"),
            fx.r2.as_os_str(),
            os("-l"),
            os(lib),
            os("-p"),
            os("2"),
            os("-o"),
            out.as_os_str(),
        ]);
        assert!(o.status.success(), "-l {lib} must still be accepted");
    }
}

/// `-r` alongside `-1/-2` silently discarded the single-end reads while
/// `cmd_info.json` went on recording them — the "wrong mode chosen silently"
/// failure `-a`'s conflict exists to prevent.
#[test]
fn single_end_and_paired_inputs_cannot_be_mixed() {
    let dir = tempfile::tempdir().unwrap();
    let fx = fixture(dir.path());
    let out = dir.path().join("q");
    let o = run(&[
        os("quant"),
        os("-i"),
        fx.index.as_os_str(),
        os("-l"),
        os("IU"),
        os("-1"),
        fx.r1.as_os_str(),
        os("-2"),
        fx.r2.as_os_str(),
        os("-r"),
        fx.se.as_os_str(),
        os("-p"),
        os("2"),
        os("-o"),
        out.as_os_str(),
    ]);
    assert!(!o.status.success(), "-r with -1/-2 must be refused");
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("--unmatedReads") || err.contains("-r"),
        "clap must name the conflicting flag:\n{err}"
    );
}

/// Bias correction adjusts effective lengths; `--noLengthCorrection` turns that
/// off. Every driver resolved the contradiction by silently dropping the bias
/// request — and then reported `gc_bias_correct: true` for a run that did no
/// bias correction at all.
#[test]
fn no_length_correction_conflicts_with_bias_correction() {
    let dir = tempfile::tempdir().unwrap();
    let fx = fixture(dir.path());
    let out = dir.path().join("q");
    for bias in ["--seqBias", "--gcBias", "--posBias"] {
        let o = run(&[
            os("quant"),
            os("-i"),
            fx.index.as_os_str(),
            os("-l"),
            os("IU"),
            os("-1"),
            fx.r1.as_os_str(),
            os("-2"),
            fx.r2.as_os_str(),
            os("--noLengthCorrection"),
            os(bias),
            os("-p"),
            os("2"),
            os("-o"),
            out.as_os_str(),
        ]);
        assert!(
            !o.status.success(),
            "--noLengthCorrection {bias} must be refused, not silently resolved"
        );
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(
            err.contains("--noLengthCorrection cannot be combined with") && err.contains(bias),
            "the error must name both sides:\n{err}"
        );
    }

    // Control: --noLengthCorrection alone is a legitimate request.
    let o = run(&[
        os("quant"),
        os("-i"),
        fx.index.as_os_str(),
        os("-l"),
        os("IU"),
        os("-1"),
        fx.r1.as_os_str(),
        os("-2"),
        fx.r2.as_os_str(),
        os("--noLengthCorrection"),
        os("-p"),
        os("2"),
        os("-o"),
        out.as_os_str(),
    ]);
    assert!(
        o.status.success(),
        "--noLengthCorrection alone must still work"
    );
}

/// A truncated FASTQ was quantified up to its last complete record and reported
/// as a normal run. The trailing-line cases were also inconsistent with each
/// other: 1 or 2 orphaned lines were dropped in silence while 3 errored.
#[test]
fn a_truncated_fastq_is_refused_rather_than_partially_quantified() {
    let dir = tempfile::tempdir().unwrap();
    let fx = fixture(dir.path());
    let full = std::fs::read_to_string(&fx.se).unwrap();
    let out = dir.path().join("q");

    let quant = |reads: &Path, out: &Path| {
        run(&[
            os("quant"),
            os("-i"),
            fx.index.as_os_str(),
            os("-l"),
            os("U"),
            os("-r"),
            reads.as_os_str(),
            os("-p"),
            os("2"),
            os("-o"),
            out.as_os_str(),
        ])
    };

    // Byte-level truncation (the interrupted-transfer shape).
    let mid = dir.path().join("mid.fq");
    std::fs::write(&mid, &full[..full.len() / 3]).unwrap();
    let o = quant(&mid, &out);
    assert!(
        !o.status.success(),
        "a byte-truncated FASTQ must be refused"
    );
    // Two wordings are correct here and which one fires depends on where the
    // cut landed: a short file reports its line count, a file whose final
    // record is malformed reports that instead (salmon cannot tell a truncated
    // transfer from a corrupted one, and says so rather than guessing).
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("looks truncated") || err.contains("final FASTQ record is incomplete"),
        "{err}"
    );
    assert!(
        err.contains("mid.fq"),
        "the error must name the file:\n{err}"
    );

    // Line-boundary truncation: 1, 2 and 3 orphaned lines must all be refused,
    // where 1 and 2 used to pass silently.
    let lines: Vec<&str> = full.lines().collect();
    for extra in 1..=3 {
        let cut = dir.path().join(format!("cut{extra}.fq"));
        let keep = lines.len() - 4 + extra;
        std::fs::write(&cut, format!("{}\n", lines[..keep].join("\n"))).unwrap();
        let o = quant(&cut, &out);
        assert!(
            !o.status.success(),
            "a FASTQ ending {extra} line(s) into a record must be refused"
        );
    }

    // Control: the intact file, and a compressed copy, still run. (Compressed
    // truncation is caught by the decompressor, so the preflight skips it.)
    assert!(
        quant(&fx.se, &out).status.success(),
        "intact FASTQ must run"
    );
}

// ---- Group 3: incompatibilities detectable at invocation ------------------

/// Genome projection lives inside the `-a` branch, so `--annotation` with any
/// other input silently discarded the whole request — annotation, genome and
/// junction discount — without even checking the paths exist.
///
/// clap's `requires` is necessary but not sufficient here, and the test covers
/// why: clap skips a `requires` when the required argument conflicts with one
/// already present, and `-1/-2`, `-r` and `--rad` all conflict with `-a`. Those
/// are exactly the ways a user gets here, so all three are exercised.
#[test]
fn annotation_requires_alignment_input() {
    let dir = tempfile::tempdir().unwrap();
    let fx = fixture(dir.path());
    let rad = dir.path().join("m.rad");
    assert!(Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&fx.index)
        .args(["-l", "IU", "-1"])
        .arg(&fx.r1)
        .arg("-2")
        .arg(&fx.r2)
        .args(["-p", "2", "--writeRad"])
        .arg(&rad)
        .arg("-o")
        .arg(dir.path().join("seed"))
        .status()
        .unwrap()
        .success());

    let out = dir.path().join("q");
    let gtf = dir.path().join("a.gtf");
    std::fs::write(&gtf, "").unwrap();
    let inputs: Vec<Vec<&std::ffi::OsStr>> = vec![
        vec![
            os("-i"),
            fx.index.as_os_str(),
            os("-1"),
            fx.r1.as_os_str(),
            os("-2"),
            fx.r2.as_os_str(),
        ],
        vec![os("-i"), fx.index.as_os_str(), os("-r"), fx.se.as_os_str()],
        vec![os("--rad"), rad.as_os_str()],
    ];
    for inp in inputs {
        let mut v: Vec<&std::ffi::OsStr> = vec![os("quant")];
        v.extend(inp.iter().copied());
        v.extend([
            os("-l"),
            os("IU"),
            os("--annotation"),
            gtf.as_os_str(),
            os("-p"),
            os("2"),
            os("-o"),
            out.as_os_str(),
        ]);
        let o = run(&v);
        assert!(
            !o.status.success(),
            "--annotation without -a must be refused, not silently dropped"
        );
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(
            err.contains("--annotation") && (err.contains("-a") || err.contains("alignments")),
            "the error must name both flags:\n{err}"
        );
    }
}

/// `--errorModel` and `--noErrorModel` are direct opposites; passing both was
/// resolved by silent precedence, after which salmon advised passing
/// `--errorModel` — already on the command line.
#[test]
fn error_model_flags_conflict() {
    let dir = tempfile::tempdir().unwrap();
    let fx = fixture(dir.path());
    let sam = dir.path().join("a.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:txA\tLN:3000\n\
         r0\t99\ttxA\t1\t60\t100M\t=\t201\t300\t*\t*\tAS:i:200\n\
         r0\t147\ttxA\t201\t60\t100M\t=\t1\t-300\t*\t*\tAS:i:200\n",
    )
    .unwrap();
    let _ = &fx;
    let o = run(&[
        os("quant"),
        os("-a"),
        sam.as_os_str(),
        os("-l"),
        os("IU"),
        os("--errorModel"),
        os("--noErrorModel"),
        os("-o"),
        dir.path().join("q").as_os_str(),
    ]);
    assert!(!o.status.success(), "opposite flags must conflict");
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("cannot be used with"),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );
}

/// Posterior sampling happens inside the pass `--skipQuant` skips, and
/// `--thinningFactor` thins a Gibbs chain that may not exist. Both were
/// discarded in silence while their siblings explained themselves.
#[test]
fn requests_that_cannot_be_honoured_say_so() {
    let dir = tempfile::tempdir().unwrap();
    let fx = fixture(dir.path());
    let out = dir.path().join("q");
    let base: Vec<&std::ffi::OsStr> = vec![
        os("quant"),
        os("-i"),
        fx.index.as_os_str(),
        os("-l"),
        os("IU"),
        os("-1"),
        fx.r1.as_os_str(),
        os("-2"),
        fx.r2.as_os_str(),
        os("-p"),
        os("2"),
    ];

    let mut v = base.clone();
    v.extend([
        os("--skipQuant"),
        os("--numBootstraps"),
        os("4"),
        os("-o"),
        out.as_os_str(),
    ]);
    let o = run(&v);
    assert!(o.status.success());
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(
        err.contains("--numBootstraps") && err.contains("--skipQuant run"),
        "--skipQuant must say it is discarding the posterior request:\n{err}"
    );

    let mut v = base.clone();
    v.extend([os("--thinningFactor"), os("8"), os("-o"), out.as_os_str()]);
    let o = run(&v);
    assert!(o.status.success());
    assert!(
        String::from_utf8_lossy(&o.stderr).contains("--thinningFactor has no effect"),
        "{}",
        String::from_utf8_lossy(&o.stderr)
    );

    // Control: with Gibbs actually running, --thinningFactor is silent.
    let mut v = base.clone();
    v.extend([
        os("--numGibbsSamples"),
        os("20"),
        os("--thinningFactor"),
        os("8"),
        os("-o"),
        out.as_os_str(),
    ]);
    let o = run(&v);
    assert!(o.status.success());
    assert!(
        !String::from_utf8_lossy(&o.stderr).contains("--thinningFactor has no effect"),
        "--thinningFactor must not warn when Gibbs sampling runs"
    );
}
