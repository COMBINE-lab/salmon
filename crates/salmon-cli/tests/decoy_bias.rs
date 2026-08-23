//! Decoys must be excluded from the reported rows *and* cost no bias-correction
//! work (#1140 follow-up). Lives in salmon-cli because it drives the real
//! binary end to end: `CARGO_BIN_EXE_salmon` is only defined for tests in the
//! crate that builds the binary.

use std::process::Command;

const SALMON: &str = env!("CARGO_BIN_EXE_salmon");

/// The end-to-end guarantee: a decoy-aware index quantified through the
/// DEFAULT (deterministic) reads path must exclude decoys from `quant.sf`,
/// conserve mass over the reported rows, and — with bias correction on — must
/// not spend any work on decoy references.
///
/// The mass assertions are the ones that would catch a regression in the bias
/// bound turning into a silent accounting error: TPM is normalized over every
/// reference and rows are filtered afterward, so if a decoy ever held
/// abundance, these sums would drift without anything else complaining.
#[test]
fn decoys_are_excluded_and_mass_is_conserved_with_bias() {
    let dir = tempfile::tempdir().unwrap();

    // gentrome: 3 transcripts + a decoy that CONTAINS a verbatim copy of txA,
    // the processed-pseudogene case decoys exist for — so the decoy is a
    // genuine competitor for txA's reads rather than inert filler.
    let mut s: u64 = 0x2545_F491_4F6C_DD1D;
    let mut base = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        b"ACGT"[(s >> 33) as usize & 3] as char
    };
    let mk = |n: usize, f: &mut dyn FnMut() -> char| (0..n).map(|_| f()).collect::<String>();
    let (txa, txb, txc) = (
        mk(3000, &mut base),
        mk(3000, &mut base),
        mk(2000, &mut base),
    );
    let decoy = format!("{}{}{}", mk(20000, &mut base), txa, mk(20000, &mut base));
    let fasta = dir.path().join("gentrome.fa");
    std::fs::write(
        &fasta,
        format!(">txA\n{txa}\n>txB\n{txb}\n>txC\n{txc}\n>decoy1\n{decoy}\n"),
    )
    .unwrap();
    let decoys = dir.path().join("decoys.txt");
    std::fs::write(&decoys, "decoy1\n").unwrap();

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
    for i in 0..600 {
        let src = [&txa, &txb, &txc][i % 3];
        let st = (i * 11) % (src.len() - 350);
        r1.push_str(&format!("@f{i}\n{}\n+\n{q}\n", &src[st..st + 100]));
        r2.push_str(&format!(
            "@f{i}\n{}\n+\n{q}\n",
            rc(&src[st + 200..st + 300])
        ));
    }
    let (p1, p2) = (dir.path().join("r1.fq"), dir.path().join("r2.fq"));
    std::fs::write(&p1, r1).unwrap();
    std::fs::write(&p2, r2).unwrap();

    let index = dir.path().join("idx");
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

    let out = dir.path().join("q");
    let run = Command::new(SALMON)
        .args(["quant", "-i"])
        .arg(&index)
        .args(["-l", "IU", "-1"])
        .arg(&p1)
        .arg("-2")
        .arg(&p2)
        .args(["--seqBias", "--gcBias", "-p", "2", "-o"])
        .arg(&out)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );

    let sf = std::fs::read_to_string(out.join("quant.sf")).unwrap();
    let rows: Vec<Vec<&str>> = sf
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    assert_eq!(rows.len(), 3, "decoy must not get a quant.sf row: {sf}");
    assert!(
        !sf.contains("decoy"),
        "no decoy may appear in quant.sf: {sf}"
    );

    let tpm: f64 = rows.iter().map(|r| r[3].parse::<f64>().unwrap()).sum();
    assert!(
        (tpm - 1e6).abs() < 1.0,
        "TPM over the reported rows must sum to 1e6, got {tpm} — mass on an \
         excluded decoy is exactly how this drifts"
    );
    let reads: f64 = rows.iter().map(|r| r[4].parse::<f64>().unwrap()).sum();
    let meta = std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap();
    let mapped: f64 = meta
        .lines()
        .find(|l| l.contains("\"num_mapped\""))
        .and_then(|l| l.rsplit(':').next())
        .map(|v| v.trim().trim_end_matches(',').parse().unwrap())
        .expect("num_mapped in meta_info.json");
    assert!(
        (reads - mapped).abs() < 1.0,
        "NumReads over reported rows ({reads}) must account for every mapped \
         fragment ({mapped})"
    );

    // The invariant guard must stay silent on a healthy run.
    let err = String::from_utf8_lossy(&run.stderr);
    assert!(
        !err.contains("excluded from quant.sf"),
        "the unreported-mass guard must not fire here: {err}"
    );
}
