//! End-to-end reads-mode quantification: build an index over several
//! transcripts, simulate paired-end reads at known abundances, quantify, and
//! check the recovered counts track the truth — for both the selective and
//! pseudoalignment paths.

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use salmon_index::{build, IndexBuildOptions};
use salmon_quant::{quantify, QuantOptions};

const READ_LEN: usize = 75;

/// Deterministic high-complexity DNA sequence (LCG over the 4 bases).
fn gen_seq(n: usize, seed: u64) -> Vec<u8> {
    const B: [u8; 4] = *b"ACGT";
    let mut x = seed;
    let mut s = Vec::with_capacity(n);
    for _ in 0..n {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        s.push(B[((x >> 33) & 3) as usize]);
    }
    s
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            o => o,
        })
        .collect()
}

struct FastqWriters {
    r1: std::io::BufWriter<std::fs::File>,
    r2: std::io::BufWriter<std::fs::File>,
}

impl FastqWriters {
    fn write_pair(&mut self, id: usize, s1: &[u8], s2: &[u8]) {
        let q1 = vec![b'I'; s1.len()];
        let q2 = vec![b'I'; s2.len()];
        writeln!(self.r1, "@r{id}/1").unwrap();
        self.r1.write_all(s1).unwrap();
        writeln!(self.r1, "\n+").unwrap();
        self.r1.write_all(&q1).unwrap();
        writeln!(self.r1).unwrap();
        writeln!(self.r2, "@r{id}/2").unwrap();
        self.r2.write_all(s2).unwrap();
        writeln!(self.r2, "\n+").unwrap();
        self.r2.write_all(&q2).unwrap();
        writeln!(self.r2).unwrap();
    }
}

/// Build a transcriptome FASTA + simulated inward paired reads. Returns the
/// fasta path, the two read files, and the true fragment counts per name.
fn simulate(dir: &Path) -> (PathBuf, PathBuf, PathBuf, HashMap<String, u64>) {
    // (name, length seed, length, true fragment count)
    let specs = [
        ("t0", 11u64, 600usize, 300u64),
        ("t1", 22, 900, 100),
        ("t2", 33, 1200, 600),
    ];

    let fasta = dir.join("txome.fa");
    let mut fa = std::fs::File::create(&fasta).unwrap();
    let mut seqs = Vec::new();
    for (name, seed, len, _) in specs {
        let seq = gen_seq(len, seed);
        writeln!(fa, ">{name}").unwrap();
        fa.write_all(&seq).unwrap();
        writeln!(fa).unwrap();
        seqs.push((name, seq));
    }
    drop(fa);

    let r1p = dir.join("reads_1.fq");
    let r2p = dir.join("reads_2.fq");
    let mut w = FastqWriters {
        r1: std::io::BufWriter::new(std::fs::File::create(&r1p).unwrap()),
        r2: std::io::BufWriter::new(std::fs::File::create(&r2p).unwrap()),
    };

    let mut truth = HashMap::new();
    let mut rng = 0x1234_5678u64;
    let mut next = || {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        rng >> 33
    };
    let mut id = 0usize;
    for ((name, seq), (_, _, len, count)) in seqs.iter().zip(specs.iter()) {
        truth.insert(name.to_string(), *count);
        for _ in 0..*count {
            let frag = 200 + (next() % 100) as usize; // 200..300
            let max_start = len - frag;
            let pos = (next() as usize) % (max_start + 1);
            let s1 = seq[pos..pos + READ_LEN].to_vec();
            let s2 = revcomp(&seq[pos + frag - READ_LEN..pos + frag]);
            w.write_pair(id, &s1, &s2);
            id += 1;
        }
    }
    w.r1.flush().unwrap();
    w.r2.flush().unwrap();
    (fasta, r1p, r2p, truth)
}

fn counts_by_name(res: &salmon_quant::QuantResult) -> HashMap<String, f64> {
    res.names
        .iter()
        .cloned()
        .zip(res.counts.iter().cloned())
        .collect()
}

#[test]
fn selective_alignment_quantification_tracks_truth() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2, truth) = simulate(tmp.path());

    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let out = tmp.path().join("quant");
    let mut opts = QuantOptions::new(idx_dir, out.clone());
    opts.mates1 = vec![r1];
    opts.mates2 = vec![r2];
    opts.lib_type = "IU".to_string();
    opts.num_threads = 1;

    let res = quantify(&opts).expect("quantify");

    let total_truth: u64 = truth.values().sum();
    assert!(
        res.num_mapped as f64 >= 0.9 * total_truth as f64,
        "only {} / {} fragments mapped",
        res.num_mapped,
        total_truth
    );

    // Estimated counts should track the true fractions (unique sequences -> the
    // EM should recover near-exact proportions).
    let counts = counts_by_name(&res);
    let total_counts: f64 = res.counts.iter().sum();
    for (name, &t) in &truth {
        let frac_true = t as f64 / total_truth as f64;
        let frac_est = counts.get(name).copied().unwrap_or(0.0) / total_counts;
        assert!(
            (frac_true - frac_est).abs() < 0.05,
            "{name}: true frac {frac_true:.3} vs estimated {frac_est:.3}"
        );
    }
    // Ordering: t2 > t0 > t1.
    assert!(
        counts["t2"] > counts["t0"] && counts["t0"] > counts["t1"],
        "{counts:?}"
    );

    // Output files exist.
    assert!(out.join("quant.sf").exists());
    assert!(out.join("aux_info").join("meta_info.json").exists());

    // quant.sf has a header + one row per transcript.
    let sf = std::fs::read_to_string(out.join("quant.sf")).unwrap();
    assert!(sf.starts_with("Name\tLength\tEffectiveLength\tTPM\tNumReads\n"));
    assert_eq!(sf.lines().count(), 1 + res.names.len());
}

/// Regression test for the stranded-library mapped-count invariant.
///
/// A fragment must be counted in `num_mapped` only if it has a strand-compatible
/// mapping that actually contributes to quantification. The mapped counter was
/// incremented before the strand-compatibility filter, so on a stranded library a
/// fragment whose every mapping was strand-incompatible (dropped, contributing no
/// count) was still counted as mapped. That inflated `num_mapped` /
/// `percent_mapped` / `num_compatible_fragments` and broke the mass-conservation
/// invariant `Σ NumReads == num_mapped` (C++ salmon counts mapped post-filter, so
/// `Σ NumReads == num_mapped` there). The simulated reads are inward FR (observed
/// `ISF`): quantified under `ISF` they map and conserve mass; quantified under the
/// opposite stranded type `ISR` every proper pair is incompatible and dropped, so
/// `num_mapped` must fall to ~0 — not stay at the fragment total.
#[test]
fn stranded_num_mapped_matches_quantified_mass() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2, truth) = simulate(tmp.path());
    let total_truth: u64 = truth.values().sum();

    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    // Mass conservation must hold for any library type: every fragment counted as
    // mapped distributes its full unit of mass across transcripts.
    let mass_conserved = |res: &salmon_quant::QuantResult| {
        let sum_counts: f64 = res.counts.iter().sum();
        let tol = 1.0 + 0.005 * res.num_mapped as f64;
        assert!(
            (res.num_mapped as f64 - sum_counts).abs() <= tol,
            "mass not conserved: num_mapped={} but Σ counts={sum_counts:.3}",
            res.num_mapped,
        );
    };

    // Compatible stranded type: reads map and mass is conserved.
    let out_isf = tmp.path().join("quant_isf");
    let mut opts = QuantOptions::new(idx_dir.clone(), out_isf);
    opts.mates1 = vec![r1.clone()];
    opts.mates2 = vec![r2.clone()];
    opts.lib_type = "ISF".to_string();
    opts.num_threads = 1;
    let res_isf = quantify(&opts).expect("quantify ISF");
    assert!(
        res_isf.num_mapped as f64 >= 0.9 * total_truth as f64,
        "ISF: only {} / {} fragments mapped",
        res_isf.num_mapped,
        total_truth
    );
    mass_conserved(&res_isf);

    // Opposite stranded type: every proper pair is strand-incompatible and
    // dropped. The mapped count must reflect post-filter assignment (~0), and mass
    // conservation must still hold.
    let out_isr = tmp.path().join("quant_isr");
    let mut opts = QuantOptions::new(idx_dir, out_isr);
    opts.mates1 = vec![r1];
    opts.mates2 = vec![r2];
    opts.lib_type = "ISR".to_string();
    opts.num_threads = 1;
    let res_isr = quantify(&opts).expect("quantify ISR");
    assert!(
        (res_isr.num_mapped as f64) < 0.1 * total_truth as f64,
        "ISR: strand-incompatible fragments counted as mapped: num_mapped={} of {}",
        res_isr.num_mapped,
        total_truth
    );
    mass_conserved(&res_isr);
}

#[test]
fn pseudoalignment_quantification_tracks_truth() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2, truth) = simulate(tmp.path());

    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let out = tmp.path().join("quant_sketch");
    let mut opts = QuantOptions::new(idx_dir, out);
    opts.mates1 = vec![r1];
    opts.mates2 = vec![r2];
    opts.lib_type = "IU".to_string();
    opts.num_threads = 1;
    opts.sketch = true; // pseudoalignment-only

    let res = quantify(&opts).expect("quantify (sketch)");
    let total_truth: u64 = truth.values().sum();
    assert!(
        res.num_mapped as f64 >= 0.9 * total_truth as f64,
        "mapped {}",
        res.num_mapped
    );

    let counts = counts_by_name(&res);
    let total_counts: f64 = res.counts.iter().sum();
    for (name, &t) in &truth {
        let frac_true = t as f64 / total_truth as f64;
        let frac_est = counts.get(name).copied().unwrap_or(0.0) / total_counts;
        assert!(
            (frac_true - frac_est).abs() < 0.05,
            "sketch {name}: true {frac_true:.3} vs est {frac_est:.3}"
        );
    }
}
