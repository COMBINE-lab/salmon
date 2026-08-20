//! End-to-end degradation normalization: simulate two libraries of the same
//! transcripts — one intact, one degraded from the 5' end — quantify both with
//! `--degCoverage`, fit the cohort, and check the degraded library is the one
//! that gets flagged.
//!
//! # Why simulate the degradation
//!
//! The degradation index is not something a run can check against itself: it is
//! defined by comparison across samples, and on real data there is no ground
//! truth to compare against. Simulating a library whose fragments are drawn
//! only from the 3' portion of each transcript makes the answer known — that
//! sample lost a known fraction of its coverage, the other lost none — which is
//! the only way to test that the whole path (accumulate during mapping, dump,
//! reload, fit) reports it.

use std::io::Write;
use std::path::{Path, PathBuf};

use salmon_degnorm::{cohort, CohortOptions, Sample};
use salmon_index::{build, IndexBuildOptions};
use salmon_quant::{quantify, QuantOptions};

const READ_LEN: usize = 75;
const FRAG_LEN: usize = 200;

/// Deterministic high-complexity DNA (LCG over the four bases), so the mapper
/// has unique sequence to work with.
fn gen_seq(n: usize, seed: u64) -> Vec<u8> {
    const B: [u8; 4] = *b"ACGT";
    let mut x = seed;
    (0..n)
        .map(|_| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            B[((x >> 33) & 3) as usize]
        })
        .collect()
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

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

const TRANSCRIPTS: [(&str, u64, usize); 3] = [("t0", 11, 1500), ("t1", 22, 2000), ("t2", 33, 1800)];

fn write_transcriptome(dir: &Path) -> (PathBuf, Vec<(String, Vec<u8>)>) {
    let fasta = dir.join("txome.fa");
    let mut fa = std::fs::File::create(&fasta).unwrap();
    let mut seqs = Vec::new();
    for (name, seed, len) in TRANSCRIPTS {
        let seq = gen_seq(len, seed);
        writeln!(fa, ">{name}").unwrap();
        fa.write_all(&seq).unwrap();
        writeln!(fa).unwrap();
        seqs.push((name.to_string(), seq));
    }
    (fasta, seqs)
}

/// Simulate one library. `keep_from` is the fraction of each transcript, from
/// the 5' end, that has been lost: `0.0` is an intact library, `0.5` one whose
/// fragments can only start in the 3' half.
fn simulate(
    dir: &Path,
    tag: &str,
    seqs: &[(String, Vec<u8>)],
    frags_per_tx: usize,
    lost_prefix: f64,
    seed: u64,
) -> (PathBuf, PathBuf) {
    let r1p = dir.join(format!("{tag}_1.fq"));
    let r2p = dir.join(format!("{tag}_2.fq"));
    let mut r1 = std::io::BufWriter::new(std::fs::File::create(&r1p).unwrap());
    let mut r2 = std::io::BufWriter::new(std::fs::File::create(&r2p).unwrap());
    let mut rng = Lcg(seed);
    let mut id = 0usize;
    for (_, seq) in seqs {
        let first = (seq.len() as f64 * lost_prefix) as usize;
        let span = seq.len() - FRAG_LEN - first;
        for _ in 0..frags_per_tx {
            let pos = first + (rng.next() as usize) % (span + 1);
            let s1 = &seq[pos..pos + READ_LEN];
            let s2 = revcomp(&seq[pos + FRAG_LEN - READ_LEN..pos + FRAG_LEN]);
            let q = vec![b'I'; READ_LEN];
            writeln!(r1, "@f{id}/1").unwrap();
            r1.write_all(s1).unwrap();
            writeln!(r1, "\n+").unwrap();
            r1.write_all(&q).unwrap();
            writeln!(r1).unwrap();
            writeln!(r2, "@f{id}/2").unwrap();
            r2.write_all(&s2).unwrap();
            writeln!(r2, "\n+").unwrap();
            r2.write_all(&q).unwrap();
            writeln!(r2).unwrap();
            id += 1;
        }
    }
    r1.flush().unwrap();
    r2.flush().unwrap();
    (r1p, r2p)
}

fn quant_with_coverage(idx: &Path, out: PathBuf, r1: PathBuf, r2: PathBuf, bins: usize) -> PathBuf {
    let mut opts = QuantOptions::new(idx.to_path_buf(), out.clone());
    opts.mates1 = vec![r1];
    opts.mates2 = vec![r2];
    opts.lib_type = "IU".to_string();
    opts.num_threads = 2;
    opts.coverage_bins = Some(bins);
    quantify(&opts).expect("quantify");
    out
}

#[test]
fn a_degraded_library_is_flagged_and_its_counts_are_restored() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, seqs) = write_transcriptome(tmp.path());

    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    // Two intact libraries and one that lost its 5' half. The degraded library
    // also yields fewer fragments, which is exactly the count bias the
    // adjustment exists to undo.
    let specs: [(&str, usize, f64, u64); 3] = [
        ("intact_a", 4000, 0.0, 0x1111),
        ("degraded", 2000, 0.5, 0x2222),
        ("intact_b", 4000, 0.0, 0x3333),
    ];
    let mut samples = Vec::new();
    for (tag, frags, lost, seed) in specs {
        let (r1, r2) = simulate(tmp.path(), tag, &seqs, frags, lost, seed);
        let out = quant_with_coverage(&idx_dir, tmp.path().join(tag), r1, r2, 50);
        assert!(
            out.join("aux_info").join("coverage.gz").exists(),
            "{tag}: --degCoverage wrote no dump"
        );
        samples.push(Sample {
            name: tag.to_string(),
            coverage: out.join("aux_info").join("coverage.gz"),
            quant: Some(out.join("quant.sf")),
        });
    }

    let opts = CohortOptions::default();
    let res = cohort::run(&samples, &opts).expect("cohort fit");
    assert_eq!(res.num_targets(), TRANSCRIPTS.len());
    assert!(
        res.fitted.iter().all(|&f| f),
        "every simulated transcript is long and deep enough to fit: {:?}",
        res.fitted
    );

    let mean = res.mean_di();
    // The degraded library lost half of every transcript; the intact ones lost
    // nothing. Only the ordering and a loose magnitude are asserted: the index
    // has a coverage-noise floor (see the `nmf` module docs), so the intact
    // samples sit near but not at zero.
    assert!(
        mean[1] > mean[0] + 0.2 && mean[1] > mean[2] + 0.2,
        "the degraded library should stand out: {mean:?}"
    );
    assert!(
        mean[1] > 0.3,
        "half the transcript is missing, index read only {:.3}",
        mean[1]
    );
    assert!(
        mean[0] < 0.25 && mean[2] < 0.25,
        "intact libraries drifted too high: {mean:?}"
    );

    // The corrected `quant.sf` files are the artefact a tximport/DESeq2
    // pipeline actually consumes, so check they are there, complete, and carry
    // the adjusted counts rather than the raw ones.
    let out = tmp.path().join("degnorm");
    cohort::write_tables(&out, &res, &opts).expect("write tables");
    cohort::write_adjusted_quants(&out, &res, &samples).expect("write quant.sf");
    for (i, s) in samples.iter().enumerate() {
        let sf = out.join("quants").join(&s.name).join("quant.sf");
        let text = std::fs::read_to_string(&sf).unwrap_or_else(|e| panic!("{}: {e}", sf.display()));
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "Name\tLength\tEffectiveLength\tTPM\tNumReads");
        assert_eq!(lines.len(), 1 + TRANSCRIPTS.len());
        let tpm: f64 = lines[1..]
            .iter()
            .map(|l| l.split('\t').nth(3).unwrap().parse::<f64>().unwrap())
            .sum();
        assert!((tpm - 1e6).abs() < 1.0, "{}: TPM sums to {tpm}", s.name);
        let first: Vec<&str> = lines[1].split('\t').collect();
        let want = res.adjusted_counts[i];
        let got: f64 = first[4].parse().unwrap();
        assert!((got - want).abs() < 0.01, "{}: {got} vs {want}", s.name);
    }

    // Per transcript, the adjustment should close most of the count gap the
    // degradation opened between the degraded library and the intact ones.
    let m = res.num_samples();
    for t in 0..res.num_targets() {
        let raw_gap = res.raw_counts[t * m + 1] / res.raw_counts[t * m];
        let adj_gap = res.adjusted_counts[t * m + 1] / res.adjusted_counts[t * m];
        assert!(
            adj_gap > raw_gap,
            "{}: adjustment moved the wrong way ({raw_gap:.3} -> {adj_gap:.3})",
            res.target_names[t]
        );
    }
}

#[test]
fn coverage_follows_the_fragments_that_were_actually_simulated() {
    // A library whose fragments only ever start in the 3' half must show
    // (near-)zero coverage over the first bins of every transcript and real
    // coverage over the last ones. This is the accumulator's own correctness,
    // independent of any model fitted on top of it.
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, seqs) = write_transcriptome(tmp.path());
    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let (r1, r2) = simulate(tmp.path(), "half", &seqs, 3000, 0.5, 0x4444);
    let out = quant_with_coverage(&idx_dir, tmp.path().join("half"), r1, r2, 20);
    let profiles =
        salmon_degnorm::CoverageProfiles::read(&out.join("aux_info").join("coverage.gz")).unwrap();
    assert_eq!(profiles.num_bins, 20);
    assert_eq!(profiles.names.len(), TRANSCRIPTS.len());
    assert!(profiles.num_mapped > 0);

    for t in 0..profiles.num_targets() {
        let row = profiles.row(t);
        // Bin 0 is entirely inside the lost 5' half; the last bin is entirely
        // inside the surviving 3' half.
        let head = row[0];
        let tail = row[profiles.num_bins - 1];
        assert!(
            head < 0.1 * tail,
            "{}: 5' bin ({head:.3}) should be far below the 3' bin ({tail:.3})",
            profiles.names[t]
        );
        assert!(tail > 1.0, "{}: no 3' coverage at all", profiles.names[t]);
    }
}
