//! End-to-end reads-mode quantification: build an index over several
//! transcripts, simulate paired-end reads at known abundances, quantify, and
//! check the recovered counts track the truth — for both the selective and
//! pseudoalignment paths.
//!
//! # Why simulate
//!
//! With real data there is no ground truth to compare against, so a test can
//! only check that nothing crashed. Simulating reads from transcripts at chosen
//! abundances makes the right answer known in advance, which is the only way to
//! test that the pipeline as a whole — index, mapping, models, EM, output —
//! actually recovers it.

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
    simulate_scaled(dir, 1)
}

/// [`simulate`] with `scale` times as many fragments per transcript. Tests that
/// depend on more than one mapping worker actually running need enough
/// fragments to fill several paraseq batches (16384 records each) — at `scale`
/// 1 the whole fixture is a single batch handled by a single worker.
fn simulate_scaled(dir: &Path, scale: u64) -> (PathBuf, PathBuf, PathBuf, HashMap<String, u64>) {
    // (name, length seed, length, true fragment count)
    let specs = [
        ("t0", 11u64, 600usize, 300 * scale),
        ("t1", 22, 900, 100 * scale),
        ("t2", 33, 1200, 600 * scale),
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

    // The one-pass path labels itself in meta_info (#1140), so the
    // online-to-deterministic transition is auditable from existing output.
    let meta: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(out.join("aux_info").join("meta_info.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(meta["inference_path"].as_str().unwrap(), "online");

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

/// `lib_format_counts.json` must report the wrong-strand fragment count.
///
/// On a stranded protocol, fragments that map but land on the strand the library
/// type does not expect are the direct signal of double-stranded input (genomic
/// DNA carry-over), because half of such a fragment population lands on the
/// expected strand by chance and is then indistinguishable from RNA. Salmon
/// already discards the wrong-strand half, but used to report
/// `compatible_fragment_ratio: 1.0` / `num_compatible_fragments == num_mapped`
/// unconditionally, so the evidence never reached the user (#1130).
///
/// The simulated reads are inward FR (observed `ISF`): under `ISF` every fragment
/// is compatible, under the opposite `ISR` every fragment is incompatible. The
/// two runs therefore pin both ends of the new field.
#[test]
fn stranded_lib_format_counts_report_incompatible_fragments() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2, truth) = simulate(tmp.path());
    let total_truth: u64 = truth.values().sum();

    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let lib_counts = |dir: &Path| -> serde_json::Value {
        let s = std::fs::read_to_string(dir.join("lib_format_counts.json")).unwrap();
        serde_json::from_str(&s).unwrap()
    };
    let u64_field = |v: &serde_json::Value, k: &str| v[k].as_u64().unwrap();

    let quant = |lib_type: &str, out: PathBuf| {
        let mut opts = QuantOptions::new(idx_dir.clone(), out);
        opts.mates1 = vec![r1.clone()];
        opts.mates2 = vec![r2.clone()];
        opts.lib_type = lib_type.to_string();
        opts.num_threads = 1;
        quantify(&opts).unwrap_or_else(|e| panic!("quantify {lib_type}: {e}"))
    };

    // Matching stranded type: essentially every mapped fragment is compatible.
    let out_isf = tmp.path().join("quant_isf");
    let res_isf = quant("ISF", out_isf.clone());
    let isf = lib_counts(&out_isf);
    assert_eq!(
        u64_field(&isf, "num_compatible_fragments"),
        res_isf.num_mapped,
        "ISF: {isf}"
    );
    assert_eq!(
        u64_field(&isf, "num_incompatible_fragments"),
        0,
        "ISF: {isf}"
    );
    assert_eq!(isf["compatible_fragment_ratio"].as_f64().unwrap(), 1.0);

    // Opposite stranded type: every fragment that mapped is on the wrong strand,
    // so the count that used to be invisible now equals (near enough) the whole
    // library, and the ratio collapses to ~0.
    let out_isr = tmp.path().join("quant_isr");
    let res_isr = quant("ISR", out_isr.clone());
    let isr = lib_counts(&out_isr);
    let incompat = u64_field(&isr, "num_incompatible_fragments");
    assert!(
        incompat as f64 >= 0.9 * total_truth as f64,
        "ISR: only {incompat} of {total_truth} fragments reported incompatible: {isr}"
    );
    assert!(
        u64_field(&isr, "num_compatible_fragments") as f64 <= 0.1 * total_truth as f64,
        "ISR: {isr}"
    );
    assert!(
        isr["compatible_fragment_ratio"].as_f64().unwrap() <= 0.1,
        "ISR: {isr}"
    );
    // The incompatible fragments were dropped (default `--incompatPrior` 0), so
    // they are not in the assigned total: the two counts answer different
    // questions and must not be conflated.
    assert_eq!(
        u64_field(&isr, "num_assigned_fragments"),
        res_isr.num_mapped,
        "ISR: {isr}"
    );
    // The observed-format histogram is measured pre-filter, so the wrong-`-l`
    // run still records the library's true (inward-FR) shape, and the derived
    // fields expose the mismatch: almost nothing is concordant with `ISR`, and
    // the strand bias sits near 1 (nearly everything on the ISF side of the
    // inward orientation).
    assert_eq!(isr["expected_format"].as_str().unwrap(), "ISR");
    assert!(
        u64_field(&isr, "ISF") as f64 >= 0.9 * total_truth as f64,
        "ISR run: {isr}"
    );
    assert!(
        u64_field(&isr, "num_frags_with_concordant_consistent_mappings") as f64
            <= 0.1 * total_truth as f64,
        "ISR run: {isr}"
    );
    // C++ convention: the bias is zeroed when nothing agreed with the expected
    // format at all (the `numAgree > 0` gate in `summarizeLibraryTypeCounts`).
    assert_eq!(
        isr["strand_mapping_bias"].as_f64().unwrap(),
        0.0,
        "ISR run: {isr}"
    );
    // …and the matching run shows the complementary picture.
    assert_eq!(isf["expected_format"].as_str().unwrap(), "ISF");
    assert!(
        u64_field(&isf, "num_frags_with_concordant_consistent_mappings") as f64
            >= 0.9 * total_truth as f64,
        "ISF run: {isf}"
    );
    assert!(
        isf["strand_mapping_bias"].as_f64().unwrap() >= 0.9,
        "ISF run: {isf}"
    );
}

/// Sketch mode must report the wrong-strand tally in `lib_format_counts.json`,
/// consistent with the strand filter it now applies (#1136).
///
/// Since the #1136 fix, a sketch proper pair carries its observed
/// `LibraryFormat`, so `is_compatible` judges it exactly as in selective
/// alignment: the tally counts every mapped fragment before the filter drops
/// the incompatible ones. A wrong `-l` therefore shows a low
/// `compatible_fragment_ratio` *and* a collapsed mapping rate, and the JSON is
/// what says why. `pseudoalignment_strand_filter_respects_library_type` pins
/// the mapping side; this test pins the reported counts against it.
#[test]
fn sketch_lib_format_counts_report_dropped_wrong_strand_fragments() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2, truth) = simulate(tmp.path());
    let total_truth: u64 = truth.values().sum();

    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let quant = |lib_type: &str, out: PathBuf| {
        let mut opts = QuantOptions::new(idx_dir.clone(), out);
        opts.mates1 = vec![r1.clone()];
        opts.mates2 = vec![r2.clone()];
        opts.lib_type = lib_type.to_string();
        opts.sketch = true;
        opts.num_threads = 1;
        quantify(&opts).unwrap_or_else(|e| panic!("sketch quantify {lib_type}: {e}"))
    };
    let lib_counts = |dir: &Path| -> serde_json::Value {
        let s = std::fs::read_to_string(dir.join("lib_format_counts.json")).unwrap();
        serde_json::from_str(&s).unwrap()
    };

    // The reads are inward FR, i.e. observed `ISF`: everything is compatible.
    let out_isf = tmp.path().join("sketch_isf");
    let res_isf = quant("ISF", out_isf.clone());
    let isf = lib_counts(&out_isf);
    assert_eq!(
        isf["num_incompatible_fragments"].as_u64().unwrap(),
        0,
        "sketch ISF: {isf}"
    );
    assert!(
        isf["compatible_fragment_ratio"].as_f64().unwrap() >= 0.999,
        "sketch ISF: {isf}"
    );
    assert_eq!(
        isf["num_assigned_fragments"].as_u64().unwrap(),
        res_isf.num_mapped,
        "sketch ISF: {isf}"
    );

    // Opposite strandedness: every pair is wrong-strand. The tally judges the
    // fragment before the filter drops it, so the count is measured even though
    // (under the default `--incompatPrior 0`) almost nothing survives to be
    // assigned.
    let out_isr = tmp.path().join("sketch_isr");
    let res_isr = quant("ISR", out_isr.clone());
    let isr = lib_counts(&out_isr);
    let incompat = isr["num_incompatible_fragments"].as_u64().unwrap();
    assert!(
        incompat as f64 >= 0.9 * total_truth as f64,
        "sketch ISR: only {incompat} of {total_truth} counted as wrong-strand: {isr}"
    );
    assert!(
        isr["compatible_fragment_ratio"].as_f64().unwrap() <= 0.1,
        "sketch ISR: {isr}"
    );
    // The dropped fragments are excluded from assignment, not silently kept:
    // the wrong-strand run maps a small fraction of what the matching run does,
    // and the JSON's assigned count agrees with the quantifier's.
    assert!(
        (res_isr.num_mapped as f64) <= 0.05 * (res_isf.num_mapped as f64),
        "sketch ISR kept {} of {} fragments a strand filter should drop",
        res_isr.num_mapped,
        res_isf.num_mapped
    );
    assert_eq!(
        isr["num_assigned_fragments"].as_u64().unwrap(),
        res_isr.num_mapped,
        "sketch ISR: {isr}"
    );
    let mass: f64 = res_isr.counts.iter().sum();
    assert!(
        (mass - res_isr.num_mapped as f64).abs() < 1.0,
        "sketch ISR mass {mass} vs num_mapped {}",
        res_isr.num_mapped
    );
    // The observed-format histogram is measured pre-filter in sketch mode too:
    // the wrong-`-l` run still records the library's true inward-FR shape (the
    // bias itself is zeroed by the C++ `numAgree > 0` convention, since nothing
    // agreed with `ISR`).
    assert!(
        isr["ISF"].as_u64().unwrap() as f64 >= 0.9 * total_truth as f64,
        "sketch ISR: {isr}"
    );
    assert!(
        isr["num_frags_with_inconsistent_or_orphan_mappings"]
            .as_u64()
            .unwrap() as f64
            >= 0.9 * total_truth as f64,
        "sketch ISR: {isr}"
    );
}

/// The `--deterministic` mapping pass must report the same wrong-strand count.
///
/// In the full `--deterministic` flow this file is later *rewritten* by the
/// phase-2 RAD requant, which measures the same tallies from the RAD's stored
/// orientations (covered by `salmon-align`'s `rad_lib_format_counts` test) —
/// so end to end the counts survive. This test pins the mapping pass's own
/// file, which is what a standalone `--writeRad --skipQuant` run ships. The
/// pass does no strand filtering of its own (incompatible placements are
/// dropped later, in the requant), but with an explicit `-l` it can still
/// count them.
#[test]
fn deterministic_mapping_pass_reports_incompatible_fragments() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2, truth) = simulate(tmp.path());
    let total_truth: u64 = truth.values().sum();

    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    // Phase 1 of `--deterministic`, as the CLI configures it: map, write the RAD,
    // skip quantification. Reads are inward FR, so `ISR` is the wrong strand.
    let out = tmp.path().join("quant_det");
    let mut opts = QuantOptions::new(idx_dir, out.clone());
    opts.mates1 = vec![r1];
    opts.mates2 = vec![r2];
    opts.lib_type = "ISR".to_string();
    opts.num_threads = 1;
    opts.write_rad = Some(out.join("intermediate_mappings.rad"));
    opts.skip_quant = true;
    opts.deterministic_fld = true;
    std::fs::create_dir_all(&out).unwrap();
    quantify(&opts).expect("deterministic mapping pass");

    let s = std::fs::read_to_string(out.join("lib_format_counts.json")).unwrap();
    let v: serde_json::Value = serde_json::from_str(&s).unwrap();
    let incompat = v["num_incompatible_fragments"].as_u64().unwrap();
    assert!(
        incompat as f64 >= 0.9 * total_truth as f64,
        "deterministic pass: only {incompat} of {total_truth} reported incompatible: {v}"
    );
    assert!(
        v["compatible_fragment_ratio"].as_f64().unwrap() <= 0.1,
        "deterministic pass: {v}"
    );
    // The deterministic pass measures the observed-format histogram too, so
    // the phase-1 file (what a standalone `--writeRad --skipQuant` run ships)
    // carries the same derived fields as every other mode.
    assert!(
        v["ISF"].as_u64().unwrap() as f64 >= 0.9 * total_truth as f64,
        "deterministic pass: {v}"
    );
    assert!(
        v["num_frags_with_inconsistent_or_orphan_mappings"]
            .as_u64()
            .unwrap() as f64
            >= 0.9 * total_truth as f64,
        "deterministic pass: {v}"
    );
}

/// Read back a salmon RAD file: detected profile, ref names, and all records.
fn read_rad(
    path: &Path,
) -> (
    salmon_rad::RadInputProfile,
    Vec<String>,
    Vec<salmon_rad::SalmonBulkRecord>,
) {
    use libradicl::chunk::Chunk;
    use libradicl::header::RadPrelude;
    use std::io::{Cursor, Read};

    let mut data = Vec::new();
    std::fs::File::open(path)
        .unwrap()
        .read_to_end(&mut data)
        .unwrap();
    let mut rc = Cursor::new(data);
    let prelude = RadPrelude::from_bytes(&mut rc).unwrap();
    let profile = salmon_rad::detect_input_profile(&prelude).unwrap();
    let _ = prelude.file_tags.parse_tags_from_bytes(&mut rc).unwrap();
    let ctx: salmon_rad::SalmonBulkContext = prelude.get_record_context().unwrap();
    let mut recs = Vec::new();
    for _ in 0..prelude.hdr.num_chunks {
        let chunk = Chunk::<salmon_rad::SalmonBulkRecord>::from_bytes(&mut rc, &ctx);
        recs.extend(chunk.reads);
    }
    (profile, prelude.hdr.ref_names, recs)
}

/// `--writeRad` produces a readable RAD file: one record per mapped fragment,
/// ref table matches quant.sf, every hit references a valid transcript, and the
/// profile (sketch vs selective-alignment, the latter carrying scores) round-
/// trips. Runs multi-threaded to exercise the concurrent chunk writer.
#[test]
/// `--writeRad` must produce a file another tool can actually read, with every
/// fragment present — a truncated or malformed RAD would only surface much later.
fn writes_rad_output_readable_and_complete() {
    use salmon_rad::{RadInputProfile, RadProfile};

    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2, _truth) = simulate(tmp.path());
    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    // (sketch?, expected profile)
    let cases = [
        (true, RadInputProfile::Salmon(RadProfile::Sketch)),
        (
            false,
            RadInputProfile::Salmon(RadProfile::SelectiveAlignment),
        ),
    ];
    for (sketch, expected_profile) in cases {
        let tag = if sketch { "sketch" } else { "sa" };
        let out = tmp.path().join(format!("quant_rad_{tag}"));
        let rad_path = tmp.path().join(format!("maps_{tag}.rad"));
        let mut opts = QuantOptions::new(idx_dir.clone(), out);
        opts.mates1 = vec![r1.clone()];
        opts.mates2 = vec![r2.clone()];
        opts.lib_type = "IU".to_string();
        opts.num_threads = 4; // exercise the concurrent writer
        opts.sketch = sketch;
        opts.write_rad = Some(rad_path.clone());

        let res = quantify(&opts).expect("quantify with --writeRad");
        let (profile, ref_names, recs) = read_rad(&rad_path);

        assert_eq!(profile, expected_profile, "{tag}: profile detection");
        // one record per mapped fragment
        assert_eq!(
            recs.len() as u64,
            res.num_mapped,
            "{tag}: RAD record count must equal num_mapped"
        );
        // RAD reference table matches quant.sf order
        assert_eq!(ref_names, res.names, "{tag}: ref names/order");
        // every hit references a valid (non-out-of-range) transcript
        for r in &recs {
            assert!(!r.hits.is_empty(), "{tag}: mapped fragment must have hits");
            for h in &r.hits {
                assert!((h.tid as usize) < res.names.len(), "{tag}: tid in range");
            }
        }
        // SA records carry a (nontrivial) score; sketch records do not persist one.
        if !sketch {
            assert!(
                recs.iter().flat_map(|r| &r.hits).any(|h| h.score != 0),
                "SA RAD should carry nonzero alignment scores"
            );
        } else {
            assert!(
                recs.iter().flat_map(|r| &r.hits).all(|h| h.score == 0),
                "sketch RAD should not persist scores"
            );
        }
    }
}

fn read_sam_records(path: &Path) -> (noodles_sam::Header, Vec<noodles_sam::alignment::RecordBuf>) {
    let file = std::io::BufReader::new(std::fs::File::open(path).unwrap());
    let mut reader = noodles_sam::io::Reader::new(file);
    let header = reader.read_header().unwrap();
    let records = reader
        .record_bufs(&header)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    (header, records)
}

fn read_bam_records(path: &Path) -> (noodles_sam::Header, Vec<noodles_sam::alignment::RecordBuf>) {
    let mut reader = noodles_bam::io::Reader::new(std::fs::File::open(path).unwrap());
    let header = reader.read_header().unwrap();
    let records = reader
        .record_bufs(&header)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    (header, records)
}

/// Total order over records that is independent of which worker emitted them.
/// Every field the two encoders could disagree on is in the key, so sorting by
/// it and comparing element-wise gives a structural diff on failure rather than
/// "these two big vectors differ".
fn record_sort_key(
    record: &noodles_sam::alignment::RecordBuf,
) -> (Vec<u8>, u16, Option<usize>, Option<usize>, u8) {
    (
        record.name().map(|n| n.to_vec()).unwrap_or_default(),
        record.flags().bits(),
        record.reference_sequence_id(),
        record.alignment_start().map(usize::from),
        record.mapping_quality().map_or(255, |q| q.get()),
    )
}

/// Quantify `reads` once per output format and return the parsed results.
fn quant_to_sam_and_bam(
    tmp: &Path,
    threads: usize,
    scale: u64,
) -> (
    (noodles_sam::Header, Vec<noodles_sam::alignment::RecordBuf>),
    (noodles_sam::Header, Vec<noodles_sam::alignment::RecordBuf>),
) {
    let (fasta, r1, r2, _truth) = simulate_scaled(tmp, scale);
    let index = tmp.join("idx");
    let mut build_opts = IndexBuildOptions::new(vec![fasta], index.clone());
    build_opts.threads = 1;
    build(&build_opts).expect("build index");

    let quant = |out: &str, set: &dyn Fn(&mut QuantOptions)| {
        let mut opts = QuantOptions::new(index.clone(), tmp.join(out));
        opts.mates1 = vec![r1.clone()];
        opts.mates2 = vec![r2.clone()];
        opts.lib_type = "IU".into();
        opts.num_threads = threads;
        set(&mut opts);
        quantify(&opts).unwrap_or_else(|e| panic!("{out} quantification: {e}"));
    };

    let sam_path = tmp.join("mappings.sam");
    quant("quant_sam", &|o: &mut QuantOptions| {
        o.write_mappings = Some(sam_path.clone())
    });
    let bam_path = tmp.join("mappings.bam");
    quant("quant_bam", &|o: &mut QuantOptions| {
        o.write_bam = Some(bam_path.clone())
    });

    (read_sam_records(&sam_path), read_bam_records(&bam_path))
}

#[test]
/// SAM and BAM are two encodings of the same records, so the two writers must
/// agree record for record; they share the field derivation precisely so that
/// they cannot drift, and this is what proves it.
fn sam_and_bam_mapping_outputs_are_semantically_equal() {
    let tmp = tempfile::tempdir().unwrap();
    let ((sam_header, mut sam_records), (bam_header, mut bam_records)) =
        quant_to_sam_and_bam(tmp.path(), 4, 1);

    // The header is built once and shared by both encoders, so it must match in
    // full — not just the reference sequences.
    assert_eq!(sam_header, bam_header);
    assert!(!bam_header.programs().as_ref().is_empty());

    // Worker completion order is intentionally unspecified, so compare the
    // record multisets. `RecordBuf: PartialEq` covers flags, coordinates, MAPQ,
    // CIGAR, mate fields, template length, sequence, qualities, and every tag.
    assert!(!bam_records.is_empty(), "no BAM records were written");
    assert_eq!(sam_records.len(), bam_records.len());
    sam_records.sort_by_key(record_sort_key);
    bam_records.sort_by_key(record_sort_key);
    for (i, (sam, bam)) in sam_records.iter().zip(&bam_records).enumerate() {
        assert_eq!(sam, bam, "record {i} differs between SAM and BAM");
    }
}

/// The one ordering guarantee the BAM path makes: every record for a fragment
/// is contiguous. Workers encode a whole fragment before the chunk-size check,
/// and chunks are written whole, so a fragment can never be split across chunks
/// and interleaved with another worker's output.
///
/// Nothing stronger is promised: the order the fragments themselves appear in
/// is whatever order the workers happen to finish chunks in, which is why
/// [`bam_content_is_independent_of_thread_count`] compares multisets and not
/// sequences.
#[test]
/// The one ordering guarantee BAM output makes: a fragment's records are
/// contiguous. Readers rely on it, and it is what the chunk-boundary logic exists
/// to preserve.
fn bam_records_for_a_fragment_are_collated() {
    let tmp = tempfile::tempdir().unwrap();
    // paraseq hands out 16384-record batches, so the fixture must be several
    // batches deep before more than one worker ever runs — at the default scale
    // the entire input is one batch and this test would pass vacuously.
    const THREADS: usize = 8;
    let (_, (_, bam_records)) = quant_to_sam_and_bam(tmp.path(), THREADS, 100);
    assert!(!bam_records.is_empty());

    let mut seen: HashMap<Vec<u8>, usize> = HashMap::new();
    let mut current: Option<Vec<u8>> = None;
    for (i, record) in bam_records.iter().enumerate() {
        let name = record.name().map(|n| n.to_vec()).unwrap_or_default();
        if current.as_ref() != Some(&name) {
            // Starting a new run: this name must not have appeared before.
            if let Some(first) = seen.insert(name.clone(), i) {
                panic!(
                    "records for {} are split: a run starts at {i} but an earlier \
                     run started at {first}",
                    String::from_utf8_lossy(&name)
                );
            }
            current = Some(name);
        }
    }
    // Sanity: the fixture really does produce multi-record fragments, otherwise
    // collation would hold vacuously.
    assert!(
        bam_records.len() > seen.len(),
        "fixture has no fragment with more than one record; collation is untested"
    );
}

/// Thread count changes the order chunks land in, but must not change what is
/// written. Comparing a 1-thread run against an 8-thread run as multisets is
/// the strongest claim the design supports — and asserting no more than that is
/// deliberate, since imposing a global record order would mean serializing the
/// writer.
#[test]
/// The records themselves must not depend on the thread count, even though their
/// order does.
fn bam_content_is_independent_of_thread_count() {
    let mut runs = Vec::new();
    for threads in [1usize, 8] {
        let tmp = tempfile::tempdir().unwrap();
        let (fasta, r1, r2, _truth) = simulate(tmp.path());
        let index = tmp.path().join("idx");
        let mut build_opts = IndexBuildOptions::new(vec![fasta], index.clone());
        build_opts.threads = 1;
        build(&build_opts).expect("build index");

        let bam_path = tmp.path().join("mappings.bam");
        let mut opts = QuantOptions::new(index, tmp.path().join("quant"));
        opts.mates1 = vec![r1];
        opts.mates2 = vec![r2];
        opts.lib_type = "IU".into();
        opts.num_threads = threads;
        opts.write_bam = Some(bam_path.clone());
        quantify(&opts).expect("BAM-output quantification");

        let (_, mut records) = read_bam_records(&bam_path);
        records.sort_by_key(record_sort_key);
        runs.push(records);
        // `tmp` must outlive the read above.
    }

    let (single, many) = (&runs[0], &runs[1]);
    assert!(!single.is_empty());
    assert_eq!(
        single.len(),
        many.len(),
        "1-thread and 8-thread runs wrote different record counts"
    );
    for (i, (a, b)) in single.iter().zip(many).enumerate() {
        assert_eq!(a, b, "record {i} differs between 1- and 8-thread runs");
    }
    // Both runs must still be well-formed BAM with names on every record.
    assert!(many.iter().all(|record| record.name().is_some()));
}

/// A BAM that could not be written must not be reported as success.
///
/// The BGZF EOF block is the last thing to land in the writer's 1 MiB
/// `BufWriter`, and `BufWriter::drop` flushes only on a best-effort basis and
/// discards the error. Without an explicit flush, a full disk silently produced
/// a truncated BAM and an exit status of 0.
///
/// `/dev/full` accepts the open and fails every write with `ENOSPC`. The
/// fixture's BAM is well under 1 MiB, so nothing reaches the device until that
/// final flush — exactly the path being guarded.
#[cfg(target_os = "linux")]
#[test]
/// A write failure on the writer thread must reach the caller. Silently
/// truncating a BAM at the end of a long run would be the worst possible outcome.
fn bam_write_failure_is_reported() {
    if std::fs::OpenOptions::new()
        .write(true)
        .open("/dev/full")
        .is_err()
    {
        eprintln!("skipping: /dev/full is not writable here");
        return;
    }

    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2, _truth) = simulate(tmp.path());
    let index = tmp.path().join("idx");
    let mut build_opts = IndexBuildOptions::new(vec![fasta], index.clone());
    build_opts.threads = 1;
    build(&build_opts).expect("build index");

    let mut opts = QuantOptions::new(index, tmp.path().join("quant"));
    opts.mates1 = vec![r1];
    opts.mates2 = vec![r2];
    opts.lib_type = "IU".into();
    opts.num_threads = 2;
    opts.write_bam = Some(PathBuf::from("/dev/full"));

    let error = quantify(&opts).expect_err("writing to /dev/full must fail");
    let chain = format!("{error:#}");
    assert!(
        chain.contains("No space left on device"),
        "expected the ENOSPC to propagate, got: {chain}"
    );
}

#[test]
/// Asking for both formats is a usage error, caught rather than half-honoured.
fn sam_and_bam_outputs_are_mutually_exclusive() {
    let tmp = tempfile::tempdir().unwrap();
    let mut opts = QuantOptions::new(tmp.path().join("idx"), tmp.path().join("quant"));
    opts.write_mappings = Some(tmp.path().join("mappings.sam"));
    opts.write_bam = Some(tmp.path().join("mappings.bam"));
    let error = quantify(&opts).unwrap_err().to_string();
    assert!(error.contains("mutually exclusive"), "{error}");
}

#[test]
/// The alignment-free path must recover the same known abundances: it is a
/// speed trade, not a different answer.
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

/// Regression test for #1115: a fragment whose two mates orphan onto *different*
/// transcripts.
///
/// When neither mate pairs concordantly to a transcript, `map_read_pair_into`
/// does not run its orphan-suppression `retain`, so mate 1 can survive as a
/// `PairedEndLeft` mapping on transcript A while mate 2 survives as a
/// `PairedEndRight` on transcript B. `record` then sees two *different*
/// `MateStatus` values in one fragment.
///
/// A `debug_assert` there used to require every surviving mapping to share one
/// status, which aborted debug builds on any index where this occurs — measured
/// at 525 of 300,000 fragments on a GRCh38 decoy index. The invariant that
/// actually holds, and the only one the code depends on, is that the mappings
/// agree on *orphan-ness*: both statuses are orphans, so the orphan counter is
/// well defined either way.
///
/// The fixture forces the case directly: each fragment's R1 is drawn from `a`
/// and its R2 from `b`, so no concordant pair exists anywhere.
///
/// In a release build the assertion is compiled out and this only checks the
/// counting; in a debug build it is the actual regression check.
#[test]
fn mates_orphaning_to_different_transcripts_are_one_orphan_fragment() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();

    // Two unrelated transcripts: no shared k-mers, so nothing pairs across them.
    let a = gen_seq(900, 101);
    let b = gen_seq(900, 202);
    let fasta = dir.join("txome.fa");
    {
        let mut fa = std::fs::File::create(&fasta).unwrap();
        for (name, seq) in [("a", &a), ("b", &b)] {
            writeln!(fa, ">{name}").unwrap();
            fa.write_all(seq).unwrap();
            writeln!(fa).unwrap();
        }
    }

    // R1 from `a`, R2 (reverse complement, as an inward pair would be) from `b`.
    // Each mate places on its own transcript; the fragment has no proper pair.
    const N: usize = 200;
    let r1p = dir.join("reads_1.fq");
    let r2p = dir.join("reads_2.fq");
    let mut w = FastqWriters {
        r1: std::io::BufWriter::new(std::fs::File::create(&r1p).unwrap()),
        r2: std::io::BufWriter::new(std::fs::File::create(&r2p).unwrap()),
    };
    for i in 0..N {
        let pos = (i * 3) % (900 - READ_LEN);
        let s1 = a[pos..pos + READ_LEN].to_vec();
        let s2 = revcomp(&b[pos..pos + READ_LEN]);
        w.write_pair(i, &s1, &s2);
    }
    w.r1.flush().unwrap();
    w.r2.flush().unwrap();
    drop(w);

    let idx_dir = dir.join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let out = dir.join("quant");
    let mut opts = QuantOptions::new(idx_dir, out);
    opts.mates1 = vec![r1p];
    opts.mates2 = vec![r2p];
    opts.lib_type = "IU".to_string();
    opts.num_threads = 1;

    // The assertion fires during `record`, so reaching a result at all is the
    // regression check in a debug build.
    let res = quantify(&opts).expect("quantify must not abort on mixed-orphan fragments");

    // Every fragment placed, and each counts once — as an orphan, not a pair.
    assert_eq!(
        res.num_mapped, N as u64,
        "expected all {N} mixed-orphan fragments to map"
    );
    let counts = counts_by_name(&res);
    let total: f64 = res.counts.iter().sum();
    assert!(
        (total - N as f64).abs() < 1.0,
        "mass must be conserved: Σ NumReads {total} vs {N} fragments ({counts:?})"
    );
    // Both transcripts carry orphan evidence, so neither should be empty.
    assert!(
        counts["a"] > 0.0 && counts["b"] > 0.0,
        "both transcripts should receive orphan mass: {counts:?}"
    );
}

/// `--writeUnmappedNames` must name exactly the fragments that did not map, and
/// must not hold them all in memory to do it.
///
/// The names used to accumulate in a `Vec<String>` — one heap `String` per
/// unmapped read, retained until the process exited. On a contaminated sample
/// that is unbounded: measured at 1620 MB peak RSS for 10M unmapped fragments,
/// against 159 MB once the names are streamed to the file per batch.
///
/// The fixture mixes reads simulated from the transcriptome with reads drawn
/// from an unrelated sequence, so the expected unmapped set is known exactly.
#[test]
fn write_unmapped_names_lists_exactly_the_unmapped_fragments() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let (fasta, r1, r2, _truth) = simulate(dir);

    // Append foreign fragments that cannot map to the indexed transcriptome.
    let foreign = gen_seq(4000, 987_654);
    const N_FOREIGN: usize = 50;
    let r1f = dir.join("mixed_1.fq");
    let r2f = dir.join("mixed_2.fq");
    std::fs::copy(&r1, &r1f).unwrap();
    std::fs::copy(&r2, &r2f).unwrap();
    let mut w = FastqWriters {
        r1: std::io::BufWriter::new(std::fs::OpenOptions::new().append(true).open(&r1f).unwrap()),
        r2: std::io::BufWriter::new(std::fs::OpenOptions::new().append(true).open(&r2f).unwrap()),
    };
    let mut expected = Vec::new();
    for i in 0..N_FOREIGN {
        let pos = i * 20;
        let s1 = foreign[pos..pos + READ_LEN].to_vec();
        let s2 = revcomp(&foreign[pos + 200..pos + 200 + READ_LEN]);
        // `write_pair` names them `r{id}`; keep clear of the simulated ids.
        let id = 1_000_000 + i;
        w.write_pair(id, &s1, &s2);
        expected.push(format!("r{id}/1"));
    }
    w.r1.flush().unwrap();
    w.r2.flush().unwrap();
    drop(w);

    let idx_dir = dir.join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let out = dir.join("quant");
    let mut opts = QuantOptions::new(idx_dir, out.clone());
    opts.mates1 = vec![r1f];
    opts.mates2 = vec![r2f];
    opts.lib_type = "IU".to_string();
    opts.num_threads = 2;
    opts.write_unmapped_names = true;
    quantify(&opts).expect("quantify");

    let txt = std::fs::read_to_string(out.join("aux_info").join("unmapped_names.txt"))
        .expect("unmapped_names.txt must exist");

    // Every line is "<name> u"; the file must be complete (final newline) so a
    // truncated flush cannot pass as a full one.
    assert!(
        txt.ends_with('\n'),
        "file must end with a newline, so a truncated write is detectable"
    );
    let mut got: Vec<&str> = txt
        .lines()
        .map(|l| {
            let (name, status) = l.rsplit_once(' ').expect("each line is \"<name> u\"");
            assert_eq!(status, "u", "unexpected status in {l:?}");
            name
        })
        .collect();
    got.sort_unstable();
    let mut want: Vec<&str> = expected.iter().map(|s| s.as_str()).collect();
    want.sort_unstable();
    assert_eq!(
        got, want,
        "unmapped_names.txt must list exactly the foreign fragments"
    );
}

/// Regression test for #1136: sketch mode must apply the same strand-
/// compatibility filter selective alignment does.
///
/// Sketch proper pairs carried `format: None`, which `is_compatible` accepts
/// rather than guessing — so on a stranded library, sketch mode silently kept
/// every wrong-strand fragment at full weight (and, because the `-l A` detector
/// samples only formatted placements, auto-detection never resolved in sketch
/// mode either). The orientation was available the whole time (`is_fw` /
/// `mate_is_fw`; the RAD writer already derives the observed format from it).
///
/// The simulated reads are inward FR (observed `ISF`): under `ISF` everything
/// must map; under the opposite `ISR` (essentially) nothing may; and `-l A`
/// must now resolve to a concrete detected type.
#[test]
fn pseudoalignment_strand_filter_respects_library_type() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2, truth) = simulate(tmp.path());
    let total_truth: u64 = truth.values().sum();

    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let quant = |lib_type: &str, out: PathBuf| {
        let mut opts = QuantOptions::new(idx_dir.clone(), out);
        opts.mates1 = vec![r1.clone()];
        opts.mates2 = vec![r2.clone()];
        opts.lib_type = lib_type.to_string();
        opts.sketch = true;
        opts.num_threads = 1;
        quantify(&opts).unwrap_or_else(|e| panic!("quantify {lib_type}: {e}"))
    };

    // Matching stranded type: everything maps, mass conserved.
    let isf = quant("ISF", tmp.path().join("sk_isf"));
    assert!(
        isf.num_mapped as f64 >= 0.9 * total_truth as f64,
        "ISF: only {} / {total_truth} mapped",
        isf.num_mapped
    );
    let sum: f64 = isf.counts.iter().sum();
    assert!(
        (isf.num_mapped as f64 - sum).abs() <= 1e-6,
        "ISF: num_mapped={} but Σ counts={sum:.6}",
        isf.num_mapped
    );

    // Opposite stranded type: every proper pair is wrong-strand and must be
    // dropped (default `--incompatPrior 0`). Before the fix this was ~100%
    // mapped — the filter never saw a sketch pair.
    let isr = quant("ISR", tmp.path().join("sk_isr"));
    assert!(
        (isr.num_mapped as f64) <= 0.05 * total_truth as f64,
        "ISR: {} of {total_truth} fragments passed a filter that should drop \
         all of them (#1136 regressed)",
        isr.num_mapped
    );
    let sum: f64 = isr.counts.iter().sum();
    assert!(
        (isr.num_mapped as f64 - sum).abs() <= 1e-6,
        "ISR: num_mapped={} but Σ counts={sum:.6}",
        isr.num_mapped
    );

    // Auto-detection now receives samples in sketch mode. The fixture is far
    // below the 50k-sample lock-in budget (`detected_library_type` stays
    // `None` for a run this small in every mode), but the end-of-run resolved
    // type is a best guess from the partial samples — which, before the fix,
    // were *zero* in sketch mode, so `-l A` always fell back to `IU`.
    let auto = quant("A", tmp.path().join("sk_auto"));
    assert_eq!(
        auto.library_type, "ISF",
        "sketch `-l A` should resolve the simulated inward-FR orientation \
         from the detector's samples, not fall back to IU"
    );
    assert!(
        auto.num_mapped as f64 >= 0.9 * total_truth as f64,
        "A: only {} / {total_truth} mapped",
        auto.num_mapped
    );
}
