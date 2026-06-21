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

#[test]
fn quant_is_byte_identical_across_thread_counts() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2, _truth) = simulate(tmp.path());

    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let run = |threads: usize, tag: &str| {
        let out = tmp.path().join(format!("quant_{tag}"));
        let mut opts = QuantOptions::new(idx_dir.clone(), out.clone());
        opts.mates1 = vec![r1.clone()];
        opts.mates2 = vec![r2.clone()];
        opts.lib_type = "IU".to_string();
        opts.num_threads = threads;
        opts.deterministic = true;
        quantify(&opts).expect("quantify");
        std::fs::read_to_string(out.join("quant.sf")).unwrap()
    };

    let p1 = run(1, "p1");
    let p4a = run(4, "p4a");
    let p4b = run(4, "p4b");

    // Multi-threaded runs must be byte-identical to each other and to the
    // single-threaded run (no thread-scheduling-dependent drift).
    assert_eq!(p4a, p4b, "two -p 4 runs differ");
    assert_eq!(p1, p4a, "-p 1 and -p 4 differ");
}

fn simulate_multimapping(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let backbone = gen_seq(1500, 99);
    let snp = |seed: u64, n: usize| {
        let mut s = backbone.clone();
        let mut x = seed;
        for _ in 0..n {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let pos = (x >> 33) as usize % s.len();
            s[pos] = match s[pos] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        s
    };
    let txs = [
        ("d0", backbone.clone()),
        ("d1", snp(1, 6)),
        ("d2", snp(2, 6)),
    ];
    let fasta = dir.join("multi.fa");
    let mut fa = std::fs::File::create(&fasta).unwrap();
    for (name, seq) in &txs {
        writeln!(fa, ">{name}").unwrap();
        fa.write_all(seq).unwrap();
        writeln!(fa).unwrap();
    }
    drop(fa);
    let r1p = dir.join("m1.fq");
    let r2p = dir.join("m2.fq");
    let mut w = FastqWriters {
        r1: std::io::BufWriter::new(std::fs::File::create(&r1p).unwrap()),
        r2: std::io::BufWriter::new(std::fs::File::create(&r2p).unwrap()),
    };
    let mut rng = 0x0BEE_F00Du64;
    let mut next = || {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        rng >> 33
    };
    // Enough fragments to span multiple paraseq batches / worker threads (so the
    // multi-threaded path is actually exercised).
    for id in 0..20000 {
        let frag = 200 + (next() % 100) as usize;
        let max_start = backbone.len() - frag;
        let pos = (next() as usize) % (max_start + 1);
        let s1 = backbone[pos..pos + READ_LEN].to_vec();
        let s2 = revcomp(&backbone[pos + frag - READ_LEN..pos + frag]);
        w.write_pair(id, &s1, &s2);
    }
    w.r1.flush().unwrap();
    w.r2.flush().unwrap();
    (fasta, r1p, r2p)
}

/// Strongest end-to-end determinism guard: compares the **full-precision** `f64`
/// abundances (`res.counts`, not the truncated `quant.sf` text) across thread
/// counts, under the hardest conditions for thread-order sensitivity — many
/// multimapping fragments (near-duplicate transcripts -> fractional eq-class
/// weights), `--seqBias --gcBias` on (online masses + bias models), and the
/// after-burn-in FLD path forced on. With stages 1 (#1027), 2a (#1028) and 2b
/// (#1030) this is bit-identical; the residual online-mass differences are far
/// too small to flip the per-fragment-seeded FLD acceptance or the binned bias
/// models.
#[test]
fn multimapping_counts_bit_identical_across_threads() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2) = simulate_multimapping(tmp.path());
    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let run = |threads: usize| -> Vec<u64> {
        let out = tmp.path().join(format!("mm_{threads}"));
        let mut opts = QuantOptions::new(idx_dir.clone(), out);
        opts.mates1 = vec![r1.clone()];
        opts.mates2 = vec![r2.clone()];
        opts.lib_type = "IU".to_string();
        opts.num_threads = threads;
        opts.seq_bias = true;
        opts.gc_bias = true;
        // force the after-burn-in FLD path on small data
        opts.num_pre_aux_model_samples = 50;
        // byte-for-byte reproducibility is opt-in (#1031 review): the two-pass,
        // key-sorted inference path is what makes this bit-identical across `-p`.
        opts.deterministic = true;
        let res = quantify(&opts).expect("quantify");
        res.counts.iter().map(|c| c.to_bits()).collect()
    };
    let p1 = run(1);
    let p8a = run(8);
    let p8b = run(8);
    assert_eq!(p8a, p8b, "two -p 8 runs differ (full precision)");
    assert_eq!(p1, p8a, "-p 1 and -p 8 differ (full precision)");
}

/// Same determinism guarantee for the single-end (unmated) path, which uses the
/// other `ParallelProcessor` impl and the single-read RAD record. Reuses the 20k
/// multimapping reads as unmated single-end input so the multi-thread store is
/// genuinely exercised.
#[test]
fn single_end_counts_bit_identical_across_threads() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, _r2) = simulate_multimapping(tmp.path());
    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let run = |threads: usize| -> Vec<u64> {
        let out = tmp.path().join(format!("se_{threads}"));
        let mut opts = QuantOptions::new(idx_dir.clone(), out);
        opts.unmated = vec![r1.clone()];
        opts.lib_type = "U".to_string();
        opts.num_threads = threads;
        opts.seq_bias = true;
        opts.num_pre_aux_model_samples = 50;
        opts.deterministic = true;
        let res = quantify(&opts).expect("quantify");
        res.counts.iter().map(|c| c.to_bits()).collect()
    };
    let p1 = run(1);
    let p8a = run(8);
    let p8b = run(8);
    assert_eq!(p8a, p8b, "two single-end -p 8 runs differ (full precision)");
    assert_eq!(p1, p8a, "single-end -p 1 and -p 8 differ (full precision)");
}

/// Determinism must hold when pass 2 streams through a multi-run external sort,
/// not just a single in-memory run. Force a tiny run size (≈20 runs over the 20k
/// multimapping fragments) so the k-way merge is exercised, and assert the result
/// is byte-identical across thread counts (and to the heavier single-run guard
/// above, which uses the same inputs).
#[test]
fn multi_run_external_sort_is_bit_identical_across_threads() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2) = simulate_multimapping(tmp.path());
    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let run = |threads: usize, run_size: Option<usize>| -> Vec<u64> {
        let out = tmp.path().join(format!("mr_{threads}_{run_size:?}"));
        let mut opts = QuantOptions::new(idx_dir.clone(), out);
        opts.mates1 = vec![r1.clone()];
        opts.mates2 = vec![r2.clone()];
        opts.lib_type = "IU".to_string();
        opts.num_threads = threads;
        opts.seq_bias = true;
        opts.gc_bias = true;
        opts.num_pre_aux_model_samples = 50;
        opts.deterministic = true;
        opts.det_run_size = run_size;
        let res = quantify(&opts).expect("quantify");
        res.counts.iter().map(|c| c.to_bits()).collect()
    };
    // ~20 runs over 20k fragments at -p 1 and -p 8, plus a single-run reference.
    let many_p1 = run(1, Some(1000));
    let many_p8 = run(8, Some(1000));
    let single = run(8, None);
    assert_eq!(
        many_p1, many_p8,
        "multi-run external sort not deterministic across -p"
    );
    assert_eq!(
        single, many_p1,
        "multi-run merge result differs from single-run result"
    );
}

/// RAD input round trip: map once with `--writeRad` to persist the mapping store,
/// then quantify from it with `--rad` (no reads). The decoupled quant must
/// reproduce the byte-identical `quant.sf`, and it must be independent of the
/// thread count used either to build the store or (irrelevant, serial) to read it.
#[test]
fn rad_input_reproduces_quant_from_written_store() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, r1, r2) = simulate_multimapping(tmp.path());
    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    let store = tmp.path().join("mappings.rad");

    // Map at -p 4, writing and keeping the RAD store.
    let out_a = tmp.path().join("a");
    let mut a = QuantOptions::new(idx_dir.clone(), out_a.clone());
    a.mates1 = vec![r1.clone()];
    a.mates2 = vec![r2.clone()];
    a.lib_type = "IU".to_string();
    a.num_threads = 4;
    a.seq_bias = true;
    a.gc_bias = true;
    a.num_pre_aux_model_samples = 50;
    a.write_rad = Some(store.clone());
    quantify(&a).expect("quantify + writeRad");
    assert!(store.exists(), "--writeRad did not keep the store");
    let sf_a = std::fs::read_to_string(out_a.join("quant.sf")).unwrap();

    // Re-quantify from the store alone (no reads), at a different thread count and
    // matching inference settings.
    let out_b = tmp.path().join("b");
    let mut b = QuantOptions::new(idx_dir.clone(), out_b.clone());
    b.lib_type = "IU".to_string();
    b.num_threads = 8;
    b.seq_bias = true;
    b.gc_bias = true;
    b.num_pre_aux_model_samples = 50;
    b.rad_input = Some(store.clone());
    quantify(&b).expect("quantify from --rad");
    let sf_b = std::fs::read_to_string(out_b.join("quant.sf")).unwrap();

    assert_eq!(
        sf_a, sf_b,
        "RAD-input quant differs from the store-producing run"
    );
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

/// Hand-write a piscem bulk/sketch-mode RAD file (`piscem map` profile) holding
/// `records`, each `(frag_map_type, [(orientation_2bit, ref_id, pos, frag_len)])`.
///
/// The prelude (header + the canonical piscem tag sections, including a
/// `ref_lengths` file tag whose *values* sit between the prelude and the chunk)
/// is written via `libradicl` so it is byte-for-byte the format piscem produces.
/// The chunk records are written by hand because `libradicl`'s record *writer*
/// encodes the orientation with a 3-bit scheme that its own reader
/// (`from_u32_paired_status`, 2 bits) does not round-trip; piscem itself writes
/// the 2-bit code, which is what we emit here so the file decodes faithfully.
fn write_piscem_rad(
    path: &Path,
    ref_names: &[&str],
    ref_lens: &[u32],
    records: &[(u8, Vec<(u32, u32, u32, u16)>)],
) {
    use libradicl::header::{RadHeader, RadPrelude};
    use libradicl::rad_types::{
        RadAtomicId, RadIntId, RadType, TagDesc, TagMap, TagSection, TagSectionLabel, TagValue,
    };

    let int_tag = |name: &str, id| TagDesc {
        name: name.to_string(),
        typeid: RadType::Int(id),
    };

    let mut file_tags = TagSection::new_with_label(TagSectionLabel::FileTags);
    file_tags.add_tag_desc(TagDesc {
        name: "ref_lengths".to_string(),
        typeid: RadType::Array(RadIntId::U32, RadAtomicId::Int(RadIntId::U32)),
    });
    let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
    read_tags.add_tag_desc(int_tag("frag_map_type", RadIntId::U8));
    let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    aln_tags.add_tag_desc(int_tag("compressed_ori_ref", RadIntId::U32));
    aln_tags.add_tag_desc(int_tag("frag_map_type", RadIntId::U32));
    aln_tags.add_tag_desc(int_tag("frag_len", RadIntId::U16));

    let prelude = RadPrelude {
        hdr: RadHeader {
            is_paired: 1,
            ref_count: ref_names.len() as u64,
            ref_names: ref_names.iter().map(|s| s.to_string()).collect(),
            num_chunks: 1,
        },
        file_tags,
        read_tags,
        aln_tags,
    };

    let mut buf: Vec<u8> = Vec::new();
    prelude.write(&mut buf).unwrap();
    // file-level tag *values*: the ref_lengths array.
    let mut file_tag_map = TagMap::with_keyset(&prelude.file_tags.tags);
    file_tag_map.add(TagValue::ArrayU32(ref_lens.to_vec()));
    file_tag_map.write_values(&mut buf).unwrap();

    // one chunk: nbytes (incl. its 8-byte header), nrec, then the records.
    let mut body: Vec<u8> = Vec::new();
    for (frag_type, alns) in records {
        body.extend_from_slice(&(alns.len() as u32).to_le_bytes());
        body.push(*frag_type);
        for (ori, ref_id, pos, flen) in alns {
            let packed: u32 = (ori << 30) | ref_id;
            body.extend_from_slice(&packed.to_le_bytes());
            body.extend_from_slice(&pos.to_le_bytes());
            body.extend_from_slice(&flen.to_le_bytes());
        }
    }
    let nbytes = (body.len() + 8) as u32;
    buf.extend_from_slice(&nbytes.to_le_bytes());
    buf.extend_from_slice(&(records.len() as u32).to_le_bytes());
    buf.extend_from_slice(&body);

    std::fs::write(path, buf).unwrap();
}

/// End-to-end `salmon quant --rad` over a piscem bulk/sketch-mode RAD: build an
/// index, synthesize a piscem RAD whose fragments map uniquely to known
/// transcripts, quantify, and confirm mass is conserved and assigned exactly.
#[test]
fn piscem_rad_input_quantifies() {
    let tmp = tempfile::tempdir().unwrap();
    let (fasta, _r1, _r2) = simulate_multimapping(tmp.path());
    let idx_dir = tmp.path().join("idx");
    let mut bopts = IndexBuildOptions::new(vec![fasta], idx_dir.clone());
    bopts.threads = 1;
    build(&bopts).expect("build index");

    // Proper-pair, forward-reverse (orientation code 2), one alignment per record,
    // each mapping a fragment uniquely to one transcript (tids 0/1/2).
    const MAPPED_PAIR: u8 = 4; // libradicl MappingType::MappedPair
    const FR: u32 = 2; // from_u32_paired_status(2, MappedPair) = ForwardReverse
    let per_tid = [120usize, 60, 20];
    let total: usize = per_tid.iter().sum();
    let mut records: Vec<(u8, Vec<(u32, u32, u32, u16)>)> = Vec::new();
    for (tid, &n) in per_tid.iter().enumerate() {
        for _ in 0..n {
            records.push((MAPPED_PAIR, vec![(FR, tid as u32, 100, 250)]));
        }
    }

    let rad = tmp.path().join("piscem.rad");
    write_piscem_rad(&rad, &["d0", "d1", "d2"], &[1500, 1500, 1500], &records);

    let out = tmp.path().join("piscem_quant");
    let mut opts = QuantOptions::new(idx_dir, out.clone());
    opts.lib_type = "IU".to_string();
    opts.num_threads = 1;
    opts.rad_input = Some(rad);
    let res = quantify(&opts).expect("quantify from piscem --rad");

    // Every fragment in the store is a mapped fragment.
    assert_eq!(res.num_mapped, total as u64, "num_mapped");
    // Read mass is conserved (unique assignment => no mass lost or invented).
    let total_counts: f64 = res.counts.iter().sum();
    assert!(
        (total_counts - total as f64).abs() < 0.5,
        "total counts {total_counts} vs {total}"
    );
    // Each fragment was assigned to exactly its transcript: the three nonzero
    // counts match the placement multiset regardless of transcript ordering.
    let mut got: Vec<f64> = res.counts.iter().copied().filter(|c| *c > 0.5).collect();
    got.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let mut want: Vec<f64> = per_tid.iter().map(|&n| n as f64).collect();
    want.sort_by(|a, b| b.partial_cmp(a).unwrap());
    assert_eq!(got.len(), want.len(), "expected exactly 3 expressed txps");
    for (g, w) in got.iter().zip(&want) {
        assert!((g - w).abs() < 0.5, "count {g} vs expected {w}");
    }
    // quant.sf was written.
    assert!(out.join("quant.sf").exists(), "quant.sf not written");
}
