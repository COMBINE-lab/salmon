//! End-to-end MEM collection against a real piscem index built by salmon-index.
//!
//! Builds a single-transcript index, then maps a forward substring and its
//! reverse complement back to it, asserting we recover a high-coverage chain on
//! the correct transcript and orientation at the expected position.

use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use salmon_core::MateStatus;
use salmon_index::{build, IndexBuildOptions, SalmonIndex};
use salmon_map::{
    align_chain, best_per_target, collect_read_mems, map_read_pair, map_read_pair_sketch,
    map_single_read, map_single_read_sketch, perfect_score, AlignConfig, MapConfig,
    MemCollectorConfig, RefProvider,
};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Adapter exposing a loaded index's reference sequences to the mapper.
struct IdxRefs<'a>(&'a SalmonIndex);
impl RefProvider for IdxRefs<'_> {
    fn num_refs(&self) -> usize {
        self.0.num_refs()
    }
    fn ref_seq(&self, tid: u32) -> &[u8] {
        self.0.ref_seq(tid)
    }
}

/// Build a single-transcript index from `seq` in a fresh temp dir; returns the
/// loaded index and the temp dir (kept alive by the caller).
fn build_index_for(seq: &[u8]) -> (SalmonIndex, tempfile::TempDir) {
    let tmp = tempfile::tempdir().unwrap();
    let fasta = write_fasta(tmp.path(), "tx0", seq);
    let out = tmp.path().join("idx");
    let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
    opts.threads = 1;
    build(&opts).expect("build index");
    (SalmonIndex::load(&out).expect("load index"), tmp)
}

/// Mutate every `step`-th base so no exact k-mer survives, while most bases
/// still match (so the read can be recovered by alignment but not by k-mer hits).
fn break_kmers(seq: &[u8], step: usize) -> Vec<u8> {
    seq.iter()
        .enumerate()
        .map(|(i, &b)| {
            if i % step == 0 {
                if b == b'A' {
                    b'C'
                } else {
                    b'A'
                }
            } else {
                b
            }
        })
        .collect()
}

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
            other => other,
        })
        .collect()
}

fn write_fasta(dir: &Path, name: &str, seq: &[u8]) -> PathBuf {
    let path = dir.join("txome.fa");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, ">{name}").unwrap();
    f.write_all(seq).unwrap();
    writeln!(f).unwrap();
    path
}

#[test]
fn maps_forward_and_revcomp_reads() {
    let tmp = tempfile::tempdir().unwrap();
    let transcript = gen_seq(300, 0xC0FFEE);
    let fasta = write_fasta(tmp.path(), "tx0", &transcript);
    let out = tmp.path().join("idx");

    let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
    opts.threads = 1;
    let info = build(&opts).expect("build index");
    assert_eq!(info.num_refs, 1);

    let idx = SalmonIndex::load(&out).expect("load index");
    let mut hs = HitSearcher::new(idx.inner());
    let cfg = MemCollectorConfig::default();

    // Forward read: bases [50, 150) of the transcript.
    let fwd = &transcript[50..150];
    let cands = best_per_target(collect_read_mems(idx.inner(), &mut hs, fwd, true, &cfg));
    assert!(!cands.is_empty(), "forward read produced no candidates");
    let best = cands
        .iter()
        .max_by_key(|c| c.chain.covered_read_bases())
        .unwrap();
    assert_eq!(best.tid, 0);
    assert!(best.is_fw, "forward substring should map forward");
    assert!(
        best.chain.covered_read_bases() >= 60,
        "low coverage: {}",
        best.chain.covered_read_bases()
    );
    // The chain should start near reference position 50.
    assert!(
        (best.chain.ref_start() - 50).abs() <= 5,
        "ref_start {} not near 50",
        best.chain.ref_start()
    );

    // Reverse-complement of the same window should map to the same place, RC.
    let rc = revcomp(fwd);
    let cands_rc = best_per_target(collect_read_mems(idx.inner(), &mut hs, &rc, true, &cfg));
    assert!(!cands_rc.is_empty(), "rc read produced no candidates");
    let best_rc = cands_rc
        .iter()
        .max_by_key(|c| c.chain.covered_read_bases())
        .unwrap();
    assert_eq!(best_rc.tid, 0);
    assert!(!best_rc.is_fw, "rc substring should map reverse");
    assert!(
        best_rc.chain.covered_read_bases() >= 60,
        "low rc coverage: {}",
        best_rc.chain.covered_read_bases()
    );
    assert!(
        (best_rc.chain.ref_start() - 50).abs() <= 5,
        "rc ref_start {} not near 50",
        best_rc.chain.ref_start()
    );
}

#[test]
fn collected_candidate_validates_by_alignment() {
    let tmp = tempfile::tempdir().unwrap();
    let transcript = gen_seq(300, 0xBEEF);
    let fasta = write_fasta(tmp.path(), "tx0", &transcript);
    let out = tmp.path().join("idx");
    let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
    opts.threads = 1;
    build(&opts).expect("build index");
    let idx = SalmonIndex::load(&out).expect("load index");
    let mut hs = HitSearcher::new(idx.inner());
    let cfg = MemCollectorConfig::default();
    let acfg = AlignConfig::default();

    // Exact substring -> best candidate aligns perfectly against the transcript.
    let read = &transcript[60..160];
    let cands = best_per_target(collect_read_mems(idx.inner(), &mut hs, read, true, &cfg));
    let best = cands
        .iter()
        .max_by_key(|c| c.chain.covered_read_bases())
        .expect("a candidate");
    let aln = align_chain(read, &transcript, &best.chain, &acfg).expect("alignment");
    assert!(aln.valid);
    assert_eq!(aln.score, perfect_score(read.len(), &acfg));

    // Introduce two mismatches -> still validates, score below perfect.
    let mut mut_read = read.to_vec();
    mut_read[20] = if mut_read[20] == b'A' { b'C' } else { b'A' };
    mut_read[70] = if mut_read[70] == b'G' { b'T' } else { b'G' };
    let cands2 = best_per_target(collect_read_mems(
        idx.inner(),
        &mut hs,
        &mut_read,
        true,
        &cfg,
    ));
    let best2 = cands2
        .iter()
        .max_by_key(|c| c.chain.covered_read_bases())
        .expect("a candidate for the mutated read");
    let aln2 = align_chain(&mut_read, &transcript, &best2.chain, &acfg).expect("alignment");
    assert!(aln2.valid);
    assert!(aln2.score < perfect_score(mut_read.len(), &acfg));
}

#[test]
fn assembled_single_read_mapper() {
    let transcript = gen_seq(600, 0xA11CE);
    let (idx, _tmp) = build_index_for(&transcript);
    let refs = IdxRefs(&idx);
    let mut hs = HitSearcher::new(idx.inner());
    let cfg = MapConfig::default();

    let read = &transcript[100..200];
    let maps = map_single_read(idx.inner(), &mut hs, &refs, read, &cfg);
    assert_eq!(maps.len(), 1, "{maps:?}");
    assert_eq!(maps[0].tid, 0);
    assert_eq!(maps[0].status, MateStatus::SingleEnd);
    assert!((maps[0].weight - 1.0).abs() < 1e-9);

    // foreign read -> unmapped
    let foreign = gen_seq(100, 0xDEAD);
    assert!(map_single_read(idx.inner(), &mut hs, &refs, &foreign, &cfg).is_empty());
}

#[test]
fn assembled_paired_read_mapper_concordant() {
    let transcript = gen_seq(600, 0xB0B);
    let (idx, _tmp) = build_index_for(&transcript);
    let refs = IdxRefs(&idx);
    let mut hs = HitSearcher::new(idx.inner());
    let cfg = MapConfig::default();

    // mate1 forward upstream, mate2 reverse-complement downstream -> inward pair
    let r1 = transcript[50..130].to_vec();
    let r2 = revcomp(&transcript[350..430]);
    let maps = map_read_pair(idx.inner(), &mut hs, &refs, &r1, &r2, &cfg);
    assert_eq!(maps.len(), 1, "{maps:?}");
    assert_eq!(maps[0].tid, 0);
    assert_eq!(maps[0].status, MateStatus::PairedEndPaired);
    // fragment spans ~ [50, 430)
    assert!(
        (maps[0].fragment_len - 380).abs() <= 5,
        "frag {}",
        maps[0].fragment_len
    );
}

#[test]
fn orphan_recovery_promotes_to_pair() {
    let transcript = gen_seq(700, 0x0CEA1);
    let (idx, _tmp) = build_index_for(&transcript);
    let refs = IdxRefs(&idx);
    let mut hs = HitSearcher::new(idx.inner());

    let r1 = transcript[60..160].to_vec(); // maps cleanly, forward
                                           // mate2: rc of a downstream window, with every 25th base broken so it has
                                           // no exact k-mers (won't be collected) but still aligns well.
    let r2 = break_kmers(&revcomp(&transcript[400..500]), 25);

    // Without recovery -> the pair degrades to a left orphan.
    let no_recover = MapConfig::default();
    let m0 = map_read_pair(idx.inner(), &mut hs, &refs, &r1, &r2, &no_recover);
    assert_eq!(m0.len(), 1, "{m0:?}");
    assert_eq!(m0[0].status, MateStatus::PairedEndLeft);

    // With recovery -> the partner is rescued by alignment, yielding a pair.
    let recover = MapConfig {
        recover_orphans: true,
        ..Default::default()
    };
    let m1 = map_read_pair(idx.inner(), &mut hs, &refs, &r1, &r2, &recover);
    assert_eq!(m1.len(), 1, "{m1:?}");
    assert_eq!(
        m1[0].status,
        MateStatus::PairedEndPaired,
        "recovery should pair"
    );
}

#[test]
fn pseudoalignment_only_mode() {
    let transcript = gen_seq(600, 0x5EED);
    let (idx, _tmp) = build_index_for(&transcript);
    let mut hs = HitSearcher::new(idx.inner());
    let mut scratch = salmon_map::SketchScratch::new(idx.inner().k());

    // single-end sketch
    let read = &transcript[100..200];
    let s = map_single_read_sketch(
        idx.inner(),
        &mut hs,
        &mut scratch,
        read,
        SkippingStrategy::Permissive,
        256,
        2500,
    );
    assert!(s.iter().any(|m| m.tid == 0), "sketch single: {s:?}");

    // paired sketch -> intersection on tid 0
    let r1 = transcript[50..130].to_vec();
    let r2 = revcomp(&transcript[350..430]);
    let sp = map_read_pair_sketch(
        idx.inner(),
        &mut hs,
        &mut scratch,
        &r1,
        &r2,
        false,
        false,
        SkippingStrategy::Permissive,
        256,
        2500,
        &idx,
        false,
    );
    assert!(
        sp.iter()
            .any(|m| m.tid == 0 && m.status == MateStatus::PairedEndPaired),
        "sketch pair: {sp:?}"
    );
}

#[test]
fn unmappable_read_yields_no_candidates() {
    let tmp = tempfile::tempdir().unwrap();
    let transcript = gen_seq(300, 0x1234);
    let fasta = write_fasta(tmp.path(), "tx0", &transcript);
    let out = tmp.path().join("idx");
    let mut opts = IndexBuildOptions::new(vec![fasta], out.clone());
    opts.threads = 1;
    build(&opts).expect("build index");
    let idx = SalmonIndex::load(&out).expect("load index");
    let mut hs = HitSearcher::new(idx.inner());

    // A different random sequence should not map.
    let foreign = gen_seq(100, 0x9999);
    let cands = collect_read_mems(
        idx.inner(),
        &mut hs,
        &foreign,
        true,
        &MemCollectorConfig::default(),
    );
    assert!(
        cands.is_empty(),
        "foreign read unexpectedly mapped: {cands:?}"
    );
}
