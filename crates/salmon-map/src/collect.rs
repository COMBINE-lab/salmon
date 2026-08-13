//! MEM collection over piscem k-mer hits.
//!
//! # Where mapping starts
//!
//! This is the first stage: turn "these k-mers of the read were found in the
//! index" into "these are the candidate placements worth considering". The index
//! answers per k-mer, and answers in terms of the de Bruijn graph's unitigs, so
//! the work here is to translate unitig hits into concrete reference positions,
//! group them per transcript, and hand each group to the chainer.
//!
//! The payoff is that most of the transcriptome disappears immediately: a
//! transcript sharing no k-mer with the read is never considered again.
//!
//! This turns piscem's raw k-mer hits into the [`Mem`] anchors that the
//! [chaining](crate::chain) DP consumes, then chains them per reference and
//! orientation to produce mapping candidates. It corresponds to the front half
//! of salmon's `MemCollector` / `findChains`.
//!
//! ## Projection
//! piscem's [`HitSearcher`] returns, per queried k-mer, a [`ProjectedHits`]
//! describing where the k-mer landed on a unitig plus the list of reference
//! occurrences of that unitig. Each occurrence decodes (via the documented
//! 4-case rule in [`ProjectedHits::decode_hit`]) to a reference position and
//! orientation. A k-mer at read offset `p` therefore yields, per reference, an
//! anchor of length `k`.
//!
//! ## Orientation frame
//! Chaining requires anchors that increase colinearly in `(read, ref)`. For
//! forward hits the read offset is used directly. For reverse-complement hits
//! the read offset is reflected to `read_len - (p + k)`, which is exactly the
//! position of that segment in the reverse-complemented read; this makes RC
//! anchors increase colinearly with the reference coordinate. The resulting
//! chains are stamped with `is_fw = false` so a later stage can map them back.

use foldhash::HashMap;

use piscem_rs::index::contig_table::EntryEncoding;
use piscem_rs::index::reference_index::ReferenceIndex;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use piscem_rs::mapping::projected_hits::ProjectedHits;
use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;
use sshash_lib::dispatch_on_k;

use crate::chain::{chain_mems, ChainConfig, MemChain};
use crate::mem::Mem;

/// Configuration for MEM collection.
#[derive(Debug, Clone)]
pub struct MemCollectorConfig {
    /// how the hit searcher advances along the read
    pub skipping: SkippingStrategy,
    /// skip k-mers whose unitig occurs in more than this many references
    /// (repetitive-hit guard; salmon's `maxOccsPerHit`)
    pub max_hit_occ: usize,
    /// chaining parameters (its `seed_len` is overridden with the index `k`)
    pub chain: ChainConfig,
    /// Consensus fraction: a per-mate target is only kept (and later aligned) if
    /// its best chain score is at least `consensus_fraction * maxChainScore`
    /// across that mate's targets. This is pufferfish's pre-alignment pruning
    /// (`consensusFraction`, salmon `1 - consensusSlack`, default 0.65) — it
    /// avoids running an SW alignment for every weakly-seeded transcript. `0.0`
    /// disables the filter (keep every target).
    pub consensus_fraction: f32,
}

impl Default for MemCollectorConfig {
    fn default() -> Self {
        Self {
            skipping: SkippingStrategy::Permissive,
            max_hit_occ: 1000,
            chain: ChainConfig::default(),
            consensus_fraction: 0.65,
        }
    }
}

/// A mapping candidate: a chain of MEMs on one reference in one orientation.
#[derive(Debug, Clone)]
pub struct MappingCandidate {
    /// reference (transcript) id
    pub tid: u32,
    /// orientation of the mapping (`true` = read maps forward to the reference)
    pub is_fw: bool,
    /// the supporting MEM chain (in the reference-forward frame)
    pub chain: MemChain,
}

/// Project raw k-mer hits onto reference anchors, grouped by `(tid, is_fw)`.
///
/// `read_len` and `k` are in bases. Hits whose unitig occurs in more than
/// `max_hit_occ` references are skipped. Pure and deterministic — the unit of
/// the collection logic that does not require a live index.
pub fn project_raw_hits(
    raw_hits: &[(i32, ProjectedHits<'_>)],
    encoding: &EntryEncoding,
    read_len: i32,
    k: i32,
    max_hit_occ: usize,
) -> HashMap<(u32, bool), Vec<Mem>> {
    let mut groups: HashMap<(u32, bool), Vec<Mem>> = HashMap::default();
    for (read_pos, phit) in raw_hits {
        if phit.num_hits() > max_hit_occ {
            continue;
        }
        for entry in phit.ref_range().iter() {
            let tid = encoding.transcript_id(entry);
            let rp = phit.decode_hit(entry, encoding);
            let read_start = if rp.is_fw {
                *read_pos
            } else {
                read_len - (*read_pos + k)
            };
            groups
                .entry((tid, rp.is_fw))
                .or_default()
                .push(Mem::new(read_start, rp.pos as i32, k));
        }
    }
    groups
}

thread_local! {
    /// Reused per-thread scratch for [`candidates_from_raw_hits`]: a flat
    /// `(tid, is_fw, mem)` buffer (sorted to group by target/orientation, instead
    /// of a fresh per-read `HashMap` of per-group `Vec`s) and a per-group `Mem`
    /// buffer handed to the chainer. Eliminates the mapper's per-read map+vec
    /// allocations.
    static PROJ_SCRATCH: std::cell::RefCell<(Vec<(u32, bool, Mem)>, Vec<Mem>)> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

/// Project then chain raw hits into mapping candidates. Pure; testable with
/// hand-built hits. Groups anchors by `(tid, is_fw)` via a sort over reused
/// scratch (cheaper than a per-read hash map for the small anchor counts of a
/// single read) and chains each group. Candidate order differs from a hash-map
/// grouping, but every downstream consumer (`best_per_target`, pairing,
/// `finalize_mappings`) keys/dedups by `tid`, so results are unchanged.
pub fn candidates_from_raw_hits(
    raw_hits: &[(i32, ProjectedHits<'_>)],
    encoding: &EntryEncoding,
    read_len: i32,
    k: i32,
    cfg: &MemCollectorConfig,
) -> Vec<MappingCandidate> {
    let mut chain_cfg = cfg.chain.clone();
    chain_cfg.seed_len = k;

    PROJ_SCRATCH.with(|cell| {
        let (flat, group) = &mut *cell.borrow_mut();
        flat.clear();
        for (read_pos, phit) in raw_hits {
            if phit.num_hits() > cfg.max_hit_occ {
                continue;
            }
            for entry in phit.ref_range().iter() {
                let tid = encoding.transcript_id(entry);
                let rp = phit.decode_hit(entry, encoding);
                let read_start = if rp.is_fw {
                    *read_pos
                } else {
                    read_len - (*read_pos + k)
                };
                flat.push((tid, rp.is_fw, Mem::new(read_start, rp.pos as i32, k)));
            }
        }
        // Group by (tid, is_fw) via a sort (stable relative anchor order within a
        // group is irrelevant — the chainer re-sorts anchors itself).
        flat.sort_unstable_by_key(|&(tid, fw, _)| (tid, fw));

        let mut candidates = Vec::new();
        let mut i = 0;
        while i < flat.len() {
            let (tid, is_fw, _) = flat[i];
            group.clear();
            let mut j = i;
            while j < flat.len() && flat[j].0 == tid && flat[j].1 == is_fw {
                group.push(flat[j].2);
                j += 1;
            }
            for chain in chain_mems(group, is_fw, &chain_cfg) {
                candidates.push(MappingCandidate { tid, is_fw, chain });
            }
            i = j;
        }
        candidates
    })
}

/// Collect mapping candidates for one read against a loaded index.
///
/// `hs` is reusable per-thread scratch (`HitSearcher::new(index)`); its hit
/// buffers are cleared on each call. `is_left` selects the left/right hit
/// buffer (use `true` for single-end / first mate).
pub fn collect_read_mems<'idx>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    read: &[u8],
    is_left: bool,
    cfg: &MemCollectorConfig,
) -> Vec<MappingCandidate> {
    let k = index.k();
    let read_len = read.len() as i32;
    if read_len < k as i32 {
        return Vec::new();
    }
    let encoding = index.encoding();

    dispatch_on_k!(k, K => {
        let mut query = PiscemStreamingQuery::<K>::new(index.dict());
        hs.get_raw_hits_sketch::<K>(read, &mut query, cfg.skipping, is_left);
        let raw = if is_left { hs.left_hits() } else { hs.right_hits() };
        candidates_from_raw_hits(raw, &encoding, read_len, k as i32, cfg)
    })
}

/// Return the best (highest read-coverage) candidate per `(tid, is_fw)`,
/// discarding weaker chains to the same target/orientation. Convenience for
/// callers that want one alignment seed per target.
pub fn best_per_target(mut candidates: Vec<MappingCandidate>) -> Vec<MappingCandidate> {
    // Highest coverage first, so the first seen per key is the best. Coverage is
    // cached on the chain, so this is a plain field read per comparison.
    candidates.sort_by(|a, b| {
        b.chain
            .covered_read_bases()
            .cmp(&a.chain.covered_read_bases())
    });
    // Keep the first candidate seen per `(tid, is_fw)`, scanning the already-kept
    // prefix rather than building a hash set.
    //
    // The set was allocated and dropped on every call to deduplicate a handful of
    // candidates; at this size a linear scan over what has been kept is cheaper
    // than hashing each key, and needs no allocation and no scratch parameter.
    // The kept elements stay in their original relative order, so the
    // highest-coverage-first ordering established above survives.
    let mut kept = 0usize;
    for i in 0..candidates.len() {
        let key = (candidates[i].tid, candidates[i].is_fw);
        if candidates[..kept].iter().all(|c| (c.tid, c.is_fw) != key) {
            candidates.swap(i, kept);
            kept += 1;
        }
    }
    candidates.truncate(kept);
    candidates
}

/// Drop candidate targets whose chain score is below `fraction * maxChainScore`
/// across `candidates` — pufferfish's consensus-fraction pruning, applied
/// per-mate before the (expensive) SW alignment step. `fraction <= 0` keeps all.
/// Mirrors `MemClusterer::findOptChain`'s `bestScore < maxChainScore * cf` gate.
pub fn consensus_filter(
    mut candidates: Vec<MappingCandidate>,
    fraction: f32,
) -> Vec<MappingCandidate> {
    if fraction <= 0.0 || candidates.len() <= 1 {
        return candidates;
    }
    let max_score = candidates
        .iter()
        .map(|c| c.chain.score)
        .fold(f32::MIN, f32::max);
    if max_score <= 0.0 {
        return candidates;
    }
    let cutoff = max_score * fraction;
    candidates.retain(|c| c.chain.score >= cutoff);
    candidates
}

#[cfg(test)]
mod tests {
    use super::*;
    use piscem_rs::index::contig_table::{ContigTable, ContigTableBuilder};
    use piscem_rs::mapping::projected_hits::ProjectedHits;

    /// One contig (id 0) occurring forward on ref 0 at position 100.
    fn fw_table() -> ContigTable {
        let mut b = ContigTableBuilder::new(1, 10_000, 1);
        b.add_occurrence(0, 0, 100, true);
        b.build()
    }

    /// One contig (id 0) occurring reverse on ref 0 at position 100.
    fn rc_table() -> ContigTable {
        let mut b = ContigTableBuilder::new(1, 10_000, 1);
        b.add_occurrence(0, 0, 100, false);
        b.build()
    }

    #[test]
    /// The basic path: k-mer hits become anchors at the right reference positions
    /// and chain into one candidate.
    fn forward_kmers_project_and_chain() {
        let table = fw_table();
        let enc = table.encoding();
        let k = 31;
        let contig_len = 200;
        // Two consecutive read k-mers (positions 0 and 10), forward on contig.
        let raw: Vec<(i32, ProjectedHits)> = [0i32, 10]
            .iter()
            .map(|&p| {
                (
                    p,
                    ProjectedHits::new(
                        0,
                        p as u32,
                        true,
                        contig_len,
                        0,
                        k as u32,
                        table.contig_entries(0),
                    ),
                )
            })
            .collect();

        let cands = candidates_from_raw_hits(&raw, &enc, 100, k, &MemCollectorConfig::default());
        assert_eq!(cands.len(), 1);
        let c = &cands[0];
        assert_eq!(c.tid, 0);
        assert!(c.is_fw);
        // contig fw -> ref pos = base_pos(100) + contig_pos; read 0 -> ref 100.
        assert_eq!(c.chain.ref_start(), 100);
        // two colinear k-mers 10 apart -> covered 10 + 31 = 41 read bases.
        assert_eq!(c.chain.covered_read_bases(), 41);
    }

    #[test]
    /// Reverse-strand hits must have their read coordinates reflected, or the
    /// anchors would not be colinear and would never chain.
    fn rc_kmers_use_reflected_read_coordinate() {
        let table = rc_table();
        let enc = table.encoding();
        let k = 31;
        let read_len = 100;
        let contig_len = 200;
        // contig stored RC on the reference, k-mer FW on the contig -> is_fw=false.
        // read_pos 0 -> read_start' = 100 - (0+31) = 69
        // read_pos 10 -> read_start' = 100 - (10+31) = 59
        let raw: Vec<(i32, ProjectedHits)> = [0i32, 10]
            .iter()
            .map(|&p| {
                (
                    p,
                    ProjectedHits::new(
                        0,
                        p as u32,
                        true,
                        contig_len,
                        0,
                        k as u32,
                        table.contig_entries(0),
                    ),
                )
            })
            .collect();

        let cands =
            candidates_from_raw_hits(&raw, &enc, read_len, k, &MemCollectorConfig::default());
        assert_eq!(cands.len(), 1);
        let c = &cands[0];
        assert_eq!(c.tid, 0);
        assert!(!c.is_fw, "contig RC, k-mer FW -> reverse mapping");
        // Both anchors chain colinearly in the reflected frame.
        assert_eq!(c.chain.mems.len(), 2);
        assert_eq!(c.chain.covered_read_bases(), 41);
    }

    #[test]
    /// A k-mer occurring in a great many places says almost nothing and would cost
    /// a candidate per occurrence, so it must be skipped.
    fn repetitive_hits_are_skipped() {
        // A contig occurring in 3 references; with max_hit_occ=2 it is skipped.
        let mut b = ContigTableBuilder::new(1, 10_000, 3);
        b.add_occurrence(0, 0, 100, true);
        b.add_occurrence(0, 1, 200, true);
        b.add_occurrence(0, 2, 300, true);
        let table = b.build();
        let enc = table.encoding();
        let raw = vec![(
            0i32,
            ProjectedHits::new(0, 0, true, 200, 0, 31, table.contig_entries(0)),
        )];

        let cfg = MemCollectorConfig {
            max_hit_occ: 2,
            ..Default::default()
        };
        assert!(candidates_from_raw_hits(&raw, &enc, 100, 31, &cfg).is_empty());

        let cfg = MemCollectorConfig {
            max_hit_occ: 3,
            ..Default::default()
        };
        assert_eq!(candidates_from_raw_hits(&raw, &enc, 100, 31, &cfg).len(), 3);
    }
}
