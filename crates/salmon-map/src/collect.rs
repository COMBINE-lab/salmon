//! MEM collection over piscem k-mer hits.
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

use std::collections::HashMap;

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
}

impl Default for MemCollectorConfig {
    fn default() -> Self {
        Self {
            skipping: SkippingStrategy::Permissive,
            max_hit_occ: 1000,
            chain: ChainConfig::default(),
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
    let mut groups: HashMap<(u32, bool), Vec<Mem>> = HashMap::new();
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

/// Project then chain raw hits into mapping candidates. Pure; testable with
/// hand-built hits.
pub fn candidates_from_raw_hits(
    raw_hits: &[(i32, ProjectedHits<'_>)],
    encoding: &EntryEncoding,
    read_len: i32,
    k: i32,
    cfg: &MemCollectorConfig,
) -> Vec<MappingCandidate> {
    let groups = project_raw_hits(raw_hits, encoding, read_len, k, cfg.max_hit_occ);
    let mut chain_cfg = cfg.chain.clone();
    chain_cfg.seed_len = k;

    let mut candidates = Vec::new();
    for ((tid, is_fw), mems) in groups {
        for chain in chain_mems(&mems, is_fw, &chain_cfg) {
            candidates.push(MappingCandidate { tid, is_fw, chain });
        }
    }
    candidates
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
    // Highest coverage first, so the first seen per key is the best.
    candidates.sort_by(|a, b| {
        b.chain
            .covered_read_bases()
            .cmp(&a.chain.covered_read_bases())
    });
    let mut seen = std::collections::HashSet::new();
    candidates.retain(|c| seen.insert((c.tid, c.is_fw)));
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
                    ProjectedHits::new(0, p as u32, true, contig_len, 0, k as u32, table.contig_entries(0)),
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
                    ProjectedHits::new(0, p as u32, true, contig_len, 0, k as u32, table.contig_entries(0)),
                )
            })
            .collect();

        let cands = candidates_from_raw_hits(&raw, &enc, read_len, k, &MemCollectorConfig::default());
        assert_eq!(cands.len(), 1);
        let c = &cands[0];
        assert_eq!(c.tid, 0);
        assert!(!c.is_fw, "contig RC, k-mer FW -> reverse mapping");
        // Both anchors chain colinearly in the reflected frame.
        assert_eq!(c.chain.mems.len(), 2);
        assert_eq!(c.chain.covered_read_bases(), 41);
    }

    #[test]
    fn repetitive_hits_are_skipped() {
        // A contig occurring in 3 references; with max_hit_occ=2 it is skipped.
        let mut b = ContigTableBuilder::new(1, 10_000, 3);
        b.add_occurrence(0, 0, 100, true);
        b.add_occurrence(0, 1, 200, true);
        b.add_occurrence(0, 2, 300, true);
        let table = b.build();
        let enc = table.encoding();
        let raw = vec![(0i32, ProjectedHits::new(0, 0, true, 200, 0, 31, table.contig_entries(0)))];

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
