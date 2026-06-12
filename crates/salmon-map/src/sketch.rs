//! Pseudoalignment-only mapping (the optional, alignment-free path).
//!
//! This is the new capability the port adds on top of salmon: quantify directly
//! from piscem's pseudoalignment, skipping MEM chaining and alignment
//! validation entirely. It calls piscem's [`map_read`] engine to collect
//! per-transcript hits and emits them as equivalence-class members with uniform
//! weight. For read pairs it uses the intersection rule (a transcript must be
//! hit by both mates), falling back to orphans when only one mate maps.

use std::collections::HashMap;

use piscem_rs::index::reference_index::ReferenceIndex;
use piscem_rs::mapping::cache::MappingCache;
use piscem_rs::mapping::engine::map_read;
use piscem_rs::mapping::filters::PoisonState;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use piscem_rs::mapping::sketch_hit_simple::SketchHitInfoSimple;
use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;
use salmon_core::MateStatus;
use sshash_lib::dispatch_on_k;

use crate::score::ScoredMapping;

/// Pseudoalign a single read; returns one uniform-weight member per hit transcript.
pub fn map_single_read_sketch<'idx>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    read: &[u8],
    strat: SkippingStrategy,
    max_hit_occ: usize,
    max_read_occ: usize,
) -> Vec<ScoredMapping> {
    let k = index.k();
    if read.len() < k {
        return Vec::new();
    }
    dispatch_on_k!(k, K => {
        let mut query = PiscemStreamingQuery::<K>::new(index.dict());
        let mut cache = MappingCache::<SketchHitInfoSimple>::new(k);
        // Honor the same repetitive-hit and multimapping caps as selective
        // alignment (--maxOccsPerHit / --maxReadOcc) rather than piscem's
        // built-in defaults (256 / 2500), so both modes filter consistently.
        cache.max_hit_occ = max_hit_occ;
        cache.max_read_occ = max_read_occ;
        let mut poison = PoisonState::new(index.poison_table());
        map_read::<K, SketchHitInfoSimple, _, _>(
            read, &mut cache, hs, &mut query, index, &mut poison, strat,
        );
        cache
            .accepted_hits
            .iter()
            .map(|h| ScoredMapping {
                tid: h.tid,
                is_fw: h.is_fw,
                status: MateStatus::SingleEnd,
                score: h.score as i32,
                fragment_len: 0,
                weight: 1.0,
                ref_pos: 0,
                fw_pos: -1,
                rc_pos: -1,
                format: None,
                r1_pos: -1,
                r2_pos: -1,
                r2_fw: false,
                r1_score: h.score as i32,
            })
            .collect::<Vec<_>>()
    })
}

/// Pseudoalign a read pair using the intersection rule.
///
/// If both mates map, transcripts hit by both are emitted as pairs. If only one
/// mate maps, its hits are emitted as orphans. An empty intersection (both map,
/// but to disjoint transcripts) yields no mapping.
pub fn map_read_pair_sketch<'idx>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    r1: &[u8],
    r2: &[u8],
    strat: SkippingStrategy,
    max_hit_occ: usize,
    max_read_occ: usize,
) -> Vec<ScoredMapping> {
    let left = map_single_read_sketch(index, hs, r1, strat, max_hit_occ, max_read_occ);
    let right = map_single_read_sketch(index, hs, r2, strat, max_hit_occ, max_read_occ);

    match (left.is_empty(), right.is_empty()) {
        (true, true) => Vec::new(),
        (true, false) => orphans(right, MateStatus::PairedEndRight),
        (false, true) => orphans(left, MateStatus::PairedEndLeft),
        (false, false) => {
            let right_score: HashMap<u32, i32> = right.iter().map(|m| (m.tid, m.score)).collect();
            left.into_iter()
                .filter_map(|m| {
                    right_score.get(&m.tid).map(|rs| ScoredMapping {
                        tid: m.tid,
                        is_fw: m.is_fw,
                        status: MateStatus::PairedEndPaired,
                        score: m.score + rs,
                        fragment_len: 0,
                        weight: 1.0,
                        ref_pos: 0,
                        fw_pos: -1,
                        rc_pos: -1,
                        format: None,
                        r1_pos: -1,
                        r2_pos: -1,
                        r2_fw: false,
                        r1_score: m.score,
                    })
                })
                .collect()
        }
    }
}

fn orphans(mut m: Vec<ScoredMapping>, status: MateStatus) -> Vec<ScoredMapping> {
    for x in &mut m {
        x.status = status;
    }
    m
}
