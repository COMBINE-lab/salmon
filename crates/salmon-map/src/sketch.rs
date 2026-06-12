//! Pseudoalignment-only mapping (the optional, alignment-free path).
//!
//! This is the new capability the port adds on top of salmon: quantify directly
//! from piscem's pseudoalignment, skipping MEM chaining and alignment
//! validation entirely. It calls piscem's [`map_read`] engine to collect
//! per-transcript hits and emits them as equivalence-class members with uniform
//! weight. For read pairs it uses the intersection rule (a transcript must be
//! hit by both mates), falling back to orphans when only one mate maps.

use piscem_rs::index::reference_index::ReferenceIndex;
use piscem_rs::mapping::cache::MappingCache;
use piscem_rs::mapping::engine::map_read;
use piscem_rs::mapping::filters::PoisonState;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use piscem_rs::mapping::hits::MappingType;
use piscem_rs::mapping::merge_pairs::merge_se_mappings;
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
/// Uses piscem's bulk pair-merge ([`merge_se_mappings`]) rather than a plain
/// per-transcript intersection, so the pairing matches piscem/salmon (and
/// kallisto) semantics:
/// - both mates accepted → concordant pairs only (opposite orientation +
///   in-range fragment length), via `merge_lists`; an empty merge is *unmapped*,
///   not orphaned;
/// - exactly one mate accepted **and the other had no matching k-mers** → that
///   mate is emitted as an orphan (C++/kallisto `check_kmers_orphans` rule —
///   `MappingCache::has_matching_kmers`);
/// - one mate accepted but the other *did* have matching k-mers (just no
///   accepted target) → unmapped (the mates conflict).
pub fn map_read_pair_sketch<'idx>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    r1: &[u8],
    r2: &[u8],
    strat: SkippingStrategy,
    max_hit_occ: usize,
    max_read_occ: usize,
) -> Vec<ScoredMapping> {
    let k = index.k();
    if r1.len() < k && r2.len() < k {
        return Vec::new();
    }
    dispatch_on_k!(k, K => {
        let mut left = MappingCache::<SketchHitInfoSimple>::new(k);
        let mut right = MappingCache::<SketchHitInfoSimple>::new(k);
        let mut out = MappingCache::<SketchHitInfoSimple>::new(k);
        for c in [&mut left, &mut right, &mut out] {
            c.max_hit_occ = max_hit_occ;
            c.max_read_occ = max_read_occ;
        }
        let mut poison = PoisonState::new(index.poison_table());
        if r1.len() >= k {
            let mut q = PiscemStreamingQuery::<K>::new(index.dict());
            map_read::<K, SketchHitInfoSimple, _, _>(r1, &mut left, hs, &mut q, index, &mut poison, strat);
        }
        if r2.len() >= k {
            let mut q = PiscemStreamingQuery::<K>::new(index.dict());
            map_read::<K, SketchHitInfoSimple, _, _>(r2, &mut right, hs, &mut q, index, &mut poison, strat);
        }
        merge_se_mappings(&mut left, &mut right, r1.len() as i32, r2.len() as i32, &mut out);

        let status = match out.map_type {
            MappingType::MappedPair => Some(MateStatus::PairedEndPaired),
            MappingType::MappedFirstOrphan => Some(MateStatus::PairedEndLeft),
            MappingType::MappedSecondOrphan => Some(MateStatus::PairedEndRight),
            _ => None,
        };
        match status {
            None => Vec::new(),
            Some(status) => out
                .accepted_hits
                .iter()
                .map(|h| ScoredMapping {
                    tid: h.tid,
                    is_fw: h.is_fw,
                    status,
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
                .collect::<Vec<_>>(),
        }
    })
}
