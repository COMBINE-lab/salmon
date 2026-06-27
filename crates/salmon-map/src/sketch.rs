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
use salmon_core::{MateStatus, RefProvider};
use sshash_lib::dispatch_on_k;

use crate::score::ScoredMapping;

/// Reusable per-thread scratch for the sketch path: the three `MappingCache`s
/// that `map_read` (left/right mate) and `merge_se_mappings` (merge output)
/// fill. Each `MappingCache::new` allocates an `AHashMap` (capacity 256) plus a
/// `Vec`; building three on every read pair is pure churn (~120 ns/pair).
/// `map_read` and `merge_se_mappings` `clear()` their caches internally, so
/// holding one `SketchScratch` per worker thread and reusing it across reads is
/// safe and removes those allocations from the hot loop — mirroring piscem's own
/// per-thread `CommonThreadState`.
pub struct SketchScratch {
    left: MappingCache<SketchHitInfoSimple>,
    right: MappingCache<SketchHitInfoSimple>,
    out: MappingCache<SketchHitInfoSimple>,
}

impl SketchScratch {
    pub fn new(k: usize) -> Self {
        Self {
            left: MappingCache::new(k),
            right: MappingCache::new(k),
            out: MappingCache::new(k),
        }
    }
}

/// Pseudoalign a single read; returns one uniform-weight member per hit transcript.
pub fn map_single_read_sketch<'idx>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    scratch: &mut SketchScratch,
    read: &[u8],
    strat: SkippingStrategy,
    max_hit_occ: usize,
    max_read_occ: usize,
) -> Vec<ScoredMapping> {
    let mut result = Vec::new();
    map_single_read_sketch_into(
        &mut result,
        index,
        hs,
        scratch,
        read,
        strat,
        max_hit_occ,
        max_read_occ,
    );
    result
}

/// [`map_single_read_sketch`] writing into a caller-provided buffer (cleared
/// first), for per-thread `Vec` reuse across reads.
#[allow(clippy::too_many_arguments)]
pub fn map_single_read_sketch_into<'idx>(
    result: &mut Vec<ScoredMapping>,
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    scratch: &mut SketchScratch,
    read: &[u8],
    strat: SkippingStrategy,
    max_hit_occ: usize,
    max_read_occ: usize,
) {
    result.clear();
    let k = index.k();
    if read.len() < k {
        return;
    }
    dispatch_on_k!(k, K => {
        let mut query = PiscemStreamingQuery::<K>::new(index.dict());
        let cache = &mut scratch.left;
        // Honor the same repetitive-hit and multimapping caps as selective
        // alignment (--maxOccsPerHit / --maxReadOcc) rather than piscem's
        // built-in defaults (256 / 2500), so both modes filter consistently.
        cache.max_hit_occ = max_hit_occ;
        cache.max_read_occ = max_read_occ;
        let mut poison = PoisonState::new(index.poison_table());
        map_read::<K, SketchHitInfoSimple, _, _>(
            read, cache, hs, &mut query, index, &mut poison, strat,
        );
        let rl = read.len() as i32;
        result.extend(cache
            .accepted_hits
            .iter()
            .map(|h| ScoredMapping {
                tid: h.tid,
                is_fw: h.is_fw,
                status: MateStatus::SingleEnd,
                score: h.score as i32,
                fragment_len: 0,
                read_len: rl,
                weight: 1.0,
                // piscem's SimpleHit.pos is the read's approximate leftmost
                // reference coordinate. ref_pos is the orientation-aware 5' end
                // (leftmost for fw, rightmost for rc) used as the bias context
                // anchor, mirroring selective alignment's single-end case.
                ref_pos: if h.is_fw { h.pos.max(0) } else { h.pos.max(0) + rl },
                fw_pos: if h.is_fw { h.pos.max(0) } else { -1 },
                rc_pos: if h.is_fw { -1 } else { h.pos.max(0) },
                format: None,
                r1_pos: h.pos.max(0),
                r2_pos: -1,
                r2_fw: false,
                r1_score: h.score as i32,
            }));
    })
}

/// Pseudoalign a read pair using the intersection rule.
///
/// Uses piscem's bulk pair-merge ([`merge_se_mappings`]) rather than a plain
/// per-transcript intersection, so the pairing matches piscem/salmon semantics:
/// - both mates accepted → concordant pairs only (opposite orientation +
///   in-range fragment length), via `merge_lists`; an empty merge is *unmapped*,
///   not orphaned;
/// - exactly one mate has an accepted target → that mate is emitted as an
///   orphan.
///
/// Orphan gate (`strict_orphan`):
/// - `false` (default): orphan whenever the *other* mate has no accepted target,
///   regardless of whether it had matching k-mers. This is closer to selective
///   alignment (which recovers orphans) and improves agreement with SA mode —
///   a read whose mate's k-mers fail consensus (e.g. an isolated repeat-element
///   k-mer) is effectively unmappable, so the confident mate should still count.
/// - `true`: only orphan when the other mate had **no matching k-mers at all**
///   (piscem's `check_kmers_orphans` / `MappingCache::has_matching_kmers` rule).
///   More conservative; mates that have k-mers but no consensus target leave the
///   pair unmapped. Selected by `--sketchStrictOrphans`.
#[allow(clippy::too_many_arguments)]
pub fn map_read_pair_sketch<'idx, R: RefProvider>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    scratch: &mut SketchScratch,
    r1: &[u8],
    r2: &[u8],
    strict_orphan: bool,
    allow_dovetail: bool,
    strat: SkippingStrategy,
    max_hit_occ: usize,
    max_read_occ: usize,
    refs: &R,
    allow_decoy_orphans: bool,
) -> Vec<ScoredMapping> {
    let mut result = Vec::new();
    map_read_pair_sketch_into(
        &mut result,
        index,
        hs,
        scratch,
        r1,
        r2,
        strict_orphan,
        allow_dovetail,
        strat,
        max_hit_occ,
        max_read_occ,
        refs,
        allow_decoy_orphans,
    );
    result
}

/// [`map_read_pair_sketch`] writing into a caller-provided buffer (cleared
/// first), so a per-thread `Vec` can be reused across reads.
#[allow(clippy::too_many_arguments)]
pub fn map_read_pair_sketch_into<'idx, R: RefProvider>(
    result: &mut Vec<ScoredMapping>,
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    scratch: &mut SketchScratch,
    r1: &[u8],
    r2: &[u8],
    strict_orphan: bool,
    allow_dovetail: bool,
    strat: SkippingStrategy,
    max_hit_occ: usize,
    max_read_occ: usize,
    refs: &R,
    allow_decoy_orphans: bool,
) {
    result.clear();
    let k = index.k();
    if r1.len() < k && r2.len() < k {
        return;
    }
    // Fragment-length lower bound for the pair merge. Default -32 (piscem's
    // small-dovetail tolerance). With `--allowDovetail`, admit genuine dovetails
    // (insert < read length): the k-mer-anchor frag_len ≈ insert − read_len, so
    // a true dovetail bottoms out at -read_len; we use that as the floor, which
    // recovers overlapping short fragments without admitting the large-negative
    // multimapping artifacts.
    let min_frag_len = if allow_dovetail {
        -(r1.len().max(r2.len()) as i32)
    } else {
        -32
    };
    let SketchScratch { left, right, out } = scratch;
    dispatch_on_k!(k, K => {
        for c in [&mut *left, &mut *right, &mut *out] {
            c.max_hit_occ = max_hit_occ;
            c.max_read_occ = max_read_occ;
            // merge_se_mappings reads the pairing policy from cache_left:
            // the dovetail floor (`--allowDovetail`) and the orphan rule
            // (relaxed by default; `--sketchStrictOrphans` restores piscem's
            // strict no-matching-k-mers gate).
            c.min_frag_len = min_frag_len;
            c.relaxed_orphans = !strict_orphan;
        }
        let mut poison = PoisonState::new(index.poison_table());
        if r1.len() >= k {
            let mut q = PiscemStreamingQuery::<K>::new(index.dict());
            map_read::<K, SketchHitInfoSimple, _, _>(r1, left, hs, &mut q, index, &mut poison, strat);
        }
        if r2.len() >= k {
            let mut q = PiscemStreamingQuery::<K>::new(index.dict());
            map_read::<K, SketchHitInfoSimple, _, _>(r2, right, hs, &mut q, index, &mut poison, strat);
        }
        merge_se_mappings(left, right, r1.len() as i32, r2.len() as i32, out);

        // Decoy-orphan rescue (`--allowDecoyOrphans`). Sketch pairing is same-tid
        // only (`merge_lists`), so a fragment whose mates map to a transcript and
        // a decoy respectively forms no concordant pair and `merge_se_mappings`
        // returns Unmapped — the SA-mode "a decoy dominates the fragment's other
        // placement" scenario, but here it is silently dropped rather than kept as
        // a transcript orphan. When the flag is set, recover the transcript mate as
        // an orphan, but only when the *other* mate's hits are entirely decoys (so
        // this never changes the discordant transcript-vs-transcript case, which
        // stays unmapped as in piscem). The retained hits are transcript-only, so
        // the downstream `filter_sketch_decoys` is a no-op on them.
        if allow_decoy_orphans && out.accepted_hits.is_empty() {
            let left_has_txp = left.accepted_hits.iter().any(|h| !refs.is_decoy(h.tid));
            let right_has_txp = right.accepted_hits.iter().any(|h| !refs.is_decoy(h.tid));
            let left_all_decoy = !left.accepted_hits.is_empty() && !left_has_txp;
            let right_all_decoy = !right.accepted_hits.is_empty() && !right_has_txp;
            if left_has_txp && right_all_decoy {
                out.accepted_hits
                    .extend(left.accepted_hits.iter().filter(|h| !refs.is_decoy(h.tid)).cloned());
                out.map_type = MappingType::MappedFirstOrphan;
            } else if right_has_txp && left_all_decoy {
                out.accepted_hits
                    .extend(right.accepted_hits.iter().filter(|h| !refs.is_decoy(h.tid)).cloned());
                out.map_type = MappingType::MappedSecondOrphan;
            }
        }

        let status = match out.map_type {
            MappingType::MappedPair => Some(MateStatus::PairedEndPaired),
            MappingType::MappedFirstOrphan => Some(MateStatus::PairedEndLeft),
            MappingType::MappedSecondOrphan => Some(MateStatus::PairedEndRight),
            _ => None,
        };
        let (r1l, r2l) = (r1.len() as i32, r2.len() as i32);
        match status {
            None => {}
            Some(status) => result.extend(
                out.accepted_hits
                .iter()
                .map(|h| {
                    // piscem's SimpleHit gives read1's leftmost in `pos` and
                    // read2's in `mate_pos` (for a pair); for an orphan only the
                    // surviving mate's `pos` is set. Derive the same position
                    // fields selective alignment fills so bias context, the
                    // ambiguous-fragment FLD term, and SAM output all work.
                    let frag = h.frag_len();
                    let (ref_pos, fw_pos, rc_pos, read_len, r1_pos, r2_pos, r2_fw) = match status
                    {
                        MateStatus::PairedEndPaired => {
                            let (p1, p2) = (h.pos.max(0), h.mate_pos.max(0));
                            let frag_start = p1.min(p2);
                            // 5' START -> fw positional model, 3' END -> rc model
                            // (matches the SA paired case).
                            (frag_start, frag_start, frag_start + frag - 1, 0, p1, p2, h.mate_is_fw)
                        }
                        // Orphan: only one mate placed; `pos` is its leftmost,
                        // `is_fw` its orientation. ref_pos is the 5' end; fw_pos/
                        // rc_pos carry the leftmost on the mate's own strand (what
                        // the ambiguous-fragment FLD term reads).
                        MateStatus::PairedEndRight => {
                            let p = h.pos.max(0);
                            let r5 = if h.is_fw { p } else { p + r2l };
                            (r5, if h.is_fw { p } else { -1 }, if h.is_fw { -1 } else { p }, r2l, -1, p, false)
                        }
                        _ => {
                            // PairedEndLeft (and any fallback): orphan = read1.
                            let p = h.pos.max(0);
                            let r5 = if h.is_fw { p } else { p + r1l };
                            (r5, if h.is_fw { p } else { -1 }, if h.is_fw { -1 } else { p }, r1l, p, -1, false)
                        }
                    };
                    ScoredMapping {
                        tid: h.tid,
                        is_fw: h.is_fw,
                        status,
                        score: h.score as i32,
                        // For concordant pairs `merge_lists` stored the template
                        // length in `fragment_length`; feed it through so the FLD
                        // is learned from the data (orphans report 0).
                        fragment_len: frag,
                        read_len,
                        weight: 1.0,
                        ref_pos,
                        fw_pos,
                        rc_pos,
                        format: None,
                        r1_pos,
                        r2_pos,
                        r2_fw,
                        r1_score: h.score as i32,
                    }
                }),
            ),
        }
    })
}
