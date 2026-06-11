//! The assembled selective-alignment mapper.
//!
//! Ties the pieces together for a read or read pair: collect MEM chains
//! ([`collect`](crate::collect)), pair them ([`pair`](crate::pair)), validate by
//! alignment ([`align`](crate::align)), optionally recover orphan mates, then
//! collapse/weight into equivalence-class members ([`score`](crate::score)).
//! The output [`ScoredMapping`]s are the read's contribution to the equivalence
//! class consumed by the inference step.

use piscem_rs::index::reference_index::ReferenceIndex;
use piscem_rs::mapping::hit_searcher::HitSearcher;
use salmon_core::{LibraryFormat, MateStatus, ReadOrientation, ReadStrandedness, ReadType};

use crate::align::{align_chain, align_in_window, revcomp, AlignConfig};
use crate::collect::{
    best_per_target, collect_read_mems, consensus_filter, MappingCandidate, MemCollectorConfig,
};
use crate::extend::{collect_read_true_unimems, collect_read_unimems};

/// How a read's seeds are turned into chaining anchors.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SeedMode {
    /// Sparse fixed-`k` k-mer anchors straight from piscem's skipping query
    /// (one length-`k` anchor per unitig transition). The default.
    #[default]
    Sparse,
    /// Reference MEMs: each seed extended against the reference transcript,
    /// crossing unitig boundaries ([`crate::extend::collect_read_unimems`]).
    RefMem,
    /// True uni-MEMs: extension clamped to each seed's unitig, reproducing
    /// pufferfish's `expandHitEfficient` ([`crate::extend::collect_read_true_unimems`]).
    UniMem,
}
use crate::pair::{join_reads_and_filter, PairingConfig};
use crate::score::{finalize_mappings_counted, RawMapping, ScoreConfig, ScoredMapping};
use salmon_core::RefProvider;

/// Full mapper configuration.
#[derive(Debug, Clone, Default)]
pub struct MapConfig {
    pub collect: MemCollectorConfig,
    pub align: AlignConfig,
    pub pair: PairingConfig,
    pub score: ScoreConfig,
    /// attempt to rescue an unmapped mate near its mapped partner
    pub recover_orphans: bool,
    /// how seeds are extended into chaining anchors (additive/experimental;
    /// see the [`extend`](crate::extend) module docs).
    pub seed_mode: SeedMode,
}

/// Per-fragment selective-alignment statistics for the meta_info counters,
/// stashed in a per-thread slot by [`map_read_pair`]/[`map_single_read`] and
/// read by the caller via [`take_last_map_stats`]. (A thread-local avoids
/// threading an out-parameter through every mapping entry point and its tests.)
#[derive(Default, Clone, Copy)]
pub struct MapStats {
    /// the fragment produced at least one chained candidate placement
    pub had_candidates: bool,
    /// dropped because its best placement was a decoy
    pub decoy_dominated: bool,
    /// mates dovetailed (each extends past the other's start)
    pub dovetail: bool,
    /// candidate alignments rejected for scoring below the acceptance threshold
    pub alns_below_threshold: u32,
}

thread_local! {
    static LAST_MAP_STATS: std::cell::Cell<MapStats> = const { std::cell::Cell::new(MapStats {
        had_candidates: false,
        decoy_dominated: false,
        dovetail: false,
        alns_below_threshold: 0,
    }) };
}

#[inline]
fn set_last_map_stats(s: MapStats) {
    LAST_MAP_STATS.with(|c| c.set(s));
}

/// Take (and reset) the [`MapStats`] of the most recent `map_*` call on this thread.
#[inline]
pub fn take_last_map_stats() -> MapStats {
    LAST_MAP_STATS.with(|c| c.take())
}

/// Collect mapping candidates for one read, using either the sparse-k-mer path
/// or the extended uni-MEM path per [`MapConfig::use_unimems`].
fn collect_candidates<'idx, R: RefProvider>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    refs: &R,
    read: &[u8],
    is_left: bool,
    cfg: &MapConfig,
) -> Vec<MappingCandidate> {
    match cfg.seed_mode {
        SeedMode::Sparse => collect_read_mems(index, hs, read, is_left, &cfg.collect),
        SeedMode::RefMem => collect_read_unimems(index, hs, refs, read, is_left, &cfg.collect),
        SeedMode::UniMem => collect_read_true_unimems(index, hs, refs, read, is_left, &cfg.collect),
    }
}

/// Map a single-end read to weighted equivalence-class members.
pub fn map_single_read<'idx, R: RefProvider>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    refs: &R,
    read: &[u8],
    cfg: &MapConfig,
) -> Vec<ScoredMapping> {
    let cands = consensus_filter(
        best_per_target(collect_candidates(index, hs, refs, read, true, cfg)),
        cfg.collect.consensus_fraction,
    );
    let had_candidates = !cands.is_empty();
    let mut below = 0u32;
    let mut raw = Vec::with_capacity(cands.len());
    for c in cands {
        if let Some(aln) = align_chain(read, refs.ref_seq(c.tid), &c.chain, &cfg.align) {
            if aln.valid {
                // single-end observed strandedness: sense if forward, else antisense
                let strand = if c.is_fw {
                    ReadStrandedness::S
                } else {
                    ReadStrandedness::A
                };
                raw.push(RawMapping {
                    tid: c.tid,
                    is_fw: c.is_fw,
                    status: MateStatus::SingleEnd,
                    score: aln.score,
                    fragment_len: 0,
                    is_decoy: refs.is_decoy(c.tid),
                    // 5' position: leftmost for a forward read, rightmost for a
                    // reverse read (orientation-aware sequence-bias context).
                    ref_pos: if c.is_fw {
                        c.chain.ref_start()
                    } else {
                        c.chain.ref_end()
                    },
                    // positional bias: the read's leftmost reference coordinate,
                    // attributed to its own strand (salmon's SINGLE_END case).
                    fw_pos: if c.is_fw { c.chain.ref_start() } else { -1 },
                    rc_pos: if c.is_fw { -1 } else { c.chain.ref_start() },
                    format: Some(LibraryFormat::new(
                        ReadType::SingleEnd,
                        ReadOrientation::None,
                        strand,
                    )),
                    r1_pos: c.chain.ref_start(),
                    r2_pos: -1,
                    r2_fw: false,
                    r1_score: aln.score,
                });
            } else {
                below += 1;
            }
        } else {
            below += 1;
        }
    }
    let (maps, decoy_dominated, below_final) = finalize_mappings_counted(raw, &cfg.score);
    set_last_map_stats(MapStats {
        had_candidates,
        decoy_dominated,
        dovetail: false,
        alns_below_threshold: below + below_final,
    });
    maps
}

/// Map a read pair to weighted equivalence-class members.
pub fn map_read_pair<'idx, R: RefProvider>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    refs: &R,
    r1: &[u8],
    r2: &[u8],
    cfg: &MapConfig,
) -> Vec<ScoredMapping> {
    let cf = cfg.collect.consensus_fraction;
    let left = consensus_filter(best_per_target(collect_candidates(index, hs, refs, r1, true, cfg)), cf);
    let right = consensus_filter(best_per_target(collect_candidates(index, hs, refs, r2, true, cfg)), cf);
    let joints = join_reads_and_filter(left, right, &cfg.pair);

    let had_candidates = !joints.is_empty();
    let mut below = 0u32;
    let mut dovetail = false;
    let mut raw = Vec::new();
    for j in joints {
        match j.status {
            MateStatus::PairedEndPaired => {
                let l = j.left.as_ref().unwrap();
                let r = j.right.as_ref().unwrap();
                // Dovetail: opposite-strand mates where the downstream (reverse)
                // mate actually starts upstream of the forward mate — they extend
                // past each other's start (salmon's `num_dovetail_fragments`).
                if l.is_fw != r.is_fw {
                    let (fw, rc) = if l.is_fw { (l, r) } else { (r, l) };
                    if rc.chain.ref_start() < fw.chain.ref_start() {
                        dovetail = true;
                    }
                }
                let refseq = refs.ref_seq(j.tid);
                let al = align_chain(r1, refseq, &l.chain, &cfg.align);
                let ar = align_chain(r2, refseq, &r.chain, &cfg.align);
                if let (Some(al), Some(ar)) = (al, ar) {
                    if al.valid && ar.valid {
                        // positional bias (salmon's PAIRED_END_PAIRED case): only
                        // opposite-strand pairs contribute. The 5' model is indexed
                        // by the fragment 5' START and the 3' model by the fragment
                        // 3' END, matching how the expected pos5/pos3 models are
                        // built (by fragStartPos as start / as end). NOTE: this uses
                        // the fragment 3' end, not the reverse mate's leftmost (which
                        // is 3'end - readLen) — fixing a coordinate mismatch between
                        // the observed and expected 3' positional models.
                        let frag_start = l.chain.ref_start().min(r.chain.ref_start());
                        let (fw_pos, rc_pos) = if l.is_fw != r.is_fw {
                            (frag_start, frag_start + j.fragment_len - 1)
                        } else {
                            (-1, -1)
                        };
                        raw.push(RawMapping {
                            tid: j.tid,
                            is_fw: l.is_fw,
                            status: MateStatus::PairedEndPaired,
                            score: al.score + ar.score,
                            fragment_len: j.fragment_len,
                            is_decoy: refs.is_decoy(j.tid),
                            ref_pos: l.chain.ref_start().min(r.chain.ref_start()),
                            fw_pos,
                            rc_pos,
                            format: Some(j.format),
                            r1_pos: l.chain.ref_start(),
                            r2_pos: r.chain.ref_start(),
                            r2_fw: r.is_fw,
                            r1_score: al.score,
                        });
                    } else if al.valid {
                        // The concordant pair was rejected because the *other* mate's
                        // alignment fell below threshold (or the pairing was an
                        // invalid orientation). Rescue the valid mate as an orphan —
                        // salmon emits an `m1`/`m2` orphan here rather than dropping
                        // the whole fragment. Forming the (failed) pair had otherwise
                        // suppressed this mate's orphan in `join_reads_and_filter`.
                        raw.push(orphan_raw(j.tid, l, al.score, true, refs.is_decoy(j.tid)));
                        below += 1;
                    } else if ar.valid {
                        raw.push(orphan_raw(j.tid, r, ar.score, false, refs.is_decoy(j.tid)));
                        below += 1;
                    } else {
                        below += 1;
                    }
                } else {
                    below += 1;
                }
            }
            MateStatus::PairedEndLeft => {
                let anchor = j.left.as_ref().unwrap();
                push_orphan_or_recovered(&mut raw, index_seq(refs, j.tid), r1, r2, anchor, true, refs, cfg, j.tid);
            }
            MateStatus::PairedEndRight => {
                let anchor = j.right.as_ref().unwrap();
                push_orphan_or_recovered(&mut raw, index_seq(refs, j.tid), r2, r1, anchor, false, refs, cfg, j.tid);
            }
            MateStatus::SingleEnd => {}
        }
    }
    // Orphans are a *fallback*: if this fragment has any concordant (proper-pair)
    // mapping, discard all orphan mappings. A lone mate matching a paralog is weak
    // evidence and would spuriously enlarge the equivalence class / leak count mass
    // to the wrong transcript; the concordant mappings are the trustworthy signal.
    // Orphans are kept only when the fragment has *no* concordant mapping at all.
    if raw.iter().any(|m| matches!(m.status, MateStatus::PairedEndPaired)) {
        raw.retain(|m| matches!(m.status, MateStatus::PairedEndPaired));
    }
    let (maps, decoy_dominated, below_final) = finalize_mappings_counted(raw, &cfg.score);
    set_last_map_stats(MapStats {
        had_candidates,
        decoy_dominated,
        dovetail,
        alns_below_threshold: below + below_final,
    });
    maps
}

#[inline]
fn index_seq<R: RefProvider>(refs: &R, tid: u32) -> &[u8] {
    refs.ref_seq(tid)
}

/// Build an orphan [`RawMapping`] for a single mate's validated chain (used when
/// a concordant pair is rejected but one mate aligns well — the salmon `m1`/`m2`
/// orphan-rescue). `is_left` selects read1 vs read2.
fn orphan_raw(
    tid: u32,
    c: &MappingCandidate,
    score: i32,
    is_left: bool,
    is_decoy: bool,
) -> RawMapping {
    let start = c.chain.ref_start();
    RawMapping {
        tid,
        is_fw: c.is_fw,
        status: if is_left {
            MateStatus::PairedEndLeft
        } else {
            MateStatus::PairedEndRight
        },
        score,
        fragment_len: 0,
        is_decoy,
        ref_pos: start,
        fw_pos: if c.is_fw { start } else { -1 },
        rc_pos: if c.is_fw { -1 } else { start },
        format: None,
        r1_pos: if is_left { start } else { -1 },
        r2_pos: if is_left { -1 } else { start },
        r2_fw: false,
        r1_score: score,
    }
}

/// Diagnostic detail for the best mapping of a single read.
#[derive(Debug, Clone)]
pub struct DebugMapping {
    pub tid: u32,
    pub is_fw: bool,
    pub ref_pos: i32,
    /// read bases covered by exact MEM seeds in the chain
    pub chain_cov: i32,
    pub read_len: usize,
    /// full-read alignment score (current backend)
    pub full_score: i32,
    /// perfect (all-match) score for this read
    pub perfect: i32,
    /// number of MEMs in the chain
    pub num_mems: usize,
}

/// Return diagnostic detail for the best (highest exact-seed coverage) mapping
/// of a single-end read, or `None` if it doesn't map. Used to characterize the
/// reads where our mapper and salmon disagree.
pub fn debug_best_mapping<'idx, R: RefProvider>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    refs: &R,
    read: &[u8],
    cfg: &MapConfig,
) -> Option<DebugMapping> {
    let cands = consensus_filter(
        best_per_target(collect_candidates(index, hs, refs, read, true, cfg)),
        cfg.collect.consensus_fraction,
    );
    let best = cands
        .iter()
        .max_by_key(|c| c.chain.covered_read_bases())?;
    let aln = align_chain(read, refs.ref_seq(best.tid), &best.chain, &cfg.align)?;
    Some(DebugMapping {
        tid: best.tid,
        is_fw: best.is_fw,
        ref_pos: best.chain.ref_start(),
        chain_cov: best.chain.covered_read_bases(),
        read_len: read.len(),
        full_score: aln.score,
        perfect: crate::align::perfect_score(read.len(), &cfg.align),
        num_mems: best.chain.mems.len(),
    })
}

/// Align the mapped (anchor) mate; if it validates, optionally try to recover
/// the partner mate near it and emit a pair, otherwise emit an orphan.
///
/// `anchor_read` is the mate that produced the chain (`anchor`), `partner_read`
/// is the unmapped mate. `anchor_is_left` records which physical mate the
/// anchor is, so the emitted orphan status is correct.
#[allow(clippy::too_many_arguments)]
fn push_orphan_or_recovered<R: RefProvider>(
    raw: &mut Vec<RawMapping>,
    refseq: &[u8],
    anchor_read: &[u8],
    partner_read: &[u8],
    anchor: &MappingCandidate,
    anchor_is_left: bool,
    refs: &R,
    cfg: &MapConfig,
    tid: u32,
) {
    let Some(anchor_aln) = align_chain(anchor_read, refseq, &anchor.chain, &cfg.align) else {
        return;
    };
    if !anchor_aln.valid {
        return;
    }
    let is_decoy = refs.is_decoy(tid);

    if cfg.recover_orphans {
        if let Some((partner_score, frag_len)) = recover_mate(partner_read, refseq, anchor, cfg) {
            raw.push(RawMapping {
                tid,
                is_fw: anchor.is_fw,
                status: MateStatus::PairedEndPaired,
                score: anchor_aln.score + partner_score,
                fragment_len: frag_len,
                is_decoy,
                ref_pos: anchor.chain.ref_start(),
                // orientation not re-derived for recovered pairs; attribute the
                // anchor's position to its own strand for positional bias.
                fw_pos: if anchor.is_fw { anchor.chain.ref_start() } else { -1 },
                rc_pos: if anchor.is_fw { -1 } else { anchor.chain.ref_start() },
                format: None, // recovered pair: orientation not re-derived
                // SAM: the partner's leftmost is estimated from the fragment length.
                r1_pos: if anchor_is_left {
                    anchor.chain.ref_start()
                } else {
                    (anchor.chain.ref_start() + frag_len - partner_read.len() as i32).max(0)
                },
                r2_pos: if anchor_is_left {
                    (anchor.chain.ref_start() + frag_len - partner_read.len() as i32).max(0)
                } else {
                    anchor.chain.ref_start()
                },
                r2_fw: !anchor.is_fw,
                r1_score: if anchor_is_left { anchor_aln.score } else { partner_score },
            });
            return;
        }
    }

    let status = if anchor_is_left {
        MateStatus::PairedEndLeft
    } else {
        MateStatus::PairedEndRight
    };
    raw.push(RawMapping {
        tid,
        is_fw: anchor.is_fw,
        status,
        score: anchor_aln.score,
        fragment_len: 0,
        is_decoy,
        ref_pos: anchor.chain.ref_start(),
        // orphan: leftmost coordinate attributed to its own strand.
        fw_pos: if anchor.is_fw { anchor.chain.ref_start() } else { -1 },
        rc_pos: if anchor.is_fw { -1 } else { anchor.chain.ref_start() },
        format: None, // orphans are not sampled for library-type detection
        r1_pos: if anchor_is_left { anchor.chain.ref_start() } else { -1 },
        r2_pos: if anchor_is_left { -1 } else { anchor.chain.ref_start() },
        r2_fw: false,
        r1_score: anchor_aln.score,
    });
}

/// Try to place the partner mate within a fragment-length window of the anchor.
/// Returns `(partner_alignment_score, fragment_length)` on success.
fn recover_mate(
    partner_read: &[u8],
    refseq: &[u8],
    anchor: &MappingCandidate,
    cfg: &MapConfig,
) -> Option<(i32, i32)> {
    let max_frag = cfg.pair.max_fragment_len;
    let reflen = refseq.len() as i32;
    let a_s = anchor.chain.ref_start();
    let a_e = anchor.chain.ref_end();

    if anchor.is_fw {
        // The partner lies downstream and maps to the reverse strand.
        let win_start = a_s;
        let win_end = (a_s + max_frag).min(reflen);
        if win_end <= win_start {
            return None;
        }
        let q = revcomp(partner_read);
        let aln = align_in_window(&q, &refseq[win_start as usize..win_end as usize], &cfg.align)?;
        if !aln.valid {
            return None;
        }
        // partner ends `end_col` bases past the anchor start -> fragment length
        Some((aln.score, aln.end_col as i32))
    } else {
        // The partner lies upstream and maps to the forward strand.
        let win_end = a_e;
        let win_start = (a_e - max_frag).max(0);
        if win_end <= win_start {
            return None;
        }
        let aln = align_in_window(partner_read, &refseq[win_start as usize..win_end as usize], &cfg.align)?;
        if !aln.valid {
            return None;
        }
        let partner_end_abs = win_start + aln.end_col as i32;
        let partner_start_abs = partner_end_abs - partner_read.len() as i32;
        Some((aln.score, (a_e - partner_start_abs).max(0)))
    }
}
