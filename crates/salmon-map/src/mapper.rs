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
use crate::collect::{best_per_target, collect_read_mems, MappingCandidate, MemCollectorConfig};
use crate::pair::{join_reads_and_filter, PairingConfig};
use crate::score::{finalize_mappings, RawMapping, ScoreConfig, ScoredMapping};
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
}

/// Map a single-end read to weighted equivalence-class members.
pub fn map_single_read<'idx, R: RefProvider>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    refs: &R,
    read: &[u8],
    cfg: &MapConfig,
) -> Vec<ScoredMapping> {
    let cands = best_per_target(collect_read_mems(index, hs, read, true, &cfg.collect));
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
                    format: Some(LibraryFormat::new(
                        ReadType::SingleEnd,
                        ReadOrientation::None,
                        strand,
                    )),
                });
            }
        }
    }
    finalize_mappings(raw, &cfg.score)
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
    let left = best_per_target(collect_read_mems(index, hs, r1, true, &cfg.collect));
    let right = best_per_target(collect_read_mems(index, hs, r2, true, &cfg.collect));
    let joints = join_reads_and_filter(left, right, &cfg.pair);

    let mut raw = Vec::new();
    for j in joints {
        match j.status {
            MateStatus::PairedEndPaired => {
                let l = j.left.as_ref().unwrap();
                let r = j.right.as_ref().unwrap();
                let refseq = refs.ref_seq(j.tid);
                let al = align_chain(r1, refseq, &l.chain, &cfg.align);
                let ar = align_chain(r2, refseq, &r.chain, &cfg.align);
                if let (Some(al), Some(ar)) = (al, ar) {
                    if al.valid && ar.valid {
                        raw.push(RawMapping {
                            tid: j.tid,
                            is_fw: l.is_fw,
                            status: MateStatus::PairedEndPaired,
                            score: al.score + ar.score,
                            fragment_len: j.fragment_len,
                            is_decoy: refs.is_decoy(j.tid),
                            ref_pos: l.chain.ref_start().min(r.chain.ref_start()),
                            format: Some(j.format),
                        });
                    }
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
    finalize_mappings(raw, &cfg.score)
}

#[inline]
fn index_seq<R: RefProvider>(refs: &R, tid: u32) -> &[u8] {
    refs.ref_seq(tid)
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
    let cands = best_per_target(collect_read_mems(index, hs, read, true, &cfg.collect));
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
                format: None, // recovered pair: orientation not re-derived
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
        format: None, // orphans are not sampled for library-type detection
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
