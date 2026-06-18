//! Paired-end joining: `joinReadsAndFilter`.
//!
//! Given the per-mate mapping candidates produced by [`collect`](crate::collect),
//! pair them into fragment mappings. Two candidates on the same reference form a
//! concordant pair when their fragment length is within bounds and (optionally)
//! they are not dovetailed; the observed [`LibraryFormat`] (orientation +
//! strandedness) is derived from the mates' strands and relative positions so a
//! [`LibraryTypeDetector`](salmon_model::LibraryTypeDetector) can sample it.
//! References that yield no concordant pair fall back to orphan mappings when
//! orphans are allowed. Mirrors pufferfish/salmon's `joinReadsAndFilter`.

use ahash::AHashSet;

use salmon_core::{LibraryFormat, MateStatus, ReadOrientation, ReadStrandedness, ReadType};

use crate::collect::MappingCandidate;

/// Pairing constraints (salmon's fragment-length / orphan / dovetail options).
#[derive(Debug, Clone)]
pub struct PairingConfig {
    /// maximum concordant fragment length
    pub max_fragment_len: i32,
    /// emit single-mate (orphan) mappings when no concordant pair exists
    pub allow_orphans: bool,
    /// permit dovetailed pairs (mates extending past each other)
    pub allow_dovetail: bool,
    /// Orphan chain sub-optimality threshold (salmon's `orphanChainSubThresh`,
    /// default there 0.95): when `> 0`, keep an orphan mapping only if its chain
    /// read-coverage is `>= orphan_chain_sub_thresh * (best orphan chain coverage)`.
    /// salmon uses this to trust the highest-chain-coverage orphan candidate and
    /// prune the rest *before* alignment (a speed heuristic that can drop a
    /// lower-coverage but better-aligning paralog). We default to **0.0 (off)** so
    /// every orphan candidate is aligned — slightly more sensitive than salmon on
    /// divergent gene-family reads. Raise toward 0.95 to mirror salmon exactly.
    pub orphan_chain_sub_thresh: f32,
    /// Post-merge concordant chain-pair sub-optimality threshold (salmon's
    /// `postMergeChainSubThresh`, default there 0.9): when `> 0`, after forming
    /// concordant pairs keep only those whose read-coverage is
    /// `>= post_merge_chain_sub_thresh * (best concordant-pair coverage)`,
    /// pruning weaker pairs before alignment. salmon thresholds on chain-pair
    /// *score*; we use read-coverage — consistent with the per-mate consensus
    /// filter and `orphan_chain_sub_thresh`. Default **0.0 (off)**: every
    /// concordant pair is aligned (current Rust behavior).
    pub post_merge_chain_sub_thresh: f32,
}

impl Default for PairingConfig {
    fn default() -> Self {
        Self {
            max_fragment_len: 1000,
            allow_orphans: true,
            allow_dovetail: false,
            orphan_chain_sub_thresh: 0.0,
            post_merge_chain_sub_thresh: 0.0,
        }
    }
}

/// A fragment mapping: a concordant pair, or a single-mate orphan.
#[derive(Debug, Clone)]
pub struct JointMapping {
    pub tid: u32,
    pub status: MateStatus,
    /// fragment length for a concordant pair; `0` for an orphan
    pub fragment_len: i32,
    /// observed library format (meaningful for concordant pairs)
    pub format: LibraryFormat,
    /// left/first-mate candidate (present unless this is a right-only orphan)
    pub left: Option<MappingCandidate>,
    /// right/second-mate candidate (present unless this is a left-only orphan)
    pub right: Option<MappingCandidate>,
}

impl JointMapping {
    /// Combined read-coverage of whichever mates are present (pairing-quality key).
    fn coverage(&self) -> i32 {
        self.left
            .as_ref()
            .map_or(0, |c| c.chain.covered_read_bases())
            + self
                .right
                .as_ref()
                .map_or(0, |c| c.chain.covered_read_bases())
    }
}

/// Observed library format for a concordant pair, where `l` is the first mate
/// and `r` the second. Reproduces salmon's orientation/strandedness encoding.
///
/// Orientation is decided from the mates' 5′ ends (the forward mate's 5′ is its
/// reference start; the reverse mate's 5′ is its reference end): inward
/// (TOWARD) when the forward mate's 5′ is upstream of the reverse mate's 5′,
/// otherwise outward (AWAY).
fn observed_format(l: &MappingCandidate, r: &MappingCandidate) -> LibraryFormat {
    let (orientation, strandedness) = if l.is_fw != r.is_fw {
        let (fw, rc) = if l.is_fw { (l, r) } else { (r, l) };
        let fw_5p = fw.chain.ref_start();
        let rc_5p = rc.chain.ref_end();
        let orientation = if fw_5p <= rc_5p {
            ReadOrientation::Toward
        } else {
            ReadOrientation::Away
        };
        // strandedness keyed on the first mate: sense (SA) if mate1 is forward.
        let strandedness = if l.is_fw {
            ReadStrandedness::SA
        } else {
            ReadStrandedness::AS
        };
        (orientation, strandedness)
    } else {
        // same strand
        let strandedness = if l.is_fw {
            ReadStrandedness::S
        } else {
            ReadStrandedness::A
        };
        (ReadOrientation::Same, strandedness)
    };
    LibraryFormat::new(ReadType::PairedEnd, orientation, strandedness)
}

/// Whether an inward (TOWARD) opposite-strand pair is dovetailed — the mates
/// overhang each other rather than nesting cleanly. Only inward pairs can
/// dovetail; outward/same configurations are left to the orientation check.
fn is_dovetailed(l: &MappingCandidate, r: &MappingCandidate, orientation: ReadOrientation) -> bool {
    if orientation != ReadOrientation::Toward || l.is_fw == r.is_fw {
        return false;
    }
    let (fw, rc) = if l.is_fw { (l, r) } else { (r, l) };
    // The forward mate should nest upstream of the reverse mate; an overhang on
    // either side is a dovetail.
    fw.chain.ref_start() > rc.chain.ref_start() || fw.chain.ref_end() > rc.chain.ref_end()
}

thread_local! {
    /// Reused per-thread `paired` set for [`join_reads_and_filter`] (the tids that
    /// formed a concordant pair). Cleared per call; reused so pairing allocates no
    /// per-read hash set. Grouping uses an in-place sort (no map), so the two
    /// per-read `AHashMap<u32, Vec<usize>>` of the old `group_by_tid` are gone.
    static PAIRED_TIDS: std::cell::RefCell<AHashSet<u32>> =
        std::cell::RefCell::new(AHashSet::new());
}

/// End (exclusive) of the contiguous same-`tid` run starting at `i` in a
/// tid-sorted candidate slice.
#[inline]
fn tid_run_end(cands: &[MappingCandidate], i: usize) -> usize {
    let tid = cands[i].tid;
    let mut e = i + 1;
    while e < cands.len() && cands[e].tid == tid {
        e += 1;
    }
    e
}

/// Pair left/right mate candidates into fragment mappings.
///
/// For each reference present on both mates, the best (highest combined
/// coverage) concordant pair within `max_fragment_len` is emitted. References
/// with no concordant pair contribute orphan mappings when `allow_orphans`.
///
/// Candidates are grouped by reference via an in-place sort + merge-join (cheap
/// struct moves, no per-read hash map or per-tid index vectors). The grouping
/// affects only the *order* of emitted records, not which targets/scores are
/// produced (downstream keys/dedups by tid).
pub fn join_reads_and_filter(
    mut left: Vec<MappingCandidate>,
    mut right: Vec<MappingCandidate>,
    cfg: &PairingConfig,
) -> Vec<JointMapping> {
    left.sort_unstable_by_key(|c| c.tid);
    right.sort_unstable_by_key(|c| c.tid);

    let mut joints = Vec::new();
    PAIRED_TIDS.with(|cell| {
        let paired = &mut *cell.borrow_mut();
        paired.clear();

        // Concordant pairs: merge-join the two tid-sorted candidate lists. For a
        // reference present on both mates we consider every (left, right) chain
        // pair (so a repeat-position chain that pairs is not lost to an earlier
        // single-best collapse) and clone only the best by combined coverage.
        let (mut i, mut j) = (0usize, 0usize);
        while i < left.len() && j < right.len() {
            let lt = left[i].tid;
            let rt = right[j].tid;
            if lt < rt {
                i = tid_run_end(&left, i);
                continue;
            }
            if rt < lt {
                j = tid_run_end(&right, j);
                continue;
            }
            let le = tid_run_end(&left, i);
            let re = tid_run_end(&right, j);
            let mut best: Option<(usize, usize, i32, LibraryFormat)> = None;
            let mut best_cov = i32::MIN;
            for li in i..le {
                for ri in j..re {
                    let l = &left[li];
                    let r = &right[ri];
                    let frag_start = l.chain.ref_start().min(r.chain.ref_start());
                    let frag_end = l.chain.ref_end().max(r.chain.ref_end());
                    let frag_len = frag_end - frag_start;
                    if frag_len <= 0 || frag_len > cfg.max_fragment_len {
                        continue;
                    }
                    let format = observed_format(l, r);
                    if !cfg.allow_dovetail && is_dovetailed(l, r, format.orientation) {
                        continue;
                    }
                    let cov = l.chain.covered_read_bases() + r.chain.covered_read_bases();
                    if best.is_none() || cov > best_cov {
                        best_cov = cov;
                        best = Some((li, ri, frag_len, format));
                    }
                }
            }
            if let Some((li, ri, frag_len, format)) = best {
                paired.insert(lt);
                joints.push(JointMapping {
                    tid: lt,
                    status: MateStatus::PairedEndPaired,
                    fragment_len: frag_len,
                    format,
                    left: Some(left[li].clone()),
                    right: Some(right[ri].clone()),
                });
            }
            i = le;
            j = re;
        }

        // Post-merge concordant chain-pair sub-optimality prune (salmon's
        // postMergeChainSubThresh; off by default). Drop concordant pairs whose
        // read-coverage is below `thresh * (best concordant coverage)`, before
        // alignment. Applied only to the concordant set (pruned tids stay in
        // `paired`, so they do not fall back to orphans).
        if cfg.post_merge_chain_sub_thresh > 0.0 && !joints.is_empty() {
            let best_cov = joints.iter().map(JointMapping::coverage).max().unwrap_or(0);
            let cutoff = (cfg.post_merge_chain_sub_thresh * best_cov as f32).ceil() as i32;
            joints.retain(|j| j.coverage() >= cutoff);
        }

        // Orphans for references that did not yield a concordant pair.
        if cfg.allow_orphans {
            // Optional orphan chain sub-optimality prune (salmon's orphanChainSubThresh;
            // off by default). When enabled, only orphan candidates whose chain
            // read-coverage is within `thresh` of the best orphan chain survive.
            let cutoff = if cfg.orphan_chain_sub_thresh > 0.0 {
                let best = left
                    .iter()
                    .chain(right.iter())
                    .filter(|c| !paired.contains(&c.tid))
                    .map(|c| c.chain.covered_read_bases())
                    .max()
                    .unwrap_or(0);
                (cfg.orphan_chain_sub_thresh * best as f32).ceil() as i32
            } else {
                0
            };
            // We pair over *all* per-target chains (no `best_per_target` pre-pass),
            // so a mate may have several chains to one reference. Emit at most one
            // orphan per (tid, is_fw): the highest-coverage unpaired chain.
            emit_best_orphans(
                &left,
                paired,
                cutoff,
                MateStatus::PairedEndLeft,
                &mut joints,
            );
            emit_best_orphans(
                &right,
                paired,
                cutoff,
                MateStatus::PairedEndRight,
                &mut joints,
            );
        }
    });

    joints
}

/// Emit one orphan per `(tid, is_fw)` for references with no concordant pair:
/// the highest-coverage unpaired chain in each orientation (≥ `cutoff`). Keeping
/// only the best per orientation avoids redundant orphan alignments when a mate
/// has several chains to one reference (repeat copies). `cands` is tid-sorted.
fn emit_best_orphans(
    cands: &[MappingCandidate],
    paired: &AHashSet<u32>,
    cutoff: i32,
    status: MateStatus,
    joints: &mut Vec<JointMapping>,
) {
    let mut i = 0;
    while i < cands.len() {
        let e = tid_run_end(cands, i);
        if !paired.contains(&cands[i].tid) {
            let (mut best_fw, mut best_rc): (Option<usize>, Option<usize>) = (None, None);
            for k in i..e {
                let cov = cands[k].chain.covered_read_bases();
                if cov < cutoff {
                    continue;
                }
                let slot = if cands[k].is_fw {
                    &mut best_fw
                } else {
                    &mut best_rc
                };
                if slot.is_none_or(|b| cov > cands[b].chain.covered_read_bases()) {
                    *slot = Some(k);
                }
            }
            for k in [best_fw, best_rc].into_iter().flatten() {
                joints.push(orphan(cands[k].clone(), status));
            }
        }
        i = e;
    }
}

fn orphan(c: MappingCandidate, status: MateStatus) -> JointMapping {
    // Orphans aren't sampled for library-type detection; give a placeholder
    // format carrying the mate's own strand.
    let strandedness = if c.is_fw {
        ReadStrandedness::S
    } else {
        ReadStrandedness::A
    };
    let format = LibraryFormat::new(ReadType::PairedEnd, ReadOrientation::None, strandedness);
    let (left, right) = match status {
        MateStatus::PairedEndLeft => (Some(c), None),
        _ => (None, Some(c)),
    };
    JointMapping {
        tid: if let Some(ref c) = left {
            c.tid
        } else {
            right.as_ref().unwrap().tid
        },
        status,
        fragment_len: 0,
        format,
        left,
        right,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chain::MemChain;
    use crate::mem::Mem;

    /// A candidate covering reference `[start, start+len)` on the given strand.
    fn cand(tid: u32, is_fw: bool, start: i32, len: i32) -> MappingCandidate {
        MappingCandidate {
            tid,
            is_fw,
            chain: MemChain::new(vec![Mem::new(0, start, len)], len as f32, is_fw),
        }
    }

    #[test]
    fn inward_pair_mate1_fw_is_isf() {
        // mate1 fw upstream, mate2 rc downstream -> ISF (Toward, SA).
        let left = vec![cand(0, true, 100, 50)];
        let right = vec![cand(0, false, 300, 50)];
        let j = join_reads_and_filter(left, right, &PairingConfig::default());
        assert_eq!(j.len(), 1);
        assert_eq!(j[0].status, MateStatus::PairedEndPaired);
        assert_eq!(j[0].tid, 0);
        assert_eq!(j[0].fragment_len, 350 - 100);
        assert_eq!(j[0].format.canonical(), "ISF");
    }

    #[test]
    fn inward_pair_mate1_rc_is_isr() {
        // mate2 fw upstream, mate1 rc downstream -> ISR (Toward, AS).
        let left = vec![cand(0, false, 300, 50)];
        let right = vec![cand(0, true, 100, 50)];
        let j = join_reads_and_filter(left, right, &PairingConfig::default());
        assert_eq!(j.len(), 1);
        assert_eq!(j[0].format.canonical(), "ISR");
        assert_eq!(j[0].fragment_len, 350 - 100);
    }

    #[test]
    fn fragment_too_long_falls_back_to_orphans() {
        let left = vec![cand(0, true, 100, 50)];
        let right = vec![cand(0, false, 5000, 50)];
        let j = join_reads_and_filter(left, right, &PairingConfig::default());
        // no concordant pair; two orphans (left + right) for tid 0
        assert_eq!(j.len(), 2);
        assert!(j.iter().all(|m| m.status.is_orphan()));
        assert!(j.iter().any(|m| m.status == MateStatus::PairedEndLeft));
        assert!(j.iter().any(|m| m.status == MateStatus::PairedEndRight));
    }

    #[test]
    fn pairs_repeat_position_chain_over_higher_coverage_nonpairing() {
        // mate1 matches tid 0 at two positions: a higher-coverage chain too far
        // from mate2 to pair, and a lower-coverage chain that pairs. Pairing over
        // ALL per-target chains must form the valid pair — the old `best_per_target`
        // pre-pass kept only the higher-coverage (non-pairing) chain and lost it.
        let left = vec![cand(0, true, 5000, 60), cand(0, true, 100, 50)];
        let right = vec![cand(0, false, 300, 50)];
        let j = join_reads_and_filter(left, right, &PairingConfig::default());
        let pairs: Vec<_> = j
            .iter()
            .filter(|m| m.status == MateStatus::PairedEndPaired)
            .collect();
        assert_eq!(pairs.len(), 1, "the valid (near) pair must be found: {j:?}");
        assert_eq!(pairs[0].tid, 0);
        // paired with the near (pairing) chain, not the far high-coverage one.
        assert_eq!(pairs[0].left.as_ref().unwrap().chain.ref_start(), 100);
        assert_eq!(pairs[0].fragment_len, 350 - 100);
    }

    #[test]
    fn orphan_when_only_one_mate_maps() {
        let left = vec![cand(7, true, 100, 50)];
        let right = vec![];
        let j = join_reads_and_filter(left, right, &PairingConfig::default());
        assert_eq!(j.len(), 1);
        assert_eq!(j[0].status, MateStatus::PairedEndLeft);
        assert_eq!(j[0].tid, 7);
        assert_eq!(j[0].fragment_len, 0);
    }

    #[test]
    fn orphans_suppressed_when_disabled() {
        let cfg = PairingConfig {
            allow_orphans: false,
            ..Default::default()
        };
        let left = vec![cand(7, true, 100, 50)];
        let j = join_reads_and_filter(left, vec![], &cfg);
        assert!(j.is_empty());
    }

    #[test]
    fn pair_preferred_over_orphans_on_same_tid() {
        // tid 0 has a concordant pair AND an extra left candidate; the pair wins
        // and no orphan is emitted for tid 0.
        let left = vec![cand(0, true, 100, 50), cand(0, true, 120, 40)];
        let right = vec![cand(0, false, 300, 50)];
        let j = join_reads_and_filter(left, right, &PairingConfig::default());
        assert_eq!(j.len(), 1);
        assert_eq!(j[0].status, MateStatus::PairedEndPaired);
    }

    #[test]
    fn post_merge_prunes_low_coverage_concordant_pair() {
        // tid 0: high-coverage concordant pair (cov 50+50=100); tid 1: low (10+10=20).
        let left = vec![cand(0, true, 100, 50), cand(1, true, 100, 10)];
        let right = vec![cand(0, false, 300, 50), cand(1, false, 300, 10)];
        let paired = |js: &[JointMapping]| {
            js.iter()
                .filter(|m| m.status == MateStatus::PairedEndPaired)
                .count()
        };
        // Off by default: both concordant pairs survive.
        let j = join_reads_and_filter(left.clone(), right.clone(), &PairingConfig::default());
        assert_eq!(paired(&j), 2);
        // thresh 0.5 -> cutoff ceil(0.5*100)=50; tid 1 (cov 20) is pruned and,
        // being already paired, does not fall back to an orphan.
        let cfg = PairingConfig {
            post_merge_chain_sub_thresh: 0.5,
            ..PairingConfig::default()
        };
        let j2 = join_reads_and_filter(left, right, &cfg);
        assert_eq!(j2.len(), 1);
        assert_eq!(j2[0].tid, 0);
        assert_eq!(paired(&j2), 1);
    }

    #[test]
    fn dovetail_filtered_by_default_but_kept_when_enabled() {
        // Inward pair (fw 5′=150 ≤ rc 5′=160) but the forward mate starts after
        // the reverse mate (150 > 100): the mates overhang -> dovetail.
        let make = || {
            (
                vec![cand(0, true, 150, 50)],  // fw mate [150,200]
                vec![cand(0, false, 100, 60)], // rc mate [100,160]
            )
        };
        // sanity: this configuration is inward and dovetailed
        let (l, r) = make();
        assert_eq!(
            observed_format(&l[0], &r[0]).orientation,
            ReadOrientation::Toward
        );
        assert!(is_dovetailed(&l[0], &r[0], ReadOrientation::Toward));

        // default: dovetail rejected -> falls back to orphans
        let (l, r) = make();
        let j = join_reads_and_filter(l, r, &PairingConfig::default());
        assert!(j.iter().all(|m| m.status.is_orphan()), "{j:?}");

        // enabled: kept as a concordant pair
        let (l, r) = make();
        let cfg = PairingConfig {
            allow_dovetail: true,
            ..Default::default()
        };
        let j = join_reads_and_filter(l, r, &cfg);
        assert_eq!(j.len(), 1);
        assert_eq!(j[0].status, MateStatus::PairedEndPaired);
    }

    #[test]
    fn normal_inward_pair_is_not_dovetailed() {
        let j = join_reads_and_filter(
            vec![cand(0, true, 100, 50)],
            vec![cand(0, false, 300, 50)],
            &PairingConfig::default(),
        );
        assert_eq!(j.len(), 1);
        assert_eq!(j[0].status, MateStatus::PairedEndPaired);
    }

    #[test]
    fn multiple_tids_each_pair() {
        let left = vec![cand(0, true, 100, 50), cand(1, true, 200, 50)];
        let right = vec![cand(0, false, 300, 50), cand(1, false, 400, 50)];
        let mut j = join_reads_and_filter(left, right, &PairingConfig::default());
        j.sort_by_key(|m| m.tid);
        assert_eq!(j.len(), 2);
        assert!(j.iter().all(|m| m.status == MateStatus::PairedEndPaired));
        assert_eq!(j[0].tid, 0);
        assert_eq!(j[1].tid, 1);
    }
}
