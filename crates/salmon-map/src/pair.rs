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
    /// Only emit orphans for a mate when the *other* mate is entirely unmapped
    /// (salmon's `--orphansRequireUnmappedMate`). Default **false**: when a pair
    /// has no concordant mapping, orphans are emitted for both mates (their
    /// union). When **true**, a read that mapped only because its mate also mapped
    /// (to a disjoint reference set) is not reported as an orphan — orphans are
    /// kept only for genuinely single-mappable fragments.
    pub orphans_require_unmapped_mate: bool,
}

impl Default for PairingConfig {
    fn default() -> Self {
        Self {
            max_fragment_len: 1000,
            allow_orphans: true,
            allow_dovetail: false,
            orphan_chain_sub_thresh: 0.0,
            post_merge_chain_sub_thresh: 0.0,
            orphans_require_unmapped_mate: false,
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

/// Whether an opposite-strand pair fails to nest cleanly — the "dovetail" test
/// salmon's `joinReadsAndFilter` applies (and discards under the default
/// no-dovetail policy). A clean inward (TOWARD) pair has the forward mate
/// nested upstream of the reverse mate (`fw.start <= rc.start` and
/// `fw.end <= rc.end`); anything else is flagged:
///   * the forward mate *ends* past the reverse mate (a true inward dovetail —
///     the mates overhang), or
///   * the forward mate *starts* past the reverse mate's start, i.e. the
///     forward mate lies (wholly or partly) downstream of the reverse mate — an
///     **outward / AWAY** pair, which is incompatible with an inward library.
///
/// Matching salmon: C++ keys on the forward mate starting after the reverse
/// mate (`fwd.start > rev.start`), which captures the outward case; we keep the
/// extra `fw.end > rc.end` term so a forward mate that merely overhangs the 3′
/// end is caught too. Same-strand pairs are not considered here. Previously this
/// returned early for non-inward orientations, so outward pairs (a fragment
/// mapping to a paralog whose homologous segments are in inverted arrangement)
/// were kept where salmon drops them.
fn is_dovetailed(l: &MappingCandidate, r: &MappingCandidate) -> bool {
    if l.is_fw == r.is_fw {
        return false;
    }
    let (fw, rc) = if l.is_fw { (l, r) } else { (r, l) };
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

/// Upper bound on concordant pairs emitted per target when several repeat loci
/// tie on combined seed coverage. Each emitted pair costs one extra pair of
/// alignments downstream; coverage ties are rare so this is almost never hit,
/// but it bounds fan-out on pathological tandem repeats.
const MAX_TIED_PAIRS_PER_TARGET: usize = 16;

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
/// For each reference present on both mates, every concordant pair within
/// `max_fragment_len` at the best combined seed coverage is emitted (pairing is
/// alignment-blind, so coverage-tied repeat loci are all forwarded and `align` +
/// `finalize` keep the highest-scoring placement per target). References with no
/// concordant pair contribute orphan mappings when `allow_orphans`.
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
            // Pass 1: the best combined seed coverage among concordant (li, ri)
            // pairs for this target.
            let mut best_cov = i32::MIN;
            for li in i..le {
                for ri in j..re {
                    let l = &left[li];
                    let r = &right[ri];
                    // A concordant pair has the mates on opposite strands (salmon's
                    // `satisfiesOri = lclust.isFw != rclust.isFw`). A same-strand
                    // pair is discordant for an inward library — it arises when a
                    // fragment maps to a paralog whose mate region is inverted, so
                    // both mates read the same strand; salmon never accepts it as
                    // concordant, and (when an opposite-strand pair exists) never
                    // falls back to the discordant round. Skip it here too.
                    if l.is_fw == r.is_fw {
                        continue;
                    }
                    let frag_start = l.chain.ref_start().min(r.chain.ref_start());
                    let frag_end = l.chain.ref_end().max(r.chain.ref_end());
                    let frag_len = frag_end - frag_start;
                    if frag_len <= 0 || frag_len > cfg.max_fragment_len {
                        continue;
                    }
                    if !cfg.allow_dovetail && is_dovetailed(l, r) {
                        continue;
                    }
                    let cov = l.chain.covered_read_bases() + r.chain.covered_read_bases();
                    if cov > best_cov {
                        best_cov = cov;
                    }
                }
            }
            // Pass 2: emit EVERY pair at the best coverage, not just the first.
            // Pairing is alignment-blind (coverage only), so when a mate matches a
            // target at several equal-coverage repeat loci the first-seen pair may
            // align worse — or be invalid (e.g. overhanging a transcript end) —
            // than another equal-coverage locus, silently dropping a co-optimal
            // paralog. Emitting all coverage-tied pairs lets `align`+`finalize`
            // (which keeps the highest-scoring placement per target) pick the
            // locus that actually aligns best. Capped to bound pathological
            // tandem-repeat fan-out; ties are rare so this adds ~no alignments.
            if best_cov > i32::MIN {
                paired.insert(lt);
                let mut emitted = 0usize;
                'pairs: for li in i..le {
                    for ri in j..re {
                        let l = &left[li];
                        let r = &right[ri];
                        // opposite-strand requirement (see pass 1)
                        if l.is_fw == r.is_fw {
                            continue;
                        }
                        let frag_start = l.chain.ref_start().min(r.chain.ref_start());
                        let frag_end = l.chain.ref_end().max(r.chain.ref_end());
                        let frag_len = frag_end - frag_start;
                        if frag_len <= 0 || frag_len > cfg.max_fragment_len {
                            continue;
                        }
                        if l.chain.covered_read_bases() + r.chain.covered_read_bases() != best_cov {
                            continue;
                        }
                        let format = observed_format(l, r);
                        if !cfg.allow_dovetail && is_dovetailed(l, r) {
                            continue;
                        }
                        joints.push(JointMapping {
                            tid: lt,
                            status: MateStatus::PairedEndPaired,
                            fragment_len: frag_len,
                            format,
                            left: Some(left[li].clone()),
                            right: Some(right[ri].clone()),
                        });
                        emitted += 1;
                        if emitted >= MAX_TIED_PAIRS_PER_TARGET {
                            break 'pairs;
                        }
                    }
                }
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
            //
            // Under `--orphansRequireUnmappedMate`, a mate's orphans are emitted
            // only when the *other* mate produced no candidates at all (it is truly
            // unmapped), so a read that mapped only alongside a mate mapping to a
            // disjoint reference set is not reported as an orphan.
            let req_unmapped = cfg.orphans_require_unmapped_mate;
            if !req_unmapped || right.is_empty() {
                emit_best_orphans(
                    &left,
                    paired,
                    cutoff,
                    MateStatus::PairedEndLeft,
                    &mut joints,
                );
            }
            if !req_unmapped || left.is_empty() {
                emit_best_orphans(
                    &right,
                    paired,
                    cutoff,
                    MateStatus::PairedEndRight,
                    &mut joints,
                );
            }
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
    fn emits_all_coverage_tied_pairs_per_target() {
        // mate1 matches tid 0 at two equal-coverage repeat loci (100 and 200),
        // both forming a valid concordant pair with mate2 at 300. Pairing is
        // alignment-blind, so it must emit BOTH tied pairs and let alignment pick
        // the locus that scores best — emitting only the first would silently drop
        // a co-optimal placement (the real-data sim.7081 bug).
        let left = vec![cand(0, true, 100, 50), cand(0, true, 200, 50)];
        let right = vec![cand(0, false, 300, 50)];
        let j = join_reads_and_filter(left, right, &PairingConfig::default());
        let pairs: Vec<_> = j
            .iter()
            .filter(|m| m.status == MateStatus::PairedEndPaired)
            .collect();
        assert_eq!(
            pairs.len(),
            2,
            "both coverage-tied loci must be emitted: {j:?}"
        );
        let starts: Vec<i32> = pairs
            .iter()
            .map(|m| m.left.as_ref().unwrap().chain.ref_start())
            .collect();
        assert!(
            starts.contains(&100) && starts.contains(&200),
            "starts={starts:?}"
        );
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
        assert!(is_dovetailed(&l[0], &r[0]));

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
    fn outward_pair_filtered_by_default_but_kept_when_enabled() {
        // OUTWARD (AWAY) pair: the forward mate lies entirely downstream of the
        // reverse mate (fw 5′=714 > rc 5′=236). This is the real-data sim.617247
        // case — a fragment that is a proper inward pair on its true transcript
        // but maps outward to a paralog whose homologous segments are in inverted
        // arrangement. salmon (and now Rust) treats this as a dovetail and drops
        // it under the default no-dovetail policy, since an outward pair is
        // incompatible with an inward library.
        let make = || {
            (
                vec![cand(0, true, 714, 50)],  // fw mate [714,764]
                vec![cand(0, false, 176, 60)], // rc mate [176,236]
            )
        };
        let (l, r) = make();
        assert_eq!(
            observed_format(&l[0], &r[0]).orientation,
            ReadOrientation::Away
        );
        assert!(is_dovetailed(&l[0], &r[0]), "outward pair must be flagged");

        // default: outward pair rejected -> falls back to orphans
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
        assert!(
            j.iter().any(|m| m.status == MateStatus::PairedEndPaired),
            "{j:?}"
        );
    }

    #[test]
    fn same_strand_pair_is_not_concordant() {
        // Both mates on the SAME strand (R1 fw, R2 fw). This is the real-data
        // sim.356845 case — a fragment that is a proper opposite-strand inward
        // pair on its true transcript but maps to a paralog whose mate region is
        // inverted, so both mates read the same strand. salmon requires opposite
        // strands for a concordant pair (`satisfiesOri = lclust.isFw != rclust.isFw`);
        // it must not be emitted as a concordant pair, regardless of allow_dovetail
        // (an orientation, not an overhang, problem). Note: auto library detection
        // leaves the downstream strand-compat filter inactive, so this geometric
        // check is the only thing that drops it — matching C++'s joinReadsAndFilter.
        let left = vec![cand(0, true, 232, 76)]; // R1 fw
        let right = vec![cand(0, true, 445, 76)]; // R2 fw (same strand)
        for allow in [false, true] {
            let cfg = PairingConfig {
                allow_dovetail: allow,
                ..Default::default()
            };
            let j = join_reads_and_filter(left.clone(), right.clone(), &cfg);
            assert!(
                !j.iter().any(|m| m.status == MateStatus::PairedEndPaired),
                "same-strand must never pair (allow_dovetail={allow}): {j:?}"
            );
        }
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
