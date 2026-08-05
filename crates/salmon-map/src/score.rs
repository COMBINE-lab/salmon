//! Scoring, deduplication, and filtering of validated mappings.
//!
//! # The last stage
//!
//! Alignment leaves a pile of scored placements, possibly several per transcript
//! and including decoys. This module reduces that to what the statistical model
//! consumes: at most one mapping per transcript, decoys removed, and each
//! survivor carrying a *soft weight*.
//!
//! The weight is the important part. A fragment compatible with several
//! transcripts should not count equally for all of them: a placement scoring far
//! below the best is much less likely to be the true origin. Weighting by score
//! difference, rather than keeping only the best, is what lets the EM apportion
//! ambiguous fragments sensibly instead of committing to a guess.
//!
//! Decoy handling is the same logic in the limit: a fragment that aligns *best*
//! to a decoy probably did not come from any transcript, so it is dropped rather
//! than forced onto whichever transcript it resembles most.
//!
//! After each mapping candidate is aligned, salmon collapses the alignments to
//! one per transcript (best score), drops the read if it aligns best to a
//! decoy, and assigns each surviving mapping a soft weight based on how far its
//! score is below the best. Those weighted transcripts become the read's
//! equivalence class. This mirrors salmon's `MappingScoreInfo` /
//! `filterAndCollectAlignments`.

use ahash::AHashMap;

use salmon_core::{LibraryFormat, MateStatus, RefProvider};

/// A validated mapping before final collapse/weighting.
#[derive(Debug, Clone, Copy)]
pub struct RawMapping {
    pub tid: u32,
    pub is_fw: bool,
    pub status: MateStatus,
    pub score: i32,
    pub fragment_len: i32,
    /// length of the mapped (anchor) read; used for the ambiguous
    /// fragment-length probability of orphans / single-end reads. `0` when not
    /// applicable.
    pub read_len: i32,
    pub is_decoy: bool,
    /// leftmost fragment position on the reference (for sequence-bias context)
    pub ref_pos: i32,
    /// leftmost reference position of the forward-strand mate/read (for
    /// positional bias); `-1` if there is no forward-strand contributor
    pub fw_pos: i32,
    /// leftmost reference position of the reverse-strand mate/read (for
    /// positional bias); `-1` if there is no reverse-strand contributor
    pub rc_pos: i32,
    /// observed library format of this mapping (for auto-detection); `None`
    /// for orphans and pseudoalignment hits
    pub format: Option<LibraryFormat>,
    /// SAM output (`--writeMappings`): leftmost reference position of read1 and
    /// read2 (`-1` if that mate is absent), and read2's strand. read1's strand is
    /// `is_fw`. Spoofed-CIGAR SAM needs per-mate positions, which `ref_pos`
    /// (a single fragment-leftmost) does not preserve.
    pub r1_pos: i32,
    pub r2_pos: i32,
    pub r2_fw: bool,
    /// read1's alignment score (for the SAM `AS` tag); read2's is `score - r1_score`.
    pub r1_score: i32,
}

/// A surviving mapping with its equivalence-class weight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScoredMapping {
    pub tid: u32,
    pub is_fw: bool,
    pub status: MateStatus,
    pub score: i32,
    pub fragment_len: i32,
    /// length of the mapped (anchor) read; used for the ambiguous
    /// fragment-length probability of orphans / single-end reads (the reverse-
    /// strand bound is `ref_pos + read_len`). `0` when not applicable.
    pub read_len: i32,
    /// equivalence-class weight (`1.0` for best-scoring; decays below)
    pub weight: f64,
    /// leftmost fragment position on the reference (for sequence-bias context)
    pub ref_pos: i32,
    /// leftmost reference position of the forward-strand mate/read (positional
    /// bias); `-1` if none
    pub fw_pos: i32,
    /// leftmost reference position of the reverse-strand mate/read (positional
    /// bias); `-1` if none
    pub rc_pos: i32,
    /// observed library format (for auto-detection), if determinable
    pub format: Option<LibraryFormat>,
    /// SAM output: per-mate leftmost positions (`-1` if absent) and read2 strand.
    pub r1_pos: i32,
    pub r2_pos: i32,
    pub r2_fw: bool,
    /// read1's alignment score (SAM `AS`); read2's is `score - r1_score`.
    pub r1_score: i32,
}

/// Filtering / weighting parameters (salmon `--scoreExp`, `--minAlnProb`,
/// `--decoyThreshold`, `--hardFilter`).
#[derive(Debug, Clone)]
pub struct ScoreConfig {
    /// soft-weight decay: `weight = exp(-score_exp * (best - score))`
    pub score_exp: f64,
    /// drop mappings whose weight falls below this
    pub min_aln_prob: f64,
    /// drop the read only when a decoy is *strictly* better than the best
    /// transcript: `best_valid < decoy_thresh * best_decoy` (so at the default
    /// `decoy_thresh = 1.0` a decoy that merely ties the transcript does **not**
    /// dominate — the transcript is kept). Matches C++ salmon's
    /// `bestScore < decoyThresh * bestDecoyScore`.
    pub decoy_thresh: f64,
    /// keep only the best-scoring mappings (each weight `1.0`)
    pub hard_filter: bool,
    /// When a fragment has a valid (non-decoy) transcript alignment but a decoy
    /// outscores it, keep the transcript alignment(s) as an orphan/mapping instead
    /// of discarding the whole fragment (`--allowDecoyOrphans`). Off by default
    /// (decoy-dominated fragments are dropped, matching salmon's default). This
    /// only affects fragments that *have* a surviving transcript candidate; it
    /// cannot recover fragments whose transcript placements were never generated
    /// (e.g. high-occurrence multimappers that only seed the genome).
    pub allow_decoy_orphans: bool,
}

impl Default for ScoreConfig {
    fn default() -> Self {
        Self {
            score_exp: 1.0,
            min_aln_prob: 1e-5,
            decoy_thresh: 1.0,
            hard_filter: false,
            allow_decoy_orphans: false,
        }
    }
}

/// Collapse to one mapping per transcript (best score), apply decoy
/// domination, and weight the survivors. Returns the read's weighted
/// equivalence-class members (sorted by `tid`), or empty if the read is
/// unmapped or decoy-dominated.
pub fn finalize_mappings(raw: Vec<RawMapping>, cfg: &ScoreConfig) -> Vec<ScoredMapping> {
    finalize_mappings_counted(raw, cfg).0
}

/// Like [`finalize_mappings`] but also reports `(decoy_dominated, num_below_thresh)`:
/// whether the fragment was dropped because its best placement was a decoy, and
/// how many surviving non-decoy alignments were filtered for weight `< min_aln_prob`.
pub fn finalize_mappings_counted(
    raw: Vec<RawMapping>,
    cfg: &ScoreConfig,
) -> (Vec<ScoredMapping>, bool, u32) {
    let mut out = Vec::new();
    let (decoy_dominated, num_below) = finalize_mappings_counted_into(&mut out, raw, cfg);
    (out, decoy_dominated, num_below)
}

/// Like [`finalize_mappings_counted`] but writes into a caller-provided buffer
/// (cleared first), so a per-thread `Vec` can be reused across reads instead of
/// allocating a fresh result each time. Returns `(decoy_dominated, num_below)`.
pub fn finalize_mappings_counted_into(
    out: &mut Vec<ScoredMapping>,
    raw: Vec<RawMapping>,
    cfg: &ScoreConfig,
) -> (bool, u32) {
    out.clear();
    if raw.is_empty() {
        return (false, 0);
    }

    // One mapping per transcript: keep the highest score.
    let mut best_per_tid: AHashMap<u32, RawMapping> = AHashMap::new();
    for m in raw {
        best_per_tid
            .entry(m.tid)
            .and_modify(|e| {
                if m.score > e.score {
                    *e = m;
                }
            })
            .or_insert(m);
    }

    let best_valid = best_per_tid
        .values()
        .filter(|m| !m.is_decoy)
        .map(|m| m.score)
        .max();
    let had_decoy = best_per_tid.values().any(|m| m.is_decoy);
    let Some(best_valid) = best_valid else {
        // only decoy mappings -> dropped, dominated by a decoy
        return (true, 0);
    };
    let best_decoy = best_per_tid
        .values()
        .filter(|m| m.is_decoy)
        .map(|m| m.score)
        .max();

    // Decoy domination (matches salmon's `bestScore < decoyThresh * bestDecoyScore`):
    // drop the read when its best transcript score falls below the decoy bar —
    // unless `--allowDecoyOrphans` asks us to keep the transcript placement(s).
    if let Some(bd) = best_decoy {
        if (best_valid as f64) < cfg.decoy_thresh * (bd as f64) && !cfg.allow_decoy_orphans {
            return (true, 0);
        }
    }
    let _ = had_decoy;

    let mut num_below = 0u32;
    for m in best_per_tid.into_values().filter(|m| !m.is_decoy) {
        let weight = if cfg.hard_filter {
            if m.score == best_valid {
                1.0
            } else {
                num_below += 1;
                continue;
            }
        } else {
            (-cfg.score_exp * (best_valid - m.score) as f64).exp()
        };
        if weight < cfg.min_aln_prob {
            num_below += 1;
            continue;
        }
        out.push(ScoredMapping {
            tid: m.tid,
            is_fw: m.is_fw,
            status: m.status,
            score: m.score,
            fragment_len: m.fragment_len,
            read_len: m.read_len,
            weight,
            ref_pos: m.ref_pos,
            fw_pos: m.fw_pos,
            rc_pos: m.rc_pos,
            format: m.format,
            r1_pos: m.r1_pos,
            r2_pos: m.r2_pos,
            r2_fw: m.r2_fw,
            r1_score: m.r1_score,
        });
    }
    out.sort_by_key(|m| m.tid);
    (false, num_below)
}

/// Apply decoy filtering to sketch-mode mappings, consistent with the decoy
/// logic [`finalize_mappings_counted`] applies on the selective-alignment path.
///
/// Sketch mappings ([`map_read_pair_sketch`](crate::map_read_pair_sketch) /
/// [`map_single_read_sketch`](crate::map_single_read_sketch)) are built directly
/// from piscem's accepted hits and carry no per-mapping `is_decoy` flag, so decoy
/// references would otherwise leak into the equivalence classes. Decoys are
/// identified here via `refs.is_decoy`, using the sketch *score* (matching-k-mer
/// coverage) as the comparison key — the same role the alignment score plays in
/// selective alignment.
///
/// Returns `true` if the fragment is decoy-dominated (in which case `maps` is
/// cleared and the fragment is unmapped); otherwise the decoy mappings are
/// removed from `maps` (so only transcript hits enter the eq-class) and `false`
/// is returned. With `cfg.allow_decoy_orphans`, a fragment that has any non-decoy
/// hit keeps its transcript hits even when a decoy scores at least as high
/// (`--allowDecoyOrphans`), mirroring the SA-mode flag.
pub fn filter_sketch_decoys<R: RefProvider>(
    maps: &mut Vec<ScoredMapping>,
    refs: &R,
    cfg: &ScoreConfig,
) -> bool {
    if maps.is_empty() {
        return false;
    }
    let best_valid = maps
        .iter()
        .filter(|m| !refs.is_decoy(m.tid))
        .map(|m| m.score)
        .max();
    let best_decoy = maps
        .iter()
        .filter(|m| refs.is_decoy(m.tid))
        .map(|m| m.score)
        .max();
    match best_valid {
        // No transcript hit at all. If any decoy was hit, the fragment is
        // decoy-dominated (its signal is best explained by the genome); either
        // way it carries no transcript placement.
        None => {
            let had_decoy = best_decoy.is_some();
            maps.clear();
            had_decoy
        }
        Some(best_valid) => {
            // Decoy domination (matches salmon's `bestScore < decoyThresh *
            // bestDecoyScore`): drop the fragment when its best transcript score
            // falls below the decoy bar, unless `--allowDecoyOrphans` keeps the
            // transcript hits.
            if let Some(bd) = best_decoy {
                if (best_valid as f64) < cfg.decoy_thresh * (bd as f64) && !cfg.allow_decoy_orphans
                {
                    maps.clear();
                    return true;
                }
            }
            // Keep only transcript mappings; decoy tids never enter the eq-class.
            maps.retain(|m| !refs.is_decoy(m.tid));
            false
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn raw(tid: u32, score: i32, is_decoy: bool) -> RawMapping {
        RawMapping {
            tid,
            is_fw: true,
            status: MateStatus::SingleEnd,
            score,
            fragment_len: 0,
            read_len: 0,
            is_decoy,
            ref_pos: 0,
            fw_pos: 0,
            rc_pos: -1,
            format: None,
            r1_pos: 0,
            r2_pos: -1,
            r2_fw: false,
            r1_score: score,
        }
    }

    #[test]
    /// Equally good placements must be weighted equally — the EM decides between
    /// them, not the mapper.
    fn equal_scores_both_weight_one() {
        let m = finalize_mappings(
            vec![raw(0, 100, false), raw(1, 100, false)],
            &ScoreConfig::default(),
        );
        assert_eq!(m.len(), 2);
        assert!(m.iter().all(|x| (x.weight - 1.0).abs() < 1e-12));
    }

    #[test]
    /// A worse placement must earn a smaller weight, and a far worse one must drop
    /// out entirely.
    fn lower_score_decays_and_can_be_filtered() {
        // score_exp 1.0: a score 10 below best -> weight exp(-10) ~ 4.5e-5
        let cfg = ScoreConfig::default();
        let m = finalize_mappings(vec![raw(0, 100, false), raw(1, 90, false)], &cfg);
        // exp(-10) ~ 4.5e-5 >= min_aln_prob 1e-5 -> kept but tiny
        let w1 = m.iter().find(|x| x.tid == 1).unwrap().weight;
        assert!(w1 < 1e-3 && w1 > 1e-5, "weight {w1}");
        // raise the floor -> the weak mapping is dropped
        let cfg2 = ScoreConfig {
            min_aln_prob: 1e-3,
            ..Default::default()
        };
        let m2 = finalize_mappings(vec![raw(0, 100, false), raw(1, 90, false)], &cfg2);
        assert_eq!(m2.len(), 1);
        assert_eq!(m2[0].tid, 0);
    }

    #[test]
    /// One mapping per transcript, the best; several would count the same evidence
    /// repeatedly.
    fn dedup_per_transcript_keeps_best() {
        let m = finalize_mappings(
            vec![raw(0, 80, false), raw(0, 100, false), raw(0, 60, false)],
            &ScoreConfig::default(),
        );
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].score, 100);
    }

    #[test]
    /// A read aligning best to a decoy probably did not come from any transcript,
    /// so it must be discarded rather than forced onto the nearest one.
    fn decoy_domination_drops_read() {
        // decoy scores higher than the best transcript -> read dropped
        let m = finalize_mappings(
            vec![raw(0, 90, false), raw(1, 95, true)],
            &ScoreConfig::default(),
        );
        assert!(m.is_empty());
    }

    #[test]
    /// A decoy that aligns worse than a transcript is not evidence against it.
    fn decoy_below_threshold_is_ignored() {
        let m = finalize_mappings(
            vec![raw(0, 100, false), raw(1, 80, true)],
            &ScoreConfig::default(),
        );
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].tid, 0);
    }

    #[test]
    /// `--hardFilter` must discard the soft weights and keep only the best
    /// placements.
    fn hard_filter_keeps_only_best() {
        let cfg = ScoreConfig {
            hard_filter: true,
            ..Default::default()
        };
        let m = finalize_mappings(
            vec![raw(0, 100, false), raw(1, 100, false), raw(2, 95, false)],
            &cfg,
        );
        assert_eq!(m.len(), 2);
        assert!(m.iter().all(|x| x.score == 100 && x.weight == 1.0));
    }

    #[test]
    /// A fragment whose only placement is a decoy counts as unmapped.
    fn only_decoy_is_unmapped() {
        let m = finalize_mappings(vec![raw(0, 100, true)], &ScoreConfig::default());
        assert!(m.is_empty());
    }

    // --- sketch-mode decoy filtering ---

    /// A `RefProvider` where every tid `>= first_decoy` is a decoy.
    struct DecoyAt(u32);
    impl RefProvider for DecoyAt {
        fn num_refs(&self) -> usize {
            u32::MAX as usize
        }
        fn ref_seq(&self, _tid: u32) -> &[u8] {
            &[]
        }
        fn is_decoy(&self, tid: u32) -> bool {
            tid >= self.0
        }
    }

    fn scored(tid: u32, score: i32) -> ScoredMapping {
        ScoredMapping {
            tid,
            is_fw: true,
            status: MateStatus::PairedEndPaired,
            score,
            fragment_len: 200,
            read_len: 0,
            weight: 1.0,
            ref_pos: 0,
            fw_pos: -1,
            rc_pos: -1,
            format: None,
            r1_pos: -1,
            r2_pos: -1,
            r2_fw: false,
            r1_score: score,
        }
    }

    #[test]
    /// The same decoy logic on the alignment-free path, which has no scores to
    /// compare and so must filter by identity.
    fn sketch_drops_decoy_tids_but_keeps_transcripts() {
        // tids 0,1 transcript; tid 2 decoy. Transcript outscores decoy.
        let refs = DecoyAt(2);
        let mut maps = vec![scored(0, 30), scored(1, 28), scored(2, 25)];
        let dominated = filter_sketch_decoys(&mut maps, &refs, &ScoreConfig::default());
        assert!(!dominated);
        assert_eq!(maps.len(), 2);
        assert!(maps.iter().all(|m| m.tid < 2), "decoy tid must be removed");
    }

    #[test]
    /// Sketch mode must also recognize a decoy-only fragment as unmapped.
    fn sketch_decoy_only_is_decoy_dominated() {
        let refs = DecoyAt(0); // everything is a decoy
        let mut maps = vec![scored(0, 30), scored(1, 25)];
        let dominated = filter_sketch_decoys(&mut maps, &refs, &ScoreConfig::default());
        assert!(dominated);
        assert!(maps.is_empty());
    }

    #[test]
    /// And a dominating decoy must still win in sketch mode.
    fn sketch_decoy_dominates_higher_score() {
        // decoy (tid 2) scores higher than the best transcript -> dropped.
        let refs = DecoyAt(2);
        let mut maps = vec![scored(0, 20), scored(2, 30)];
        let dominated = filter_sketch_decoys(&mut maps, &refs, &ScoreConfig::default());
        assert!(dominated);
        assert!(maps.is_empty());
    }

    #[test]
    /// With decoy orphans allowed, a genuine transcript placement must survive
    /// alongside the decoy rather than being suppressed by it.
    fn sketch_allow_decoy_orphans_keeps_transcript() {
        // Same as above but --allowDecoyOrphans keeps the transcript hit even
        // though the decoy scores higher (the mate maps to the decoy).
        let refs = DecoyAt(2);
        let cfg = ScoreConfig {
            allow_decoy_orphans: true,
            ..Default::default()
        };
        let mut maps = vec![scored(0, 20), scored(2, 30)];
        let dominated = filter_sketch_decoys(&mut maps, &refs, &cfg);
        assert!(!dominated);
        assert_eq!(maps.len(), 1);
        assert_eq!(maps[0].tid, 0);
    }
}
