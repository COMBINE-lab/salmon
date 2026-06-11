//! Scoring, deduplication, and filtering of validated mappings.
//!
//! After each mapping candidate is aligned, salmon collapses the alignments to
//! one per transcript (best score), drops the read if it aligns best to a
//! decoy, and assigns each surviving mapping a soft weight based on how far its
//! score is below the best. Those weighted transcripts become the read's
//! equivalence class. This mirrors salmon's `MappingScoreInfo` /
//! `filterAndCollectAlignments`.

use ahash::AHashMap;

use salmon_core::{LibraryFormat, MateStatus};

/// A validated mapping before final collapse/weighting.
#[derive(Debug, Clone, Copy)]
pub struct RawMapping {
    pub tid: u32,
    pub is_fw: bool,
    pub status: MateStatus,
    pub score: i32,
    pub fragment_len: i32,
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
    /// drop the read if `best_decoy >= decoy_thresh * best_valid`
    pub decoy_thresh: f64,
    /// keep only the best-scoring mappings (each weight `1.0`)
    pub hard_filter: bool,
}

impl Default for ScoreConfig {
    fn default() -> Self {
        Self {
            score_exp: 1.0,
            min_aln_prob: 1e-5,
            decoy_thresh: 1.0,
            hard_filter: false,
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
    if raw.is_empty() {
        return (Vec::new(), false, 0);
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
        return (Vec::new(), true, 0);
    };
    let best_decoy = best_per_tid
        .values()
        .filter(|m| m.is_decoy)
        .map(|m| m.score)
        .max();

    // Decoy domination (matches salmon's `bestScore < decoyThresh * bestDecoyScore`):
    // drop the read when its best transcript score falls below the decoy bar.
    if let Some(bd) = best_decoy {
        if (best_valid as f64) < cfg.decoy_thresh * (bd as f64) {
            return (Vec::new(), true, 0);
        }
    }
    let _ = had_decoy;

    let mut num_below = 0u32;
    let mut out: Vec<ScoredMapping> = best_per_tid
        .into_values()
        .filter(|m| !m.is_decoy)
        .filter_map(|m| {
            let weight = if cfg.hard_filter {
                if m.score == best_valid {
                    1.0
                } else {
                    num_below += 1;
                    return None;
                }
            } else {
                (-cfg.score_exp * (best_valid - m.score) as f64).exp()
            };
            if weight < cfg.min_aln_prob {
                num_below += 1;
                return None;
            }
            Some(ScoredMapping {
                tid: m.tid,
                is_fw: m.is_fw,
                status: m.status,
                score: m.score,
                fragment_len: m.fragment_len,
                weight,
                ref_pos: m.ref_pos,
                fw_pos: m.fw_pos,
                rc_pos: m.rc_pos,
                format: m.format,
                r1_pos: m.r1_pos,
                r2_pos: m.r2_pos,
                r2_fw: m.r2_fw,
                r1_score: m.r1_score,
            })
        })
        .collect();
    out.sort_by_key(|m| m.tid);
    (out, false, num_below)
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
    fn equal_scores_both_weight_one() {
        let m = finalize_mappings(vec![raw(0, 100, false), raw(1, 100, false)], &ScoreConfig::default());
        assert_eq!(m.len(), 2);
        assert!(m.iter().all(|x| (x.weight - 1.0).abs() < 1e-12));
    }

    #[test]
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
    fn dedup_per_transcript_keeps_best() {
        let m = finalize_mappings(
            vec![raw(0, 80, false), raw(0, 100, false), raw(0, 60, false)],
            &ScoreConfig::default(),
        );
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].score, 100);
    }

    #[test]
    fn decoy_domination_drops_read() {
        // decoy scores higher than the best transcript -> read dropped
        let m = finalize_mappings(vec![raw(0, 90, false), raw(1, 95, true)], &ScoreConfig::default());
        assert!(m.is_empty());
    }

    #[test]
    fn decoy_below_threshold_is_ignored() {
        let m = finalize_mappings(vec![raw(0, 100, false), raw(1, 80, true)], &ScoreConfig::default());
        assert_eq!(m.len(), 1);
        assert_eq!(m[0].tid, 0);
    }

    #[test]
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
    fn only_decoy_is_unmapped() {
        let m = finalize_mappings(vec![raw(0, 100, true)], &ScoreConfig::default());
        assert!(m.is_empty());
    }
}
