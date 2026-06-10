//! Colinear MEM chaining.
//!
//! Given a set of [`Mem`] anchors that all lie on one reference in one
//! orientation, find high-scoring colinear chains. This is the
//! minimap2-style dynamic program that pufferfish/salmon use to assemble
//! uni-MEMs into mapping candidates before alignment validation.
//!
//! ## Algorithm
//! Anchors are sorted by `(ref_start, read_start)`. A DP computes, for each
//! anchor `i`, the best chain score ending at `i`:
//!
//! ```text
//! f[i] = max( len_i,  max_{j precedes i} f[j] + gain(j,i) - gap_cost(j,i) )
//! ```
//!
//! where a valid predecessor `j` is strictly before `i` on both the read and
//! the reference and within the gap bandwidth; `gain` is the number of newly
//! matched bases (`min(read_gap, ref_gap, len_i)`); and `gap_cost` is
//! minimap2's affine-plus-log penalty on the difference between the read and
//! reference gaps. Chains are then recovered by backtracking from the highest
//! peaks, and reported if their score is within `chain_subopt_thresh` of the
//! best chain.
//!
//! Exact numerical parity with pufferfish's chainer is validated empirically
//! against `sample_data` in a later phase; the tunables here map onto salmon's
//! `consensusSlack` / `pre`/`postMergeChainSubThresh` options.

use crate::mem::{covered_read_bases, Mem};

/// Tuning parameters for [`chain_mems`].
#[derive(Debug, Clone)]
pub struct ChainConfig {
    /// Maximum allowed gap (in bases) between consecutive anchors, on either
    /// the read or the reference axis. Bounds the chaining bandwidth.
    pub max_gap: i32,
    /// Average seed length used in the gap cost (typically the index k).
    pub seed_len: i32,
    /// Report chains whose score is at least `best_score * chain_subopt_thresh`.
    /// `1.0` keeps only optimal-scoring chains; smaller values keep more.
    pub chain_subopt_thresh: f32,
    /// Maximum number of predecessors to examine per anchor (lookback bound).
    /// `0` means examine all earlier anchors.
    pub max_lookback: usize,
}

impl Default for ChainConfig {
    fn default() -> Self {
        Self {
            max_gap: 5000,
            seed_len: 31,
            chain_subopt_thresh: 0.8,
            max_lookback: 0,
        }
    }
}

/// A colinear chain of MEMs on a single reference / orientation.
#[derive(Debug, Clone)]
pub struct MemChain {
    /// the anchors in the chain, in colinear order
    pub mems: Vec<Mem>,
    /// DP score of the chain
    pub score: f32,
    /// orientation of the underlying mapping (stamped by the caller)
    pub is_fw: bool,
}

impl MemChain {
    /// First matched read position in the chain.
    pub fn read_start(&self) -> i32 {
        self.mems.first().map(|m| m.read_start).unwrap_or(0)
    }

    /// One past the last matched read position in the chain.
    pub fn read_end(&self) -> i32 {
        self.mems.iter().map(|m| m.read_end()).max().unwrap_or(0)
    }

    /// First matched reference position in the chain.
    pub fn ref_start(&self) -> i32 {
        self.mems.first().map(|m| m.ref_start).unwrap_or(0)
    }

    /// One past the last matched reference position in the chain.
    pub fn ref_end(&self) -> i32 {
        self.mems.iter().map(|m| m.ref_end()).max().unwrap_or(0)
    }

    /// Number of read bases covered (overlapping anchors counted once). This is
    /// the quantity compared against the read length for the consensus filter.
    pub fn covered_read_bases(&self) -> i32 {
        let idx: Vec<usize> = (0..self.mems.len()).collect();
        covered_read_bases(&self.mems, &idx)
    }
}

/// minimap2-style gap cost: 0 when the read and reference gaps match, otherwise
/// an affine term proportional to the gap plus a logarithmic term.
#[inline]
fn gap_cost(gap: i32, seed_len: i32) -> f32 {
    if gap == 0 {
        0.0
    } else {
        0.01 * seed_len as f32 * gap as f32 + 0.5 * (gap as f32).log2()
    }
}

/// Find colinear chains among `mems` (all on one reference / orientation).
///
/// Returns chains sorted by descending score, filtered to those within
/// `cfg.chain_subopt_thresh` of the best. `is_fw` stamps the orientation onto
/// each returned chain. An empty input yields no chains.
pub fn chain_mems(mems: &[Mem], is_fw: bool, cfg: &ChainConfig) -> Vec<MemChain> {
    if mems.is_empty() {
        return Vec::new();
    }

    // Sort anchors colinearly; keep an index permutation so we can report the
    // original Mem values.
    let mut order: Vec<usize> = (0..mems.len()).collect();
    order.sort_by(|&a, &b| {
        (mems[a].ref_start, mems[a].read_start).cmp(&(mems[b].ref_start, mems[b].read_start))
    });
    let sorted: Vec<Mem> = order.iter().map(|&i| mems[i]).collect();
    let n = sorted.len();

    let mut f = vec![0f32; n]; // best score ending at i
    let mut p = vec![usize::MAX; n]; // predecessor (MAX = none)

    for i in 0..n {
        f[i] = sorted[i].len as f32;
        let lb = if cfg.max_lookback == 0 {
            0
        } else {
            i.saturating_sub(cfg.max_lookback)
        };
        for j in (lb..i).rev() {
            let dr = sorted[i].ref_start - sorted[j].ref_start;
            let dq = sorted[i].read_start - sorted[j].read_start;
            // strictly colinear and increasing on both axes
            if dr <= 0 || dq <= 0 {
                continue;
            }
            if dr > cfg.max_gap || dq > cfg.max_gap {
                continue;
            }
            let gap = (dr - dq).abs();
            // newly matched bases contributed by anchor i
            let gain = dq.min(dr).min(sorted[i].len) as f32;
            let sc = f[j] + gain - gap_cost(gap, cfg.seed_len);
            if sc > f[i] {
                f[i] = sc;
                p[i] = j;
            }
        }
    }

    // Recover chains: repeatedly take the highest-scoring unused peak and
    // backtrack until reaching an unused-chain boundary.
    let mut peaks: Vec<usize> = (0..n).collect();
    peaks.sort_by(|&a, &b| f[b].total_cmp(&f[a]));
    let mut used = vec![false; n];
    let mut chains: Vec<MemChain> = Vec::new();

    for &peak in &peaks {
        if used[peak] {
            continue;
        }
        let score = f[peak];
        let mut idx = peak;
        let mut chain_idx = Vec::new();
        loop {
            if used[idx] {
                break;
            }
            used[idx] = true;
            chain_idx.push(idx);
            if p[idx] == usize::MAX {
                break;
            }
            idx = p[idx];
        }
        chain_idx.reverse(); // now in colinear order
        let chain_mems: Vec<Mem> = chain_idx.iter().map(|&k| sorted[k]).collect();
        chains.push(MemChain {
            mems: chain_mems,
            score,
            is_fw,
        });
    }

    // Filter by sub-optimality threshold and sort by descending score.
    let best = chains.iter().map(|c| c.score).fold(f32::MIN, f32::max);
    let cutoff = best * cfg.chain_subopt_thresh;
    chains.retain(|c| c.score >= cutoff);
    chains.sort_by(|a, b| b.score.total_cmp(&a.score));
    chains
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg() -> ChainConfig {
        ChainConfig {
            max_gap: 5000,
            seed_len: 10,
            chain_subopt_thresh: 0.8,
            max_lookback: 0,
        }
    }

    #[test]
    fn empty_input() {
        assert!(chain_mems(&[], true, &cfg()).is_empty());
    }

    #[test]
    fn single_mem_one_chain() {
        let mems = [Mem::new(0, 100, 20)];
        let chains = chain_mems(&mems, true, &cfg());
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].mems.len(), 1);
        assert_eq!(chains[0].score, 20.0);
        assert!(chains[0].is_fw);
    }

    #[test]
    fn colinear_mems_merge_into_one_chain() {
        // Two gap-free anchors on the same diagonal -> single chain.
        let mems = [Mem::new(0, 100, 15), Mem::new(20, 120, 15)];
        let chains = chain_mems(&mems, true, &cfg());
        assert_eq!(chains.len(), 1, "expected a single chain, got {chains:?}");
        assert_eq!(chains[0].mems.len(), 2);
        // diagonal is consistent (gap == 0), so no gap penalty: 15 + min(20,20,15)=15 => 30
        assert!((chains[0].score - 30.0).abs() < 1e-3, "score={}", chains[0].score);
        assert_eq!(chains[0].covered_read_bases(), 30);
    }

    #[test]
    fn anti_colinear_mems_split() {
        // Second anchor is earlier on the read but later on the ref -> not
        // colinear, so two separate single-anchor chains.
        let mems = [Mem::new(50, 100, 15), Mem::new(0, 200, 15)];
        let chains = chain_mems(&mems, true, &cfg());
        assert_eq!(chains.len(), 2);
        assert!(chains.iter().all(|c| c.mems.len() == 1));
    }

    #[test]
    fn gap_penalty_reduces_score() {
        // Same anchors, but a large indel between them (read gap != ref gap).
        let mems = [Mem::new(0, 100, 15), Mem::new(20, 220, 15)];
        let chains = chain_mems(&mems, true, &cfg());
        // ref gap=120, read gap=20 -> gap=100; still within max_gap on read but
        // dr=120 < max_gap so chainable, with a penalty.
        assert_eq!(chains.len(), 1);
        let two = &chains[0];
        assert_eq!(two.mems.len(), 2);
        // score is the colinear score minus the gap cost, so < 30
        assert!(two.score < 30.0, "score={} should be penalized", two.score);
    }

    #[test]
    fn gap_beyond_bandwidth_does_not_chain() {
        let mut c = cfg();
        c.max_gap = 50;
        // ref gap 500 exceeds bandwidth
        let mems = [Mem::new(0, 100, 15), Mem::new(20, 600, 15)];
        let chains = chain_mems(&mems, true, &c);
        assert_eq!(chains.len(), 2, "should not chain across a 500bp gap");
    }

    #[test]
    fn subopt_threshold_filters_weak_chains() {
        // One strong colinear pair plus one isolated short anchor far away.
        let mems = [
            Mem::new(0, 100, 30),
            Mem::new(40, 140, 30),    // chains with the first -> strong
            Mem::new(5, 5000, 5),     // isolated, weak
        ];
        let mut c = cfg();
        c.chain_subopt_thresh = 0.5;
        let chains = chain_mems(&mems, true, &c);
        // strong chain ~60; weak chain =5; 5 < 0.5*60 -> filtered out
        assert_eq!(chains.len(), 1, "weak chain should be filtered: {chains:?}");
        assert!(chains[0].score >= 60.0);
    }

    #[test]
    fn chains_sorted_by_descending_score() {
        let mems = [
            Mem::new(0, 100, 30),
            Mem::new(40, 140, 30), // strong chain
            Mem::new(0, 3000, 25), // medium isolated
        ];
        let mut c = cfg();
        c.chain_subopt_thresh = 0.1;
        let chains = chain_mems(&mems, true, &c);
        assert!(chains.len() >= 2);
        for w in chains.windows(2) {
            assert!(w[0].score >= w[1].score);
        }
    }
}
