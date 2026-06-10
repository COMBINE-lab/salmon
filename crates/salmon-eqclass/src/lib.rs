//! Equivalence classes for salmon-style quantification.
//!
//! A fragment that is compatible with a set of transcripts contributes to an
//! *equivalence class* keyed by that transcript set. Counts and per-transcript
//! auxiliary weights are aggregated per class, then handed to the EM/VBEM
//! optimizer. This mirrors salmon's `EquivalenceClassBuilder`,
//! `TranscriptGroup`, and `TGValue` (`include/.../quant/EquivalenceClassBuilder.hpp`,
//! `model/TranscriptGroup.hpp`).
//!
//! Building is concurrent (worker threads insert as they map reads); after all
//! reads are processed the builder is [`finish`](EquivalenceClassBuilder::finish)ed
//! into a flat [`CollapsedEqClasses`] for inference.

use std::hash::{Hash, Hasher};
use xxhash_rust::xxh3::Xxh3;

/// The label of an equivalence class: the sorted set of transcript ids a
/// fragment is compatible with, optionally refined by per-transcript
/// range-factorization bins (see [`range_factorize_bins`]). Two fragments with
/// the same transcript set but different weight-bin patterns form distinct
/// classes, which sharpens the inference for multimappers (salmon's
/// `--rangeFactorizationBins`).
#[derive(Debug, Clone)]
pub struct TranscriptGroup {
    /// transcript ids, ascending and deduplicated
    pub txps: Vec<u32>,
    /// range-factorization bins, aligned to `txps`; empty when factorization is
    /// off. Participates in identity/hashing but **not** in the EM (which only
    /// uses `txps` + weights).
    pub bins: Vec<u32>,
    /// precomputed hash of the full label (`txps` + `bins`)
    hash: u64,
    /// EM may invalidate a degenerate class
    pub valid: bool,
}

impl TranscriptGroup {
    /// Build a group from transcript ids. The ids are sorted and deduplicated;
    /// the hash is computed once. No range factorization.
    pub fn new(mut txps: Vec<u32>) -> Self {
        txps.sort_unstable();
        txps.dedup();
        let hash = Self::hash_label(&txps, &[]);
        Self {
            txps,
            bins: Vec::new(),
            hash,
            valid: true,
        }
    }

    /// Build a group from already-sorted, deduplicated ids (no re-sorting), no
    /// range factorization.
    pub fn from_sorted(txps: Vec<u32>) -> Self {
        debug_assert!(txps.windows(2).all(|w| w[0] < w[1]), "txps must be strictly increasing");
        let hash = Self::hash_label(&txps, &[]);
        Self {
            txps,
            bins: Vec::new(),
            hash,
            valid: true,
        }
    }

    /// Build a range-factorized group: sorted ids plus per-transcript bins
    /// (same length, aligned to `txps`).
    pub fn with_bins(txps: Vec<u32>, bins: Vec<u32>) -> Self {
        debug_assert!(txps.windows(2).all(|w| w[0] < w[1]), "txps must be strictly increasing");
        debug_assert_eq!(txps.len(), bins.len(), "bins must align with txps");
        let hash = Self::hash_label(&txps, &bins);
        Self {
            txps,
            bins,
            hash,
            valid: true,
        }
    }

    fn hash_label(txps: &[u32], bins: &[u32]) -> u64 {
        // Hash the raw little-endian bytes of both the ids and the bins so that
        // groups differing only in their bin pattern hash differently.
        let mut h = Xxh3::new();
        h.update(bytemuck_cast(txps));
        h.update(bytemuck_cast(bins));
        h.digest()
    }

    pub fn len(&self) -> usize {
        self.txps.len()
    }
    pub fn is_empty(&self) -> bool {
        self.txps.is_empty()
    }
}

/// Compute range-factorization bins for one fragment's per-transcript
/// conditional probabilities. `weights` must be normalized to sum to 1.
///
/// Mirrors salmon: `rangeCount = floor(sqrt(n)) + bins`, then each weight maps
/// to `floor(weight * rangeCount)`. With `bins == 0` factorization is disabled
/// and an empty vector is returned.
pub fn range_factorize_bins(weights: &[f64], bins: u32) -> Vec<u32> {
    if bins == 0 {
        return Vec::new();
    }
    let n = weights.len();
    let range_count = (n as f64).sqrt() as i64 + bins as i64;
    weights
        .iter()
        .map(|&w| ((w * range_count as f64) as i64).max(0) as u32)
        .collect()
}

/// Reinterpret a `&[u32]` as `&[u8]` without an extra dependency.
fn bytemuck_cast(s: &[u32]) -> &[u8] {
    // SAFETY: u32 has no padding and any bit pattern is valid for u8; the
    // resulting slice covers exactly the same bytes and lives as long as `s`.
    #[allow(unsafe_code)]
    unsafe {
        std::slice::from_raw_parts(s.as_ptr() as *const u8, std::mem::size_of_val(s))
    }
}

impl PartialEq for TranscriptGroup {
    fn eq(&self, other: &Self) -> bool {
        self.hash == other.hash && self.txps == other.txps && self.bins == other.bins
    }
}
impl Eq for TranscriptGroup {}

impl Hash for TranscriptGroup {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.hash);
    }
}

/// The aggregated value for an equivalence class.
#[derive(Debug, Clone)]
pub struct TGValue {
    /// per-transcript auxiliary weights (bias / positional / orientation), in
    /// the same order as the group's `txps`
    pub weights: Vec<f64>,
    /// `weights[i] / effLen(txps[i])`, filled by [`CollapsedEqClasses::update_eff_lengths`]
    pub combined_weights: Vec<f64>,
    /// number of fragments assigned to this class
    pub count: u64,
}

impl TGValue {
    pub fn new(weights: Vec<f64>, count: u64) -> Self {
        let combined_weights = weights.clone();
        Self {
            weights,
            combined_weights,
            count,
        }
    }

    /// Accumulate another observation's weights and count into this class.
    fn accumulate(&mut self, weights: &[f64], count: u64) {
        debug_assert_eq!(self.weights.len(), weights.len());
        for (w, &add) in self.weights.iter_mut().zip(weights) {
            *w += add;
        }
        self.count += count;
    }
}

/// Concurrent equivalence-class builder backed by `scc::HashMap`.
///
/// Worker threads call [`add_group`](Self::add_group) with `&self` while
/// mapping reads; insertion is lock-free per bucket.
pub struct EquivalenceClassBuilder {
    map: scc::HashMap<TranscriptGroup, TGValue>,
}

impl Default for EquivalenceClassBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl EquivalenceClassBuilder {
    pub fn new() -> Self {
        Self {
            map: scc::HashMap::new(),
        }
    }

    /// Add one fragment's worth of evidence: the transcript group, its
    /// per-transcript weights, and a count (usually 1). Thread-safe.
    pub fn add_group(&self, group: TranscriptGroup, weights: Vec<f64>, count: u64) {
        debug_assert_eq!(group.txps.len(), weights.len());
        self.map
            .entry(group)
            .and_modify(|v| v.accumulate(&weights, count))
            .or_insert_with(|| TGValue::new(weights, count));
    }

    /// Number of distinct equivalence classes accumulated so far.
    pub fn len(&self) -> usize {
        self.map.len()
    }
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Finalize into a flat, index-stable collection for inference. Classes are
    /// sorted by their transcript label for determinism.
    pub fn finish(self) -> CollapsedEqClasses {
        let mut classes: Vec<(TranscriptGroup, TGValue)> = Vec::with_capacity(self.map.len());
        self.map.scan(|k, v| classes.push((k.clone(), v.clone())));
        classes.sort_by(|a, b| a.0.txps.cmp(&b.0.txps));
        let total_count = classes.iter().map(|(_, v)| v.count).sum();
        CollapsedEqClasses {
            classes,
            total_count,
        }
    }
}

/// A finalized, flat set of equivalence classes ready for the EM optimizer.
#[derive(Debug, Default)]
pub struct CollapsedEqClasses {
    pub classes: Vec<(TranscriptGroup, TGValue)>,
    pub total_count: u64,
}

impl CollapsedEqClasses {
    pub fn len(&self) -> usize {
        self.classes.len()
    }
    pub fn is_empty(&self) -> bool {
        self.classes.is_empty()
    }

    /// Recompute `combined_weights[i] = weights[i] / effLen(txps[i])` for every
    /// class, mirroring salmon's `updateEqClassWeights`. Called whenever
    /// effective lengths change (initial setup and after bias retraining).
    ///
    /// `eff_lengths` is indexed by transcript id and must be in linear (not
    /// log) space.
    pub fn update_eff_lengths(&mut self, eff_lengths: &[f64]) {
        for (group, value) in &mut self.classes {
            let mut sum = 0.0;
            for (i, &tid) in group.txps.iter().enumerate() {
                let el = eff_lengths[tid as usize].max(1.0);
                let w = value.weights[i] / el;
                value.combined_weights[i] = w;
                sum += w;
            }
            // Normalize per class to sum to 1 (does not affect the EM fixpoint,
            // but matches salmon and keeps the weights well-scaled).
            if sum > 0.0 {
                for w in &mut value.combined_weights {
                    *w /= sum;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn group_normalizes_and_hashes_consistently() {
        let a = TranscriptGroup::new(vec![3, 1, 2, 1]);
        assert_eq!(a.txps, vec![1, 2, 3]);
        let b = TranscriptGroup::new(vec![1, 2, 3]);
        assert_eq!(a, b);
        // hashes equal for equal groups
        let mut ha = std::collections::hash_map::DefaultHasher::new();
        let mut hb = std::collections::hash_map::DefaultHasher::new();
        a.hash(&mut ha);
        b.hash(&mut hb);
        assert_eq!(std::hash::Hasher::finish(&ha), std::hash::Hasher::finish(&hb));
    }

    #[test]
    fn builder_aggregates_counts_and_weights() {
        let b = EquivalenceClassBuilder::new();
        b.add_group(TranscriptGroup::new(vec![1, 2]), vec![0.5, 0.5], 1);
        b.add_group(TranscriptGroup::new(vec![2, 1]), vec![0.5, 0.5], 1);
        b.add_group(TranscriptGroup::new(vec![3]), vec![1.0], 1);
        assert_eq!(b.len(), 2);

        let collapsed = b.finish();
        assert_eq!(collapsed.total_count, 3);
        // sorted: [1,2] then [3]
        let (g0, v0) = &collapsed.classes[0];
        assert_eq!(g0.txps, vec![1, 2]);
        assert_eq!(v0.count, 2);
        assert_eq!(v0.weights, vec![1.0, 1.0]); // accumulated
        let (g1, v1) = &collapsed.classes[1];
        assert_eq!(g1.txps, vec![3]);
        assert_eq!(v1.count, 1);
    }

    #[test]
    fn range_factorization_disabled_returns_empty() {
        assert!(range_factorize_bins(&[0.5, 0.5], 0).is_empty());
    }

    #[test]
    fn range_factorization_bins_match_salmon_formula() {
        // n=2 -> range_count = floor(sqrt(2)) + 4 = 1 + 4 = 5
        // weights 0.8 / 0.2 -> floor(0.8*5)=4, floor(0.2*5)=1
        let bins = range_factorize_bins(&[0.8, 0.2], 4);
        assert_eq!(bins, vec![4, 1]);
    }

    #[test]
    fn bins_split_classes_with_same_transcripts() {
        // Same transcript set {0,1}, but different weight shapes -> two classes.
        let b = EquivalenceClassBuilder::new();
        let txps = vec![0u32, 1];
        let g1 = TranscriptGroup::with_bins(txps.clone(), range_factorize_bins(&[0.9, 0.1], 4));
        let g2 = TranscriptGroup::with_bins(txps.clone(), range_factorize_bins(&[0.5, 0.5], 4));
        assert_ne!(g1, g2, "different weight shapes must differ");
        b.add_group(g1, vec![0.9, 0.1], 1);
        b.add_group(g2, vec![0.5, 0.5], 1);
        // a third fragment matching g1's shape merges with it
        let g3 = TranscriptGroup::with_bins(txps, range_factorize_bins(&[0.92, 0.08], 4));
        b.add_group(g3, vec![0.92, 0.08], 1);
        assert_eq!(b.len(), 2, "expected 2 factorized classes");
    }

    #[test]
    fn combined_weights_divide_by_eff_length_and_normalize() {
        let b = EquivalenceClassBuilder::new();
        b.add_group(TranscriptGroup::new(vec![0, 1]), vec![1.0, 1.0], 4);
        let mut collapsed = b.finish();
        // eff lengths 100 and 300 -> raw combined [1/100, 1/300], normalized.
        collapsed.update_eff_lengths(&[100.0, 300.0]);
        let cw = &collapsed.classes[0].1.combined_weights;
        let s: f64 = cw.iter().sum();
        assert!((s - 1.0).abs() < 1e-12);
        // ratio should be 3:1 in favor of the shorter transcript
        assert!((cw[0] / cw[1] - 3.0).abs() < 1e-9);
    }
}
