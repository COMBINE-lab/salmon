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
        debug_assert!(
            txps.windows(2).all(|w| w[0] < w[1]),
            "txps must be strictly increasing"
        );
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
        debug_assert!(
            txps.windows(2).all(|w| w[0] < w[1]),
            "txps must be strictly increasing"
        );
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

/// A passthrough `Hasher` for keys that already carry a precomputed, well-mixed
/// hash. [`TranscriptGroup::hash`] writes a single `u64` (its xxh3 digest), so
/// running it through the default SipHash again is pure overhead — profiling on
/// human-scale data showed SipHash (`DefaultHasher::write`) eating ~9% of the
/// quant runtime. This hasher simply returns the `u64` that was written.
#[derive(Default, Clone, Copy)]
pub struct IdentityHasher(u64);

impl Hasher for IdentityHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.0 = i;
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        // Fallback for any non-`u64` writes (none on the hot path): fold the
        // bytes in so distinct inputs still differ.
        for &b in bytes {
            self.0 = self.0.rotate_left(8) ^ u64::from(b);
        }
    }
}

/// `BuildHasher` producing [`IdentityHasher`]s.
#[derive(Default, Clone, Copy)]
pub struct IdentityBuildHasher;

impl std::hash::BuildHasher for IdentityBuildHasher {
    type Hasher = IdentityHasher;
    #[inline]
    fn build_hasher(&self) -> IdentityHasher {
        IdentityHasher(0)
    }
}

/// Scale for fixed-point accumulation of per-class weights. Each fragment's
/// per-transcript weight is a conditional probability in `[0, 1]`; scaling by
/// 2^40 and summing as integers makes the per-class total **independent of the
/// order** fragments were added in. f64 addition is not associative, so the
/// arrival order — which varies with thread count and (in RAD mode) chunk order
/// — would otherwise perturb the sum and, through range-factorization bin
/// boundaries, the equivalence-class set itself. 2^40 (~1.1e12) gives ~1e-12
/// weight resolution, finer than range-factorization needs, and u128
/// accumulation cannot overflow for any realistic fragment count.
const WEIGHT_SCALE: f64 = (1u64 << 40) as f64;

#[inline]
fn weight_to_fixed(w: f64) -> u128 {
    (w.max(0.0) * WEIGHT_SCALE).round() as u128
}

/// The aggregated value for an equivalence class.
#[derive(Debug, Clone)]
pub struct TGValue {
    /// order-independent fixed-point accumulator of the summed per-transcript
    /// weights (integer addition is associative; see [`WEIGHT_SCALE`]).
    /// Materialized into `weights` at [`EquivalenceClassBuilder::finish`] and
    /// then emptied.
    acc: Vec<u128>,
    /// per-transcript auxiliary weights (bias / positional / orientation), in the
    /// same order as the group's `txps`. Materialized from `acc` at finish.
    pub weights: Vec<f64>,
    /// `weights[i] / effLen(txps[i])`, filled by [`CollapsedEqClasses::update_eff_lengths`]
    pub combined_weights: Vec<f64>,
    /// number of fragments assigned to this class
    pub count: u64,
}

impl TGValue {
    pub fn new(weights: Vec<f64>, count: u64) -> Self {
        let acc = weights.iter().map(|&w| weight_to_fixed(w)).collect();
        Self {
            acc,
            weights: Vec::new(),
            combined_weights: Vec::new(),
            count,
        }
    }

    /// Accumulate another observation's weights and count into this class.
    fn accumulate(&mut self, weights: &[f64], count: u64) {
        debug_assert_eq!(self.acc.len(), weights.len());
        for (a, &add) in self.acc.iter_mut().zip(weights) {
            *a += weight_to_fixed(add);
        }
        self.count += count;
    }

    /// Convert the fixed-point accumulator into the f64 `weights` (and seed
    /// `combined_weights`), then drop the accumulator. Called once at finish.
    fn materialize(&mut self) {
        self.weights = self.acc.iter().map(|&a| a as f64 / WEIGHT_SCALE).collect();
        self.combined_weights = self.weights.clone();
        self.acc = Vec::new();
    }
}

/// Concurrent equivalence-class builder backed by `scc::HashMap`.
///
/// Worker threads call [`add_group`](Self::add_group) with `&self` while
/// mapping reads; insertion is lock-free per bucket.
pub struct EquivalenceClassBuilder {
    map: dashmap::DashMap<TranscriptGroup, TGValue, IdentityBuildHasher>,
}

impl Default for EquivalenceClassBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl EquivalenceClassBuilder {
    pub fn new() -> Self {
        Self {
            map: dashmap::DashMap::with_hasher(IdentityBuildHasher),
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
        self.map
            .iter()
            .for_each(|e| classes.push((e.key().clone(), e.value().clone())));
        classes.sort_by(|a, b| a.0.txps.cmp(&b.0.txps));
        // Materialize f64 weights from each class's order-independent fixed-point
        // accumulator (the sort above already makes the class *order* stable).
        for (_, v) in &mut classes {
            v.materialize();
        }
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

/// `fmt_id` sentinel for a placement with no observed paired format (orphan /
/// single-end); compatibility is then decided by strand + status.
pub const NAIVE_NO_FMT: u8 = u8::MAX;

/// One placement of a NAIVE fragment signature: a transcript plus enough
/// orientation/strand info to drop it later if it is library-incompatible.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct NaivePlacement {
    pub tid: u32,
    /// observed paired-format id ([`salmon_core::LibraryFormat::format_id`]), or
    /// [`NAIVE_NO_FMT`] for orphan/single-end (filtered by `is_fw`/`status`).
    pub fmt_id: u8,
    pub is_fw: bool,
    pub status: salmon_core::MateStatus,
}

/// Concurrent accumulator of orientation-tagged **naive** (weightless) fragment
/// signatures, for the rough seed EM that produces bias-weighting abundances when
/// none are baked.
///
/// Naive classes carry no conditional-probability weights (the FLD isn't known
/// when they are built) and no strand filter yet — instead each placement keeps
/// its observed orientation, so that once the library type is determined (after
/// the FLD pass) library-incompatible placements can be **dropped** in
/// [`finish`](Self::finish) (the rough seed treats `incompat_prior` as 0) before
/// forming the uniform eq-classes the rough EM consumes. The output is built via
/// a normal [`EquivalenceClassBuilder`], so it is deterministic (sorted finish +
/// integer counts) regardless of this accumulator's hashing/insertion order.
#[derive(Default)]
pub struct NaiveEqBuilder {
    sigs: dashmap::DashMap<Vec<NaivePlacement>, u64>,
}

impl NaiveEqBuilder {
    pub fn new() -> Self {
        Self {
            sigs: dashmap::DashMap::new(),
        }
    }

    /// Record one fragment's placements (the full compatible set, orientation
    /// tagged). Thread-safe; order-independent (counts are commutative).
    pub fn add(&self, mut placements: Vec<NaivePlacement>) {
        if placements.is_empty() {
            return;
        }
        // Canonicalize the signature so identical fragments aggregate.
        placements.sort_unstable_by_key(|p| (p.tid, p.fmt_id, p.is_fw));
        *self.sigs.entry(placements).or_insert(0) += 1;
    }

    /// Build uniform-weight equivalence classes, dropping placements that are
    /// incompatible with `expected` (`None` ⇒ keep all, e.g. unresolved `-l A`).
    /// Incompatible placements are removed (not down-weighted): for the rough seed
    /// `incompat_prior` is treated as 0.
    pub fn finish(&self, expected: Option<salmon_core::LibraryFormat>) -> CollapsedEqClasses {
        use salmon_core::LibraryFormat;
        let out = EquivalenceClassBuilder::new();
        for e in self.sigs.iter() {
            let count = *e.value();
            let mut tids: Vec<u32> = e
                .key()
                .iter()
                .filter(|p| {
                    expected.is_none_or(|exp| {
                        let obs = (p.fmt_id != NAIVE_NO_FMT)
                            .then(|| LibraryFormat::from_format_id(p.fmt_id));
                        salmon_core::is_compatible(exp, obs, p.is_fw, p.status)
                    })
                })
                .map(|p| p.tid)
                .collect();
            tids.sort_unstable();
            tids.dedup();
            let n = tids.len();
            if n > 0 {
                let w = vec![1.0 / n as f64; n];
                out.add_group(TranscriptGroup::from_sorted(tids), w, count);
            }
        }
        out.finish()
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
        assert_eq!(
            std::hash::Hasher::finish(&ha),
            std::hash::Hasher::finish(&hb)
        );
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
