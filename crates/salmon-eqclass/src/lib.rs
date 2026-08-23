//! Equivalence classes for salmon-style quantification.
//!
//! # What this crate is for
//!
//! RNA-seq gives us short DNA fragments ("reads") sequenced from a soup of RNA
//! molecules ("transcripts"). We want to know how much of each transcript was
//! in the soup. The hard part is that transcripts share sequence: one fragment
//! is often compatible with many transcripts, so we cannot simply count "this
//! fragment came from transcript X".
//!
//! The standard trick is to group fragments that carry the *same information*.
//! A fragment compatible with transcripts {3, 7, 12} tells the statistical
//! model exactly the same thing as another fragment compatible with {3, 7, 12}.
//! So instead of storing millions of fragments we store one *equivalence class*
//! per distinct compatible-transcript set, together with:
//!
//! * how many fragments landed in that class (a count), and
//! * a per-transcript weight capturing how *likely* each transcript in the set
//!   is to have produced such a fragment (from fragment length, sequence bias,
//!   position bias, orientation).
//!
//! The EM/VBEM optimizer downstream then iterates over a few hundred thousand
//! classes rather than hundreds of millions of fragments. This is the single
//! biggest reason salmon is fast.
//!
//! This mirrors salmon's C++ `EquivalenceClassBuilder`, `TranscriptGroup`, and
//! `TGValue` (`include/.../quant/EquivalenceClassBuilder.hpp`,
//! `model/TranscriptGroup.hpp`).
//!
//! # Lifecycle
//!
//! 1. Many worker threads map reads in parallel and call
//!    [`EquivalenceClassBuilder::add_group`] concurrently. The builder is a
//!    concurrent hash map, so no thread has to wait for a global lock.
//! 2. When every read has been processed,
//!    [`finish`](EquivalenceClassBuilder::finish) converts the map into a flat,
//!    sorted [`CollapsedEqClasses`] vector. Flat + sorted matters: the EM
//!    indexes classes by position, and sorting makes those positions the same
//!    on every run regardless of thread scheduling.
//! 3. [`CollapsedEqClasses::update_eff_lengths`] folds transcript effective
//!    lengths into the weights, and the optimizer takes it from there.

// `Hash`/`Hasher` are the Rust standard-library traits ("interfaces") that let a
// type be used as a hash-map key. We implement them by hand below so that a
// transcript group hashes by a digest we precompute once.
use std::hash::{Hash, Hasher};
// xxh3 is a very fast non-cryptographic hash. We only need speed and good bit
// mixing here, not any security property.
use xxhash_rust::xxh3::Xxh3;

/// The label ("key") of an equivalence class.
///
/// It is the sorted set of transcript ids a fragment is compatible with,
/// optionally refined by per-transcript range-factorization bins (see
/// [`range_factorize_bins`]).
///
/// **Why the bins are part of the key.** Two fragments can be compatible with
/// exactly the same transcripts and yet be very different evidence: one may be
/// 90% likely to come from transcript A, the other 50/50. Collapsing them into
/// one class would throw that away. Range factorization discretizes each
/// fragment's probability vector into coarse bins and makes the bin pattern
/// part of the class label, so those two fragments form two distinct classes.
/// This sharpens inference for multi-mapping reads (salmon's
/// `--rangeFactorizationBins`).
#[derive(Debug, Clone)]
pub struct TranscriptGroup {
    /// Transcript ids, ascending and deduplicated. Sorting is what makes the
    /// label canonical: the sets {7, 3} and {3, 7} must be the *same* class, so
    /// we always store them in the same order.
    pub txps: Vec<u32>,
    /// Range-factorization bins, aligned to `txps` (`bins[i]` describes
    /// `txps[i]`); empty when factorization is off. Participates in
    /// identity/hashing but **not** in the EM, which only uses `txps` plus the
    /// weights stored in [`TGValue`].
    pub bins: Vec<u32>,
    /// Precomputed hash of the full label (`txps` + `bins`).
    ///
    /// Computed once at construction rather than on every hash-map lookup: a
    /// group can hold hundreds of ids, and rehashing that vector on each of
    /// hundreds of millions of lookups would dominate the runtime.
    hash: u64,
}

impl TranscriptGroup {
    /// Build a group from arbitrary transcript ids.
    ///
    /// The ids are sorted and deduplicated here, so callers may pass them in
    /// any order and with repeats. No range factorization: `bins` stays empty.
    pub fn new(mut txps: Vec<u32>) -> Self {
        // `sort_unstable` is the faster sort; "unstable" only means equal
        // elements may be reordered, which is irrelevant for plain integers.
        txps.sort_unstable();
        // `dedup` removes *consecutive* duplicates, which after sorting means
        // all duplicates.
        txps.dedup();
        let hash = Self::hash_label(&txps, &[]);
        Self {
            txps,
            bins: Vec::new(),
            hash,
        }
    }

    /// Build a group from ids the caller has already sorted and deduplicated,
    /// skipping the sort. No range factorization.
    ///
    /// The mapping code usually produces ids in ascending order already, and
    /// this is on the hottest path in the program, so the redundant sort is
    /// worth avoiding.
    pub fn from_sorted(txps: Vec<u32>) -> Self {
        // `debug_assert!` checks the caller kept their promise, but only in
        // debug/test builds — it compiles to nothing in a release build, so the
        // check costs release users nothing. `windows(2)` walks adjacent pairs.
        debug_assert!(
            txps.windows(2).all(|w| w[0] < w[1]),
            "txps must be strictly increasing"
        );
        let hash = Self::hash_label(&txps, &[]);
        Self {
            txps,
            bins: Vec::new(),
            hash,
        }
    }

    /// Build a range-factorized group: sorted ids plus per-transcript bins of
    /// the same length, aligned to `txps`.
    pub fn with_bins(txps: Vec<u32>, bins: Vec<u32>) -> Self {
        debug_assert!(
            txps.windows(2).all(|w| w[0] < w[1]),
            "txps must be strictly increasing"
        );
        debug_assert_eq!(txps.len(), bins.len(), "bins must align with txps");
        let hash = Self::hash_label(&txps, &bins);
        Self { txps, bins, hash }
    }

    /// Hash the full label into a single 64-bit digest.
    fn hash_label(txps: &[u32], bins: &[u32]) -> u64 {
        // Hash the raw little-endian bytes of both the ids and the bins so that
        // groups differing only in their bin pattern hash differently. Feeding
        // both arrays into one hasher (rather than combining two digests) means
        // the bin pattern is mixed through the whole digest.
        let mut h = Xxh3::new();
        h.update(bytemuck_cast(txps));
        h.update(bytemuck_cast(bins));
        h.digest()
    }

    /// How many transcripts are in this class label.
    pub fn len(&self) -> usize {
        self.txps.len()
    }
    /// Whether the label is empty (no compatible transcript). Rust's linter
    /// asks for this whenever a type has `len`.
    pub fn is_empty(&self) -> bool {
        self.txps.is_empty()
    }
}

/// Compute range-factorization bins for one fragment's per-transcript
/// conditional probabilities. `weights` must be normalized to sum to 1.
///
/// **What it does.** Turns a vector of probabilities like `[0.8, 0.2]` into a
/// vector of small integers like `[4, 1]`. Fragments whose probability vectors
/// fall in the same bins are treated as the same evidence and merge into one
/// equivalence class; fragments that differ enough to land in different bins
/// stay separate.
///
/// **How.** Mirrors salmon exactly: `rangeCount = floor(sqrt(n)) + bins`, then
/// each weight maps to `floor(weight * rangeCount)`. The `sqrt(n)` term makes
/// the resolution grow with the size of the transcript set, so large ambiguous
/// sets are not all crushed into the same handful of bins.
///
/// **Why bin at all** instead of keeping exact probabilities: exact vectors are
/// almost never equal, so nothing would ever collapse and the class count would
/// approach the fragment count, defeating the purpose.
///
/// With `bins == 0` factorization is disabled and an empty vector is returned.
pub fn range_factorize_bins(weights: &[f64], bins: u32) -> Vec<u32> {
    if bins == 0 {
        return Vec::new();
    }
    let n = weights.len();
    // `as i64` on a float truncates toward zero, which for a non-negative
    // sqrt is the same as floor.
    let range_count = (n as f64).sqrt() as i64 + bins as i64;
    weights
        .iter()
        // `.max(0)` guards against a tiny negative weight arriving from
        // floating-point noise; a negative bin index would be meaningless.
        .map(|&w| ((w * range_count as f64) as i64).max(0) as u32)
        .collect()
}

/// Reinterpret a `&[u32]` as `&[u8]` (view the same memory as raw bytes)
/// without pulling in an extra dependency.
///
/// The hasher wants bytes; our data is 32-bit integers. Copying them into a
/// byte buffer first would allocate on the hottest path, so we hand the hasher
/// a view of the same memory instead.
fn bytemuck_cast(s: &[u32]) -> &[u8] {
    // SAFETY: u32 has no padding and any bit pattern is valid for u8; the
    // resulting slice covers exactly the same bytes and lives as long as `s`.
    #[allow(unsafe_code)]
    unsafe {
        std::slice::from_raw_parts(s.as_ptr() as *const u8, std::mem::size_of_val(s))
    }
}

// Two groups are equal when their labels are equal. We compare the cheap
// precomputed hash first: unequal groups almost always differ there, so the
// expensive vector comparisons are skipped. The vector comparisons are still
// required for correctness, because two different labels can in principle
// produce the same 64-bit hash.
impl PartialEq for TranscriptGroup {
    fn eq(&self, other: &Self) -> bool {
        self.hash == other.hash && self.txps == other.txps && self.bins == other.bins
    }
}
// `Eq` adds no methods; it is a promise that `eq` is a proper equality
// (reflexive, symmetric, transitive), which hash maps rely on.
impl Eq for TranscriptGroup {}

// Hashing just replays the digest we computed at construction time.
impl Hash for TranscriptGroup {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.hash);
    }
}

/// A passthrough `Hasher` for keys that already carry a precomputed, well-mixed
/// hash.
///
/// **Why this exists.** Rust's default hash map runs every key through SipHash,
/// a hash chosen to resist denial-of-service attacks from hostile input. Our
/// keys are transcript ids from our own index, not hostile input, and
/// [`TranscriptGroup::hash`] already writes a single well-mixed `u64` (its xxh3
/// digest). Re-hashing that digest is pure overhead: profiling on human-scale
/// data showed SipHash (`DefaultHasher::write`) eating ~9% of the quant
/// runtime. This hasher simply returns the `u64` that was written.
#[derive(Default, Clone, Copy)]
pub struct IdentityHasher(u64);

impl Hasher for IdentityHasher {
    /// Return the digest we were handed.
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    /// The only path we actually use: store the digest verbatim.
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.0 = i;
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        // Fallback for any non-`u64` writes (none on the hot path): fold the
        // bytes in so distinct inputs still differ. Rotating before XOR-ing
        // keeps byte order significant, so "ab" and "ba" do not collide.
        for &b in bytes {
            self.0 = self.0.rotate_left(8) ^ u64::from(b);
        }
    }
}

/// `BuildHasher` producing [`IdentityHasher`]s.
///
/// A hash map needs a *factory* for hashers (one fresh hasher per key), not a
/// hasher; this zero-sized type is that factory.
#[derive(Default, Clone, Copy)]
pub struct IdentityBuildHasher;

impl std::hash::BuildHasher for IdentityBuildHasher {
    type Hasher = IdentityHasher;
    #[inline]
    fn build_hasher(&self) -> IdentityHasher {
        IdentityHasher(0)
    }
}

/// Scale for fixed-point accumulation of per-class weights.
///
/// **The problem.** Each fragment's per-transcript weight is a conditional
/// probability in `[0, 1]`, and a class sums the weights of every fragment that
/// landed in it. Floating-point addition is *not* associative: `(a + b) + c` can
/// differ in the last bits from `a + (b + c)`. Fragments arrive in an order that
/// depends on thread scheduling and (in RAD mode) chunk order, so a plain `f64`
/// sum would differ slightly between runs. Worse, those last bits can push a
/// weight across a range-factorization bin boundary and change the
/// equivalence-class *set* itself, so results would not reproduce.
///
/// **The fix.** Multiply each weight by 2^40, round to an integer, and sum as
/// integers. Integer addition *is* associative, so the total is identical
/// whatever the arrival order. 2^40 (~1.1e12) gives ~1e-12 weight resolution,
/// finer than range factorization needs, and `u128` accumulation cannot
/// overflow for any realistic fragment count.
const WEIGHT_SCALE: f64 = (1u64 << 40) as f64;

/// Convert one probability into the fixed-point representation described on
/// [`WEIGHT_SCALE`]. `.max(0.0)` clamps floating-point noise below zero.
#[inline]
fn weight_to_fixed(w: f64) -> u128 {
    (w.max(0.0) * WEIGHT_SCALE).round() as u128
}

/// The aggregated value for an equivalence class: how many fragments landed in
/// it and how much each transcript in the label is credited.
#[derive(Debug, Clone)]
pub struct TGValue {
    /// Order-independent fixed-point accumulator of the summed per-transcript
    /// weights (integer addition is associative; see [`WEIGHT_SCALE`]).
    /// Materialized into `weights` at [`EquivalenceClassBuilder::finish`] and
    /// then emptied to release the memory.
    acc: Vec<u128>,
    /// Per-transcript auxiliary weights (bias / positional / orientation), in
    /// the same order as the group's `txps`. Materialized from `acc` at finish,
    /// so this is empty while the builder is still accumulating.
    pub weights: Vec<f64>,
    /// `weights[i] / effLen(txps[i])`, filled by
    /// [`CollapsedEqClasses::update_eff_lengths`]. This is what the EM actually
    /// reads: dividing by effective length turns "how likely is this fragment
    /// given the transcript" into "how much abundance does it imply".
    pub combined_weights: Vec<f64>,
    /// Number of fragments assigned to this class.
    pub count: u64,
}

impl TGValue {
    /// Seed a brand-new class from its first observation.
    pub fn new(weights: Vec<f64>, count: u64) -> Self {
        let acc = weights.iter().map(|&w| weight_to_fixed(w)).collect();
        Self {
            acc,
            // `weights`/`combined_weights` stay empty until `materialize`; only
            // `acc` is live during accumulation.
            weights: Vec::new(),
            combined_weights: Vec::new(),
            count,
        }
    }

    /// Accumulate another observation's weights and count into this class.
    fn accumulate(&mut self, weights: &[f64], count: u64) {
        debug_assert_eq!(self.acc.len(), weights.len());
        // `zip` walks the accumulator and the incoming weights in lockstep;
        // both are indexed by position in the class label.
        for (a, &add) in self.acc.iter_mut().zip(weights) {
            *a += weight_to_fixed(add);
        }
        self.count += count;
    }

    /// Convert the fixed-point accumulator into the f64 `weights` (and seed
    /// `combined_weights` with the same values), then drop the accumulator.
    /// Called exactly once, at finish, when no more fragments can arrive.
    fn materialize(&mut self) {
        self.weights = self.acc.iter().map(|&a| a as f64 / WEIGHT_SCALE).collect();
        self.combined_weights = self.weights.clone();
        // Replacing with an empty vector frees the u128 buffer, which is twice
        // the size of the f64 one we just built.
        self.acc = Vec::new();
    }
}

/// Concurrent equivalence-class builder backed by a sharded concurrent hash map.
///
/// Worker threads call [`add_group`](Self::add_group) with `&self` (a shared,
/// non-exclusive reference) while mapping reads. `dashmap` splits the table into
/// many independently locked shards, so two threads touching different classes
/// never block each other.
pub struct EquivalenceClassBuilder {
    map: dashmap::DashMap<TranscriptGroup, TGValue, IdentityBuildHasher>,
}

// `Default` lets callers write `EquivalenceClassBuilder::default()`; Rust's
// linter asks for it whenever a type has a no-argument `new`.
impl Default for EquivalenceClassBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl EquivalenceClassBuilder {
    /// Create an empty builder using the passthrough hasher.
    pub fn new() -> Self {
        Self {
            map: dashmap::DashMap::with_hasher(IdentityBuildHasher),
        }
    }

    /// Add one fragment's worth of evidence: the transcript group, its
    /// per-transcript weights, and a count (usually 1). Thread-safe.
    pub fn add_group(&self, group: TranscriptGroup, weights: Vec<f64>, count: u64) {
        debug_assert_eq!(group.txps.len(), weights.len());
        // One map operation does both cases: if the class already exists add to
        // it (`and_modify`), otherwise create it (`or_insert_with`). Doing it in
        // one call holds the shard lock once and avoids the race where two
        // threads both see "absent" and both insert.
        self.map
            .entry(group)
            .and_modify(|v| v.accumulate(&weights, count))
            .or_insert_with(|| TGValue::new(weights, count));
    }

    /// Number of distinct equivalence classes accumulated so far.
    pub fn len(&self) -> usize {
        self.map.len()
    }
    /// Whether no class has been recorded yet.
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Finalize into a flat, index-stable collection for inference.
    ///
    /// Takes `self` by value (consuming the builder) because nothing may be
    /// added afterwards. Classes are sorted by their transcript label so that
    /// class *indices* — which the EM uses everywhere — are identical on every
    /// run, no matter what order the threads happened to insert them in.
    pub fn finish(self) -> CollapsedEqClasses {
        // Move the entries out of the concurrent map rather than cloning them.
        // `self` is owned here, so the map is dropped either way; cloning would
        // deep-copy every `txps`/`bins`/`acc` vector and hold two full copies of
        // the largest structure in the run at once.
        //
        // Reserve first: `dashmap`'s owning iterator implements neither
        // `size_hint` nor `ExactSizeIterator`, so `collect` would report `(0,
        // None)` and grow this vector from empty by doubling — a realloc and a
        // memcpy of every element at each step, plus a transient ~1.5x peak at
        // the last one. At 144 bytes per entry that is the same cost this change
        // removes from `PackedEqClasses::from_collapsed`. `len` has to be read
        // before `into_iter` consumes the map.
        //
        // Iteration order of a hash map is arbitrary, hence the sort below.
        let num_classes = self.map.len();
        let mut classes: Vec<(TranscriptGroup, TGValue)> = Vec::with_capacity(num_classes);
        classes.extend(self.map);
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
///
/// "Collapsed" is salmon's term for the post-aggregation form: millions of
/// fragments have been collapsed into this comparatively short list.
#[derive(Debug, Default)]
pub struct CollapsedEqClasses {
    /// One entry per class: its label and its aggregated value. The position of
    /// an entry in this vector is the class id used by the optimizer.
    pub classes: Vec<(TranscriptGroup, TGValue)>,
    /// Total number of fragments across all classes.
    pub total_count: u64,
}

impl CollapsedEqClasses {
    /// Number of classes.
    pub fn len(&self) -> usize {
        self.classes.len()
    }
    /// Whether there are no classes at all (nothing mapped).
    pub fn is_empty(&self) -> bool {
        self.classes.is_empty()
    }

    /// Recompute `combined_weights[i] = weights[i] / effLen(txps[i])` for every
    /// class, mirroring salmon's `updateEqClassWeights`.
    ///
    /// **Why divide by effective length.** A long transcript offers more places
    /// for a fragment to start, so at equal molar abundance it yields more
    /// fragments. Dividing the compatibility weight by the transcript's
    /// effective length converts "fragments observed" into "molecules implied",
    /// which is the quantity the EM estimates.
    ///
    /// Called whenever effective lengths change: at initial setup, and again
    /// after the bias models are retrained and shift the effective lengths.
    ///
    /// `eff_lengths` is indexed by transcript id and must be in linear (not
    /// log) space.
    pub fn update_eff_lengths(&mut self, eff_lengths: &[f64]) {
        for (group, value) in &mut self.classes {
            let mut sum = 0.0;
            // `enumerate` gives the position `i` alongside the transcript id, so
            // we can index the parallel `weights`/`combined_weights` arrays.
            for (i, &tid) in group.txps.iter().enumerate() {
                // Clamp at 1.0: a transcript shorter than the typical fragment
                // can get an effective length at or below zero, and dividing by
                // that would produce infinities.
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
///
/// A sentinel value is used instead of `Option<u8>` because these placements are
/// stored in the millions and the extra byte per placement would matter.
pub const NAIVE_NO_FMT: u8 = u8::MAX;

/// Write `aux_info/eq_classes.txt.gz` in salmon's format: the transcript count,
/// the equivalence-class count, all transcript names (one per line, in the index
/// order the classes reference), then one line per class.
///
/// Without `with_weights`, classes are collapsed by transcript set (`groupSize`,
/// tids, summed count), matching salmon's `--dumpEq`; with `with_weights` each
/// class is emitted as-is with its per-transcript combined weights before the
/// count (`--dumpEqWeights`).
///
/// This is the intermediate the EM actually consumes, so dumping it is how to
/// see what salmon inferred *before* the statistics: which transcript sets the
/// fragments were compatible with, and how many fell into each.
///
/// It lives here, beside [`CollapsedEqClasses`], because both the read-mapping
/// path and the RAD path have to produce the same bytes for the same classes.
/// When it lived in one of them, the other could not dump at all
/// (COMBINE-lab/salmon#1140).
pub fn write_eq_classes(
    dir: &std::path::Path,
    names: &[String],
    collapsed: &CollapsedEqClasses,
    with_weights: bool,
) -> std::io::Result<()> {
    use std::collections::BTreeMap;
    use std::io::Write as _;
    let num_txps = names.len();
    std::fs::create_dir_all(dir.join("aux_info"))?;
    let f = std::fs::File::create(dir.join("aux_info").join("eq_classes.txt.gz"))?;
    // One `write!` per field per class is a lot of small writes; buffer them
    // before they reach the encoder.
    let f = std::io::BufWriter::new(f);
    let mut w = flate2::write::GzEncoder::new(f, flate2::Compression::new(6));

    if with_weights {
        writeln!(w, "{num_txps}")?;
        writeln!(w, "{}", collapsed.classes.len())?;
        for name in names {
            writeln!(w, "{name}")?;
        }
        for (group, value) in &collapsed.classes {
            write!(w, "{}", group.txps.len())?;
            for &tid in &group.txps {
                write!(w, "\t{tid}")?;
            }
            for wt in &value.combined_weights {
                write!(w, "\t{wt}")?;
            }
            writeln!(w, "\t{}", value.count)?;
        }
    } else {
        // Collapse range-factorized sub-classes that share a transcript set.
        let mut merged: BTreeMap<&Vec<u32>, u64> = BTreeMap::new();
        for (group, value) in &collapsed.classes {
            *merged.entry(&group.txps).or_insert(0) += value.count;
        }
        writeln!(w, "{num_txps}")?;
        writeln!(w, "{}", merged.len())?;
        for name in names {
            writeln!(w, "{name}")?;
        }
        for (txps, count) in &merged {
            write!(w, "{}", txps.len())?;
            for &tid in *txps {
                write!(w, "\t{tid}")?;
            }
            writeln!(w, "\t{count}")?;
        }
    }
    w.finish()?;
    Ok(())
}

/// One placement of a NAIVE fragment signature: a transcript plus enough
/// orientation/strand info to drop it later if it is library-incompatible.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct NaivePlacement {
    /// Transcript this fragment could have come from.
    pub tid: u32,
    /// Observed paired-format id ([`salmon_core::LibraryFormat::format_id`]), or
    /// [`NAIVE_NO_FMT`] for orphan/single-end (filtered by `is_fw`/`status`).
    pub fmt_id: u8,
    /// Whether the fragment aligned to the forward strand of the transcript.
    pub is_fw: bool,
    /// Which mate(s) of the pair were placed (both, left only, right only).
    pub status: salmon_core::MateStatus,
}

/// Concurrent accumulator of orientation-tagged **naive** (weightless) fragment
/// signatures, for the rough seed EM that produces bias-weighting abundances
/// when none are baked.
///
/// **Why a second, cruder accumulator exists.** Bias correction needs to know
/// roughly how abundant each transcript is *before* it can weight the bias
/// observations. But the real equivalence classes cannot be built yet, because
/// their weights depend on the fragment-length distribution, which is still
/// being learned. So salmon first builds these naive classes: no weights, no
/// strand filter, just "which transcripts was this fragment compatible with",
/// runs a rough EM on them, and uses the result as bias weights.
///
/// Each placement keeps its observed orientation so that once the library type
/// is determined (after the FLD pass) library-incompatible placements can be
/// **dropped** in [`finish`](Self::finish) — the rough seed treats
/// `incompat_prior` as 0 — before forming the uniform eq-classes the rough EM
/// consumes. The output is built via a normal [`EquivalenceClassBuilder`], so it
/// is deterministic (sorted finish + integer counts) regardless of this
/// accumulator's hashing/insertion order.
#[derive(Default)]
pub struct NaiveEqBuilder {
    /// Maps a canonical placement signature to how many fragments produced it.
    sigs: dashmap::DashMap<Vec<NaivePlacement>, u64>,
}

impl NaiveEqBuilder {
    /// Create an empty accumulator.
    pub fn new() -> Self {
        Self {
            sigs: dashmap::DashMap::new(),
        }
    }

    /// Record one fragment's placements (the full compatible set, orientation
    /// tagged). Thread-safe; order-independent, because only counts are
    /// aggregated and integer addition commutes.
    pub fn add(&self, mut placements: Vec<NaivePlacement>) {
        if placements.is_empty() {
            return;
        }
        // Canonicalize the signature so identical fragments aggregate: without
        // this, the same evidence arriving in a different order would hash to a
        // different key and split into two entries.
        placements.sort_unstable_by_key(|p| (p.tid, p.fmt_id, p.is_fw));
        *self.sigs.entry(placements).or_insert(0) += 1;
    }

    /// Build uniform-weight equivalence classes, dropping placements that are
    /// incompatible with `expected` (`None` ⇒ keep all, e.g. unresolved `-l A`).
    ///
    /// Incompatible placements are removed rather than down-weighted: for the
    /// rough seed, `incompat_prior` is treated as 0.
    pub fn finish(&self, expected: Option<salmon_core::LibraryFormat>) -> CollapsedEqClasses {
        use salmon_core::LibraryFormat;
        let out = EquivalenceClassBuilder::new();
        for e in self.sigs.iter() {
            let count = *e.value();
            let mut tids: Vec<u32> = e
                .key()
                .iter()
                .filter(|p| {
                    // `is_none_or` means "no expectation set, or the expectation
                    // is satisfied" — with `-l A` unresolved we keep everything.
                    expected.is_none_or(|exp| {
                        // Orphans/single-end carry no pair orientation, so their
                        // observed format is `None` and compatibility falls back
                        // to strand + mate status.
                        let obs = (p.fmt_id != NAIVE_NO_FMT)
                            .then(|| LibraryFormat::from_format_id(p.fmt_id));
                        salmon_core::is_compatible(exp, obs, p.is_fw, p.status)
                    })
                })
                .map(|p| p.tid)
                .collect();
            // One transcript can appear in several surviving placements (e.g.
            // both orientations); collapse to a set.
            tids.sort_unstable();
            tids.dedup();
            let n = tids.len();
            if n > 0 {
                // "Naive" = every compatible transcript is equally credited.
                let w = vec![1.0 / n as f64; n];
                out.add_group(TranscriptGroup::from_sorted(tids), w, count);
            }
        }
        out.finish()
    }
}

#[cfg(test)]
mod tests {

    /// The `--dumpEq` file is what a downstream tool reads instead of re-deriving
    /// the classes, so its two shapes have to stay exactly what salmon promises:
    /// collapsed by transcript set without weights, verbatim with them.
    #[test]
    fn eq_class_dump_collapses_only_without_weights() {
        let dir = tempfile::tempdir().unwrap();
        let names = vec!["txA".to_string(), "txB".to_string(), "txC".to_string()];
        // Two classes over the same transcript set (as range factorization
        // produces) plus one over another set.
        let classes = vec![
            (
                TranscriptGroup::new(vec![0, 1]),
                TGValue {
                    acc: Vec::new(),
                    weights: vec![0.5, 0.5],
                    combined_weights: vec![0.25, 0.75],
                    count: 3,
                },
            ),
            (
                TranscriptGroup::new(vec![0, 1]),
                TGValue {
                    acc: Vec::new(),
                    weights: vec![0.5, 0.5],
                    combined_weights: vec![0.6, 0.4],
                    count: 4,
                },
            ),
            (
                TranscriptGroup::new(vec![2]),
                TGValue {
                    acc: Vec::new(),
                    weights: vec![1.0],
                    combined_weights: vec![1.0],
                    count: 5,
                },
            ),
        ];
        let collapsed = CollapsedEqClasses {
            classes,
            total_count: 12,
        };

        let read = |with_weights: bool| -> Vec<String> {
            write_eq_classes(dir.path(), &names, &collapsed, with_weights).unwrap();
            let f =
                std::fs::File::open(dir.path().join("aux_info").join("eq_classes.txt.gz")).unwrap();
            let mut text = String::new();
            std::io::Read::read_to_string(&mut flate2::read::GzDecoder::new(f), &mut text).unwrap();
            text.lines().map(str::to_string).collect()
        };

        let plain = read(false);
        assert_eq!(plain[0], "3", "transcript count");
        assert_eq!(
            plain[1], "2",
            "the two classes on {{txA,txB}} collapse into one"
        );
        assert_eq!(&plain[2..5], &["txA", "txB", "txC"]);
        // Collapsed: group size, tids, summed count. No weights.
        assert_eq!(plain[5], "2\t0\t1\t7");
        assert_eq!(plain[6], "1\t2\t5");

        let weighted = read(true);
        assert_eq!(weighted[1], "3", "with weights nothing is collapsed");
        assert_eq!(weighted[5], "2\t0\t1\t0.25\t0.75\t3");
        assert_eq!(weighted[6], "2\t0\t1\t0.6\t0.4\t4");
        assert_eq!(weighted[7], "1\t2\t1\t5");
    }
    use super::*;

    /// A group must be order-insensitive and duplicate-insensitive: the same
    /// transcript set built two different ways has to be one class, not two.
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

    /// Repeated observations of the same label must merge into one class with
    /// summed counts and summed weights.
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

    /// `--rangeFactorizationBins 0` must switch factorization off entirely.
    #[test]
    fn range_factorization_disabled_returns_empty() {
        assert!(range_factorize_bins(&[0.5, 0.5], 0).is_empty());
    }

    /// Pin the exact binning arithmetic to salmon's, since a different rounding
    /// rule would silently produce a different class set.
    #[test]
    fn range_factorization_bins_match_salmon_formula() {
        // n=2 -> range_count = floor(sqrt(2)) + 4 = 1 + 4 = 5
        // weights 0.8 / 0.2 -> floor(0.8*5)=4, floor(0.2*5)=1
        let bins = range_factorize_bins(&[0.8, 0.2], 4);
        assert_eq!(bins, vec![4, 1]);
    }

    /// The point of range factorization: same transcripts, different evidence
    /// shape, therefore different classes — while similar shapes still merge.
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

    /// Effective-length division must favour the shorter transcript in exact
    /// inverse proportion, and the per-class weights must end up normalized.
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
