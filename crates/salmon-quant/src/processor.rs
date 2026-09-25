//! paraseq parallel processor that maps reads and accumulates equivalence
//! classes, fragment lengths, observed library formats, and (for `--seqBias`)
//! observed sequence-bias contexts.
//!
//! # This is the hot loop
//!
//! Everything in this file runs once per sequenced fragment — hundreds of
//! millions of times in a real run — so it is where the run's time goes and why
//! several structures elsewhere in the codebase exist in the shape they do.
//!
//! For each fragment the work is: place the read(s) on the transcriptome, score
//! and filter the placements, then fold the result into whatever the run needs —
//! an equivalence class, the fragment-length distribution, the bias models, the
//! library-type tally, and optionally a SAM/BAM/RAD record.
//!
//! # The sharing model
//!
//! `paraseq` (the FASTQ reader) drives this: it hands each worker thread batches
//! of records and a *clone* of the processor. So the type is split in two.
//!
//! Shared `&'a` references hold the index and thread-safe accumulators; the
//! per-thread scratch (a [`HitSearcher`] and, when bias-correcting, observed
//! sequence-bias models) is rebuilt on `Clone` (paraseq clones the processor
//! per worker) and merged into shared state in `on_thread_complete`.
//!
//! Accumulating into per-thread models and merging once at the end — rather than
//! into a shared model per fragment — is what keeps the workers from contending,
//! and (because the merges are integer sums) is also what makes the models
//! independent of how the fragments happened to be split across threads.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

use paraseq::fastx::RefRecord;
use paraseq::parallel::{PairedParallelProcessor, ParallelProcessor};
use paraseq::Record;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};

use salmon_core::math::{log_add, LOG_0, LOG_1};
use salmon_core::{is_compatible, observed_paired_format, LibraryFormat, MateStatus};
use salmon_eqclass::{
    range_factorize_bins, EquivalenceClassBuilder, NaiveEqBuilder, NaivePlacement, TranscriptGroup,
    NAIVE_NO_FMT,
};
use salmon_index::SalmonIndex;
use salmon_infer::OnlineInference;
use salmon_map::{
    map_read_pair_into, map_read_pair_sketch_into, map_single_read_into,
    map_single_read_sketch_into, MapConfig, ScoredMapping,
};
use salmon_model::seqbias::{SBModel, CONTEXT_LEFT, CONTEXT_LENGTH, CONTEXT_RIGHT};
use salmon_model::{
    gc_desc, FragmentLengthDistribution, GcFragModel, LibraryTypeDetector, SimplePosBias,
    NUM_LENGTH_CLASSES,
};

/// Shared, thread-safe state (all `Copy` references).
///
/// Everything here is either read-only (the index, the configuration) or
/// internally synchronized (atomics, mutexes, concurrent maps), so every worker
/// can hold a copy. Being `Copy` means passing it around costs nothing.
///
/// Most of the fields are flags gating optional work. They are resolved once,
/// before mapping starts, so the per-fragment path tests a boolean rather than
/// re-deriving what this run needs.
#[derive(Clone, Copy)]
pub(crate) struct Shared<'a> {
    /// Thread-time the mapping threads spend mapping, for the thread broker.
    ///
    /// Lives on `Shared` because it is the one thing already handed to every
    /// processor clone, so metering needs no extra plumbing. Published every few
    /// hundred fragments rather than once per batch: a batch takes far longer
    /// than one sampling window, and a per-batch counter would leave windows
    /// reading zero work -- reporting maximum starvation for a thread mapping
    /// flat out. See `thread_broker::BusyMeter`.
    pub busy: &'a thread_broker::BusyMeter,
    pub salmon: &'a SalmonIndex,
    pub eq: &'a EquivalenceClassBuilder,
    pub fld: &'a FragmentLengthDistribution,
    /// library-type detector. `quantify()` always constructs one — not just
    /// under `-l A`: with an explicit `-l` it still observes the orientation
    /// distribution (sampled pre-strand-filter) to power the
    /// `library_type_mismatch` diagnostic, without ever overriding the choice.
    pub detector: Option<&'a LibraryTypeDetector>,
    pub map_cfg: &'a MapConfig,
    pub sketch: bool,
    /// sketch mode: only orphan a pair when the other mate had no matching
    /// k-mers (strict), instead of the default empty-accepted-target rule
    /// (`--sketchStrictOrphans`)
    pub sketch_strict_orphan: bool,
    /// fragments processed without the auxiliary (FLD) model before it is folded
    /// into the online posterior (salmon's `--numPreAuxModelSamples`)
    pub pre_burnin: u64,
    /// discard a fragment mapping to more than this many places (`--maxReadOcc`)
    pub max_read_occ: usize,
    pub skip: SkippingStrategy,
    /// range-factorization bins (0 = disabled)
    pub range_factorization_bins: u32,
    /// expected library format for strand-compatibility filtering; `None` when
    /// auto-detecting (filtering only applies with an explicit `-l`)
    pub expected_format: Option<LibraryFormat>,
    /// linear-space weight multiplier for incompatible mappings
    pub incompat_prior: f64,
    /// drop incompatible mappings entirely (set when `incompat_prior <= 0`)
    pub ignore_incompat: bool,
    /// collect observed sequence-bias contexts (`--seqBias`)
    pub collect_seqbias: bool,
    /// shared observed (fw, rc) sequence-bias models to merge per-thread results into
    pub seqbias_obs: Option<&'a Mutex<(SBModel, SBModel)>>,
    /// collect observed fragment-GC contexts (`--gcBias`)
    pub collect_gcbias: bool,
    /// GC bias model bin counts (`--conditionalGCBins` × `--numGCBins`)
    pub cond_gc_bins: usize,
    pub gc_bins: usize,
    /// per-transcript cumulative-GC source (dense or rank), when GC-correcting
    pub gc_store: Option<salmon_model::GcStore<'a>>,
    /// shared observed fragment-GC model to merge per-thread results into
    pub gcbias_obs: Option<&'a Mutex<GcFragModel>>,
    /// collect observed positional bias (`--posBias`)
    pub collect_posbias: bool,
    /// per-transcript length-class index (for positional bias), when pos-correcting
    pub length_class: Option<&'a [usize]>,
    /// shared observed (5', 3') positional-bias models (per length class) to merge into
    pub posbias_obs: Option<&'a Mutex<(Vec<SimplePosBias>, Vec<SimplePosBias>)>>,
    /// shared per-transcript coverage accumulator (`--degCoverage`). Unlike the
    /// bias models this is not per-thread: a per-worker copy would be one
    /// `bins × transcripts` array *per thread*, so workers update it in place
    /// through relaxed atomics instead (see
    /// [`salmon_degnorm::CoverageAccumulator`]).
    pub coverage: Option<&'a salmon_degnorm::CoverageAccumulator>,
    /// online (dual-phase) inference state; when present, observed bias models
    /// are collected with abundance-aware posteriors instead of score-only weights
    pub online: Option<&'a OnlineInference>,
    /// whether the library is paired (for the online orphan penalty)
    pub paired_lib: bool,
    /// model the fragment-length probability of orphan / single-end mappings via
    /// the bounded-CMF "ambiguous" weight (salmon default). When `false`
    /// (`--noSingleFragProb`) orphans fall back to a flat penalty / weight-1.
    pub model_single_frag_prob: bool,
    /// disable the fragment-length distribution in the per-fragment assignment
    /// probability entirely (`--noFragLengthDist`).
    pub no_frag_length_dist: bool,
    /// assemble equivalence classes during the mapping pass (consumed by the EM
    /// or an eq-class dump): `!skip_quant || dump_eq || dump_eq_weights`. When
    /// `false` (a map-only `--writeRad` run) [`record`] still counts fragments and
    /// trains the FLD but builds no classes. Orthogonal to range-factorization
    /// (driven by `range_factorization_bins`) and to [`Self::use_online`].
    pub build_eq: bool,
    /// compute the abundance-aware online posterior and train the FLD by online
    /// acceptance sampling. When `false` the FLD is trained deterministically from
    /// uniquely-mapped fragments (order-independent; matches the RAD reader's
    /// `derive_fld`) and bias weights fall back to score-only. False whenever the
    /// run does not quantify (`skip_quant`), where the online estimate has no
    /// consumer.
    pub use_online: bool,
    pub num_processed: &'a AtomicU64,
    pub num_mapped: &'a AtomicU64,
    /// mapped fragments whose representative mapping is an orphan (only one mate
    /// of a paired-end fragment mapped); feeds `meta_info.json`'s `num_orphan`.
    /// (`lib_format_counts.json`'s
    /// `num_frags_with_inconsistent_or_orphan_mappings` is NOT built from this —
    /// it is derived purely from the observed-format histogram below.)
    pub num_orphan: &'a AtomicU64,
    /// fragments with at least one mapping that is strand-compatible with the
    /// expected library format, and fragments with mappings of which *none* is.
    /// The second is the direct, free measure of how much of the library sits on
    /// the wrong strand. On a stranded protocol the leading cause is
    /// double-stranded (genomic DNA) input, half of which lands on the expected
    /// strand by chance and cannot be told apart from RNA (#1130). Both feed
    /// `lib_format_counts.json` (`num_compatible_fragments`,
    /// `num_incompatible_fragments`, `compatible_fragment_ratio`).
    ///
    /// Only fragments filtered against a *known* expected format are counted, so
    /// under `-l A` the fragments consumed before the detector locks in are in
    /// neither tally.
    pub num_frags_compat: &'a AtomicU64,
    pub num_frags_incompat: &'a AtomicU64,
    /// per-observed-format fragment counts (`counts[format_id]`), the raw
    /// histogram behind `lib_format_counts.json`'s per-format keys and its
    /// `num_frags_with_concordant_consistent_mappings` /
    /// `strand_mapping_bias` derivations. Tallied for every mapped fragment
    /// regardless of the expected format — the histogram records what was
    /// *observed*, so it is independent of what was declared or detected.
    pub lib_format_counts: &'a [AtomicU64; salmon_core::NUM_LIB_FORMATS],
    /// selective-alignment meta counters (salmon's `aux_info/meta_info.json`):
    /// decoy-dominated fragments, dovetail fragments, fragments that had
    /// candidates but no surviving mapping, and (for mapped fragments) candidate
    /// alignments dropped below the score threshold.
    pub num_decoy: &'a AtomicU64,
    pub num_dovetail: &'a AtomicU64,
    pub num_frags_filtered_vm: &'a AtomicU64,
    pub num_below_threshold_vm: &'a AtomicU64,
    /// when set, the open `unmapped_names.txt` writer (`--writeUnmappedNames`).
    /// Workers stream into it per batch rather than accumulating every name for
    /// the whole run; see [`QuantProcessor::unmapped`].
    pub unmapped_names: Option<&'a Mutex<UnmappedWriter>>,
    /// when set, write per-mapping SAM records (`--writeMappings`)
    pub sam: Option<&'a crate::sam::SamWriter>,
    /// when set, encode per-mapping BAM records into reusable worker chunks
    pub bam: Option<&'a crate::bam::BamOutput>,
    /// when set, write per-fragment mappings to a RAD file (`--writeRad`)
    pub rad: Option<&'a salmon_rad::RadOutputWriter>,
    /// `--deterministic` mode: collect an order-independent fragment-length
    /// distribution (+ library-format tally) from uniquely-mapped proper pairs
    /// during the mapping pass, instead of training the online (log-space) FLD or
    /// using the prefix library-type detector. Built into a baked FLD at end of
    /// pass. When `Some`, the per-fragment work is the cheap [`record_discrete`].
    pub discrete_fld: Option<&'a salmon_model::DiscreteFld>,
    /// `--deterministic` with bias: also build *naive* (uniform-weight, no-FLD)
    /// equivalence classes during the mapping pass, so a rough EM at end of mapping
    /// produces the `initial_abundances` baked into the RAD to seed the requant's
    /// bias model (so the requant is a single fused RAD read). Set only alongside
    /// [`Self::discrete_fld`]. Orientation-tagged so incompatible placements can be
    /// dropped once the library type is resolved at end of mapping.
    pub naive_eq: Option<&'a NaiveEqBuilder>,
}

/// Per-thread mapping processor: the shared state plus this worker's private
/// scratch.
///
/// Every `Option` field is `None` unless the corresponding feature is on, so a
/// plain run allocates none of them. The private buffers exist so the hot path
/// touches no shared memory at all; they are merged or flushed at batch and
/// thread boundaries.
/// Per-worker tallies for the counters that are only ever *written* on the hot
/// path, folded into the shared atomics once in `on_thread_complete`.
///
/// Each of these used to be a `fetch_add` on an `AtomicU64` that every worker
/// shares. Relaxed ordering makes an individual increment cheap, but the cost
/// was never the ordering — it was the cache line: every thread writing the same
/// word serialises behind a coherence transfer, so the tax grows with `-p` and
/// is invisible at low thread counts. Accumulating privately and merging once
/// per thread is the same pattern the bias models already use.
///
/// `num_processed` is deliberately NOT here: it is also *read* on the hot path,
/// to gate `use_aux`, so making it thread-local would change which fragments see
/// the auxiliary model and therefore the output. That gate is its own question
/// (see the `use_aux` issue); this change is output-preserving.
#[derive(Default)]
pub(crate) struct Counters {
    pub orphan: u64,
    /// see [`Shared::num_frags_compat`]
    pub frags_compat: u64,
    pub frags_incompat: u64,
    /// see [`Shared::lib_format_counts`]
    pub lib_fmt: salmon_core::LibFormatCountsArray,
    pub decoy: u64,
    pub dovetail: u64,
    pub frags_filtered_vm: u64,
    pub below_threshold_vm: u64,
}

impl Counters {
    /// Add this worker's tallies to the shared atomics. One `fetch_add` per
    /// counter per thread, instead of one per counter per fragment.
    fn merge_into(&self, sh: &Shared) {
        let add = |a: &AtomicU64, v: u64| {
            if v != 0 {
                a.fetch_add(v, Ordering::Relaxed);
            }
        };
        add(sh.num_orphan, self.orphan);
        add(sh.num_frags_compat, self.frags_compat);
        add(sh.num_frags_incompat, self.frags_incompat);
        for (a, &v) in sh.lib_format_counts.iter().zip(self.lib_fmt.iter()) {
            add(a, v);
        }
        add(sh.num_decoy, self.decoy);
        add(sh.num_dovetail, self.dovetail);
        add(sh.num_frags_filtered_vm, self.frags_filtered_vm);
        add(sh.num_below_threshold_vm, self.below_threshold_vm);
    }
}

/// Index into a fragment's `placements` slice.
///
/// The record path carries two different kinds of position: an index into
/// `placements` (a fragment's candidate mappings) and a position within
/// `compat` (the strand-compatible subset). They coincide only when nothing is
/// filtered out — which is to say on unstranded libraries but not stranded ones
/// — so confusing them mis-weights bias on stranded libraries specifically,
/// with no panic and no test failure on unstranded data. Giving each its own
/// type makes that mistake a compile error.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(transparent)]
pub(crate) struct MapIdx(usize);

impl MapIdx {
    /// The placement this refers to.
    ///
    /// A method rather than `placements[i]`, so a map index cannot silently be
    /// used on a buffer that is aligned to `compat` instead.
    #[inline]
    fn get(self, placements: &[ScoredMapping]) -> &ScoredMapping {
        &placements[self.0]
    }
}

/// Position within `compat`, and so within every buffer aligned to it.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(transparent)]
pub(crate) struct CompatIdx(usize);

/// A per-fragment buffer holding exactly one entry per `compat` element.
///
/// Indexable only by [`CompatIdx`], so indexing it with a [`MapIdx`] does not
/// compile. The alignment was previously enforced by comments alone.
///
/// `#[repr(transparent)]` with `#[inline]` accessors, so this compiles to the
/// same code as indexing the `Vec` directly.
#[derive(Default)]
#[repr(transparent)]
pub(crate) struct CompatAligned<T>(Vec<T>);

impl<T> CompatAligned<T> {
    #[inline]
    fn clear(&mut self) {
        self.0.clear();
    }
    #[inline]
    fn extend<I: IntoIterator<Item = T>>(&mut self, it: I) {
        self.0.extend(it);
    }
    #[inline]
    fn iter(&self) -> std::slice::Iter<'_, T> {
        self.0.iter()
    }
}

impl<T> From<Vec<T>> for CompatAligned<T> {
    #[inline]
    fn from(v: Vec<T>) -> Self {
        Self(v)
    }
}

impl<T> std::ops::Index<CompatIdx> for CompatAligned<T> {
    type Output = T;
    #[inline]
    fn index(&self, i: CompatIdx) -> &T {
        &self.0[i.0]
    }
}

/// Reusable per-fragment working buffers for [`record`].
///
/// `record` runs once per fragment and used to allocate five throwaway `Vec`s
/// every time (the compatible-mapping list, the online-posterior inputs, the
/// bias weights, the per-mapping fragment-length terms, and the tid/weight
/// pairs). Only the two vectors that are *moved* into the equivalence-class
/// builder genuinely have to be fresh; the rest are scratch and live here for
/// the lifetime of the worker.
#[derive(Default)]
pub(crate) struct RecordScratch {
    /// strand-compatible mappings and their (possibly penalised) weights
    compat: Vec<(MapIdx, f64)>,
    /// `(tid, log-weight)` inputs to the online posterior
    online_in: Vec<(u32, f64)>,
    /// abundance-aware (or score-only) per-mapping bias-collection weights
    bias_w: CompatAligned<f64>,
    /// per-mapping fragment-length log-probability
    log_fps: CompatAligned<f64>,
    /// `(tid, conditional probability)` before duplicate-tid merging
    pairs: Vec<(u32, f64)>,
    /// reusable buffers the mapper fills per fragment: raw mappings, the mate
    /// reverse complement used in orphan recovery, and the dedup table for
    /// fragments large enough to take the hashed path
    map: salmon_map::MapScratch,
}

pub(crate) struct QuantProcessor<'a> {
    pub shared: Shared<'a>,
    /// per-thread counter tallies, merged into the shared atomics at thread end
    pub counters: Counters,
    /// per-thread reusable buffers for the per-fragment `record` path
    pub scratch: RecordScratch,
    pub hs: Option<HitSearcher<'a>>,
    /// per-thread reusable sketch caches (avoids per-read MappingCache allocs)
    pub sketch_scratch: Option<salmon_map::SketchScratch>,
    /// per-thread observed (fw, rc) sequence-bias models (when collecting)
    pub seqbias: Option<(SBModel, SBModel)>,
    /// per-thread observed fragment-GC model (when collecting)
    pub gcbias: Option<GcFragModel>,
    /// per-thread observed (5', 3') positional-bias models, per length class
    pub posbias: Option<(Vec<SimplePosBias>, Vec<SimplePosBias>)>,
    /// per-thread unmapped-fragment name lines, flushed to the shared writer per
    /// batch (like `sam_buf`) rather than held for the whole run. This used to be
    /// a `Vec<String>` that grew until the process exited — one heap `String` per
    /// unmapped read, retained to the end. On a heavily contaminated sample that
    /// is tens of GB of small strings for output that could have been written
    /// as it was produced.
    pub unmapped: String,
    /// per-thread SAM record buffer (flushed to the shared writer per batch)
    pub sam_buf: String,
    /// per-thread raw BAM record chunk; allocated only for `--writeBam`. Holds
    /// its own borrow of the writer, so "have a chunk ⇒ have somewhere to send
    /// it" is enforced by the type rather than re-checked at each use.
    pub bam_scratch: Option<crate::bam::BamScratch<'a>>,
    /// per-thread RAD chunk buffer (flushed as one chunk per batch); `None`
    /// unless `--writeRad` is set
    pub rad_buf: Option<salmon_rad::FragmentChunkBuf>,
}

impl<'a> QuantProcessor<'a> {
    pub fn new(shared: Shared<'a>) -> Self {
        let seqbias = shared
            .collect_seqbias
            .then(|| (SBModel::new(), SBModel::new()));
        let gcbias = shared
            .collect_gcbias
            .then(|| GcFragModel::new(shared.cond_gc_bins, shared.gc_bins));
        let posbias = shared.collect_posbias.then(|| {
            (
                (0..NUM_LENGTH_CLASSES)
                    .map(|_| SimplePosBias::default())
                    .collect(),
                (0..NUM_LENGTH_CLASSES)
                    .map(|_| SimplePosBias::default())
                    .collect(),
            )
        });
        Self {
            shared,
            counters: Counters::default(),
            scratch: RecordScratch::default(),
            hs: None,
            sketch_scratch: None,
            seqbias,
            gcbias,
            posbias,
            unmapped: String::new(),
            sam_buf: String::new(),
            bam_scratch: shared.bam.map(|output| output.scratch()),
            rad_buf: shared.rad.map(|rad| {
                salmon_rad::FragmentChunkBuf::with_capacity_codec(64 * 1024, rad.codec())
            }),
        }
    }
}

// paraseq clones the processor once per worker thread. Cloning the *scratch*
// would be wrong (two threads would share a buffer) and wasteful, so this
// deliberately ignores it and builds fresh state from the shared references.
impl Clone for QuantProcessor<'_> {
    fn clone(&self) -> Self {
        // Per-thread scratch is rebuilt fresh in each worker.
        QuantProcessor::new(self.shared)
    }
}

/// salmon's `LOG_EPSILON = log(0.375e-10)` — the small log-probability assigned
/// to implausible placements (out-of-range fragment lengths, unexpected orphans).
///
/// Small but finite: a placement that looks impossible is heavily penalized
/// rather than discarded outright, so a fragment whose every placement is
/// implausible still contributes somewhere instead of vanishing.
const LOG_EPSILON: f64 = -23.998_158_637_57; // (0.375e-10f64).ln()
/// salmon's `numPreBurninFrags`: fold the fragment-length term into the online
/// posterior only after this many fragments have been assigned.
pub(crate) const NUM_PRE_BURNIN: u64 = 5000;

/// Per-mapping fragment-length log-probability (`logFragProb`) — the single term
/// shared by the abundance-aware online posterior and the persistent eq-class
/// weight. C++ salmon computes this once as part of `aln.logProb` and reuses it
/// for both; the Rust port previously duplicated it inconsistently (online used a
/// flat `LOG_EPSILON` for orphans, the eq-class used `LOG_1`), so this unifies
/// them.
///
/// * proper pair → length-conditioned `pmf(flen) − cmf(txpLen)` (salmon's
///   `lenProb − refLengthCM`); `cmf(txpLen)` ≈ 0 for transcripts longer than the
///   FLD support, so this only differs from the bare PMF for short transcripts.
/// * orphan / single-end → the bounded-CMF "ambiguous" weight
///   ([`ambig_frag_log_prob`](salmon_model::ambig_frag_log_prob)), which is
///   already conditioned on the transcript length.
///
/// `pmf`/`cmf` are the per-fragment FLD snapshots (empty before the auxiliary
/// model is trained). With `model_single_frag_prob` off (`--noSingleFragProb`)
/// orphans fall back to a flat `LOG_EPSILON` (paired library) / `LOG_1`
/// (single-end), matching salmon. Returns `LOG_1` (= 0, unmodelled) before the
/// model is trained or when no position is available (sketch mode).
#[allow(clippy::too_many_arguments)]
fn frag_log_prob(
    m: &ScoredMapping,
    txp_len: i32,
    use_aux: bool,
    model_single_frag_prob: bool,
    paired_lib: bool,
    no_frag_length_dist: bool,
    pmf: &[f64],
    cmf: &[f64],
) -> f64 {
    // --noFragLengthDist: do not consider the fragment-length distribution in the
    // per-fragment assignment probability (salmon's flag). Applies to both the
    // proper-pair PMF term and the orphan/single-end ambiguous term.
    if no_frag_length_dist {
        return LOG_1;
    }
    if m.status == MateStatus::PairedEndPaired && m.fragment_len > 0 {
        if use_aux && !pmf.is_empty() {
            let flen = (m.fragment_len as usize).min(pmf.len() - 1);
            // length-conditioning normalizer log P(L ≤ txpLen); ~0 for long txps.
            let cond = cmf[(txp_len.max(0) as usize).min(cmf.len().saturating_sub(1))];
            pmf[flen] - cond
        } else {
            LOG_1
        }
    } else {
        // orphan (PairedEndLeft / PairedEndRight) or single-end read
        let pos = if m.is_fw { m.fw_pos } else { m.rc_pos };
        if model_single_frag_prob && use_aux && pos >= 0 {
            salmon_model::ambig_frag_log_prob(cmf, m.is_fw, pos, m.read_len, txp_len)
        } else if paired_lib
            && matches!(
                m.status,
                MateStatus::PairedEndLeft | MateStatus::PairedEndRight
            )
        {
            LOG_EPSILON // unexpected orphan, fragment length unmodelled
        } else {
            LOG_1
        }
    }
}

thread_local! {
    /// Per-thread PRNG state for stochastic FLD-sample acceptance (mirrors
    /// salmon's per-thread RNG used when training the fragment-length model).
    static FLD_RNG: std::cell::Cell<u64> = const { std::cell::Cell::new(0x2545_F491_4F6C_DD1D) };
}

/// Draw a pseudo-random value in `[0, 1)` from the per-thread xorshift state.
///
/// Used only for acceptance sampling when training the fragment-length
/// distribution online: rather than folding every fragment into the model (which
/// would over-weight the early, poorly-mapped ones), a fragment is accepted with
/// some probability. Per-thread state means no contention; the resulting
/// thread-order dependence is exactly why `--deterministic` uses the integer
/// histogram instead.
fn fld_rng_u01() -> f64 {
    FLD_RNG.with(|s| {
        let mut x = s.get();
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        s.set(x);
        ((x >> 11) as f64) * (1.0 / ((1u64 << 53) as f64))
    })
}

/// Collect the orientation-aware 5'/3' sequence-bias contexts for one mapping,
/// weighted by `weight` (the fragment-transcript posterior). A forward read's 5'
/// context feeds the forward model; a reverse read's 5' context
/// (reverse-complemented) feeds the RC model. For a paired fragment the opposite
/// end feeds the other model.
fn collect_context(
    salmon: &SalmonIndex,
    m: &ScoredMapping,
    weight: f64,
    obs: &mut (SBModel, SBModel),
) {
    let seq = salmon.ref_seq(m.tid);
    let n = seq.len() as i32;
    let (cl, cr, k) = (
        CONTEXT_LEFT as i32,
        CONTEXT_RIGHT as i32,
        CONTEXT_LENGTH as i32,
    );

    let add_fwd = |obs: &mut SBModel, five_prime: i32| {
        let s = five_prime - cl;
        if s >= 0 && s + k <= n {
            obs.add_context(&seq[s as usize..(s + k) as usize], false, weight);
        }
    };
    let add_rev = |obs: &mut SBModel, five_prime: i32| {
        let s = five_prime - cr;
        if s >= 0 && s + k <= n {
            obs.add_context(&seq[s as usize..(s + k) as usize], true, weight);
        }
    };

    if m.is_fw {
        add_fwd(&mut obs.0, m.ref_pos); // 5' -> forward model
        if m.fragment_len > 0 {
            add_rev(&mut obs.1, m.ref_pos + m.fragment_len - 1); // 3' -> RC model
        }
    } else {
        add_rev(&mut obs.1, m.ref_pos); // reverse read's 5' -> RC model
        if m.fragment_len > 0 {
            add_fwd(&mut obs.0, m.ref_pos - m.fragment_len + 1); // 3' -> forward model
        }
    }
}

/// Collect observed fragment-GC contexts for one fragment's compatible paired
/// mappings, each weighted by `bias_w[k]` (the abundance-aware posterior under
/// online inference, else the normalized aux weight). Mirrors salmon's
/// per-alignment `observedGCMass.inc(gcDesc(start, stop), logProb)`.
fn collect_gc(
    sh: &Shared,
    placements: &[ScoredMapping],
    compat: &[(MapIdx, f64)],
    bias_w: &CompatAligned<f64>,
    gc: &mut GcFragModel,
) {
    let Some(store) = sh.gc_store else { return };
    for (k, &(i, _)) in compat.iter().enumerate().map(|(k, e)| (CompatIdx(k), e)) {
        let m = i.get(placements);
        if m.status != MateStatus::PairedEndPaired || m.fragment_len <= 0 {
            continue;
        }
        let start = m.ref_pos;
        let stop = m.ref_pos + m.fragment_len - 1;
        let view = store.view(m.tid as usize);
        let ref_len = view.ref_len() as i32;
        if start < 0 || stop >= ref_len {
            continue;
        }
        if let Some((ff, cf)) = gc_desc(&view, start, stop) {
            // `bias_w` is aligned to `compat`, so it is indexed by the position
            // in `compat` (`k`), not by the mapping index (`i`).
            gc.inc(ff, cf, bias_w[k]);
        }
    }
}

/// Collect observed positional-bias mass for one fragment's compatible mappings,
/// each weighted by `bias_w[k]` (log space). Mirrors salmon's
/// `observedPosBias{Fwd,RC}[lengthClass].addMass(pos, refLen, logProb)`.
fn collect_pos(
    sh: &Shared,
    placements: &[ScoredMapping],
    compat: &[(MapIdx, f64)],
    bias_w: &CompatAligned<f64>,
    pos: &mut (Vec<SimplePosBias>, Vec<SimplePosBias>),
) {
    let Some(length_class) = sh.length_class else {
        return;
    };
    for (k, &(i, _)) in compat.iter().enumerate().map(|(k, e)| (CompatIdx(k), e)) {
        let m = i.get(placements);
        // `bias_w` is aligned to `compat`, indexed by `k`, not by `i`.
        if bias_w[k] <= 0.0 {
            continue;
        }
        let lc = length_class[m.tid as usize];
        let ref_len = sh.salmon.ref_len(m.tid as usize) as i32;
        // add_mass now takes LINEAR mass (accumulated in fixed-point); pass the
        // posterior weight directly rather than its log.
        let mass = bias_w[k];
        if m.fw_pos >= 0 {
            pos.0[lc].add_mass(m.fw_pos, ref_len, mass);
        }
        if m.rc_pos >= 0 {
            pos.1[lc].add_mass(m.rc_pos, ref_len, mass);
        }
    }
}

/// Tally a mapped fragment into the per-worker observed-format histogram: one
/// count per distinct observed format among its placements (ports salmon's
/// per-fragment `libTypeCountsPerFrag[i] > 0` accumulation). Pairs use their
/// recorded format (falling back to the mate-strand derivation, which is how
/// the RAD writer and `DiscreteFld` classify the same placement); orphans and
/// single-end mates are a single-end observation on their strand.
fn tally_observed_formats(
    placements: &[ScoredMapping],
    lib_fmt: &mut salmon_core::LibFormatCountsArray,
) {
    let mut mask: u16 = 0;
    for m in placements {
        let f = match (m.format, m.status) {
            (Some(f), _) => f,
            (None, MateStatus::PairedEndPaired) => observed_paired_format(m.is_fw, m.r2_fw),
            (None, _) => salmon_core::observed_single_format(m.is_fw),
        };
        mask |= 1 << f.format_id();
    }
    let mut bits = mask;
    while bits != 0 {
        let i = bits.trailing_zeros() as usize;
        lib_fmt[i] += 1;
        bits &= bits - 1;
    }
}

/// Lay one fragment's placements onto the per-transcript coverage curves
/// (`--degCoverage`), each weighted by `bias_w[k]`.
///
/// The span is the fragment's footprint on the transcript: the inferred insert
/// for a proper pair, the read itself for an orphan or a single-end library.
/// That is the footprint degradation removes, so it is what the downstream
/// model has to see — the *start* alone (what positional bias records) would
/// leave the interior of every fragment invisible.
///
/// Unlike the bias models this is not stopped at the end of the training
/// window: a coverage curve is a summary of the whole library, and truncating
/// it after the first few million fragments would describe the head of the
/// FASTQ rather than the sample. The weights themselves do change at that
/// point — `bias_w` is the abundance-aware online posterior inside the window
/// and the score-only normalized weight after it — which matters only for
/// multi-mapping fragments, and only in how they are split between the
/// transcripts they are compatible with.
fn collect_coverage(
    placements: &[ScoredMapping],
    compat: &[(MapIdx, f64)],
    bias_w: &CompatAligned<f64>,
    cov: &salmon_degnorm::CoverageAccumulator,
) {
    for (k, &(i, _)) in compat.iter().enumerate().map(|(k, e)| (CompatIdx(k), e)) {
        let m = i.get(placements);
        // `bias_w` is aligned to `compat`, indexed by `k`, not by `i`.
        let mass = bias_w[k];
        if mass <= 0.0 {
            continue;
        }
        let span = if m.status == MateStatus::PairedEndPaired && m.fragment_len > 0 {
            m.fragment_len
        } else {
            m.read_len
        };
        cov.add_fragment(m.tid as usize, m.ref_pos, span, mass);
    }
}

/// Cheap per-fragment work for `--deterministic`'s mapping pass: count the
/// fragment and, for a uniquely-mapped proper pair, record its length + observed
/// format into the order-independent [`salmon_model::DiscreteFld`]. No online
/// inference, eq-class assembly, bias collection, or prefix library-type
/// detection — the FLD and library format are derived deterministically from the
/// accumulator at end of pass and baked into the RAD.
fn record_discrete(
    sh: &Shared,
    placements: &[ScoredMapping],
    acc: &salmon_model::DiscreteFld,
    counters: &mut Counters,
) {
    sh.num_processed.fetch_add(1, Ordering::Relaxed);
    if placements.is_empty() {
        return;
    }
    // Stays shared, unlike the other tallies: `record` reads this running total
    // to gate the auxiliary models (see `use_aux`), so a per-worker count merged
    // at thread end would leave it at zero for the whole run.
    sh.num_mapped.fetch_add(1, Ordering::Relaxed);
    // Same orphan definition `record` uses: a fragment's surviving mappings share
    // one status, so the first is representative. Counting it here is what lets
    // `--deterministic` report the concordant/orphan split in
    // `lib_format_counts.json` instead of calling every mapped fragment
    // concordant.
    if matches!(
        placements[0].status,
        MateStatus::PairedEndLeft | MateStatus::PairedEndRight
    ) {
        counters.orphan += 1;
    }
    tally_observed_formats(placements, &mut counters.lib_fmt);
    // Library-format accounting (see `record`). The deterministic pass does no
    // strand filtering (incompatible placements are dropped later, when the RAD
    // is quantified), but with an explicit `-l` the expected format is already
    // known here, so the same wrong-strand tally is free. Under `-l A` the format
    // is only resolved at end of pass, so nothing is counted.
    if let Some(exp) = sh.expected_format {
        if placements
            .iter()
            .any(|m| is_compatible(exp, m.format, m.is_fw, m.status))
        {
            counters.frags_compat += 1;
        } else {
            counters.frags_incompat += 1;
        }
    }
    // A uniquely-mapped proper pair unambiguously implies its fragment length and
    // orientation; its observed format also feeds order-independent `-l A`
    // detection. Uses the raw mapping (like the RAD records written this pass, and
    // like the reader's `derive_fld`), so a baked FLD matches a re-derived one.
    if placements.len() == 1 {
        let m = &placements[0];
        if m.status == MateStatus::PairedEndPaired && m.fragment_len > 0 {
            // Orientation from the mate strands rather than `m.format`, exactly
            // as `build_rad_record` derives the RAD hit — the two must classify
            // a placement identically for a baked FLD to match a re-derived one.
            acc.add(m.fragment_len as usize, m.is_fw, m.r2_fw);
        }
    }
    // Naive (kallisto-style) equivalence class for the rough seed EM: the
    // compatible transcripts, orientation-tagged so library-incompatible
    // placements can be dropped once the library type is resolved at end of
    // mapping. No FLD/score weighting (the FLD isn't known yet). The observed
    // paired format is derived from the mate strands (as `build_rad_record` /
    // `DiscreteFld` do), keeping every deterministic classifier identical.
    if let Some(nb) = sh.naive_eq {
        let sig: Vec<NaivePlacement> = placements
            .iter()
            .map(|m| {
                let fmt_id = if m.status == MateStatus::PairedEndPaired {
                    observed_paired_format(m.is_fw, m.r2_fw).format_id()
                } else {
                    NAIVE_NO_FMT
                };
                NaivePlacement {
                    tid: m.tid,
                    fmt_id,
                    is_fw: m.is_fw,
                    status: m.status,
                }
            })
            .collect();
        nb.add(sig);
    }
}

/// Record one fragment's weighted mappings into the shared accumulators.
///
/// The per-fragment quant work decomposes into independent axes, each gated by a
/// flag on [`Shared`] so a run does exactly what it needs and no more:
///   * [`Shared::use_online`] — compute the abundance-aware online posterior
///     (advancing the online masses) and train the fragment-length distribution
///     by online acceptance sampling. When `false` the FLD is instead trained
///     *deterministically* from uniquely-mapped concordant pairs (the
///     order-independent, commutative-histogram estimate the RAD reader's
///     `derive_fld` also uses) and bias weights fall back to score-only.
///   * [`Shared::build_eq`] — assemble the equivalence class and add it to the
///     shared builder (consumed by the EM or an eq-class dump). Whether those
///     classes are *range-factorized* is a further, orthogonal choice driven by
///     `range_factorization_bins` (0 ⇒ basic classes, e.g. for a rough warm-up
///     EM). When `false` (a map-only `--writeRad` run) no eq-class is built.
///
/// Counting, library detection, and FLD training always run; `log_fm` is the
/// current minibatch's forgetting mass (used only when `use_online`).
///
/// `fld_pmf`/`fld_cmf` are the batch's fragment-length snapshot pair, taken once
/// by the caller (see the batch loops). They are passed in rather than read here
/// because the snapshot only changes at `refresh_online`, which runs at the batch
/// boundary — so re-reading it per fragment would re-acquire two `RwLock`s and
/// bump two `Arc` refcounts, on cache lines every worker shares, to obtain a value
/// that provably cannot have changed.
#[allow(clippy::too_many_arguments)]
fn record(
    sh: &Shared,
    placements: &[ScoredMapping],
    log_fm: f64,
    fld_pmf: &[f64],
    fld_cmf: &[f64],
    counters: &mut Counters,
    scratch: &mut RecordScratch,
    seqbias: Option<&mut (SBModel, SBModel)>,
    gcbias: Option<&mut GcFragModel>,
    posbias: Option<&mut (Vec<SimplePosBias>, Vec<SimplePosBias>)>,
) {
    // Still a shared atomic: unlike the other counters this one is *read* below,
    // to gate `use_aux`. See `Counters`.
    sh.num_processed.fetch_add(1, Ordering::Relaxed);
    if placements.is_empty() {
        return;
    }
    // Independent `&mut` borrows of the reusable buffers.
    let RecordScratch {
        compat,
        online_in,
        bias_w,
        log_fps,
        pairs,
        // Belongs to the mapping step, not to `record`.
        map: _,
    } = scratch;

    // Sample the observed library format for auto-detection (auto mode only).
    // Sampled from the raw mappings, before the strand-compatibility filter
    // below, so detection sees the true observed orientation/strand.
    if let Some(det) = sh.detector {
        if det.is_active() {
            if let Some(m) = placements
                .iter()
                .filter(|m| m.format.is_some())
                .max_by(|a, b| a.weight.total_cmp(&b.weight))
            {
                det.add_sample(m.format.unwrap());
            }
        }
    }

    // The observed-format histogram, from the raw mappings for the same reason
    // as the detector sample above: a fragment whose every placement is
    // wrong-strand is dropped by the filter below, and its observation is
    // precisely what the histogram exists to report.
    tally_observed_formats(placements, &mut counters.lib_fmt);

    // Strand-compatibility filtering against the expected format: the explicit
    // `-l` type, or — in auto (`-l A`) mode — the format the detector locks in
    // once it has seen enough of the read prefix (`resolved_format`). Until the
    // prefix is consumed this is `None` and nothing is filtered; afterwards the
    // detected type is applied for the rest of the sample (salmon's behavior).
    let expected = sh
        .expected_format
        .or_else(|| sh.detector.and_then(|d| d.resolved_format()));
    // Indices into `placements` rather than references, so the buffer can outlive the
    // call and be reused by the next fragment.
    compat.clear();
    // Tallied here rather than derived from `compat.len()` afterwards, because
    // with `--incompatPrior > 0` an incompatible mapping is kept (down-weighted)
    // and is indistinguishable from a compatible one once it is in the buffer.
    let mut num_compat_placements = 0usize;
    compat.extend(
        placements
            .iter()
            .enumerate()
            .filter_map(|(i, m)| match expected {
                Some(exp) => {
                    if is_compatible(exp, m.format, m.is_fw, m.status) {
                        num_compat_placements += 1;
                        Some((MapIdx(i), m.weight))
                    } else if sh.ignore_incompat {
                        None
                    } else {
                        Some((MapIdx(i), m.weight * sh.incompat_prior))
                    }
                }
                None => Some((MapIdx(i), m.weight)),
            }),
    );
    // Library-format accounting, independent of what is quantified below: a
    // fragment that mapped somewhere but nowhere on the expected strand is the
    // wrong-strand evidence `lib_format_counts.json` reports. It is counted under
    // every `--incompatPrior`, including the default 0 where the fragment is then
    // dropped, and only when there is an expected format to compare against.
    if expected.is_some() {
        if num_compat_placements > 0 {
            counters.frags_compat += 1;
        } else {
            counters.frags_incompat += 1;
        }
    }
    if compat.is_empty() {
        return; // no compatible mapping -> fragment is unassigned (not mapped)
    }

    // The fragment has at least one strand-compatible mapping, so it contributes
    // mass and counts as mapped. Count it *after* the compatibility filter above:
    // counting on any mapping (before the filter) over-reported `num_mapped` /
    // `percent_mapped` / `num_compatible_fragments` on stranded libraries, because
    // a fragment whose every mapping was strand-incompatible is dropped here yet
    // quantifies as nothing — which also broke the invariant `Σ NumReads ==
    // num_mapped`. C++ salmon likewise counts a fragment as mapped only once it
    // has a compatible mapping.
    // `fetch_add` returns the count *before* this fragment, which is what the
    // burn-in gate below needs: the online gate reads `num_assigned()` before its
    // own increment, and C++ evaluates `useAuxParams` before counting the current
    // read. Reading the post-increment value would gate this fragment on a total
    // that includes itself, differing from both by one.
    //
    // This one counter therefore stays shared while the rest became per-worker:
    // it is *read* during the run, not merely accumulated, so a thread-local
    // tally merged at thread end would gate every fragment on zero.
    let assigned_before_this = sh.num_mapped.fetch_add(1, Ordering::Relaxed);
    // Classify the fragment as orphan when only one mate of a pair was placed
    // (PairedEndLeft / PairedEndRight). Single-end libraries never count as
    // orphan here.
    //
    // Reading `compat[0]` is sound because every surviving mapping of a fragment
    // agrees on *orphan-ness*, which is the only thing this test distinguishes:
    //   * selective alignment: `map_read_pair_into` drops all orphans as soon as
    //     one non-decoy concordant pair survives, and `finalize_mappings_counted_into`
    //     never emits decoy mappings — so the set is either all-paired or
    //     all-orphan;
    //   * sketch: `map_read_pair_sketch_into` derives one `status` from the
    //     fragment's `map_type` and stamps it on every accepted hit.
    //
    // "All-orphan" does NOT mean one shared status: when neither mate pairs
    // concordantly to a transcript, mate 1 can orphan onto transcript A while
    // mate 2 orphans onto transcript B, leaving `PairedEndLeft` and
    // `PairedEndRight` side by side. Measured on a GRCh38 decoy index, 525 of
    // 300,000 fragments (0.19%) do exactly that, and no fragment ever mixes
    // `PairedEndPaired` or `SingleEnd` with anything else. Both orphan statuses
    // answer this test the same way, so which one sorts first cannot change the
    // count — but asserting a single shared status was wrong, and aborted debug
    // builds on any decoy-bearing index (#1115).
    //
    // Note `compat` is ordered by transcript id (finalize sorts by `tid`), not by
    // score, so this must not be read as "the best mapping" — it is "the status
    // this fragment mapped with", which is well defined.
    let is_orphan = |i: MapIdx| {
        matches!(
            i.get(placements).status,
            MateStatus::PairedEndLeft | MateStatus::PairedEndRight
        )
    };
    debug_assert!(
        compat
            .iter()
            .all(|&(i, _)| is_orphan(i) == is_orphan(compat[0].0)),
        "a fragment's surviving mappings must agree on orphan-ness"
    );
    if is_orphan(compat[0].0) {
        counters.orphan += 1;
    }

    // `fld_pmf`/`fld_cmf` are ONE immutable FLD snapshot pair (PMF + CMF), taken
    // once per batch by the caller and shared here between the online posterior
    // and the persistent eq-class weight. Indexing the same snapshot for every
    // mapping guarantees that transcripts sharing a fragment length (in
    // particular exact-duplicate transcripts) receive an identical
    // fragment-length probability even if another thread refreshes the shared
    // snapshot meanwhile. Reading the live FLD per mapping (racing atomic loads
    // on a concurrently-updated distribution) returned slightly different values
    // for the same length and let the VBEM α<1 prior break duplicate symmetry.
    // Mirrors C++ salmon's per-worker `LogCMFCache`. Empty until the first online
    // refresh (the pre-burn-in window, where the FLD term is not folded in
    // anyway).

    // Per-fragment bias-collection weights. With online (dual-phase) inference
    // these are abundance-aware posteriors `softmax(mass_t + log w_t)`, which
    // also advances the online masses by `logForgettingMass + log(posterior)`;
    // without it they fall back to the normalized aux (score) weights.
    let collecting =
        seqbias.is_some() || gcbias.is_some() || posbias.is_some() || sh.coverage.is_some();
    // Abundance-aware online posterior, computed once per fragment within the
    // model-training window (`o.collecting()`), advancing the online masses.
    // Mirrors salmon's `aln.logProb = transcriptLogCount + auxProb + startPosProb`
    // normalized over the fragment's mappings, where
    //   logFragCov   = ln(score weight)
    //   startPosProb = proper pair -> -ln(refLen - flen + 1) (flen<=refLen, else
    //                  LOG_EPSILON); otherwise -ln(refLen)
    //   logFragProb  = the shared `frag_log_prob` term (proper-pair conditioned
    //                  PMF, or the orphan / single-end ambiguous-length weight).
    // Used for abundance-aware bias collection and abundance-aware FLD training.
    let online_post: Option<CompatAligned<f64>> = sh
        .online
        .filter(|o| sh.use_online && o.collecting())
        .map(|online| {
            let use_aux = online.num_assigned() >= sh.pre_burnin;
            online_in.clear();
            online_in.extend(compat.iter().map(|&(i, w)| {
                let m = i.get(placements);
                {
                    let ref_len = sh.salmon.ref_len(m.tid as usize);
                    let rl = ref_len.max(1) as f64;
                    let proper = m.status == MateStatus::PairedEndPaired && m.fragment_len > 0;
                    let flen = m.fragment_len as f64;
                    let start_pos_prob = if proper {
                        if flen <= rl {
                            -((rl - flen + 1.0).ln())
                        } else {
                            LOG_EPSILON
                        }
                    } else {
                        -(rl.ln())
                    };
                    let log_frag_prob = frag_log_prob(
                        m,
                        ref_len as i32,
                        use_aux,
                        sh.model_single_frag_prob,
                        sh.paired_lib,
                        sh.no_frag_length_dist,
                        fld_pmf,
                        fld_cmf,
                    );
                    let log_cov = if w > 0.0 { w.ln() } else { f64::NEG_INFINITY };
                    (m.tid, log_cov + start_pos_prob + log_frag_prob)
                }
            }));
            CompatAligned::from(online.assign_fragment(online_in, log_fm))
        });
    bias_w.clear();
    if collecting {
        if let Some(post) = &online_post {
            bias_w.extend(post.iter().copied());
        } else {
            let wsum: f64 = compat.iter().map(|&(_, w)| w).sum();
            bias_w.extend(
                compat
                    .iter()
                    .map(|&(_, w)| if wsum > 0.0 { w / wsum } else { 0.0 }),
            );
        }
    }
    // After burn-in salmon freezes model collection (still advances masses).
    let collect_now = collecting && sh.online.is_none_or(|o| o.collecting());

    if collect_now {
        if let Some(obs) = seqbias {
            for (k, &(i, _)) in compat.iter().enumerate().map(|(k, e)| (CompatIdx(k), e)) {
                collect_context(sh.salmon, i.get(placements), bias_w[k], obs);
            }
        }
        if let Some(gc) = gcbias {
            collect_gc(sh, placements, compat, bias_w, gc);
        }
        if let Some(pos) = posbias {
            collect_pos(sh, placements, compat, bias_w, pos);
        }
    }
    // Deliberately outside `collect_now`: see `collect_coverage`.
    if let Some(cov) = sh.coverage {
        collect_coverage(placements, compat, bias_w, cov);
    }

    // Train the fragment-length distribution from this fragment. The FLD is
    // always needed — to quantify now, to set effective lengths, or to bake into
    // a `--writeRad` header — so this runs regardless of `build_eq`. Two methods,
    // selected by `use_online`:
    if let Some(post) = &online_post {
        // Online (abundance-aware) acceptance: accept each concordant compatible
        // pair's fragment length with probability = its online posterior
        // (salmon's `if (r < exp(aln.logProb)) fragLengthDist.addVal(...)`, where
        // aln.logProb includes transcriptLogCount). For reads shared between
        // near-duplicates this preferentially samples the dominant transcript's
        // implied length, concentrating the FLD as salmon does (vs adding every
        // best pair at full weight, which overdisperses it). Frozen after the
        // training window (`online_post` is `None`).
        for (k, &(i, _)) in compat.iter().enumerate().map(|(k, e)| (CompatIdx(k), e)) {
            let m = i.get(placements);
            let conc = m.status == MateStatus::PairedEndPaired
                && m.fragment_len > 0
                && expected.is_none_or(|exp| is_compatible(exp, m.format, m.is_fw, m.status));
            if conc && fld_rng_u01() < post[k] {
                sh.fld.add_val(m.fragment_len as usize, 0.0);
            }
        }
    } else if !sh.use_online {
        // Deterministic: a *uniquely*-mapped concordant proper pair implies its
        // fragment length exactly. Adding only unique fragments at full weight is
        // order-independent (a commutative histogram) and matches the RAD reader's
        // `derive_fld`, so a baked FLD reproduces on requant. (When `use_online`
        // is true but `online_post` is `None`, the training window has closed and
        // the FLD is intentionally frozen — so this branch is skipped.)
        if compat.len() == 1 {
            let m = &compat[0].0.get(placements);
            if m.status == MateStatus::PairedEndPaired
                && m.fragment_len > 0
                && expected.is_none_or(|exp| is_compatible(exp, m.format, m.is_fw, m.status))
            {
                sh.fld.add_val(m.fragment_len as usize, 0.0);
            }
        }
    }

    // Equivalence-class assembly is pure quantification machinery (consumed by
    // the EM or an eq-class dump). A map-only `--writeRad` run produces the RAD
    // above and stops here. Range-factorization of the weights below is a further,
    // orthogonal choice driven by `range_factorization_bins` (0 ⇒ basic classes).
    if !sh.build_eq {
        return;
    }

    // Fold the fragment-length probability into the equivalence-class weight.
    // salmon's eq-class auxProb is `logFragProb + logFragCov + logAlignCompatProb`
    // (SalmonQuantify.cpp) -- i.e. it includes the FLD term but NOT the
    // start-position term (that enters only the abundance/`aln.logProb`; the
    // effective-length factor is applied separately in `update_eff_lengths`).
    // The port previously built the eq-class weight from the bare coverage
    // weight, flattening it for paralogs/isoforms whose implied insert size
    // differs and coarsening the range-factorized classes. The FLD term is only
    // applied once the auxiliary model is trained (after `pre_burnin` fragments),
    // matching salmon's `numPreAuxModelSamples` gating.
    // Gate on *assigned* fragments, matching the online path above
    // (`online.num_assigned()`) and C++'s `useAuxParams`, which tests
    // `localNumAssignedFragments + numAssignedFragments` against
    // `numPreBurninFrags`. This previously read `num_processed`, which counts
    // every fragment reaching `record` — including ones with no mappings at all,
    // and ones whose mappings were all library-incompatible — so the gate opened
    // earlier than salmon intends, by a margin that grows as the mapping rate
    // falls. See #1089.
    let use_aux = assigned_before_this >= sh.pre_burnin;
    // Per-mapping FLD log-probability via the shared `frag_log_prob` term — the
    // same one the online posterior uses, reading the same per-fragment snapshot
    // pair captured above (so the two paths share one consistent FLD view, as
    // C++ does via a single `aln.logProb`). Proper pairs get the length-
    // conditioned PMF; orphans / single-end reads get the bounded-CMF ambiguous
    // weight. `LOG_1` (= 0) before the auxiliary model is trained. The snapshot
    // is primed from the prior before the run starts and refreshed at every
    // mini-batch boundary, so it is never empty here.
    log_fps.clear();
    log_fps.extend(compat.iter().map(|&(i, _)| {
        let m = i.get(placements);
        let ref_len = sh.salmon.ref_len(m.tid as usize) as i32;
        frag_log_prob(
            m,
            ref_len,
            use_aux,
            sh.model_single_frag_prob,
            sh.paired_lib,
            sh.no_frag_length_dist,
            fld_pmf,
            fld_cmf,
        )
    }));
    // Per-fragment conditional probabilities, normalized in *log* space — the
    // log weight of each mapping is `ln(score weight) + logFragProb`, and we
    // subtract the log-sum-exp over the fragment's mappings (C++'s
    // `exp(auxProb - auxDenom)`). This is the same normalization salmon's
    // C++ does and is mathematically identical to the linear `w/Σw`, but doing it
    // in log space keeps it well-defined when every FLD weight underflows: a
    // fragment whose implied lengths all have ~0 FLD probability (logFragProb at
    // the no-mass sentinel) would, in linear space, give all-zero weights — the
    // `wsum > 0` normalization below then leaves them zero, the eq-class denom is
    // 0, and the VBEM silently drops the fragment's count (lost mapped mass; the
    // EM's degenerate-class branch is a no-op, as in C++). In log space the same
    // fragment yields well-defined relative weights, so no mass is lost.
    let log_denom = compat
        .iter()
        .zip(log_fps.iter())
        .fold(LOG_0, |acc, (&(_, w), &lfp)| log_add(acc, w.ln() + lfp));
    pairs.clear();
    pairs.extend(
        compat
            .iter()
            .zip(log_fps.iter())
            .map(|(&(i, w), &lfp)| (i.get(placements).tid, (w.ln() + lfp - log_denom).exp())),
    );

    // Build the equivalence class: sorted transcript ids + weights, combining
    // duplicate ids by SUMMING their conditional probabilities. A fragment that
    // maps to the same transcript more than once (e.g. both orientations of an
    // unstranded library, or two positions) contributes P = Σ of the placements,
    // not just the first — matching the RAD path's per-tid `logsumexp`.
    pairs.sort_by_key(|p| p.0);
    pairs.dedup_by(|a, b| {
        if a.0 == b.0 {
            b.1 += a.1;
            true
        } else {
            false
        }
    });
    let tids: Vec<u32> = pairs.iter().map(|p| p.0).collect();
    let mut weights: Vec<f64> = pairs.iter().map(|p| p.1).collect();

    // Normalize to per-fragment conditional probabilities (sum to 1).
    let wsum: f64 = weights.iter().sum();
    if wsum > 0.0 {
        for w in &mut weights {
            *w /= wsum;
        }
    }

    let group = if sh.range_factorization_bins > 0 {
        let bins = range_factorize_bins(&weights, sh.range_factorization_bins);
        TranscriptGroup::with_bins(tids, bins)
    } else {
        TranscriptGroup::from_sorted(tids)
    };
    sh.eq.add_group(group, weights, 1);
}

/// Fold the most recent fragment's selective-alignment [`MapStats`] into the
/// shared meta counters. Call once per mapped/attempted fragment on the SA path.
fn accumulate_vm_stats(maps_empty: bool, counters: &mut Counters) {
    let s = salmon_map::take_last_map_stats();
    if maps_empty {
        if s.decoy_dominated {
            counters.decoy += 1;
        } else if s.had_candidates {
            counters.frags_filtered_vm += 1;
        }
    } else if s.alns_below_threshold > 0 {
        counters.below_threshold_vm += s.alns_below_threshold as u64;
    }
    if s.dovetail {
        counters.dovetail += 1;
    }
}

/// Merge this thread's observed bias models (sequence and GC) into the shared
/// accumulators.
/// The open `unmapped_names.txt`, shared by every worker.
///
/// Workers buffer their own lines and hand over a block per batch, so the lock
/// is taken once per batch rather than once per unmapped read, and no name is
/// held in memory beyond the batch that produced it.
pub struct UnmappedWriter {
    w: std::io::BufWriter<std::fs::File>,
}

impl UnmappedWriter {
    /// Create (or truncate) `unmapped_names.txt` at `path`.
    pub fn create(path: &std::path::Path) -> std::io::Result<Self> {
        Ok(Self {
            // 1 MiB: this is a bulk sequential text stream, like the SAM writer.
            w: std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(path)?),
        })
    }

    /// Append one worker's buffered lines.
    fn write_block(&mut self, block: &str) -> std::io::Result<()> {
        use std::io::Write as _;
        self.w.write_all(block.as_bytes())
    }

    /// Flush buffered bytes to the file. Called once, after every worker has
    /// finished, so a truncated file cannot be mistaken for a complete one.
    pub fn finish(mut self) -> std::io::Result<()> {
        use std::io::Write as _;
        self.w.flush()
    }
}

/// Append one `"<name> u"` line to this worker's buffer.
///
/// The name is the id up to the first whitespace (the conventional fragment
/// name), lossily decoded. Written straight into the buffer, so an unmapped read
/// costs no heap allocation of its own — the old path built a `String` per read
/// and pushed it into a `Vec` that lived for the whole run.
fn push_unmapped_name(buf: &mut String, id: &[u8]) {
    let end = id
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(id.len());
    buf.push_str(&String::from_utf8_lossy(&id[..end]));
    buf.push_str(" u\n");
}

/// Write this worker's buffered unmapped names to the shared file and clear the
/// buffer. Called per batch and again when the thread finishes, so the names
/// leave memory as they are produced instead of accumulating for the whole run.
/// Fold this worker's private counter tallies into the shared atomics.
fn merge_counters(proc: &QuantProcessor) {
    proc.counters.merge_into(&proc.shared);
}

fn flush_unmapped(proc: &mut QuantProcessor) -> std::io::Result<()> {
    if let Some(shared) = proc.shared.unmapped_names {
        if !proc.unmapped.is_empty() {
            shared.lock().unwrap().write_block(&proc.unmapped)?;
            proc.unmapped.clear();
        }
    }
    Ok(())
}

/// Merge this worker's private bias models into the shared ones.
///
/// Called once per thread at the end of its work, not per fragment: the mutexes
/// below would be ruinously contended otherwise. The merges are integer sums, so
/// the result does not depend on the order threads finish in.
fn merge_bias(proc: &QuantProcessor) {
    if let (Some(local), Some(shared_mtx)) = (proc.seqbias.as_ref(), proc.shared.seqbias_obs) {
        let mut g = shared_mtx.lock().unwrap();
        g.0.combine_counts(&local.0);
        g.1.combine_counts(&local.1);
    }
    if let (Some(local), Some(shared_mtx)) = (proc.gcbias.as_ref(), proc.shared.gcbias_obs) {
        let mut g = shared_mtx.lock().unwrap();
        g.combine_counts(local);
    }
    if let (Some(local), Some(shared_mtx)) = (proc.posbias.as_ref(), proc.shared.posbias_obs) {
        let mut g = shared_mtx.lock().unwrap();
        for (a, b) in g.0.iter_mut().zip(&local.0) {
            a.combine(b);
        }
        for (a, b) in g.1.iter_mut().zip(&local.1) {
            a.combine(b);
        }
    }
}

// The paired-end entry point. `PairedParallelProcessor` is paraseq's interface
// for "given a batch of mate pairs, do something with them"; paraseq owns the
// reading and the threading, this owns the per-fragment work.
impl<'a, 'r> PairedParallelProcessor<RefRecord<'r>> for QuantProcessor<'a> {
    /// Map and record one worker's batch of mate pairs.
    fn process_record_pair_batch(
        &mut self,
        pairs: impl Iterator<Item = (RefRecord<'r>, RefRecord<'r>)>,
    ) -> paraseq::Result<()> {
        // Destructure `self` into its fields up front: the loop below needs
        // several of them mutably at once, which the borrow checker only permits
        // on distinct fields, not through repeated `self.` accesses.
        let mut busy = self.shared.busy.timer();
        let QuantProcessor {
            shared,
            counters,
            scratch,
            hs,
            sketch_scratch,
            seqbias,
            gcbias,
            posbias,
            unmapped,
            sam_buf,
            bam_scratch,
            rad_buf,
        } = self;
        let sh = *shared;
        let idx = sh.salmon.inner();
        // The seed searcher is created lazily on this worker's first batch, so
        // its (substantial) state is never allocated for threads paraseq does not
        // end up using.
        if hs.is_none() {
            *hs = Some(HitSearcher::new(idx));
        }
        let hs = hs.as_mut().unwrap();
        if sh.sketch && sketch_scratch.is_none() {
            *sketch_scratch = Some(salmon_map::SketchScratch::new(idx.k()));
        }
        // one forgetting-mass timestep per minibatch (online inference)
        let log_fm = if sh.use_online {
            sh.online.map_or(0.0, |o| o.next_log_fm())
        } else {
            0.0
        };
        // One FLD snapshot pair for the whole batch. The shared snapshot is only
        // replaced by `refresh_online`, which runs at the batch boundary, so a
        // per-fragment re-read would take two `RwLock`s and bump two `Arc`
        // refcounts — on cache lines shared by every worker — to fetch a value
        // that cannot have changed. Taking it here also makes the per-fragment
        // view *more* stable than before: previously a concurrent refresh could
        // land between two fragments of the same batch.
        let fld_pmf = sh.fld.online_snapshot();
        let fld_cmf = sh.fld.online_cmf_snapshot();
        // Reused across the batch's reads (the `*_into` calls clear it), so the
        // per-fragment result Vec isn't reallocated per read.
        let mut placements: Vec<ScoredMapping> = Vec::new();
        for (r1, r2) in pairs {
            // Per fragment, not per batch: a batch runs far longer than one
            // broker sampling window, so a per-batch counter would leave whole
            // windows reading zero work and report a thread mapping flat out as
            // fully starved. `BusyMeter` publishes on its own cadence.
            busy.tick();
            let s1 = r1.seq();
            let s2 = r2.seq();
            if sh.sketch {
                map_read_pair_sketch_into(
                    &mut placements,
                    idx,
                    hs,
                    sketch_scratch.as_mut().unwrap(),
                    s1.as_ref(),
                    s2.as_ref(),
                    sh.sketch_strict_orphan,
                    sh.map_cfg.pair.allow_dovetail,
                    sh.skip,
                    sh.map_cfg.collect.max_hit_occ,
                    sh.max_read_occ,
                    sh.salmon,
                    sh.map_cfg.score.allow_decoy_orphans,
                );
            } else {
                map_read_pair_into(
                    &mut placements,
                    idx,
                    hs,
                    sh.salmon,
                    s1.as_ref(),
                    s2.as_ref(),
                    sh.map_cfg,
                    &mut scratch.map,
                );
            }
            // Sketch mappings carry no per-hit decoy flag and bypass the
            // selective-alignment finalize, so decoys would otherwise leak into
            // the eq-classes. Apply the same decoy policy here (drop decoy tids,
            // decoy-domination, --allowDecoyOrphans) and account decoy-dominated
            // fragments. SA mode already handled decoys inside finalize.
            if sh.sketch && sh.salmon.info().num_decoys > 0 {
                let decoy_dominated =
                    salmon_map::filter_sketch_decoys(&mut placements, sh.salmon, &sh.map_cfg.score);
                if decoy_dominated {
                    counters.decoy += 1;
                }
            }
            // A fragment mapping to too many places is discarded (salmon's
            // tooManyHits / maxReadOccs): treat as unmapped everywhere below.
            if placements.len() > sh.max_read_occ {
                placements.clear();
            }
            if placements.is_empty() && sh.unmapped_names.is_some() {
                push_unmapped_name(unmapped, r1.id());
            }
            if !sh.sketch {
                accumulate_vm_stats(placements.is_empty(), counters);
            }
            if sh.sam.is_some() && !placements.is_empty() {
                crate::sam::write_fragment(
                    sam_buf,
                    sh.salmon,
                    r1.id(),
                    s1.as_ref(),
                    Some((r2.id(), s2.as_ref())),
                    &placements,
                );
            }
            if let Some(scratch) = bam_scratch.as_mut() {
                if !placements.is_empty() {
                    scratch.write_fragment(
                        sh.salmon,
                        r1.id(),
                        s1.as_ref(),
                        Some((r2.id(), s2.as_ref())),
                        &placements,
                    )?;
                }
            }
            if let (Some(rad), Some(buf)) = (sh.rad, rad_buf.as_mut()) {
                if !placements.is_empty() {
                    let rec = build_rad_record(&placements);
                    buf.write(&rec, rad.context()).map_err(|e| {
                        paraseq::Error::from(std::io::Error::other(format!(
                            "serializing RAD record: {e:#}"
                        )))
                    })?;
                }
            }
            if let Some(acc) = sh.discrete_fld {
                record_discrete(&sh, &placements, acc, counters);
            } else {
                record(
                    &sh,
                    &placements,
                    log_fm,
                    &fld_pmf,
                    &fld_cmf,
                    counters,
                    scratch,
                    seqbias.as_mut(),
                    gcbias.as_mut(),
                    posbias.as_mut(),
                );
            }
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        // Refresh the FLD online PMF snapshot at the mini-batch boundary so the
        // next batch's per-fragment reads are stable (see `record`). No-op once the
        // final FLD cache is taken; amortized over the whole batch.
        self.shared.fld.refresh_online();
        flush_unmapped(self)?;
        flush_sam(self)?;
        flush_bam(self)?;
        flush_rad(self)?;
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        merge_counters(self);
        merge_bias(self);
        flush_unmapped(self)?;
        flush_sam(self)?;
        flush_bam(self)?;
        flush_rad(self)?;
        Ok(())
    }
}

// The single-end entry point: the same work with one read per fragment, so no
// pairing, no fragment length, and every placement is `SingleEnd`.
impl<'a, 'r> ParallelProcessor<RefRecord<'r>> for QuantProcessor<'a> {
    /// Map and record one worker's batch of single-end reads.
    fn process_record_batch(
        &mut self,
        records: impl Iterator<Item = RefRecord<'r>>,
    ) -> paraseq::Result<()> {
        let mut busy = self.shared.busy.timer();
        let QuantProcessor {
            shared,
            counters,
            scratch,
            hs,
            sketch_scratch,
            seqbias,
            gcbias,
            posbias,
            unmapped,
            sam_buf,
            bam_scratch,
            rad_buf,
        } = self;
        let sh = *shared;
        let idx = sh.salmon.inner();
        if hs.is_none() {
            *hs = Some(HitSearcher::new(idx));
        }
        let hs = hs.as_mut().unwrap();
        if sh.sketch && sketch_scratch.is_none() {
            *sketch_scratch = Some(salmon_map::SketchScratch::new(idx.k()));
        }
        let log_fm = if sh.use_online {
            sh.online.map_or(0.0, |o| o.next_log_fm())
        } else {
            0.0
        };
        // One FLD snapshot pair for the whole batch (see the paired-end loop).
        let fld_pmf = sh.fld.online_snapshot();
        let fld_cmf = sh.fld.online_cmf_snapshot();
        // Reused across the batch's reads (the `*_into` calls clear it).
        let mut placements: Vec<ScoredMapping> = Vec::new();
        for rec in records {
            busy.tick(); // see the paired-end loop
            let s = rec.seq();
            if sh.sketch {
                map_single_read_sketch_into(
                    &mut placements,
                    idx,
                    hs,
                    sketch_scratch.as_mut().unwrap(),
                    s.as_ref(),
                    sh.skip,
                    sh.map_cfg.collect.max_hit_occ,
                    sh.max_read_occ,
                );
            } else {
                map_single_read_into(
                    &mut placements,
                    idx,
                    hs,
                    sh.salmon,
                    s.as_ref(),
                    sh.map_cfg,
                    &mut scratch.map,
                );
            }
            // Sketch decoy policy (see the paired-end branch): drop decoy tids /
            // decoy-dominated fragments that SA mode handles inside finalize.
            if sh.sketch && sh.salmon.info().num_decoys > 0 {
                let decoy_dominated =
                    salmon_map::filter_sketch_decoys(&mut placements, sh.salmon, &sh.map_cfg.score);
                if decoy_dominated {
                    counters.decoy += 1;
                }
            }
            if placements.len() > sh.max_read_occ {
                placements.clear();
            }
            if placements.is_empty() && sh.unmapped_names.is_some() {
                push_unmapped_name(unmapped, rec.id());
            }
            if !sh.sketch {
                accumulate_vm_stats(placements.is_empty(), counters);
            }
            if sh.sam.is_some() && !placements.is_empty() {
                crate::sam::write_fragment(
                    sam_buf,
                    sh.salmon,
                    rec.id(),
                    s.as_ref(),
                    None,
                    &placements,
                );
            }
            if let Some(scratch) = bam_scratch.as_mut() {
                if !placements.is_empty() {
                    scratch.write_fragment(sh.salmon, rec.id(), s.as_ref(), None, &placements)?;
                }
            }
            if let (Some(rad), Some(buf)) = (sh.rad, rad_buf.as_mut()) {
                if !placements.is_empty() {
                    let radrec = build_rad_record(&placements);
                    buf.write(&radrec, rad.context()).map_err(|e| {
                        paraseq::Error::from(std::io::Error::other(format!(
                            "serializing RAD record: {e:#}"
                        )))
                    })?;
                }
            }
            if let Some(acc) = sh.discrete_fld {
                record_discrete(&sh, &placements, acc, counters);
            } else {
                record(
                    &sh,
                    &placements,
                    log_fm,
                    &fld_pmf,
                    &fld_cmf,
                    counters,
                    scratch,
                    seqbias.as_mut(),
                    gcbias.as_mut(),
                    posbias.as_mut(),
                );
            }
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        // Refresh the FLD online PMF snapshot at the mini-batch boundary so the
        // next batch's per-fragment reads are stable (see `record`). No-op once the
        // final FLD cache is taken; amortized over the whole batch.
        self.shared.fld.refresh_online();
        flush_unmapped(self)?;
        flush_sam(self)?;
        flush_bam(self)?;
        flush_rad(self)?;
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        merge_counters(self);
        merge_bias(self);
        flush_unmapped(self)?;
        flush_sam(self)?;
        flush_bam(self)?;
        flush_rad(self)?;
        Ok(())
    }
}

/// Flush this thread's accumulated SAM buffer to the shared writer (under lock).
/// Hand this worker's accumulated SAM text to the shared writer and clear the
/// buffer. Called at batch boundaries, so the lock is taken once per batch
/// rather than once per record.
fn flush_sam(proc: &mut QuantProcessor) -> paraseq::Result<()> {
    if let Some(sw) = proc.shared.sam {
        if !proc.sam_buf.is_empty() {
            sw.write_block(&proc.sam_buf)?;
            proc.sam_buf.clear();
        }
    }
    Ok(())
}

/// Hand this thread's accumulated BAM chunk to the writer thread. No-op unless
/// `--writeBam` is set; the scratch carries its own writer borrow.
/// Same for BAM: hand the encoded chunk to the writer thread's queue.
fn flush_bam(proc: &mut QuantProcessor) -> paraseq::Result<()> {
    if let Some(scratch) = proc.bam_scratch.as_mut() {
        scratch.flush()?;
    }
    Ok(())
}

/// Build a salmon RAD record for one fragment's finalized, decoy-filtered
/// mappings. The read-name hash uses the trimmed id (stable regardless of FASTQ
/// comment fields). `frag_len`/`mate_fw` follow piscem's bulk convention: a real
/// fragment length and the mate strand for proper pairs, the unpaired sentinel
/// otherwise. `score` is recorded for the selective-alignment profile and
/// ignored when writing the sketch profile.
fn build_rad_record(placements: &[ScoredMapping]) -> salmon_rad::SalmonBulkRecord {
    let frag_type = salmon_rad::frag_map_type::fragment_level(placements.iter().map(|m| m.status));
    let hits = placements
        .iter()
        .map(|m| {
            let paired = m.status == MateStatus::PairedEndPaired;
            // Shared bulk frag_len/pos convention (see RadHit::for_placement):
            // proper pair → real fragment_len; orphan/SE → the mapped mate's
            // read_len in the slot, so the reader recovers the bounded-CMF orphan
            // weight rather than a flat penalty.
            salmon_rad::RadHit::for_placement(
                m.tid,
                m.is_fw,
                m.r2_fw,
                m.ref_pos,
                paired,
                m.fragment_len,
                m.read_len,
                m.score,
            )
        })
        .collect();
    salmon_rad::SalmonBulkRecord::new(frag_type, hits)
}

/// Flush this thread's accumulated RAD records as one chunk to the shared writer.
///
/// Errors here are fatal and must propagate: this is the only point at which a
/// RAD write actually reaches the filesystem, so discarding the result means a
/// full disk (`ENOSPC`) or an I/O error silently yields a truncated RAD that the
/// run then reports as successful. See salmon#1105.
fn flush_rad(proc: &mut QuantProcessor) -> paraseq::Result<()> {
    if let (Some(rad), Some(buf)) = (proc.shared.rad, proc.rad_buf.as_mut()) {
        if buf.nrec() > 0 {
            let bytes = buf
                .take_bytes()
                .map_err(|e| paraseq::Error::from(std::io::Error::other(e.to_string())))?;
            rad.append_chunk_bytes(&bytes)
                .map_err(|e| paraseq::Error::from(std::io::Error::other(format!("{e:#}"))))?;
        }
    }
    Ok(())
}
