//! Parallel RAD-format input quantification (`salmon quant --rad`).
//!
//! Quantifies from a RAD file of mappings instead of FASTQ/BAM, taking a path
//! analogous to alignment mode: each fragment's placements are turned into a
//! weighted equivalence class, then the shared inference (EM/VBEM + bootstrap/
//! Gibbs) and output (`write_outputs`) run exactly as for BAM input.
//!
//! Accepts:
//!   * **piscem `map-bulk`** RAD (`bulk_with_pos`) — sketch-style, no scores, all
//!     placements equally-best (uniform eq-class weights);
//!   * **salmon-generated** RAD — sketch (uniform) or selective-alignment (per-hit
//!     alignment scores → soft-weighted like internal selective alignments).
//!
//! The RAD chunks are already collated by read, so we read them in parallel with
//! libradicl's [`ParallelRadReader`] (a producer feeding a bounded work queue)
//! and a pool of worker threads, mirroring piscem-infer's `process_parallel`.
//!
//! Quantification is **deterministic** and order-independent: there is no online
//! (dual-phase) phase. Instead the fragment-length distribution is *fixed before*
//! equivalence-class assembly, so every fragment's weight is a pure function of
//! its placements and that FLD, and the eq-class builder (sorted `finish`) +
//! range-factorization + EM are all order-independent. Two-pass when the FLD must
//! be derived, one-pass when it is already in the file:
//!
//!   * **FLD baked in the RAD header** (salmon RAD, written by `--writeRad`): load
//!     it directly — **one pass** (just build the eq-classes).
//!   * **No baked FLD** (piscem RAD, or salmon RAD without one): derive it in a
//!     first pass from **uniquely-mapped** fragments (single placement, proper
//!     pair) bucketed by orientation, then decide the library orientation from
//!     global counts (under `-A`) and build a fixed FLD — **two passes**.
//!
//! Replacing the online phase also removed its ~6× cost in RAD mode (atomic-mass
//! contention had nothing to hide behind without a mapping step).
//!
//! v1 scope: no sequence/GC/positional bias correction (RAD carries no read
//! sequences; bias would need `-t`, the expected-model machinery, and the baked
//! `initial_abundances`) and no decoy filtering for piscem RAD (salmon RAD already
//! excludes decoys at write time). These are tracked follow-ups.

use std::fs::File;
use std::io::BufReader;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{Context, Result};
use libradicl::header::RadPrelude;
use libradicl::rad_types::{MappedFragmentOrientation, TagMap};
use libradicl::readers::{ParallelRadReader, EMPTY_METACHUNK_CALLBACK};
use libradicl::record::{MappedRecord, PiscemBulkReadRecord, RecordContext};

use salmon_core::{
    is_compatible, LibraryFormat, MateStatus, PhaseTimer, ReadOrientation, ReadStrandedness,
    ReadType,
};
use salmon_eqclass::{
    range_factorize_bins, EquivalenceClassBuilder, NaiveEqBuilder, NaivePlacement, TranscriptGroup,
    NAIVE_NO_FMT,
};
use salmon_model::{
    DiscreteFld, FragLengthSource, FragmentLengthDistribution, GcFragModel, GcStore, SBModel,
    SimplePosBias,
};
use salmon_rad::{RadInputProfile, SalmonBulkRecord};

use crate::{
    asctime_now, logsumexp, AlignQuantOptions, AlignQuantResult, ExplicitFldArgs, FldPolicy,
};

/// Floor for abundances feeding the bias-collection posterior (avoid `ln(0)`).
const BIAS_MIN_ALPHA: f64 = 1e-8;

// ---------------------------------------------------------------------------
// Placement adapter
// ---------------------------------------------------------------------------

/// One placement parsed from a RAD record.
struct RadPlacement {
    tid: u32,
    /// alignment score (meaningful only for the salmon selective-alignment profile)
    score: i32,
    is_fw: bool,
    mate_fw: bool,
    status: MateStatus,
    /// leftmost reference position of the mapped read (used to bound the orphan
    /// ambiguous fragment-length probability).
    pos: u32,
    /// fragment length for proper pairs; for orphan/single-end placements this slot
    /// instead carries the mapped mate's read length (salmon RAD) or
    /// [`salmon_rad::FRAG_LEN_UNPAIRED`] when unavailable (piscem RAD).
    frag_len: u16,
}

/// A RAD record that can be turned into placements.
trait RadFragment {
    fn placements(&self) -> Vec<RadPlacement>;
}

fn mapping_type_to_status(frag_type: u8) -> MateStatus {
    use salmon_rad::frag_map_type::*;
    match frag_type {
        MAPPED_PAIR => MateStatus::PairedEndPaired,
        LEFT_ORPHAN => MateStatus::PairedEndLeft,
        RIGHT_ORPHAN => MateStatus::PairedEndRight,
        _ => MateStatus::SingleEnd,
    }
}

/// piscem `MappedFragmentOrientation` → `(read1_fw, mate_fw)`.
fn orientation_to_strands(o: MappedFragmentOrientation) -> (bool, bool) {
    match o {
        MappedFragmentOrientation::Forward => (true, false),
        MappedFragmentOrientation::Reverse => (false, false),
        MappedFragmentOrientation::ForwardReverse => (true, false),
        MappedFragmentOrientation::ReverseForward => (false, true),
        MappedFragmentOrientation::ForwardForward => (true, true),
        MappedFragmentOrientation::ReverseReverse => (false, false),
        MappedFragmentOrientation::Unknown => (true, false),
    }
}

impl RadFragment for PiscemBulkReadRecord {
    fn placements(&self) -> Vec<RadPlacement> {
        let status = mapping_type_to_status(self.frag_type);
        (0..self.refs.len())
            .map(|i| {
                let (is_fw, mate_fw) = orientation_to_strands(self.dirs[i]);
                RadPlacement {
                    tid: self.refs[i],
                    score: 0,
                    is_fw,
                    mate_fw,
                    status,
                    pos: self.positions[i],
                    frag_len: self.frag_lengths[i],
                }
            })
            .collect()
    }
}

impl RadFragment for SalmonBulkRecord {
    fn placements(&self) -> Vec<RadPlacement> {
        let status = mapping_type_to_status(self.frag_type);
        self.hits
            .iter()
            .map(|h| RadPlacement {
                tid: h.tid,
                score: h.score,
                is_fw: h.is_fw,
                mate_fw: h.mate_fw,
                status,
                pos: h.pos,
                frag_len: h.frag_len,
            })
            .collect()
    }
}

/// Best-effort observed library format for a RAD placement (proper pairs only).
/// Mirrors `frag_format` for the inward case; orientation defaults to `Toward`
/// (RAD bulk records carry only the leftmost fragment position, so outward pairs
/// cannot be distinguished — acceptable for the common inward libraries).
fn rad_frag_format(p: &RadPlacement) -> (Option<LibraryFormat>, bool, MateStatus) {
    if p.status == MateStatus::PairedEndPaired {
        let (orientation, strandedness) = if p.is_fw != p.mate_fw {
            (
                ReadOrientation::Toward,
                if p.is_fw {
                    ReadStrandedness::SA
                } else {
                    ReadStrandedness::AS
                },
            )
        } else {
            (
                ReadOrientation::Same,
                if p.is_fw {
                    ReadStrandedness::S
                } else {
                    ReadStrandedness::A
                },
            )
        };
        (
            Some(LibraryFormat::new(
                ReadType::PairedEnd,
                orientation,
                strandedness,
            )),
            p.is_fw,
            MateStatus::PairedEndPaired,
        )
    } else {
        (None, p.is_fw, p.status)
    }
}

// ---------------------------------------------------------------------------
// Per-fragment equivalence-class assembly
// ---------------------------------------------------------------------------

struct FragCfg<'a> {
    /// salmon selective-alignment profile → soft-weight by score; else uniform
    scored: bool,
    score_exp: f64,
    /// how to turn a hit's stored score into the log-weight basis (see
    /// [`salmon_rad::SCORE_KIND_TAG`]): AS soft-weight vs. quantized log-LR.
    score_kind: u8,
    /// fixed log-PMF of the fragment-length distribution, indexed by raw length.
    pmf: &'a [f64],
    /// fixed log-CMF of the FLD, indexed by raw length (also passed whole to the
    /// orphan ambiguous-length term).
    cmf: &'a [f64],
    lengths: &'a [u32],
    expected_format: Option<LibraryFormat>,
    ignore_incompat: bool,
    incompat_prior: f64,
    paired_lib: bool,
    range_factorization_bins: u32,
    eq_builder: &'a EquivalenceClassBuilder,
}

/// The eq-class log-weight basis for one hit, given the file's score
/// interpretation. Uniform (`0`) when the profile is unscored; for a raw AS
/// score, the salmon soft-weight `−scoreExp·(best−score)` (best hit → `0`, worse
/// hits penalized); for a quantized log-LR ([`salmon_rad::SCORE_KIND_LOGWEIGHT`]),
/// the log-weight itself (`score / scale`), used absolutely like the online error
/// model's `Σ(fg−bg)` (the eq-class softmax handles the relative scaling).
#[inline]
fn score_weight_basis(scored: bool, score_kind: u8, score_exp: f64, best: i32, score: i32) -> f64 {
    if !scored {
        0.0
    } else if score_kind == salmon_rad::SCORE_KIND_LOGWEIGHT {
        score as f64 / salmon_rad::SCORE_LOG_SCALE
    } else {
        -score_exp * (best - score) as f64
    }
}

/// Weight one fragment's placements and add the equivalence class. Returns
/// `true` if the fragment was assigned (had a surviving placement). A trimmed
/// analog of [`crate::process_fragment`] with no error model, bias, or online
/// phase: the FLD (`cfg.pmf`/`cfg.cmf`) is fixed for the whole run, so the weight
/// is a deterministic function of the placements alone — duplicate/paralogous
/// transcripts sharing a fragment length get identical weights, and the result is
/// independent of fragment processing order.
fn process_rad_fragment(placements: &[RadPlacement], cfg: &FragCfg) -> bool {
    if placements.is_empty() {
        return false;
    }
    let best = placements.iter().map(|p| p.score).max().unwrap_or(0);
    let at = |snap: &[f64], i: usize| -> f64 {
        if snap.is_empty() {
            0.0
        } else {
            snap[i.min(snap.len() - 1)]
        }
    };

    let mut sp_tid: Vec<u32> = Vec::with_capacity(placements.len());
    let mut sp_eq: Vec<f64> = Vec::with_capacity(placements.len());
    for p in placements {
        // eq-class log-weight basis: soft-weight by score (SA) or uniform (sketch).
        let basis = score_weight_basis(cfg.scored, cfg.score_kind, cfg.score_exp, best, p.score);
        let rl = cfg.lengths.get(p.tid as usize).copied().unwrap_or(0) as i32;
        let proper = p.status == MateStatus::PairedEndPaired
            && p.frag_len != salmon_rad::FRAG_LEN_UNPAIRED
            && p.frag_len > 0;
        let log_frag_prob = if proper {
            // length-conditioned proper-pair probability. The fixed PMF is
            // normalized by `log P(L ≤ txpLen)` (the CMF at the transcript length)
            // so short transcripts — where only the left tail of insert sizes can
            // fit — are not under-weighted; ~0 for transcripts longer than the FLD
            // support. Matches the direct/reads path the round-trip is compared to.
            let flen = (p.frag_len as i32).min(rl.max(1));
            at(cfg.pmf, flen as usize) - at(cfg.cmf, rl.max(1) as usize)
        } else if cfg.paired_lib {
            // Orphan / single-end in a paired library: bounded-CMF ambiguous
            // fragment-length probability (salmon's `getAmbigFragLengthProb`), the
            // same orphan weight the one-pass path uses, NOT a flat penalty. The
            // mapped mate's position + orientation + transcript length bound the
            // possible fragment lengths; a forward read is exact from `txp_len −
            // pos`, a reverse read uses the mate's read length carried in the
            // orphan `frag_len` slot (salmon RAD) or `0` (piscem RAD ⇒ forward
            // bound only).
            let read_len = if p.frag_len != salmon_rad::FRAG_LEN_UNPAIRED {
                p.frag_len as i32
            } else {
                0
            };
            salmon_model::ambig_frag_log_prob(cfg.cmf, p.is_fw, p.pos as i32, read_len, rl)
        } else {
            0.0
        };
        let mut aux = basis + log_frag_prob;
        if let Some(exp) = cfg.expected_format {
            let (obs, is_fw, status) = rad_frag_format(p);
            if !is_compatible(exp, obs, is_fw, status) {
                if cfg.ignore_incompat {
                    continue;
                }
                aux += cfg.incompat_prior.ln();
            }
        }
        sp_tid.push(p.tid);
        sp_eq.push(aux);
    }
    if sp_tid.is_empty() {
        return false;
    }

    // Aggregate by distinct transcript id (sorted), softmax over the logsumexp'd
    // per-transcript weights — identical to alignment mode's eq-class formation.
    let mut agg: std::collections::BTreeMap<u32, Vec<usize>> = std::collections::BTreeMap::new();
    for (k, &t) in sp_tid.iter().enumerate() {
        agg.entry(t).or_default().push(k);
    }
    let tids: Vec<u32> = agg.keys().cloned().collect();
    let eq_log: Vec<f64> = agg
        .values()
        .map(|ks| logsumexp(&ks.iter().map(|&k| sp_eq[k]).collect::<Vec<_>>()))
        .collect();
    let maxe = eq_log.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut weights: Vec<f64> = eq_log.iter().map(|&l| (l - maxe).exp()).collect();
    let wsum: f64 = weights.iter().sum();
    if wsum > 0.0 {
        for w in &mut weights {
            *w /= wsum;
        }
    }

    let group = if cfg.range_factorization_bins > 0 {
        let bins = range_factorize_bins(&weights, cfg.range_factorization_bins);
        TranscriptGroup::with_bins(tids, bins)
    } else {
        TranscriptGroup::from_sorted(tids)
    };
    cfg.eq_builder.add_group(group, weights, 1);
    true
}

// ---------------------------------------------------------------------------
// RAD file helpers
// ---------------------------------------------------------------------------

/// Open the file and read its prelude + file-tag values, returning the reader
/// positioned at the first chunk.
fn open_prelude(path: &Path) -> Result<(BufReader<File>, RadPrelude, TagMap)> {
    let mut reader = BufReader::new(
        File::open(path).with_context(|| format!("opening RAD file {}", path.display()))?,
    );
    let prelude = RadPrelude::from_bytes(&mut reader).context("parsing RAD prelude")?;
    let file_tag_map = prelude
        .file_tags
        .parse_tags_from_bytes(&mut reader)
        .context("parsing RAD file tags")?;
    Ok((reader, prelude, file_tag_map))
}

/// Reference lengths from the RAD `ref_lengths` file tag (one per reference).
fn ref_lengths_from_tags(prelude: &RadPrelude, file_tag_map: &TagMap) -> Result<Vec<u32>> {
    use libradicl::rad_types::TagValue;
    let v = match file_tag_map.get("ref_lengths") {
        Some(TagValue::ArrayU32(v)) => v.clone(),
        _ => anyhow::bail!(
            "RAD file is missing a `ref_lengths` (Array<u32>) file tag; \
             cannot compute effective lengths"
        ),
    };
    debug_assert_eq!(v.len(), prelude.hdr.ref_names.len());
    Ok(v)
}

/// Warn that a baked fragment-length distribution has made explicitly-supplied
/// `--fldMean`/`--fldSD`/`--fldMax` ineffective (issue #1062).
///
/// Accepting those flags and then provably ignoring them is how a
/// fragment-length sensitivity analysis silently becomes a no-op: every run
/// produces byte-identical output, and nothing in the logs or `meta_info.json`
/// says why. Naming the ignored flags and the remedy costs nothing and makes
/// that failure impossible to miss.
///
/// Stays quiet when the user supplied none of them, so an ordinary requant —
/// which inherits their defaults — is not nagged about flags it never used.
fn warn_baked_fld_overrides(
    explicit: &ExplicitFldArgs,
    baked: &FragmentLengthDistribution,
    source: FragLengthSource,
) {
    if !explicit.any() {
        return;
    }
    // `derive` is only a real alternative for a paired-end RAD; single-end data
    // has no fragment lengths to re-derive from, so offering it there would send
    // the user straight back to the prior branch by a longer route.
    let (provenance, remedy) = if source == FragLengthSource::RadBakedPrior {
        (
            "That distribution was baked by a single-end run, so it is itself an \
             --fldMean/--fldSD prior rather than observed fragment lengths.",
            "Pass `--fldPolicy prior` to use your values instead.",
        )
    } else {
        (
            "That distribution was observed by the run that wrote the RAD.",
            "Pass `--fldPolicy prior` to use your values instead, or \
             `--fldPolicy derive` to re-derive the distribution from this RAD's \
             own fragment lengths.",
        )
    };
    tracing::warn!(
        "{} had no effect: this RAD carries a baked fragment-length distribution \
         (mean={:.3}, sd={:.3}), which takes precedence. {} {}",
        explicit.names(),
        baked.mean(),
        baked.sd(),
        provenance,
        remedy
    );
}

/// Load a fragment-length distribution baked into the RAD header by salmon's
/// writer, if present (the `baked_flags` byte has the FLD bit set and a non-empty
/// `frag_length_dist` log-PMF tag exists).
///
/// Returns `None` only for RAD salmon did not write — piscem or third-party.
/// salmon's writer bakes the FLD unconditionally (`--skipQuant` suppresses the
/// *abundances*, not the FLD), so every salmon RAD takes this path and the
/// caller's `--fldMean`/`--fldSD` go unused unless `--fldPolicy` says otherwise.
///
/// For a single-end write run the baked distribution is that run's own prior,
/// since no fragment lengths existed to observe; the caller distinguishes this
/// via the baked library format.
fn load_baked_fld(file_tag_map: &TagMap) -> Option<FragmentLengthDistribution> {
    use libradicl::rad_types::TagValue;
    let flags = match file_tag_map.get(salmon_rad::BAKED_FLAGS_TAG) {
        Some(TagValue::U8(f)) => *f,
        _ => return None,
    };
    if flags & salmon_rad::BAKED_FLD == 0 {
        return None;
    }
    match file_tag_map.get(salmon_rad::FRAG_LENGTH_DIST_TAG) {
        Some(TagValue::ArrayF64(pmf)) if !pmf.is_empty() => {
            Some(FragmentLengthDistribution::from_log_pmf(pmf))
        }
        _ => None,
    }
}

/// Load the initial per-reference abundance estimates baked by salmon's writer,
/// if present (the `baked_flags` byte has the abundance bit set). Used as a prior
/// for abundance-aware bias-model collection so a salmon RAD written with a quant
/// can do bias-aware requant in fewer passes. `None` ⇒ derive abundances by EM.
fn load_baked_abundances(file_tag_map: &TagMap) -> Option<Vec<f64>> {
    use libradicl::rad_types::TagValue;
    let flags = match file_tag_map.get(salmon_rad::BAKED_FLAGS_TAG) {
        Some(TagValue::U8(f)) => *f,
        _ => return None,
    };
    if flags & salmon_rad::BAKED_ABUND == 0 {
        return None;
    }
    match file_tag_map.get(salmon_rad::INITIAL_ABUNDANCES_TAG) {
        Some(TagValue::ArrayF64(a)) if !a.is_empty() => Some(a.clone()),
        _ => None,
    }
}

/// Load the resolved library format baked by salmon's writer, if present (the
/// `baked_flags` byte has the libfmt bit set). Used under `-l A` so the reader
/// can apply the same concordance filtering the write run used, without
/// re-inferring the type. `None` ⇒ auto-detect from the RAD (see [`derive_fld`]).
fn load_baked_library_format(file_tag_map: &TagMap) -> Option<LibraryFormat> {
    use libradicl::rad_types::TagValue;
    let flags = match file_tag_map.get(salmon_rad::BAKED_FLAGS_TAG) {
        Some(TagValue::U8(f)) => *f,
        _ => return None,
    };
    if flags & salmon_rad::BAKED_LIBFMT == 0 {
        return None;
    }
    match file_tag_map.get(salmon_rad::LIBRARY_FORMAT_TAG) {
        Some(TagValue::U8(id)) if *id <= LibraryFormat::MAX_FORMAT_ID => {
            Some(LibraryFormat::from_format_id(*id))
        }
        _ => None,
    }
}

/// How to interpret the per-hit `alignment_score`, from the writer's baked
/// [`salmon_rad::SCORE_KIND_TAG`]. Defaults to [`salmon_rad::SCORE_KIND_AS`] (raw
/// AS, soft-weighted) when the tag or its baked-flag bit is absent — so every
/// pre-existing salmon RAD and all piscem RAD read back unchanged.
fn load_baked_score_kind(file_tag_map: &TagMap) -> u8 {
    use libradicl::rad_types::TagValue;
    match file_tag_map.get(salmon_rad::BAKED_FLAGS_TAG) {
        Some(TagValue::U8(f)) if f & salmon_rad::BAKED_SCORE_KIND != 0 => {}
        _ => return salmon_rad::SCORE_KIND_AS,
    }
    match file_tag_map.get(salmon_rad::SCORE_KIND_TAG) {
        Some(TagValue::U8(k)) => *k,
        _ => salmon_rad::SCORE_KIND_AS,
    }
}

// ---------------------------------------------------------------------------
// Parallel weighting pass
// ---------------------------------------------------------------------------

/// Run the parallel RAD pass: `ParallelRadReader` feeds MetaChunks to a worker
/// pool that weights each fragment into the shared equivalence-class builder.
#[allow(clippy::too_many_arguments)]
fn run_rad_pass<R>(
    path: &Path,
    nthreads: usize,
    cfg: &FragCfg,
    num_processed: &AtomicU64,
    num_mapped: &AtomicU64,
) -> Result<()>
where
    R: MappedRecord + RadFragment + Send,
    R::ParsingContext: RecordContext + Clone + Send + Sync,
{
    let (reader, prelude, file_tag_map) = open_prelude(path)?;
    let nz = NonZeroUsize::new(nthreads.max(1)).unwrap();
    let mut prr =
        ParallelRadReader::<R, _>::from_prelude_and_file_tag_map(reader, prelude, file_tag_map, nz);
    let queue = prr.get_queue();
    let done = prr.is_done();

    std::thread::scope(|s| -> Result<()> {
        for _ in 0..nthreads.max(1) {
            let queue = queue.clone();
            let done = done.clone();
            s.spawn(move || {
                // Set once `done` is observed; see the drain comment below.
                let mut draining = false;
                loop {
                    if let Some(mc) = queue.pop() {
                        for chunk in mc.iter() {
                            for rec in &chunk.reads {
                                let assigned = process_rad_fragment(&rec.placements(), cfg);
                                num_processed.fetch_add(1, Ordering::Relaxed);
                                if assigned {
                                    num_mapped.fetch_add(1, Ordering::Relaxed);
                                }
                            }
                        }
                    } else if draining {
                        break;
                    } else if done.load(Ordering::Acquire) {
                        // `done` observed after an empty pop. The producer pushes
                        // every meta-chunk before storing the flag, so anything it
                        // enqueued between our failed pop and that store is still
                        // waiting. Sweep once more before exiting.
                        draining = true;
                    } else {
                        std::hint::spin_loop();
                    }
                }
            });
        }
        // Producer: fill the queue until the file is exhausted (blocks).
        prr.start_chunk_parsing(EMPTY_METACHUNK_CALLBACK)
            .context("parsing RAD chunks")
    })
}

// ---------------------------------------------------------------------------
// Deterministic FLD derivation (order-independent, unique fragments)
// ---------------------------------------------------------------------------

/// Derive a FIXED fragment-length distribution from the RAD file in one
/// order-independent pass, using only **uniquely-mapped proper pairs** (a single
/// placement). Lengths are accumulated into two orientation buckets
/// (opposite-strand vs same-strand) via the FLD's commutative `add_val`, so the
/// result is independent of chunk/worker order. The bucket with more observations
/// defines the library's concordant orientation — robust under `-A` and to a
/// mis-specified `-l` — and its FLD is returned cached. Mirrors piscem-infer's
/// unique-fragment FLD, made order-independent by using *all* unique fragments
/// plus a global-count orientation choice instead of a file-order sample quota.
fn derive_fld<R>(
    path: &Path,
    nthreads: usize,
    fld_max: usize,
    fld_mean: f64,
    fld_sd: f64,
    naive_eq: Option<&NaiveEqBuilder>,
) -> Result<(FragmentLengthDistribution, Option<LibraryFormat>)>
where
    R: MappedRecord + RadFragment + Send,
    R::ParsingContext: RecordContext + Clone + Send + Sync,
{
    // Tally uniquely-mapped proper pairs into an order-independent integer
    // accumulator (see [`DiscreteFld`]); the FLD + format are built once, after the
    // pass, deterministically. This is the SAME accumulator the deterministic
    // mapping pass uses, so a derived FLD and a baked one agree.
    let acc = DiscreteFld::new(fld_max);

    let (reader, prelude, file_tag_map) = open_prelude(path)?;
    let nz = NonZeroUsize::new(nthreads.max(1)).unwrap();
    let mut prr =
        ParallelRadReader::<R, _>::from_prelude_and_file_tag_map(reader, prelude, file_tag_map, nz);
    let queue = prr.get_queue();
    let done = prr.is_done();
    std::thread::scope(|s| -> Result<()> {
        for _ in 0..nthreads.max(1) {
            let queue = queue.clone();
            let done = done.clone();
            let acc = &acc;
            s.spawn(move || {
                // Set once `done` is observed; see the drain comment below.
                let mut draining = false;
                loop {
                    if let Some(mc) = queue.pop() {
                        for chunk in mc.iter() {
                            for rec in &chunk.reads {
                                let pls = rec.placements();
                                // Naive (kallisto-style) eq for the rough seed EM: the
                                // compatible transcripts, orientation-tagged (so library-
                                // incompatible placements can be dropped once the lib
                                // type is known), no FLD/score weighting. Built in the
                                // same read as the FLD counts; order-independent.
                                if let Some(nb) = naive_eq {
                                    let sig: Vec<NaivePlacement> = pls
                                        .iter()
                                        .map(|p| {
                                            let (obs, is_fw, status) = rad_frag_format(p);
                                            NaivePlacement {
                                                tid: p.tid,
                                                fmt_id: obs
                                                    .map(|f| f.format_id())
                                                    .unwrap_or(NAIVE_NO_FMT),
                                                is_fw,
                                                status,
                                            }
                                        })
                                        .collect();
                                    nb.add(sig);
                                }
                                // "unique" = exactly one placement; only proper pairs
                                // carry a meaningful fragment length.
                                if pls.len() != 1 {
                                    continue;
                                }
                                let p = &pls[0];
                                if p.status != MateStatus::PairedEndPaired
                                    || p.frag_len == salmon_rad::FRAG_LEN_UNPAIRED
                                    || p.frag_len == 0
                                {
                                    continue;
                                }
                                acc.add(p.frag_len as usize, p.is_fw, p.mate_fw);
                            }
                        }
                    } else if draining {
                        break;
                    } else if done.load(Ordering::Acquire) {
                        // `done` observed after an empty pop. The producer pushes
                        // every meta-chunk before storing the flag, so anything it
                        // enqueued between our failed pop and that store is still
                        // waiting. Sweep once more before exiting.
                        draining = true;
                    } else {
                        std::hint::spin_loop();
                    }
                }
            });
        }
        prr.start_chunk_parsing(EMPTY_METACHUNK_CALLBACK)
            .context("scanning RAD chunks for FLD estimation")
    })?;

    let (fld, detected) = acc.finish(fld_mean, fld_sd);
    tracing::info!(
        "derived FLD from {} unique proper pairs{}",
        acc.count(),
        match detected {
            Some(f) => format!(" (auto-detected {})", f.canonical()),
            None => String::new(),
        }
    );
    Ok((fld, detected))
}

// ---------------------------------------------------------------------------
// Bias-model collection (abundance-aware, from a fixed abundance vector)
// ---------------------------------------------------------------------------

/// Per-thread accumulator of sequence / GC / positional bias observations.
/// Mirrors alignment mode's `Local`; merged across workers with the models'
/// `combine_counts` / `combine`.
struct BiasLocal {
    seq_obs: Option<(SBModel, SBModel)>,
    gc_obs: Option<GcFragModel>,
    pos_obs: Option<(Vec<SimplePosBias>, Vec<SimplePosBias>)>,
}

impl BiasLocal {
    fn new(seq: bool, gc: bool, pos: bool, cond_gc_bins: usize, gc_bins: usize) -> Self {
        let mk = || {
            (0..salmon_model::NUM_LENGTH_CLASSES)
                .map(|_| SimplePosBias::default())
                .collect::<Vec<_>>()
        };
        Self {
            seq_obs: seq.then(|| (SBModel::new(), SBModel::new())),
            gc_obs: gc.then(|| GcFragModel::new(cond_gc_bins, gc_bins)),
            pos_obs: pos.then(|| (mk(), mk())),
        }
    }

    fn merge(&mut self, other: &BiasLocal) {
        if let (Some(a), Some(b)) = (self.seq_obs.as_mut(), other.seq_obs.as_ref()) {
            a.0.combine_counts(&b.0);
            a.1.combine_counts(&b.1);
        }
        if let (Some(a), Some(b)) = (self.gc_obs.as_mut(), other.gc_obs.as_ref()) {
            a.combine_counts(b);
        }
        if let (Some(a), Some(b)) = (self.pos_obs.as_mut(), other.pos_obs.as_ref()) {
            for (x, y) in a.0.iter_mut().zip(&b.0) {
                x.combine(y);
            }
            for (x, y) in a.1.iter_mut().zip(&b.1) {
                x.combine(y);
            }
        }
    }
}

/// Read-only context for the bias-collection pass: the same per-placement
/// weighting inputs as [`FragCfg`] plus the FIXED abundance vector and the
/// reference-derived inputs (sequence bytes, GC store, length classes).
struct BiasCfg<'a> {
    scored: bool,
    score_exp: f64,
    /// see [`FragCfg::score_kind`].
    score_kind: u8,
    pmf: &'a [f64],
    cmf: &'a [f64],
    lengths: &'a [u32],
    expected_format: Option<LibraryFormat>,
    ignore_incompat: bool,
    incompat_prior: f64,
    paired_lib: bool,
    /// fixed per-reference abundances (baked or first-EM) for the posterior.
    abundances: &'a [f64],
    ref_bytes: &'a [Vec<u8>],
    gc_store: &'a GcStore<'a>,
    length_class: Option<&'a [usize]>,
}

/// Collect one fragment's sequence/GC/positional bias contributions, weighted by
/// an abundance-aware posterior derived from the FIXED abundance vector (the RAD
/// analog of `process_fragment`'s online posterior). The reference context comes
/// from `ref_bytes` at the fragment's positions — RAD never needs the read bases.
fn collect_bias_fragment(placements: &[RadPlacement], cfg: &BiasCfg, local: &mut BiasLocal) {
    use salmon_model::seqbias::{CONTEXT_LEFT, CONTEXT_LENGTH, CONTEXT_RIGHT};
    if placements.is_empty() {
        return;
    }
    let best = placements.iter().map(|p| p.score).max().unwrap_or(0);
    let at = |snap: &[f64], i: usize| -> f64 {
        if snap.is_empty() {
            0.0
        } else {
            snap[i.min(snap.len() - 1)]
        }
    };

    let mut sp_tid: Vec<u32> = Vec::with_capacity(placements.len());
    // abundance-independent online log-weight (aux + start-position term).
    let mut sp_online: Vec<f64> = Vec::with_capacity(placements.len());
    // (frag_start, frag_end, proper, fwd_5'_pos, rev_5'_pos)
    let mut sp_geom: Vec<(usize, usize, bool, Option<usize>, Option<usize>)> =
        Vec::with_capacity(placements.len());
    for p in placements {
        let basis = score_weight_basis(cfg.scored, cfg.score_kind, cfg.score_exp, best, p.score);
        let rl = cfg.lengths.get(p.tid as usize).copied().unwrap_or(0) as i32;
        let proper = p.status == MateStatus::PairedEndPaired
            && p.frag_len != salmon_rad::FRAG_LEN_UNPAIRED
            && p.frag_len > 0;
        let flen = (p.frag_len as i32).min(rl.max(1));
        let log_frag_prob = if proper {
            at(cfg.pmf, flen as usize) - at(cfg.cmf, rl.max(1) as usize)
        } else if cfg.paired_lib {
            let read_len = if p.frag_len != salmon_rad::FRAG_LEN_UNPAIRED {
                p.frag_len as i32
            } else {
                0
            };
            salmon_model::ambig_frag_log_prob(cfg.cmf, p.is_fw, p.pos as i32, read_len, rl)
        } else {
            0.0
        };
        let start_pos = if proper && flen <= rl {
            -(((rl - flen + 1) as f64).ln())
        } else {
            -((rl.max(1) as f64).ln())
        };
        let mut aux = basis + log_frag_prob;
        if let Some(exp) = cfg.expected_format {
            let (obs, is_fw, status) = rad_frag_format(p);
            if !is_compatible(exp, obs, is_fw, status) {
                if cfg.ignore_incompat {
                    continue;
                }
                aux += cfg.incompat_prior.ln();
            }
        }
        // Fragment geometry on the transcript: for a proper pair the forward
        // mate's 5' end is the fragment's left and the reverse mate's 5' end is the
        // fragment's right; an orphan contributes only its own mate's 5'.
        let fs = p.pos as usize;
        let (fe, fwd_five, rev_five) = if proper {
            let fe = fs + p.frag_len as usize - 1;
            (fe, Some(fs), Some(fe))
        } else if p.is_fw {
            (fs, Some(fs), None)
        } else {
            let read_len = if p.frag_len != salmon_rad::FRAG_LEN_UNPAIRED {
                p.frag_len as usize
            } else {
                0
            };
            (fs, None, Some(fs + read_len.saturating_sub(1)))
        };
        sp_tid.push(p.tid);
        sp_online.push(aux + start_pos);
        sp_geom.push((fs, fe, proper, fwd_five, rev_five));
    }
    if sp_tid.is_empty() {
        return;
    }

    // Aggregate by transcript and form the abundance-aware posterior:
    // softmax over tids of `ln(alpha_tid) + online_log_tid` (the RAD stand-in for
    // `OnlineInference::assign_fragment`, using the fixed abundance vector).
    let mut agg: std::collections::BTreeMap<u32, Vec<usize>> = std::collections::BTreeMap::new();
    for (k, &t) in sp_tid.iter().enumerate() {
        agg.entry(t).or_default().push(k);
    }
    let tids: Vec<u32> = agg.keys().cloned().collect();
    let online_log: Vec<f64> = agg
        .values()
        .map(|ks| logsumexp(&ks.iter().map(|&k| sp_online[k]).collect::<Vec<_>>()))
        .collect();
    let unnorm: Vec<f64> = tids
        .iter()
        .zip(&online_log)
        .map(|(&t, &ol)| {
            cfg.abundances
                .get(t as usize)
                .copied()
                .unwrap_or(0.0)
                .max(BIAS_MIN_ALPHA)
                .ln()
                + ol
        })
        .collect();
    let denom = logsumexp(&unnorm);
    if !denom.is_finite() {
        return;
    }
    let post: Vec<f64> = unnorm.iter().map(|&u| (u - denom).exp()).collect();

    for (ti, (tid, ks)) in agg.iter().enumerate() {
        let p_tid = post[ti];
        if p_tid <= 0.0 {
            continue;
        }
        // Sequence/GC need the reference bytes; positional bias needs only the
        // transcript length, so it works even without `-t` (empty `ref_bytes`).
        let refseq = cfg
            .ref_bytes
            .get(*tid as usize)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
        let rl = if refseq.is_empty() {
            cfg.lengths.get(*tid as usize).copied().unwrap_or(0) as usize
        } else {
            refseq.len()
        };
        if rl == 0 {
            continue;
        }
        for &k in ks {
            // split the transcript posterior across its placements
            let p = p_tid * (sp_online[k] - online_log[ti]).exp();
            if p <= 0.0 {
                continue;
            }
            let (fs, fe, proper, fwd_five, rev_five) = sp_geom[k];
            if let Some(obs) = local.seq_obs.as_mut() {
                if !refseq.is_empty() {
                    if let Some(five) = fwd_five {
                        let s = five as i32 - CONTEXT_LEFT as i32;
                        if s >= 0 && (s as usize + CONTEXT_LENGTH) <= rl {
                            obs.0.add_context(
                                &refseq[s as usize..s as usize + CONTEXT_LENGTH],
                                false,
                                p,
                            );
                        }
                    }
                    if let Some(five) = rev_five {
                        let s = five as i32 - CONTEXT_RIGHT as i32;
                        if s >= 0 && (s as usize + CONTEXT_LENGTH) <= rl {
                            obs.1.add_context(
                                &refseq[s as usize..s as usize + CONTEXT_LENGTH],
                                true,
                                p,
                            );
                        }
                    }
                }
            }
            if let (Some(gc), true) = (local.gc_obs.as_mut(), proper && fe < rl) {
                let view = cfg.gc_store.view(*tid as usize);
                if let Some((ff, cf)) = salmon_model::gc_desc(&view, fs as i32, fe as i32) {
                    gc.inc(ff, cf, p);
                }
            }
            if let Some(pos) = local.pos_obs.as_mut() {
                let lc = cfg.length_class.expect("pos bias needs length classes")[*tid as usize];
                if let Some(five) = fwd_five {
                    pos.0[lc].add_mass(five as i32, rl as i32, p);
                }
                if let Some(five) = rev_five {
                    pos.1[lc].add_mass(five as i32, rl as i32, p);
                }
            }
        }
    }
}

/// Run the parallel bias-collection pass over the RAD file, returning the merged
/// sequence / GC / positional observations.
#[allow(clippy::type_complexity)]
fn run_bias_pass<R>(
    path: &Path,
    nthreads: usize,
    cfg: &BiasCfg,
    seq: bool,
    gc: bool,
    pos: bool,
    cond_gc_bins: usize,
    gc_bins: usize,
) -> Result<(
    Option<(SBModel, SBModel)>,
    Option<GcFragModel>,
    Option<(Vec<SimplePosBias>, Vec<SimplePosBias>)>,
)>
where
    R: MappedRecord + RadFragment + Send,
    R::ParsingContext: RecordContext + Clone + Send + Sync,
{
    let (reader, prelude, file_tag_map) = open_prelude(path)?;
    let nz = NonZeroUsize::new(nthreads.max(1)).unwrap();
    let mut prr =
        ParallelRadReader::<R, _>::from_prelude_and_file_tag_map(reader, prelude, file_tag_map, nz);
    let queue = prr.get_queue();
    let done = prr.is_done();
    let merged = std::sync::Mutex::new(BiasLocal::new(seq, gc, pos, cond_gc_bins, gc_bins));

    std::thread::scope(|s| -> Result<()> {
        for _ in 0..nthreads.max(1) {
            let queue = queue.clone();
            let done = done.clone();
            let merged = &merged;
            s.spawn(move || {
                let mut local = BiasLocal::new(seq, gc, pos, cond_gc_bins, gc_bins);
                // Set once `done` is observed; see the drain comment below.
                let mut draining = false;
                loop {
                    if let Some(mc) = queue.pop() {
                        for chunk in mc.iter() {
                            for rec in &chunk.reads {
                                collect_bias_fragment(&rec.placements(), cfg, &mut local);
                            }
                        }
                    } else if draining {
                        break;
                    } else if done.load(Ordering::Acquire) {
                        // `done` observed after an empty pop. The producer
                        // pushes every meta-chunk before storing the flag, so
                        // anything enqueued between our failed pop and that
                        // store is still waiting. Sweep once more before
                        // exiting.
                        draining = true;
                    } else {
                        std::hint::spin_loop();
                    }
                }
                merged.lock().unwrap().merge(&local);
            });
        }
        prr.start_chunk_parsing(EMPTY_METACHUNK_CALLBACK)
            .context("parsing RAD chunks for bias collection")
    })?;

    let merged = merged.into_inner().unwrap();
    Ok((merged.seq_obs, merged.gc_obs, merged.pos_obs))
}

/// Fused pass-2: build the (range-factorized) equivalence classes **and** collect
/// bias observations in a SINGLE RAD read, weighting bias by an abundance vector
/// known up front (baked, or a rough EM from a prepare pass). Reuses the same
/// per-fragment [`process_rad_fragment`] and [`collect_bias_fragment`] routines as
/// the separate passes (re-deriving a few per-placement scalars twice — cheap CPU
/// vs. a second RAD read). Lets the bias path avoid both the standalone first EM
/// and a second RAD read (see the pass matrix in [`quantify_rad`]).
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn run_rad_eq_and_bias_pass<R>(
    path: &Path,
    nthreads: usize,
    cfg: &FragCfg,
    bias_cfg: &BiasCfg,
    seq: bool,
    gc: bool,
    pos: bool,
    cond_gc_bins: usize,
    gc_bins: usize,
    num_processed: &AtomicU64,
    num_mapped: &AtomicU64,
) -> Result<(
    Option<(SBModel, SBModel)>,
    Option<GcFragModel>,
    Option<(Vec<SimplePosBias>, Vec<SimplePosBias>)>,
)>
where
    R: MappedRecord + RadFragment + Send,
    R::ParsingContext: RecordContext + Clone + Send + Sync,
{
    let (reader, prelude, file_tag_map) = open_prelude(path)?;
    let nz = NonZeroUsize::new(nthreads.max(1)).unwrap();
    let mut prr =
        ParallelRadReader::<R, _>::from_prelude_and_file_tag_map(reader, prelude, file_tag_map, nz);
    let queue = prr.get_queue();
    let done = prr.is_done();
    let merged = std::sync::Mutex::new(BiasLocal::new(seq, gc, pos, cond_gc_bins, gc_bins));

    std::thread::scope(|s| -> Result<()> {
        for _ in 0..nthreads.max(1) {
            let queue = queue.clone();
            let done = done.clone();
            let merged = &merged;
            s.spawn(move || {
                let mut local = BiasLocal::new(seq, gc, pos, cond_gc_bins, gc_bins);
                // Set once `done` is observed; see the drain comment below.
                let mut draining = false;
                loop {
                    if let Some(mc) = queue.pop() {
                        for chunk in mc.iter() {
                            for rec in &chunk.reads {
                                let pls = rec.placements();
                                let assigned = process_rad_fragment(&pls, cfg);
                                num_processed.fetch_add(1, Ordering::Relaxed);
                                if assigned {
                                    num_mapped.fetch_add(1, Ordering::Relaxed);
                                }
                                collect_bias_fragment(&pls, bias_cfg, &mut local);
                            }
                        }
                    } else if draining {
                        break;
                    } else if done.load(Ordering::Acquire) {
                        // `done` observed after an empty pop. The producer
                        // pushes every meta-chunk before storing the flag, so
                        // anything enqueued between our failed pop and that
                        // store is still waiting. Sweep once more before
                        // exiting.
                        draining = true;
                    } else {
                        std::hint::spin_loop();
                    }
                }
                merged.lock().unwrap().merge(&local);
            });
        }
        prr.start_chunk_parsing(EMPTY_METACHUNK_CALLBACK)
            .context("parsing RAD chunks (fused eq + bias)")
    })?;

    let merged = merged.into_inner().unwrap();
    Ok((merged.seq_obs, merged.gc_obs, merged.pos_obs))
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Quantify from a RAD file (`salmon quant --rad`).
pub fn quantify_rad(opts: &AlignQuantOptions, rad_path: &Path) -> Result<AlignQuantResult> {
    let start_time = asctime_now();
    let run_timer = std::time::Instant::now();
    let mut timer = PhaseTimer::new();

    // Bias correction is computed from the REFERENCE sequence at each fragment's
    // RAD position (never the read bases). seq/GC need the reference bases;
    // positional needs only lengths, but the shared correction tail folds all
    // three over the reference, so `-t` (the transcriptome FASTA) is required for
    // any bias model — as in reads/alignment mode. (Decoy filtering for piscem RAD
    // without `-i` is still a separate follow-up.)
    let bias_on = opts.seq_bias || opts.gc_bias || opts.pos_bias;
    anyhow::ensure!(
        !bias_on || opts.ref_seqs.is_some() || opts.transcripts.is_some(),
        "--seqBias/--gcBias/--posBias with --rad require -t/--targets (the transcriptome FASTA)"
    );

    // Header: profile, reference names + lengths (from the RAD file tags).
    let (_reader, prelude, file_tag_map) = open_prelude(rad_path)?;
    let profile = salmon_rad::detect_input_profile(&prelude)?;
    let names: Vec<String> = prelude.hdr.ref_names.clone();
    let lengths = ref_lengths_from_tags(&prelude, &file_tag_map)?;
    let num_refs = names.len();
    anyhow::ensure!(num_refs > 0, "RAD file has no reference sequences");
    anyhow::ensure!(
        lengths.len() == num_refs,
        "RAD ref_lengths ({}) does not match ref_names ({num_refs})",
        lengths.len()
    );
    let scored = matches!(
        profile,
        RadInputProfile::Salmon(salmon_rad::RadProfile::SelectiveAlignment)
    );
    tracing::info!(
        "quantifying from RAD ({}), {num_refs} references",
        match profile {
            RadInputProfile::PiscemBulk => "piscem bulk",
            RadInputProfile::Salmon(salmon_rad::RadProfile::Sketch) => "salmon sketch",
            RadInputProfile::Salmon(salmon_rad::RadProfile::SelectiveAlignment) => "salmon SA",
        }
    );

    // Library-type handling (same as alignment mode).
    let paired_lib = !matches!(opts.lib_type.as_str(), "U" | "SF" | "SR" | "S");
    let auto_lib = opts.lib_type.as_str() == "A";
    let mut expected_format = match opts.lib_type.as_str() {
        "A" => None,
        s => LibraryFormat::parse(s).ok(),
    };
    // Under `-l A` we resolve the format from (in priority order) a value baked by
    // the writer (salmon RAD — it already detected the type) or, failing that,
    // auto-detection during FLD derivation below; then apply concordance filtering
    // exactly like an explicit `-l` (reads mode does the same via its detector).
    let baked_libfmt = load_baked_library_format(&file_tag_map);
    // Score interpretation: raw AS (default) or a quantized error-model log-LR.
    let score_kind = load_baked_score_kind(&file_tag_map);
    let ignore_incompat = opts.incompat_prior <= 0.0;
    let nthreads = rayon::current_num_threads().max(1);

    // Fixed fragment-length distribution, acquired BEFORE eq-class assembly so the
    // per-fragment weights (and hence the whole pass) are deterministic.
    //
    // Priority under the default `--fldPolicy baked`:
    //   1. a RAD that baked its FLD into the header — use it directly (one pass,
    //      exact parity with the run that wrote the file). NB salmon's writer
    //      bakes the FLD *unconditionally*, `--skipQuant` included, so this is
    //      every salmon-written RAD; only piscem/third-party RAD reach step 2;
    //   2. otherwise derive it in one order-independent pass from uniquely-mapped
    //      proper pairs;
    //   3. single-end libraries have no fragment length — use the prior.
    //
    // `--fldPolicy derive` skips step 1, and `--fldPolicy prior` skips steps 1
    // and 2, which is what makes `--fldMean`/`--fldSD` meaningful on a RAD that
    // carries a baked distribution.
    let fld_baked = match opts.fld_policy {
        FldPolicy::Baked => load_baked_fld(&file_tag_map),
        FldPolicy::Derive | FldPolicy::Prior => None,
    };
    let derive_fld_pass = fld_baked.is_none() && paired_lib && opts.fld_policy != FldPolicy::Prior;
    // When the FLD must be derived (piscem / unbaked / `--fldPolicy derive`) AND we
    // need bias correction, build *naive* eq-classes in that same prepare read, so a
    // rough EM can produce the up-front abundances that let pass-2 fuse — keeping
    // this case to 2 RAD reads (prepare + fused quant) rather than 3.
    let prepare_eqb = (bias_on && !opts.skip_quant && derive_fld_pass).then(NaiveEqBuilder::new);

    let mut detected_fmt: Option<LibraryFormat> = None;
    let (fld, fld_source) = if let Some(baked) = fld_baked {
        // A baked FLD written by a single-end run is not observed data: that run
        // had no fragment lengths to measure, so what it baked is its *own*
        // `--fldMean`/`--fldSD` prior. Distinguish the two, because it changes
        // what the number means and how safe it is to override.
        let baked_single_end = baked_libfmt
            .map(|f| f.read_type == ReadType::SingleEnd)
            .unwrap_or(!paired_lib);
        let source = if baked_single_end {
            FragLengthSource::RadBakedPrior
        } else {
            FragLengthSource::RadBaked
        };
        tracing::info!("using fragment-length distribution baked in the RAD header");
        warn_baked_fld_overrides(&opts.explicit_fld_args, &baked, source);
        (baked, source)
    } else if derive_fld_pass {
        let (f, detected) = match profile {
            RadInputProfile::PiscemBulk => derive_fld::<PiscemBulkReadRecord>(
                rad_path,
                nthreads,
                opts.fld_max,
                opts.fld_mean,
                opts.fld_sd,
                prepare_eqb.as_ref(),
            )?,
            RadInputProfile::Salmon(_) => derive_fld::<SalmonBulkRecord>(
                rad_path,
                nthreads,
                opts.fld_max,
                opts.fld_mean,
                opts.fld_sd,
                prepare_eqb.as_ref(),
            )?,
        };
        detected_fmt = detected;
        (f, FragLengthSource::RadDerived)
    } else {
        if opts.fld_policy == FldPolicy::Prior && paired_lib {
            tracing::info!(
                "--fldPolicy prior: using the --fldMean/--fldSD prior alone \
                 (mean={:.3}, sd={:.3}); no fragment lengths were read from the RAD",
                opts.fld_mean,
                opts.fld_sd
            );
        }
        let mut f = FragmentLengthDistribution::new(
            1.0,
            opts.fld_max,
            opts.fld_mean,
            opts.fld_sd,
            4,
            0.5,
            1,
        );
        f.cache();
        (f, FragLengthSource::Prior)
    };
    // Resolve the `-l A` format: baked (authoritative) else auto-detected.
    if auto_lib {
        expected_format = baked_libfmt.or(detected_fmt);
        match expected_format {
            Some(ef) => tracing::info!(
                "library format under `-l A`: {} ({})",
                ef.canonical(),
                if baked_libfmt.is_some() {
                    "baked in RAD header"
                } else {
                    "auto-detected from RAD"
                }
            ),
            None => {
                tracing::info!("library format under `-l A` undetermined; no concordance filtering")
            }
        }
    }
    // Materialize the fixed log-PMF / log-CMF once as flat slices (indexed by raw
    // length), shared read-only across every fragment in the eq-class pass. The
    // PMF *is* the cached normalized log-PMF; the CMF is its cumulative companion.
    let pmf: Vec<f64> = fld.log_pmf().to_vec();
    let cmf: Vec<f64> = (0..pmf.len()).map(|l| fld.cmf(l)).collect();

    // Base effective lengths from the fixed FLD (used by the rough seed EM below
    // and by the final inference); bias correction mutates them in place later.
    let cond_means = fld.conditional_means();
    let mut eff_lengths = vec![0f64; num_refs];
    for (tid, &len) in lengths.iter().enumerate() {
        eff_lengths[tid] = salmon_model::smoothed_effective_length(&cond_means, len as usize);
    }
    // Rough seed abundances from the naive eq-classes gathered during the prepare
    // read (piscem/unbaked + bias) → the up-front abundances that let pass-2 fuse.
    let prepared_abund: Option<Vec<f64>> = prepare_eqb.map(|eqb| {
        // Drop library-incompatible placements (incompat_prior = 0) now that the
        // lib type is resolved, then build uniform eq-classes for the rough EM.
        let mut c = eqb.finish(expected_format);
        c.update_eff_lengths(&eff_lengths);
        let mut ro = opts.em.clone();
        ro.min_iter = opts.bias_seed_em_iters;
        ro.max_iter = opts.bias_seed_em_iters;
        ro.min_alpha = 0.0;
        salmon_infer::optimize(&c, num_refs, &ro, Some(&eff_lengths)).alphas
    });

    // Reference-derived inputs for bias correction (built only when a bias model
    // is requested; mirrors alignment mode). Seq/GC need the transcript bytes; the
    // GC store is a rank bitvector over the concatenated references; positional
    // bias needs only the per-transcript length classes.
    let ref_bytes: Vec<Vec<u8>> = if !bias_on {
        Vec::new()
    } else if let Some(rs) = &opts.ref_seqs {
        // Sequences supplied directly (e.g. by `--deterministic`, from the index).
        // They are in transcript-id order, which is the RAD ref-name order.
        anyhow::ensure!(
            rs.len() == names.len(),
            "supplied ref_seqs ({}) do not match the RAD reference count ({})",
            rs.len(),
            names.len()
        );
        rs.clone()
    } else {
        crate::load_ref_bytes(opts.transcripts.as_ref().unwrap(), &names)?
    };
    let (gc_rank, gc_offsets): (Option<salmon_model::GcRank>, Vec<u64>) = if opts.gc_bias {
        let mut concat: Vec<u8> = Vec::new();
        let mut offs: Vec<u64> = Vec::with_capacity(ref_bytes.len() + 1);
        offs.push(0);
        for s in &ref_bytes {
            concat.extend_from_slice(s);
            offs.push(concat.len() as u64);
        }
        (Some(salmon_model::GcRank::new(&concat)), offs)
    } else {
        (None, Vec::new())
    };
    let gc_store = match &gc_rank {
        Some(r) => salmon_model::GcStore::Rank {
            rank: r,
            offsets: &gc_offsets,
        },
        None => salmon_model::GcStore::Dense(&[]),
    };
    let length_quantiles: Option<Vec<u32>> = opts.pos_bias.then(|| {
        salmon_model::compute_length_quantiles(&lengths, salmon_model::NUM_LENGTH_CLASSES)
    });
    let length_class: Option<Vec<usize>> = length_quantiles.as_ref().map(|q| {
        lengths
            .iter()
            .map(|&l| salmon_model::length_class_index(q, l))
            .collect()
    });

    let eq_builder = EquivalenceClassBuilder::new();
    let num_processed = AtomicU64::new(0);
    let num_mapped = AtomicU64::new(0);

    // When the abundances that weight bias collection are known up FRONT (baked by
    // the write run / `--deterministic`), eq-building and bias collection can share
    // ONE RAD read and the standalone first EM is unnecessary. The no-baked-abundance
    // case keeps the separate-pass route below.
    let upfront_abund: Option<Vec<f64>> = if bias_on && !opts.skip_quant {
        // Baked abundances (salmon RAD / `--deterministic`) take precedence; else
        // the rough seed from the prepare read (piscem / unbaked-FLD).
        load_baked_abundances(&file_tag_map).or(prepared_abund)
    } else {
        None
    };

    let fused_bias_obs = {
        let cfg = FragCfg {
            scored,
            score_exp: opts.score_exp,
            score_kind,
            pmf: &pmf,
            cmf: &cmf,
            lengths: &lengths,
            expected_format,
            ignore_incompat,
            incompat_prior: opts.incompat_prior,
            paired_lib,
            range_factorization_bins: opts.range_factorization_bins,
            eq_builder: &eq_builder,
        };
        if let Some(abund) = upfront_abund.as_ref() {
            // Fused pass-2: eq-build + bias collection in one read.
            let bias_cfg = BiasCfg {
                scored,
                score_exp: opts.score_exp,
                score_kind,
                pmf: &pmf,
                cmf: &cmf,
                lengths: &lengths,
                expected_format,
                ignore_incompat,
                incompat_prior: opts.incompat_prior,
                paired_lib,
                abundances: abund,
                ref_bytes: &ref_bytes,
                gc_store: &gc_store,
                length_class: length_class.as_deref(),
            };
            let obs = match profile {
                RadInputProfile::PiscemBulk => run_rad_eq_and_bias_pass::<PiscemBulkReadRecord>(
                    rad_path,
                    nthreads,
                    &cfg,
                    &bias_cfg,
                    opts.seq_bias,
                    opts.gc_bias,
                    opts.pos_bias,
                    opts.cond_gc_bins,
                    opts.gc_bins,
                    &num_processed,
                    &num_mapped,
                )?,
                RadInputProfile::Salmon(_) => run_rad_eq_and_bias_pass::<SalmonBulkRecord>(
                    rad_path,
                    nthreads,
                    &cfg,
                    &bias_cfg,
                    opts.seq_bias,
                    opts.gc_bias,
                    opts.pos_bias,
                    opts.cond_gc_bins,
                    opts.gc_bins,
                    &num_processed,
                    &num_mapped,
                )?,
            };
            Some(obs)
        } else {
            match profile {
                RadInputProfile::PiscemBulk => run_rad_pass::<PiscemBulkReadRecord>(
                    rad_path,
                    nthreads,
                    &cfg,
                    &num_processed,
                    &num_mapped,
                )?,
                RadInputProfile::Salmon(_) => run_rad_pass::<SalmonBulkRecord>(
                    rad_path,
                    nthreads,
                    &cfg,
                    &num_processed,
                    &num_mapped,
                )?,
            }
            None
        }
    };
    let num_processed = num_processed.into_inner();
    let num_mapped = num_mapped.into_inner();
    tracing::info!("mapped {num_mapped} / {num_processed} fragments from RAD");

    // ---- inference ---------------------------------------------------------
    // Effective lengths were computed up front (from the fixed FLD) so the rough
    // seed EM could use them; they are reused here and corrected by bias below.
    let mut collapsed = eq_builder.finish();
    collapsed.update_eff_lengths(&eff_lengths);
    let num_eq_classes = collapsed.len();
    // RAD read + equivalence-class build (the input phase for RAD/deterministic/
    // alignment input; analogous to "mapping" in the reads driver).
    timer.mark("rad_read");

    // Deterministic EM from a uniform start. With a fixed FLD the collapsed
    // classes are independent of fragment/worker order, so the EM fixed point is
    // reproducible; there is no online warm-start (it was proven not to affect the
    // converged estimate, and dropping it removes the only order-dependent step).
    //
    // Pass structure for bias correction (see the matrix in the module docs):
    //  * `fused_bias_obs` is `Some` ⇒ bias was already collected alongside
    //    eq-building, weighted by up-front (baked) abundances. No standalone first
    //    EM — just correct effective lengths and run ONE final EM.
    //  * otherwise ⇒ first EM, then (if bias) a SECOND RAD read for bias
    //    collection (weighted by that EM or a baked / diagnostic prior) + re-EM.
    let mut bias_dump = salmon_model::dumps::BiasDump::default();
    let em = if opts.skip_quant {
        salmon_infer::EmResult {
            alphas: vec![0.0; num_refs],
            iters: 0,
            converged: true,
            dropped_mass: 0.0,
        }
    } else if let Some((seq_obs, gc_obs, pos_obs)) = fused_bias_obs {
        let abund = upfront_abund
            .as_ref()
            .expect("fused bias observations imply up-front abundances");
        bias_dump = crate::apply_bias_correction(
            num_refs,
            &ref_bytes,
            &gc_store,
            &lengths,
            length_class.as_deref(),
            length_quantiles.as_deref(),
            &fld,
            abund,
            &mut eff_lengths,
            seq_obs,
            gc_obs,
            pos_obs,
            opts.seq_bias,
            opts.cond_gc_bins,
            opts.gc_bins,
            opts.bias_speed_samp,
            opts.no_bias_length_threshold,
        );
        collapsed.update_eff_lengths(&eff_lengths);
        salmon_infer::optimize(&collapsed, num_refs, &opts.em, Some(&eff_lengths))
    } else {
        let mut em = salmon_infer::optimize(&collapsed, num_refs, &opts.em, Some(&eff_lengths));
        if bias_on {
            // Abundances weighting bias collection: baked prior if present, else
            // this run's (converged) first EM. (Reached only when no up-front
            // abundance was available, e.g. baked FLD without baked abundances or
            // an unpaired library; the common path fuses via `upfront_abund`.)
            let bias_abund: Vec<f64> =
                load_baked_abundances(&file_tag_map).unwrap_or_else(|| em.alphas.clone());
            let bias_cfg = BiasCfg {
                scored,
                score_exp: opts.score_exp,
                score_kind,
                pmf: &pmf,
                cmf: &cmf,
                lengths: &lengths,
                expected_format,
                ignore_incompat,
                incompat_prior: opts.incompat_prior,
                paired_lib,
                abundances: &bias_abund,
                ref_bytes: &ref_bytes,
                gc_store: &gc_store,
                length_class: length_class.as_deref(),
            };
            let (seq_obs, gc_obs, pos_obs) = match profile {
                RadInputProfile::PiscemBulk => run_bias_pass::<PiscemBulkReadRecord>(
                    rad_path,
                    nthreads,
                    &bias_cfg,
                    opts.seq_bias,
                    opts.gc_bias,
                    opts.pos_bias,
                    opts.cond_gc_bins,
                    opts.gc_bins,
                )?,
                RadInputProfile::Salmon(_) => run_bias_pass::<SalmonBulkRecord>(
                    rad_path,
                    nthreads,
                    &bias_cfg,
                    opts.seq_bias,
                    opts.gc_bias,
                    opts.pos_bias,
                    opts.cond_gc_bins,
                    opts.gc_bins,
                )?,
            };
            bias_dump = crate::apply_bias_correction(
                num_refs,
                &ref_bytes,
                &gc_store,
                &lengths,
                length_class.as_deref(),
                length_quantiles.as_deref(),
                &fld,
                &bias_abund,
                &mut eff_lengths,
                seq_obs,
                gc_obs,
                pos_obs,
                opts.seq_bias,
                opts.cond_gc_bins,
                opts.gc_bins,
                opts.bias_speed_samp,
                opts.no_bias_length_threshold,
            );
            collapsed.update_eff_lengths(&eff_lengths);
            em = salmon_infer::optimize(&collapsed, num_refs, &opts.em, Some(&eff_lengths));
        }
        em
    };
    let inference_truncated_mass = em.dropped_mass;
    let (em_iters, em_converged) = (em.iters, em.converged);
    let counts = em.alphas;
    // Bias correction + EM/VBEM point estimate (matches the reads driver's phase).
    timer.mark("em_bias");

    let rates: Vec<f64> = (0..num_refs)
        .map(|i| {
            if eff_lengths[i] > 0.0 {
                counts[i] / eff_lengths[i]
            } else {
                0.0
            }
        })
        .collect();
    let rate_sum: f64 = rates.iter().sum();
    let tpm: Vec<f64> = rates
        .iter()
        .map(|r| {
            if rate_sum > 0.0 {
                r / rate_sum * 1e6
            } else {
                0.0
            }
        })
        .collect();

    let length_classes =
        salmon_model::compute_length_quantiles(&lengths, salmon_model::NUM_LENGTH_CLASSES);

    let packed = salmon_infer::PackedEqClasses::from_collapsed(&collapsed, num_refs);
    let ambig = salmon_infer::ambiguity_counts(&packed);
    let bootstraps: Vec<Vec<f64>> = if opts.skip_quant {
        Vec::new()
    } else if opts.num_bootstraps > 0 {
        salmon_infer::bootstrap(&packed, &opts.em, opts.num_bootstraps, 0x5A13_0000)
    } else if opts.num_gibbs_samples > 0 {
        let prior = if opts.em.use_vbem {
            opts.em.vb_prior.max(1.0)
        } else {
            1e-3
        };
        let gopts = salmon_infer::GibbsOptions {
            num_samples: opts.num_gibbs_samples,
            thinning: opts.thinning_factor,
            prior,
            per_transcript_prior: true,
        };
        salmon_infer::gibbs_sample(&packed, &eff_lengths, &counts, &gopts, 0x6217_0000)
    } else {
        Vec::new()
    };
    // Posterior sampling (bootstrap / Gibbs); empty when neither is requested.
    timer.mark("posterior");

    // Run diagnostics from end-of-run aggregates (no per-fragment cost): the
    // observed (baked/auto-detected) library format powers the strandedness
    // sanity check, gated on `-l A` so an explicit `-l` is never contradicted.
    let detected_library_type = detected_fmt
        .or(baked_libfmt)
        .map(|f| f.canonical().to_string());
    let auto_detect = LibraryFormat::is_auto(&opts.lib_type);
    let requested_lib = LibraryFormat::parse(&opts.lib_type)
        .map(|f| f.canonical().to_string())
        .unwrap_or_else(|_| opts.lib_type.clone());
    let diagnostics = salmon_core::input_diagnostics(
        num_processed,
        num_mapped,
        auto_detect,
        &requested_lib,
        detected_library_type.as_deref(),
    );
    for d in &diagnostics {
        if d.severity == "error" {
            tracing::error!("{}", d.message);
        } else {
            tracing::warn!("{}", d.message);
        }
    }

    let result = AlignQuantResult {
        names,
        lengths,
        eff_lengths,
        tpm,
        counts,
        num_processed,
        num_mapped,
        inference_truncated_mass,
        num_eq_classes,
        frag_len_mean: fld.mean(),
        frag_len_sd: fld.sd(),
        frag_len_source: fld_source,
        length_classes,
        frag_len_dist: fld.log_pmf().iter().map(|lp| lp.exp()).collect(),
        start_time,
        bias_dump,
        ambig,
        bootstraps,
        em_iters,
        em_converged,
        detected_library_type,
        total_seconds: run_timer.elapsed().as_secs_f64(),
        peak_rss_kb: salmon_core::peak_rss_kb(),
        diagnostics,
    };
    crate::write_outputs(opts, &result)?;
    timer.mark("output");
    Ok(result)
}
