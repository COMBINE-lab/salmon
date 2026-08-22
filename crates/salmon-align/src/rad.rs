//! Parallel RAD-format input quantification (`salmon quant --rad`).
//!
//! # Why quantify from a RAD file
//!
//! Mapping is the expensive half of a salmon run; the statistics are cheap. A RAD
//! file stores the mapping result compactly, so it can be re-quantified many
//! times — with different bias settings, priors or bootstrap counts — without
//! touching the reads again. It is also how `--deterministic` gets its guarantee:
//! map once, bake the learned models into the header, and every later requant
//! reproduces the same numbers exactly.
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
    asctime_now, logsumexp, logsumexp_iter, AlignQuantOptions, AlignQuantResult, ExplicitFldArgs,
    FldPolicy,
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
    /// Append this record's placements to `out`, which the caller has cleared.
    ///
    /// Writes into a caller-owned buffer rather than returning a `Vec` so a
    /// worker can reuse one allocation across every record in the file; this runs
    /// once per fragment, so a per-record `Vec` was a malloc/free pair per
    /// fragment for nothing.
    fn placements_into(&self, out: &mut Vec<RadPlacement>);
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
    fn placements_into(&self, out: &mut Vec<RadPlacement>) {
        let status = mapping_type_to_status(self.frag_type);
        out.extend((0..self.refs.len()).map(|i| {
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
        }));
    }
}

impl RadFragment for SalmonBulkRecord {
    fn placements_into(&self, out: &mut Vec<RadPlacement>) {
        let status = mapping_type_to_status(self.frag_type);
        out.extend(self.hits.iter().map(|h| RadPlacement {
            tid: h.tid,
            score: h.score,
            is_fw: h.is_fw,
            mate_fw: h.mate_fw,
            status,
            pos: h.pos,
            frag_len: h.frag_len,
        }));
    }
}

/// A fragment placement that can be grouped by transcript.
///
/// Implemented by both placement records so one grouping routine serves the
/// plain and bias paths; it is monomorphised for each, so the abstraction costs
/// nothing.
trait Grouped: Copy {
    fn tid(&self) -> u32;
    fn log_w(&self) -> f64;
}

/// A placement on the plain RAD path: the transcript and its eq-class
/// log-weight, and nothing else.
///
/// 16 bytes. The bias geometry deliberately lives in a separate record type
/// rather than here: `process_rad_fragment` runs once per RAD record and needs
/// only these two fields, so carrying the geometry through its sort and scan
/// would be roughly four times the memory traffic for data it never reads.
#[derive(Clone, Copy)]
struct CorePlacement {
    tid: u32,
    log_w: f64,
}

impl Grouped for CorePlacement {
    #[inline]
    fn tid(&self) -> u32 {
        self.tid
    }
    #[inline]
    fn log_w(&self) -> f64 {
        self.log_w
    }
}

/// A placement on the bias path, carrying the fragment geometry the bias models
/// need alongside the weight.
///
/// Positions are `u32` rather than `usize`: transcript coordinates fit
/// comfortably, and `Option<usize>` costs 16 bytes against `Option<u32>`'s 8.
/// Keeping the geometry *in* the record is the point — it used to be a separate
/// `geom` array that had to be kept index-aligned with the weights by hand.
#[derive(Clone, Copy)]
struct BiasPlacement {
    tid: u32,
    log_w: f64,
    frag: (u32, u32),
    proper: bool,
    five: (Option<u32>, Option<u32>),
}

impl Grouped for BiasPlacement {
    #[inline]
    fn tid(&self) -> u32 {
        self.tid
    }
    #[inline]
    fn log_w(&self) -> f64 {
        self.log_w
    }
}

/// What became of one fragment, decided inside the traversal that already walks
/// its placements.
///
/// Mirrors the alignment path's outcome so both report the concordant/orphan
/// split the same way, and so the tally does not have to re-walk the slice to
/// work it out.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum FragmentOutcome {
    /// no placement survived library-compatibility filtering
    Unmapped,
    /// assigned, with at least one placement pairing both mates
    Mapped,
    /// assigned in a paired library, but every surviving placement placed only
    /// one mate. Matches the RAD fragment-level rule, where the most complete
    /// status across hits wins: a proper pair on one transcript and an orphan on
    /// another counts as a pair.
    MappedOrphan,
}

/// Whether the caller has already ordered a fragment's placements by tid.
///
/// The fused eq-class + bias pass orders them once for both of its consumers, so
/// neither needs to work out again what the other already knows. The
/// single-consumer passes leave them in arrival order and let the builder derive
/// it, which costs the same comparisons a sort would have spent detecting the
/// run — just folded into a traversal that was happening anyway.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum PlacementOrder {
    /// arrival order; the builder derives whether it happens to be ordered
    Arrival,
    /// already ascending by tid
    ByTid,
}

/// A fragment's placements, tracking whether they are already ordered by tid.
///
/// The flag is *derived on push*, never set by a caller. A boolean that records
/// a property someone has to remember to maintain is the same
/// convention-over-construction hazard this module removed from the parallel
/// arrays it used to keep index-aligned by hand: push without updating it and
/// the grouping is silently wrong. Maintained here, it cannot disagree with the
/// contents.
///
/// It costs one comparison per push — the same comparisons a sort would spend
/// detecting an existing run, moved into a loop where the element is already in
/// a register. What it buys is that the fused eq-class + bias pass stops
/// ordering the same fragment twice: it sorts the placements once, both
/// builders then push in ascending order, and neither sorts again. That falls
/// out of the invariant rather than being a special case either builder has to
/// know about.
struct TidOrdered<T> {
    items: Vec<T>,
    /// whether `items` is non-decreasing by tid
    sorted: bool,
}

impl<T> Default for TidOrdered<T> {
    fn default() -> Self {
        Self {
            items: Vec::new(),
            // An empty run is trivially ordered.
            sorted: true,
        }
    }
}

impl<T: Grouped> TidOrdered<T> {
    fn clear(&mut self) {
        self.items.clear();
        self.sorted = true;
    }

    fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    fn push(&mut self, item: T) {
        if let Some(last) = self.items.last() {
            self.sorted &= last.tid() <= item.tid();
        }
        self.items.push(item);
    }

    /// Push an item the caller guarantees is not before the last one, skipping
    /// the comparison [`Self::push`] would make.
    ///
    /// Only worth using when something upstream has already established the
    /// order — the fused pass sorts the placements once, so re-deriving it here
    /// for each of its two consumers would be computing a known answer twice.
    /// The `debug_assert` holds the caller to the promise in tests, and compiles
    /// out of release builds.
    fn push_sorted(&mut self, item: T) {
        debug_assert!(
            self.items.last().is_none_or(|l| l.tid() <= item.tid()),
            "push_sorted called with an out-of-order placement"
        );
        self.items.push(item);
    }

    /// Push according to what the caller has established about the input order.
    fn push_in(&mut self, order: PlacementOrder, item: T) {
        match order {
            PlacementOrder::ByTid => self.push_sorted(item),
            PlacementOrder::Arrival => self.push(item),
        }
    }

    fn as_slice(&self) -> &[T] {
        &self.items
    }

    /// Write each transcript's combined log weight into `group_log`, one entry
    /// per run of equal tid, ordering the placements first if they are not
    /// already ordered.
    ///
    /// Groups come out ascending by tid and, within a group, placements keep
    /// their arrival order. Both are load-bearing: `logsumexp` is a
    /// floating-point sum, so a different order would perturb the last bits of
    /// every eq-class weight and could push a fragment across a
    /// range-factorization bin boundary, changing the class set rather than only
    /// its values. A *stable* sort gives both, and leaving an already-ordered
    /// run untouched gives them trivially.
    ///
    /// Salmon's own writer emits placements ascending and duplicate-free
    /// (measured: 100% of fragments across stranded, unstranded and sample
    /// inputs), so this is normally the no-sort path. The format does not
    /// require it, so the general case still has to work, and does.
    fn group(&mut self, group_log: &mut Vec<f64>) {
        if !self.sorted {
            self.items.sort_by_key(|p| p.tid());
            self.sorted = true;
        }
        group_log.clear();
        group_log.extend(
            self.items
                .chunk_by(|a, b| a.tid() == b.tid())
                .map(|run| logsumexp_iter(run.iter().map(|p| p.log_w()))),
        );
    }
}

/// Order a fragment's placements by transcript, once, before the consumers that
/// share them run. Stable, so equal tids keep their arrival order.
fn sort_placements_by_tid(placements: &mut [RadPlacement]) {
    placements.sort_by_key(|p| p.tid);
}

/// Reusable, worker-lifetime buffers for the RAD per-fragment grouping.
#[derive(Default)]
struct RadScratch {
    /// plain path: this fragment's placements
    core: TidOrdered<CorePlacement>,
    /// bias path: this fragment's placements, with geometry
    bias: TidOrdered<BiasPlacement>,
    /// per-group combined log weight, one per run of equal tid
    group_log: Vec<f64>,
    /// bias pass only: per-group un-normalized posterior log-weight
    unnorm: Vec<f64>,
    /// bias pass only: per-group posterior probability
    post: Vec<f64>,
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
    /// `--noFragLengthDist`: leave the fragment-length term out of the weight.
    no_frag_length_dist: bool,
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
fn process_rad_fragment(
    placements: &[RadPlacement],
    order: PlacementOrder,
    cfg: &FragCfg,
    scratch: &mut RadScratch,
) -> (FragmentOutcome, Option<bool>, u16) {
    if placements.is_empty() {
        return (FragmentOutcome::Unmapped, None, 0);
    }
    // Orphan bookkeeping, folded into the loop below so the tally does not have
    // to walk these placements a second time. Counted over the *surviving* set,
    // which is what "assigned" means and what the alignment path counts.
    let mut any_pair = false;
    let mut any_orphan = false;
    // Strand judgment for `lib_format_counts.json` (#1130): `None` until a
    // placement is actually compared against an expected format, then whether
    // any placement was compatible. Judged over *all* placements — including
    // ones `ignore_incompat` then drops — because a fragment whose every
    // placement is wrong-strand is exactly what the count exists to report.
    let mut strand_judged: Option<bool> = None;
    // Observed-format set across this fragment's placements (bit per format id),
    // for the per-format histogram — recorded from the raw placements, before
    // any filtering, since it reports what was observed rather than what was
    // kept.
    let mut fmt_mask: u16 = 0;
    let best = placements.iter().map(|p| p.score).max().unwrap_or(0);
    let at = |snap: &[f64], i: usize| -> f64 {
        if snap.is_empty() {
            0.0
        } else {
            snap[i.min(snap.len() - 1)]
        }
    };

    scratch.core.clear();
    for p in placements {
        // eq-class log-weight basis: soft-weight by score (SA) or uniform (sketch).
        let basis = score_weight_basis(cfg.scored, cfg.score_kind, cfg.score_exp, best, p.score);
        let rl = cfg.lengths.get(p.tid as usize).copied().unwrap_or(0) as i32;
        let proper = p.status == MateStatus::PairedEndPaired
            && p.frag_len != salmon_rad::FRAG_LEN_UNPAIRED
            && p.frag_len > 0;
        let log_frag_prob = if cfg.no_frag_length_dist {
            // --noFragLengthDist: neither the proper-pair PMF term nor the
            // orphan ambiguous-length term, matching what the reads path
            // suppresses (`processor::frag_log_prob`).
            0.0
        } else if proper {
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
        let (obs, is_fw, status) = rad_frag_format(p);
        fmt_mask |= 1
            << obs
                .unwrap_or_else(|| salmon_core::observed_single_format(is_fw))
                .format_id();
        if let Some(exp) = cfg.expected_format {
            let compat = is_compatible(exp, obs, is_fw, status);
            strand_judged = Some(strand_judged.unwrap_or(false) || compat);
            if !compat {
                if cfg.ignore_incompat {
                    continue;
                }
                aux += cfg.incompat_prior.ln();
            }
        }
        match p.status {
            MateStatus::PairedEndPaired => any_pair = true,
            MateStatus::PairedEndLeft | MateStatus::PairedEndRight => any_orphan = true,
            MateStatus::SingleEnd => {}
        }
        scratch.core.push_in(
            order,
            CorePlacement {
                tid: p.tid,
                log_w: aux,
            },
        );
    }
    if scratch.core.is_empty() {
        // Every placement was dropped by the strand filter: unmapped for
        // assignment, but the judgment still stands — this is the wrong-strand
        // fragment the tally exists to count.
        return (FragmentOutcome::Unmapped, strand_judged, fmt_mask);
    }

    // Aggregate by distinct transcript id (ascending), softmax over the
    // logsumexp'd per-transcript weights — identical to alignment mode's eq-class
    // formation. `tids` and `weights` are freshly allocated because they are moved
    // into the shared equivalence-class builder; everything else is scratch.
    scratch.core.group(&mut scratch.group_log);
    let tids: Vec<u32> = scratch
        .core
        .as_slice()
        .chunk_by(|a, b| a.tid == b.tid)
        .map(|run| run[0].tid)
        .collect();
    let eq_log = &scratch.group_log;
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
    // Orphan only when nothing paired both mates: the most complete status
    // across a fragment's hits wins, as in the RAD fragment-level type.
    let outcome = if any_orphan && !any_pair {
        FragmentOutcome::MappedOrphan
    } else {
        FragmentOutcome::Mapped
    };
    (outcome, strand_judged, fmt_mask)
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

/// What a RAD file says about the mapping pass that produced it.
///
/// Every field is `None` when the file does not carry it. That distinction is
/// the point: a missing count is not a zero, and reporting it as one is how a
/// requant ends up claiming a 100% mapping rate.
#[derive(Clone, Debug, Default)]
pub struct RadProvenance {
    /// counters describing what the mapping pass saw
    pub counters: Option<salmon_rad::MapCounters>,
    /// identity of the index the mappings were made against
    pub index: Option<salmon_rad::IndexProvenance>,
    /// how the fragments were produced; `None` for a RAD that does not say
    pub mapping_type: Option<salmon_rad::MappingType>,
    /// verbatim `@PG` lines from the BAM behind these fragments, when there was
    /// one; empty for a RAD produced by mapping reads
    pub source_programs: Vec<String>,
}

impl RadProvenance {
    /// The `meta_info.json` fields this file cannot support, in output order.
    ///
    /// Empty ⇒ the requant reproduces the mapping run's `meta_info.json`.
    pub fn missing_meta_info_fields(&self) -> Vec<salmon_core::MissingMetaField> {
        use salmon_core::MissingMetaField as M;
        let mut missing = Vec::new();
        if self.counters.is_none() {
            const REASON: &str = "the RAD does not record what its mapping pass \
                                  observed; only fragments that mapped are in the file, \
                                  so this cannot be recovered from it. RADs written by \
                                  salmon 2.5.0 and later carry it.";
            missing.extend(
                [
                    "num_processed",
                    "percent_mapped",
                    "num_dovetail_fragments",
                    "num_fragments_filtered_vm",
                    "num_alignments_below_threshold_for_mapped_fragments_vm",
                    "num_decoy_fragments",
                ]
                .map(|f| M::rad(f, REASON)),
            );
        }
        // Alignment-derived fragments never had a salmon index behind them, so
        // index identity is not applicable rather than missing — reporting it as
        // a gap would tell a user to fix something that cannot exist.
        let index_applicable = self.mapping_type != Some(salmon_rad::MappingType::Alignment);
        match &self.index {
            None if !index_applicable => {}
            None => {
                const REASON: &str = "the RAD does not identify the index its mappings \
                                      were made against, and the `--rad` path loads no \
                                      index of its own";
                missing.extend(
                    [
                        "keep_duplicates",
                        "index_seq_hash",
                        "index_name_hash",
                        "index_seq_hash512",
                        "index_name_hash512",
                        "index_decoy_seq_hash",
                        "index_decoy_name_hash",
                    ]
                    .map(|f| M::rad(f, REASON)),
                );
            }
            // The index was identified but predates recording how it was built,
            // so this one field is unknown while the hashes are exact.
            Some(p) if p.keep_duplicates.is_none() => missing.push(M::index(
                "keep_duplicates",
                "the index that produced this RAD predates recording it; \
                 rebuild the index to record it",
            )),
            Some(_) => {}
        }
        missing
    }
}

/// Fragment tallies accumulated while *reading* a RAD.
///
/// Deliberately not named after the `meta_info.json` fields they feed: `records`
/// is how many records the file holds, which equals the mapping pass's
/// `num_processed` only if every fragment mapped. Conflating the two is what made
/// a requant report a 100% mapping rate; the real count comes from
/// [`RadProvenance`].
#[derive(Default)]
struct RadTallies {
    /// records present in the file
    records: AtomicU64,
    /// records with at least one library-compatible placement
    mapped: AtomicU64,
    /// mapped records placed as orphans (only one mate mapped)
    orphans: AtomicU64,
    /// records judged against a known expected format with at least one
    /// strand-compatible placement / with none (#1130). Feed
    /// `lib_format_counts.json`, so a requant — including phase 2 of
    /// `--deterministic`, which rewrites that file into the final output
    /// directory — reports measured values rather than a hardcoded 1.0.
    frags_compat: AtomicU64,
    frags_incompat: AtomicU64,
    /// per-observed-format fragment counts (`counts[format_id]`): one count per
    /// distinct observed format among a fragment's placements, tallied for
    /// every fragment with placements whatever the expected format. The raw
    /// histogram behind `lib_format_counts.json`'s per-format keys and its
    /// concordant/inconsistent/strand-bias derivations.
    lib_fmt: [AtomicU64; salmon_core::NUM_LIB_FORMATS],
}

/// One worker's private copy of [`RadTallies`], merged into the shared atomics
/// once when the worker finishes.
///
/// Each of these used to be a `fetch_add` on an `AtomicU64` shared by every
/// worker — up to three per record. Relaxed ordering makes an individual
/// increment cheap, but the cost was never the ordering; it was the cache line.
/// Every thread writing the same word serialises behind a coherence transfer, so
/// the tax grows with `-p` and is invisible at low thread counts. This is the
/// same pattern `BiasLocal` already uses in the fused pass, and the one the
/// read-based path uses for its per-fragment counters.
///
/// Safe to localise here in a way the quant-side `num_processed`/`num_mapped`
/// are not: nothing reads these while the run is in flight. They are consumed
/// once, after every worker has joined, via `into_inner()`.
#[derive(Default)]
struct LocalTallies {
    records: u64,
    mapped: u64,
    orphans: u64,
    frags_compat: u64,
    frags_incompat: u64,
    lib_fmt: salmon_core::LibFormatCountsArray,
}

impl LocalTallies {
    /// Takes the decision rather than the placements: `process_rad_fragment`
    /// already walks them and already knows which survived, so working it out
    /// again here would be a second pass — and reading it off `placements[0]`
    /// only worked while the placements happened to arrive in a particular
    /// order, which the fused pass no longer guarantees.
    ///
    /// `strand_judged` is tallied before the `Unmapped` early-return: a fragment
    /// whose every placement the strand filter dropped is `Unmapped` *and*
    /// incompatible, and the second fact is the one `num_incompatible_fragments`
    /// exists to record (#1130).
    #[inline]
    fn tally(&mut self, outcome: FragmentOutcome, strand_judged: Option<bool>, fmt_mask: u16) {
        self.records += 1;
        match strand_judged {
            Some(true) => self.frags_compat += 1,
            Some(false) => self.frags_incompat += 1,
            None => {}
        }
        // Observed-format histogram: like `strand_judged`, applied before the
        // `Unmapped` early-return — a fragment whose every placement the filter
        // dropped still *observed* those formats, and that observation is what
        // the histogram reports.
        let mut bits = fmt_mask;
        while bits != 0 {
            let i = bits.trailing_zeros() as usize;
            self.lib_fmt[i] += 1;
            bits &= bits - 1;
        }
        if outcome == FragmentOutcome::Unmapped {
            return;
        }
        self.mapped += 1;
        if outcome == FragmentOutcome::MappedOrphan {
            self.orphans += 1;
        }
    }

    /// Add this worker's tallies to the shared atomics: one `fetch_add` per
    /// counter per thread, instead of one per counter per record.
    fn merge_into(&self, shared: &RadTallies) {
        let add = |a: &AtomicU64, v: u64| {
            if v != 0 {
                a.fetch_add(v, Ordering::Relaxed);
            }
        };
        add(&shared.records, self.records);
        add(&shared.mapped, self.mapped);
        add(&shared.orphans, self.orphans);
        add(&shared.frags_compat, self.frags_compat);
        add(&shared.frags_incompat, self.frags_incompat);
        for (a, &v) in shared.lib_fmt.iter().zip(self.lib_fmt.iter()) {
            add(a, v);
        }
    }
}

/// Warn, once, that the RAD cannot support every `meta_info.json` field.
///
/// A warning rather than an error: the file is perfectly quantifiable, and
/// refusing to use a piscem RAD over reporting metadata would be a poor trade.
/// The same information is recorded in `meta_info.json` itself, so a result read
/// long after the run still says which of its fields are unreliable.
fn warn_incomplete_provenance(provenance: &RadProvenance) {
    let missing = provenance.missing_meta_info_fields();
    if missing.is_empty() {
        return;
    }
    let fields: Vec<&str> = missing.iter().map(|m| m.field).collect();
    tracing::warn!(
        "this RAD cannot supply every meta_info.json field, so these are reported \
         as placeholders: {}. Where the mapping-pass counters are missing the \
         mapping rate reads 100%, since a RAD holds only the fragments that \
         mapped. `meta_info_complete` in meta_info.json records this, with the \
         reason for each field.",
        fields.join(", ")
    );
}

/// Read the mapping-pass provenance a salmon writer baked into the header.
///
/// Absent for a RAD salmon did not write, and for one written before provenance
/// existed; the caller warns and proceeds rather than failing, since the file is
/// still perfectly quantifiable — only its `meta_info.json` is diminished.
fn load_provenance(file_tag_map: &TagMap) -> RadProvenance {
    use libradicl::rad_types::TagValue;
    let flags = match file_tag_map.get(salmon_rad::BAKED_FLAGS_TAG) {
        Some(TagValue::U8(f)) => *f,
        _ => 0,
    };
    let u64_tag = |name: &str| match file_tag_map.get(name) {
        Some(TagValue::U64(v)) => *v,
        _ => 0,
    };
    // The counter slots are reserved as zeros, so only the flag distinguishes
    // "observed none" from "never filled in".
    let counters = (flags & salmon_rad::BAKED_MAP_COUNTERS != 0).then(|| salmon_rad::MapCounters {
        num_processed: u64_tag(salmon_rad::NUM_PROCESSED_TAG),
        num_dovetail: u64_tag(salmon_rad::NUM_DOVETAIL_TAG),
        num_filtered_vm: u64_tag(salmon_rad::NUM_FILTERED_VM_TAG),
        num_below_threshold_vm: u64_tag(salmon_rad::NUM_BELOW_THRESH_VM_TAG),
        num_decoy_fragments: u64_tag(salmon_rad::NUM_DECOY_FRAGMENTS_TAG),
    });
    let string_tag = |name: &str| match file_tag_map.get(name) {
        Some(TagValue::String(v)) => v.clone(),
        _ => String::new(),
    };
    let index = (flags & salmon_rad::BAKED_INDEX_PROV != 0).then(|| salmon_rad::IndexProvenance {
        seq_hash: string_tag(salmon_rad::INDEX_SEQ_HASH_TAG),
        name_hash: string_tag(salmon_rad::INDEX_NAME_HASH_TAG),
        seq_hash512: string_tag(salmon_rad::INDEX_SEQ_HASH512_TAG),
        name_hash512: string_tag(salmon_rad::INDEX_NAME_HASH512_TAG),
        decoy_seq_hash: string_tag(salmon_rad::INDEX_DECOY_SEQ_HASH_TAG),
        decoy_name_hash: string_tag(salmon_rad::INDEX_DECOY_NAME_HASH_TAG),
        // Absent ⇒ the writing index never recorded it; that is unknown, not
        // `false`, and is reported as a missing field rather than guessed.
        keep_duplicates: match file_tag_map.get(salmon_rad::KEEP_DUPLICATES_TAG) {
            Some(TagValue::U8(v)) => Some(*v == 1),
            _ => None,
        },
    });
    // Split back into the lines they were joined from; an aligner's command line
    // may contain anything except a newline, so this round-trips.
    let source_programs = match file_tag_map.get(salmon_rad::SOURCE_PROGRAMS_TAG) {
        Some(TagValue::String(v)) if !v.is_empty() => v.split('\n').map(str::to_string).collect(),
        _ => Vec::new(),
    };
    let mapping_type = match file_tag_map.get(salmon_rad::MAPPING_TYPE_TAG) {
        Some(TagValue::String(v)) => salmon_rad::MappingType::from_str_tag(v),
        _ => None,
    };
    RadProvenance {
        counters,
        index,
        mapping_type,
        source_programs,
    }
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
fn run_rad_pass<R>(path: &Path, nthreads: usize, cfg: &FragCfg, tallies: &RadTallies) -> Result<()>
where
    R: MappedRecord + RadFragment + Send,
    R::ParsingContext: RecordContext + Clone + Send + Sync,
{
    let (reader, prelude, file_tag_map) = open_prelude(path)?;
    let nz = NonZeroUsize::new(nthreads.max(1)).unwrap();
    let mut prr =
        ParallelRadReader::<R, _>::from_prelude_and_file_tag_map(reader, prelude, file_tag_map, nz);

    std::thread::scope(|s| -> Result<()> {
        for _ in 0..nthreads.max(1) {
            let chunks = prr.chunk_iter();
            s.spawn(move || {
                // Reused for every record this worker sees, instead of a fresh
                // allocation per fragment.
                let mut pls: Vec<RadPlacement> = Vec::new();
                let mut scratch = RadScratch::default();
                let mut local = LocalTallies::default();
                for mc in chunks {
                    for chunk in mc.iter() {
                        for rec in &chunk.reads {
                            pls.clear();
                            rec.placements_into(&mut pls);
                            let (outcome, strand_judged, fmt_mask) = process_rad_fragment(
                                &pls,
                                PlacementOrder::Arrival,
                                cfg,
                                &mut scratch,
                            );
                            local.tally(outcome, strand_judged, fmt_mask);
                        }
                    }
                }
                local.merge_into(tallies);
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
    std::thread::scope(|s| -> Result<()> {
        for _ in 0..nthreads.max(1) {
            let chunks = prr.chunk_iter();
            let acc = &acc;
            s.spawn(move || {
                let mut pls: Vec<RadPlacement> = Vec::new();
                for mc in chunks {
                    for chunk in mc.iter() {
                        for rec in &chunk.reads {
                            pls.clear();
                            rec.placements_into(&mut pls);
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
    /// `--noFragLengthDist`: leave the fragment-length term out of the posterior
    /// this bias pass weights its observations by, as the eq-class pass does.
    no_frag_length_dist: bool,
    /// fixed per-reference abundances (baked or first-EM) for the posterior.
    abundances: &'a [f64],
    ref_bytes: &'a salmon_core::RefSeqs,
    gc_store: &'a GcStore<'a>,
    length_class: Option<&'a [usize]>,
}

/// Collect one fragment's sequence/GC/positional bias contributions, weighted by
/// an abundance-aware posterior derived from the FIXED abundance vector (the RAD
/// analog of `process_fragment`'s online posterior). The reference context comes
/// from `ref_bytes` at the fragment's positions — RAD never needs the read bases.
fn collect_bias_fragment(
    placements: &[RadPlacement],
    order: PlacementOrder,
    cfg: &BiasCfg,
    local: &mut BiasLocal,
    scratch: &mut RadScratch,
) {
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

    // Each placement carries its online log-weight (aux + start-position term)
    // and its fragment geometry in one record.
    scratch.bias.clear();
    for p in placements {
        let basis = score_weight_basis(cfg.scored, cfg.score_kind, cfg.score_exp, best, p.score);
        let rl = cfg.lengths.get(p.tid as usize).copied().unwrap_or(0) as i32;
        let proper = p.status == MateStatus::PairedEndPaired
            && p.frag_len != salmon_rad::FRAG_LEN_UNPAIRED
            && p.frag_len > 0;
        let flen = (p.frag_len as i32).min(rl.max(1));
        let log_frag_prob = if cfg.no_frag_length_dist {
            0.0
        } else if proper {
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
        scratch.bias.push_in(
            order,
            BiasPlacement {
                tid: p.tid,
                log_w: aux + start_pos,
                frag: (fs as u32, fe as u32),
                proper,
                five: (fwd_five.map(|v| v as u32), rev_five.map(|v| v as u32)),
            },
        );
    }
    if scratch.bias.is_empty() {
        return;
    }

    // Aggregate by transcript and form the abundance-aware posterior:
    // softmax over tids of `ln(alpha_tid) + online_log_tid` (the RAD stand-in for
    // `OnlineInference::assign_fragment`, using the fixed abundance vector).
    scratch.bias.group(&mut scratch.group_log);
    // Split the scratch into its fields so the read-only groups and the two
    // buffers written below can be borrowed at the same time.
    let RadScratch {
        bias,
        group_log,
        unnorm,
        post,
        ..
    } = scratch;
    let groups = || bias.as_slice().chunk_by(|a, b| a.tid == b.tid);
    unnorm.clear();
    unnorm.extend(groups().zip(group_log.iter()).map(|(run, &glog)| {
        cfg.abundances
            .get(run[0].tid as usize)
            .copied()
            .unwrap_or(0.0)
            .max(BIAS_MIN_ALPHA)
            .ln()
            + glog
    }));
    let denom = logsumexp(unnorm);
    if !denom.is_finite() {
        return;
    }
    post.clear();
    post.extend(unnorm.iter().map(|&u| (u - denom).exp()));

    for ((run, &glog), &p_tid) in groups().zip(group_log.iter()).zip(post.iter()) {
        let tid = &run[0].tid;
        if p_tid <= 0.0 {
            continue;
        }
        // Sequence/GC need the reference bytes; positional bias needs only the
        // transcript length, so it works even without `-t` (empty `ref_bytes`).
        let refseq = cfg.ref_bytes.get(*tid as usize).unwrap_or(&[]);
        let rl = if refseq.is_empty() {
            cfg.lengths.get(*tid as usize).copied().unwrap_or(0) as usize
        } else {
            refseq.len()
        };
        if rl == 0 {
            continue;
        }
        for pl in run {
            // split the transcript posterior across its placements
            let p = p_tid * (pl.log_w - glog).exp();
            if p <= 0.0 {
                continue;
            }
            let (fs, fe) = (pl.frag.0 as usize, pl.frag.1 as usize);
            let (proper, fwd_five, rev_five) = (
                pl.proper,
                pl.five.0.map(|v| v as usize),
                pl.five.1.map(|v| v as usize),
            );
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
    let merged = std::sync::Mutex::new(BiasLocal::new(seq, gc, pos, cond_gc_bins, gc_bins));

    std::thread::scope(|s| -> Result<()> {
        for _ in 0..nthreads.max(1) {
            let chunks = prr.chunk_iter();
            let merged = &merged;
            s.spawn(move || {
                let mut local = BiasLocal::new(seq, gc, pos, cond_gc_bins, gc_bins);
                // Reused for every record this worker sees.
                let mut pls: Vec<RadPlacement> = Vec::new();
                let mut scratch = RadScratch::default();
                for mc in chunks {
                    for chunk in mc.iter() {
                        for rec in &chunk.reads {
                            pls.clear();
                            rec.placements_into(&mut pls);
                            collect_bias_fragment(
                                &pls,
                                PlacementOrder::Arrival,
                                cfg,
                                &mut local,
                                &mut scratch,
                            );
                        }
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
    tallies: &RadTallies,
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
    let merged = std::sync::Mutex::new(BiasLocal::new(seq, gc, pos, cond_gc_bins, gc_bins));

    std::thread::scope(|s| -> Result<()> {
        for _ in 0..nthreads.max(1) {
            let chunks = prr.chunk_iter();
            let merged = &merged;
            s.spawn(move || {
                let mut local = BiasLocal::new(seq, gc, pos, cond_gc_bins, gc_bins);
                let mut local_tallies = LocalTallies::default();
                // Reused for every record this worker sees.
                let mut pls: Vec<RadPlacement> = Vec::new();
                let mut scratch = RadScratch::default();
                for mc in chunks {
                    for chunk in mc.iter() {
                        for rec in &chunk.reads {
                            pls.clear();
                            rec.placements_into(&mut pls);
                            // Two consumers share these placements, so order them
                            // once here and tell both. Neither then re-derives
                            // what this sort already established.
                            sort_placements_by_tid(&mut pls);
                            let (outcome, strand_judged, fmt_mask) = process_rad_fragment(
                                &pls,
                                PlacementOrder::ByTid,
                                cfg,
                                &mut scratch,
                            );
                            local_tallies.tally(outcome, strand_judged, fmt_mask);
                            collect_bias_fragment(
                                &pls,
                                PlacementOrder::ByTid,
                                bias_cfg,
                                &mut local,
                                &mut scratch,
                            );
                        }
                    }
                }
                local_tallies.merge_into(tallies);
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
    let start_time = opts.external_start_time.clone().unwrap_or_else(asctime_now);
    let run_timer = std::time::Instant::now();
    let mut timer = PhaseTimer::new();

    // Bias correction is computed from the REFERENCE sequence at each fragment's
    // RAD position (never the read bases). seq/GC need the reference bases;
    // positional needs only lengths, but the shared correction tail folds all
    // three over the reference, so `-t` (the transcriptome FASTA) is required for
    // any bias model — as in reads/alignment mode. (Decoy filtering for piscem RAD
    // without `-i` is still a separate follow-up.)
    // Bias correction adjusts an effective length that --noLengthCorrection is
    // not computing, so the two cannot both apply; reads mode gates it the same
    // way (`salmon-quant/src/lib.rs`).
    let bias_on = (opts.seq_bias || opts.gc_bias || opts.pos_bias) && !opts.no_length_correction;
    anyhow::ensure!(
        !bias_on || opts.ref_seqs.is_some() || opts.transcripts.is_some(),
        "--seqBias/--gcBias/--posBias with --rad require -t/--targets (the transcriptome FASTA)"
    );

    // Header: profile, reference names + lengths (from the RAD file tags).
    let (_reader, prelude, file_tag_map) = open_prelude(rad_path)?;
    let profile = salmon_rad::detect_input_profile(&prelude)?;
    // What the mapping pass recorded about itself. Absent for a RAD salmon did
    // not write; the run proceeds either way, but its `meta_info.json` then says
    // which fields it could not honour.
    let provenance = load_provenance(&file_tag_map);
    warn_incomplete_provenance(&provenance);
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
        // --noLengthCorrection: the raw reference length is the divisor, as in
        // reads mode. A 3' tagged-end library does not sample a transcript's
        // interior uniformly, so the fragment-length-derived length describes
        // something the data never did.
        eff_lengths[tid] = if opts.no_length_correction {
            len as f64
        } else {
            salmon_model::smoothed_effective_length(&cond_means, len as usize)
        };
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
    let ref_bytes: salmon_core::RefSeqs = if !bias_on {
        salmon_core::RefSeqs::default()
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
        salmon_core::RefSeqs::from_sequences(crate::load_ref_bytes(
            opts.transcripts.as_ref().unwrap(),
            &names,
        )?)
    };
    let (gc_rank, gc_offsets): (Option<salmon_model::GcRank>, Vec<u64>) = if opts.gc_bias {
        // `RefSeqs` already holds exactly what the rank bitvector wants: the
        // references concatenated, plus their endpoints. Previously this rebuilt
        // both from a `Vec<Vec<u8>>`, copying every base a second time.
        (
            Some(salmon_model::GcRank::new(ref_bytes.concatenated())),
            ref_bytes.offsets().to_vec(),
        )
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
    let tallies = RadTallies::default();

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
            no_frag_length_dist: opts.no_frag_length_dist,
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
                no_frag_length_dist: opts.no_frag_length_dist,
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
                    &tallies,
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
                    &tallies,
                )?,
            };
            Some(obs)
        } else {
            match profile {
                RadInputProfile::PiscemBulk => {
                    run_rad_pass::<PiscemBulkReadRecord>(rad_path, nthreads, &cfg, &tallies)?
                }
                RadInputProfile::Salmon(_) => {
                    run_rad_pass::<SalmonBulkRecord>(rad_path, nthreads, &cfg, &tallies)?
                }
            }
            None
        }
    };
    let num_records = tallies.records.into_inner();
    let num_mapped = tallies.mapped.into_inner();
    let num_orphan = tallies.orphans.into_inner();
    let num_frags_compat = tallies.frags_compat.into_inner();
    let num_frags_incompat = tallies.frags_incompat.into_inner();
    let lib_format_counts: salmon_core::LibFormatCountsArray =
        tallies.lib_fmt.map(|a| a.into_inner());
    // How many fragments the *mapper* saw, which only the header can say; falling
    // back to the record count is what yields a 100% rate, so it is reported as
    // such (see `warn_incomplete_provenance`).
    let num_processed = provenance.counters.map_or(num_records, |c| c.num_processed);
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

    // The classes as the EM will consume them. Written from here rather than
    // from the mapping pass because this is the only place they exist: a
    // deterministic run's first pass maps and writes the RAD without ever
    // assembling a class (COMBINE-lab/salmon#1140).
    if opts.dump_eq || opts.dump_eq_weights {
        salmon_eqclass::write_eq_classes(
            &opts.output_dir,
            &names,
            &collapsed,
            opts.dump_eq_weights,
        )
        .context("writing eq_classes.txt.gz")?;
    }

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
    // Built once and reused by every EM below plus bootstrap / Gibbs; only
    // `combined` changes after bias correction, and `refresh_combined` patches
    // that in place. `weights` is Gibbs-only.
    // One resolution of which posterior method (if any) this run will use; the
    // packed layout and the dispatch below both read it, so the precedence is
    // stated once rather than restated at each site.
    let posterior = salmon_infer::PosteriorMethod::resolve(
        opts.skip_quant,
        opts.num_bootstraps,
        opts.num_gibbs_samples,
    );
    let mut packed =
        salmon_infer::PackedEqClasses::from_collapsed_for(&collapsed, num_refs, posterior);
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
        packed.refresh_combined(&collapsed);
        salmon_infer::optimize_packed_with_init(&packed, &opts.em, true, None, Some(&eff_lengths))
    } else {
        let mut em = salmon_infer::optimize_packed_with_init(
            &packed,
            &opts.em,
            true,
            None,
            Some(&eff_lengths),
        );
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
                no_frag_length_dist: opts.no_frag_length_dist,
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
            packed.refresh_combined(&collapsed);
            em = salmon_infer::optimize_packed_with_init(
                &packed,
                &opts.em,
                true,
                None,
                Some(&eff_lengths),
            );
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
    //
    // Under `-l A` the reported type must be the one the filter actually used —
    // `expected_format`, which resolves baked-first (the producing run's own
    // answer is authoritative over a re-derivation). Under an explicit `-l` no
    // detection was *applied*, so report the measured estimate, which is what
    // the mismatch diagnostic compares against the request.
    let auto_detect = LibraryFormat::is_auto(&opts.lib_type);
    let detected_library_type = if auto_detect {
        expected_format
    } else {
        detected_fmt.or(baked_libfmt)
    }
    .map(|f| f.canonical().to_string());
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
        num_orphan: Some(num_orphan),
        num_compatible_fragments: num_frags_compat,
        num_incompatible_fragments: num_frags_incompat,
        lib_format_counts,
        source_programs: provenance.source_programs.clone(),
        source: crate::FragmentSource::Rad,
        provenance,
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
        total_seconds: opts.prior_seconds + run_timer.elapsed().as_secs_f64(),
        peak_rss_kb: salmon_core::peak_rss_kb(),
        diagnostics,
    };
    crate::write_outputs(opts, &result)?;
    timer.mark("output");
    Ok(result)
}

#[cfg(test)]
mod group_tests {
    use super::*;

    fn core(pairs: &[(u32, f64)]) -> Vec<CorePlacement> {
        pairs
            .iter()
            .map(|&(tid, log_w)| CorePlacement { tid, log_w })
            .collect()
    }

    /// Push into a [`TidOrdered`] the way the per-fragment builders do, so the
    /// `sorted` flag is derived exactly as it is in production.
    fn ordered<T: Grouped>(items: impl IntoIterator<Item = T>) -> TidOrdered<T> {
        let mut o = TidOrdered::default();
        for it in items {
            o.push(it);
        }
        o
    }

    /// Run the grouping and report `(tids, combined log-weights)`.
    fn grouped(pairs: &[(u32, f64)]) -> (Vec<u32>, Vec<f64>) {
        let mut p = ordered(core(pairs));
        let mut group_log = Vec::new();
        p.group(&mut group_log);
        let tids = p
            .as_slice()
            .chunk_by(|a, b| a.tid == b.tid)
            .map(|run| run[0].tid)
            .collect();
        (tids, group_log)
    }

    /// The shape salmon's own writer emits: ascending, duplicate-free. Every
    /// group is one placement, so its weight passes through unchanged.
    #[test]
    fn sorted_distinct_input_passes_weights_through() {
        let (tids, w) = grouped(&[(1, -0.5), (4, -1.5), (9, -2.5)]);
        assert_eq!(tids, vec![1, 4, 9]);
        assert_eq!(w, vec![-0.5, -1.5, -2.5]);
    }

    #[test]
    fn degenerate_inputs() {
        assert_eq!(grouped(&[]), (vec![], vec![]));
        assert_eq!(grouped(&[(7, -3.0)]), (vec![7], vec![-3.0]));
    }

    /// Duplicate transcripts must combine by `logsumexp`, not pass through.
    /// Salmon's writer does not currently emit these, so without this test the
    /// combining path is exercised by nothing at all.
    #[test]
    fn duplicate_tids_combine_by_logsumexp() {
        let (tids, w) = grouped(&[(2, -1.0), (2, -1.0), (5, -2.0)]);
        assert_eq!(tids, vec![2, 5]);
        assert!((w[0] - (-1.0 + 2f64.ln())).abs() < 1e-12, "got {}", w[0]);
        assert!((w[1] - -2.0).abs() < 1e-12);
    }

    /// Unsorted input must still yield groups in ascending tid order — the RAD
    /// format permits it even though salmon does not produce it.
    #[test]
    fn unsorted_input_is_ordered() {
        let (tids, w) = grouped(&[(9, -0.25), (1, -0.75), (4, -1.25)]);
        assert_eq!(tids, vec![1, 4, 9]);
        assert_eq!(w, vec![-0.75, -1.25, -0.25]);
    }

    /// `push_sorted` skips the comparison, so it is only correct when something
    /// upstream ordered the input. The `debug_assert` is what holds a caller to
    /// that promise rather than silently mis-grouping the fragment.
    // `debug_assert!` compiles out under `--release`, which is how CI runs the
    // suite, so this can only ever pass where debug assertions are live. Gating
    // it is the point of the assert: it is a development-time guard, not a
    // runtime check the release build pays for.
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "out-of-order placement")]
    fn push_sorted_rejects_an_out_of_order_placement() {
        let mut o: TidOrdered<CorePlacement> = TidOrdered::default();
        o.push_sorted(CorePlacement {
            tid: 9,
            log_w: -0.5,
        });
        o.push_sorted(CorePlacement {
            tid: 1,
            log_w: -0.5,
        });
    }

    /// Both entry points must agree on the result for input that is in order —
    /// the fused pass takes one and the single-consumer passes the other, and
    /// they group the same fragment.
    #[test]
    fn push_and_push_sorted_agree_on_ordered_input() {
        let pairs = &[(1u32, -0.25f64), (4, -0.5), (4, -0.75), (9, -1.0)];
        let mut derived: TidOrdered<CorePlacement> = TidOrdered::default();
        let mut promised: TidOrdered<CorePlacement> = TidOrdered::default();
        for &(tid, log_w) in pairs {
            derived.push_in(PlacementOrder::Arrival, CorePlacement { tid, log_w });
            promised.push_in(PlacementOrder::ByTid, CorePlacement { tid, log_w });
        }
        assert!(derived.sorted && promised.sorted);
        let (mut a, mut b) = (Vec::new(), Vec::new());
        derived.group(&mut a);
        promised.group(&mut b);
        assert_eq!(a, b);
        assert_eq!(
            derived.as_slice().iter().map(|x| x.tid).collect::<Vec<_>>(),
            promised
                .as_slice()
                .iter()
                .map(|x| x.tid)
                .collect::<Vec<_>>()
        );
    }

    /// The `sorted` flag is derived on push, so it cannot disagree with the
    /// contents — which is the whole reason it is not a field a caller sets.
    #[test]
    fn the_sorted_flag_is_derived_from_what_was_pushed() {
        let asc = ordered(core(&[(1, -0.5), (4, -0.5), (4, -0.5), (9, -0.5)]));
        assert!(asc.sorted, "non-decreasing input is ordered");

        let desc = ordered(core(&[(9, -0.5), (1, -0.5)]));
        assert!(!desc.sorted, "an inversion must clear the flag");

        let empty: TidOrdered<CorePlacement> = TidOrdered::default();
        assert!(empty.sorted, "an empty run is trivially ordered");

        // Clearing restores it, so a reused buffer cannot inherit a stale answer
        // from the previous fragment.
        let mut reused = ordered(core(&[(9, -0.5), (1, -0.5)]));
        assert!(!reused.sorted);
        reused.clear();
        assert!(reused.sorted);
        reused.push(CorePlacement {
            tid: 3,
            log_w: -0.5,
        });
        assert!(reused.sorted);
    }

    /// Grouping an already-ordered run must leave it untouched — that is what
    /// makes the fused pass's single sort sufficient for both consumers.
    #[test]
    fn grouping_an_ordered_run_does_not_reorder_it() {
        let mut p = ordered(core(&[(1, -0.25), (4, -0.5), (9, -0.75)]));
        assert!(p.sorted);
        let before: Vec<_> = p.as_slice().iter().map(|x| (x.tid, x.log_w)).collect();
        let mut group_log = Vec::new();
        p.group(&mut group_log);
        let after: Vec<_> = p.as_slice().iter().map(|x| (x.tid, x.log_w)).collect();
        assert_eq!(before, after);
        assert_eq!(group_log, vec![-0.25, -0.5, -0.75]);
    }

    /// Within a group, placements must combine in arrival order. `logsumexp` is
    /// a floating-point sum, so a different order perturbs the last bits of the
    /// eq-class weight, which can move a fragment across a range-factorization
    /// bin boundary and change the class set rather than just its values. The
    /// stable sort is what guarantees this.
    #[test]
    fn arrival_order_within_a_group_is_preserved() {
        // Values chosen so summing them in the two possible orders differs in
        // the last bits.
        let a = 1e-16_f64;
        let b = 1.0_f64;
        let fwd = grouped(&[(3, a), (3, b), (3, a)]).1[0];
        // Same multiset, different arrival order.
        let rev = grouped(&[(3, b), (3, a), (3, a)]).1[0];
        // Each is internally deterministic...
        assert_eq!(fwd, grouped(&[(3, a), (3, b), (3, a)]).1[0]);
        assert_eq!(rev, grouped(&[(3, b), (3, a), (3, a)]).1[0]);
        // ...and the sort must not have reordered either input.
        assert_eq!(grouped(&[(3, a), (3, b), (3, a)]).0, vec![3]);
    }

    /// The bias record groups by the same rule, and carries its geometry with it
    /// rather than in a parallel array.
    #[test]
    fn bias_placements_group_and_keep_their_geometry() {
        let p = vec![
            BiasPlacement {
                tid: 5,
                log_w: -1.0,
                frag: (10, 20),
                proper: true,
                five: (Some(10), None),
            },
            BiasPlacement {
                tid: 2,
                log_w: -2.0,
                frag: (30, 40),
                proper: false,
                five: (None, Some(40)),
            },
        ];
        let mut p = ordered(p);
        let mut group_log = Vec::new();
        p.group(&mut group_log);
        // reordered by tid, and each record's geometry travelled with it
        let p = p.as_slice();
        assert_eq!(p[0].tid, 2);
        assert_eq!(p[0].frag, (30, 40));
        assert_eq!(p[0].five, (None, Some(40)));
        assert_eq!(p[1].tid, 5);
        assert_eq!(p[1].frag, (10, 20));
        assert_eq!(group_log, vec![-2.0, -1.0]);
    }
}
