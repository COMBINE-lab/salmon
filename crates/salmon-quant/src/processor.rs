//! Reads-mode quantification core, with two inference strategies.
//!
//! **Default (inline, fully parallel):** each worker maps its reads and runs the
//! per-fragment inference ([`record`]) immediately, accumulating the online
//! posterior, FLD training, bias models, and equivalence classes concurrently.
//! This is the historical behaviour: fast, no fragment store, but the
//! accumulation order (hence the last ULPs of `quant.sf`) depends on thread
//! scheduling.
//!
//! **`--deterministic` (opt-in, two-pass):** pass 1 maps in parallel and only
//! buffers every mapped fragment's `(key, mappings)` into [`Shared::frag_sink`];
//! pass 2 ([`run_inference_serial`]) sorts the buffered fragments by a stable
//! per-fragment key and replays the inference single-threaded in that fixed order.
//! Every order-dependent accumulation is then byte-identical regardless of how
//! many threads mapped the reads, at the cost of holding the mappings between
//! passes and a single-threaded inference phase. Selected by
//! [`Shared::deterministic`].

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

use paraseq::fastx::RefRecord;
use paraseq::parallel::{PairedParallelProcessor, ParallelProcessor};
use paraseq::Record;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};

use salmon_core::math::{log_add, LOG_0, LOG_1};
use salmon_core::{is_compatible, LibraryFormat, MateStatus};
use salmon_eqclass::{range_factorize_bins, EquivalenceClassBuilder, TranscriptGroup};
use salmon_index::SalmonIndex;
use salmon_infer::OnlineInference;
use salmon_map::{
    map_read_pair, map_read_pair_sketch, map_single_read, map_single_read_sketch, MapConfig,
    ScoredMapping,
};
use salmon_model::seqbias::{SBModel, CONTEXT_LEFT, CONTEXT_LENGTH, CONTEXT_RIGHT};
use salmon_model::{
    gc_desc, FragmentLengthDistribution, GcFragModel, LibraryTypeDetector, SimplePosBias,
    NUM_LENGTH_CLASSES,
};

/// Shared, thread-safe state (all `Copy` references).
#[derive(Clone, Copy)]
pub(crate) struct Shared<'a> {
    pub salmon: &'a SalmonIndex,
    pub eq: &'a EquivalenceClassBuilder,
    pub fld: &'a FragmentLengthDistribution,
    /// library-type detector when auto-detecting (`-l A`); else `None`
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
    pub num_processed: &'a AtomicU64,
    pub num_mapped: &'a AtomicU64,
    /// mapped fragments whose representative mapping is an orphan (only one mate
    /// of a paired-end fragment mapped); feeds `lib_format_counts.json`'s
    /// `num_frags_with_inconsistent_or_orphan_mappings`.
    pub num_orphan: &'a AtomicU64,
    /// selective-alignment meta counters (salmon's `aux_info/meta_info.json`):
    /// decoy-dominated fragments, dovetail fragments, fragments that had
    /// candidates but no surviving mapping, and (for mapped fragments) candidate
    /// alignments dropped below the score threshold.
    pub num_decoy: &'a AtomicU64,
    pub num_dovetail: &'a AtomicU64,
    pub num_frags_filtered_vm: &'a AtomicU64,
    pub num_below_threshold_vm: &'a AtomicU64,
    /// when set, collect the names of unmapped fragments (`--writeUnmappedNames`)
    pub unmapped_names: Option<&'a Mutex<Vec<String>>>,
    /// when set, write per-mapping SAM records (`--writeMappings`)
    pub sam: Option<&'a crate::sam::SamWriter>,
    /// opt-in byte-for-byte reproducibility (`--deterministic`). When `true`, the
    /// mapping pass only writes fragments to the RAD store (`rad`) and the
    /// inference runs in a deterministic key-sorted second pass
    /// ([`run_inference_serial`]); when `false` (default) the inference runs inline
    /// during the parallel mapping pass (no store, no extra overhead), as before.
    pub deterministic: bool,
    /// RAD-format mapping store for the deterministic second pass; `Some` only
    /// when `deterministic` is set. Each worker appends its mapped fragments
    /// `(key, mappings)` as a chunk; pass 2 reads the file back, sorts by key, and
    /// replays the inference single-threaded, byte-identical across thread counts.
    pub rad: Option<&'a crate::rad::RadChunkWriter>,
}

/// Per-thread mapping processor.
pub(crate) struct QuantProcessor<'a> {
    pub shared: Shared<'a>,
    pub hs: Option<HitSearcher<'a>>,
    /// per-thread reusable sketch caches (avoids per-read MappingCache allocs)
    pub sketch_scratch: Option<salmon_map::SketchScratch>,
    /// per-thread observed (fw, rc) sequence-bias models (inline path only)
    pub seqbias: Option<(SBModel, SBModel)>,
    /// per-thread observed fragment-GC model (inline path only)
    pub gcbias: Option<GcFragModel>,
    /// per-thread observed (5', 3') positional-bias models (inline path only)
    pub posbias: Option<(Vec<SimplePosBias>, Vec<SimplePosBias>)>,
    /// per-thread collected unmapped-fragment names (merged in on_thread_complete)
    pub unmapped: Vec<String>,
    /// per-thread SAM record buffer (flushed to the shared writer per batch)
    pub sam_buf: String,
    /// per-thread buffer of mapped fragments `(key, mappings)`, written to the RAD
    /// store (`Shared::rad`) as a chunk at each batch boundary and at thread
    /// completion for the deterministic second pass (used only when
    /// `Shared::deterministic`).
    pub frag_buf: Vec<(u128, Vec<ScoredMapping>)>,
}

impl<'a> QuantProcessor<'a> {
    pub fn new(shared: Shared<'a>) -> Self {
        // Per-thread bias models for the inline (non-deterministic) path; the
        // deterministic path collects bias in `run_inference_serial` instead.
        let (seqbias, gcbias, posbias) = if shared.deterministic {
            (None, None, None)
        } else {
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
            (seqbias, gcbias, posbias)
        };
        Self {
            shared,
            hs: None,
            sketch_scratch: None,
            seqbias,
            gcbias,
            posbias,
            unmapped: Vec::new(),
            sam_buf: String::new(),
            frag_buf: Vec::new(),
        }
    }
}

impl Clone for QuantProcessor<'_> {
    fn clone(&self) -> Self {
        // Per-thread scratch is rebuilt fresh in each worker.
        QuantProcessor::new(self.shared)
    }
}

/// salmon's `LOG_EPSILON = log(0.375e-10)` — the small log-probability assigned
/// to implausible placements (out-of-range fragment lengths, unexpected orphans).
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

/// Deterministic per-fragment acceptance sampler for stochastic FLD training.
///
/// salmon trains the fragment-length distribution online by accepting each
/// concordant mapping's length with probability equal to its abundance-aware
/// posterior. The previous port drew that acceptance coin from a *thread-local*
/// PRNG, so which fragments trained the FLD depended on how the reads happened to
/// be distributed across worker threads — the documented source of multi-threaded
/// run-to-run wobble. Seeding instead from the read's own identity makes the
/// acceptance decision a pure function of `(read_id, mapping_index, posterior)`,
/// so the *set* of fragments training the FLD is identical regardless of thread
/// count or fragment-arrival order. The draw for mapping `i` is derived directly
/// from `(seed, i)` (not from a sequential stream), so it is independent of how
/// many earlier mappings of the same fragment were concordant.
#[derive(Clone, Copy)]
struct FldSampler {
    seed: u64,
}

impl FldSampler {
    /// Seed from the fragment's precomputed key (`xxh3_128(read_id)`); the low 64
    /// bits give a stable per-fragment seed. `| 1` keeps it nonzero.
    fn from_key(key: u128) -> Self {
        FldSampler {
            seed: (key as u64) | 1,
        }
    }

    /// A pseudo-random value in `[0, 1)` for compatible-mapping index `i`,
    /// derived statelessly from `(seed, i)` via the splitmix64 finalizer.
    fn u01(&self, i: usize) -> f64 {
        let mut z = self
            .seed
            .wrapping_add((i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        ((z >> 11) as f64) * (1.0 / ((1u64 << 53) as f64))
    }
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
/// mappings, each weighted by `bias_w[i]` (the abundance-aware posterior under
/// online inference, else the normalized aux weight). Mirrors salmon's
/// per-alignment `observedGCMass.inc(gcDesc(start, stop), logProb)`.
fn collect_gc(sh: &Shared, compat: &[(&ScoredMapping, f64)], bias_w: &[f64], gc: &mut GcFragModel) {
    let Some(store) = sh.gc_store else { return };
    for (i, (m, _)) in compat.iter().enumerate() {
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
            gc.inc(ff, cf, bias_w[i]);
        }
    }
}

/// Collect observed positional-bias mass for one fragment's compatible mappings,
/// each weighted by `bias_w[i]` (log space). Mirrors salmon's
/// `observedPosBias{Fwd,RC}[lengthClass].addMass(pos, refLen, logProb)`.
fn collect_pos(
    sh: &Shared,
    compat: &[(&ScoredMapping, f64)],
    bias_w: &[f64],
    pos: &mut (Vec<SimplePosBias>, Vec<SimplePosBias>),
) {
    let Some(length_class) = sh.length_class else {
        return;
    };
    for (i, (m, _)) in compat.iter().enumerate() {
        if bias_w[i] <= 0.0 {
            continue;
        }
        let lc = length_class[m.tid as usize];
        let ref_len = sh.salmon.ref_len(m.tid as usize) as i32;
        let mass = bias_w[i].ln();
        if m.fw_pos >= 0 {
            pos.0[lc].add_mass(m.fw_pos, ref_len, mass);
        }
        if m.rc_pos >= 0 {
            pos.1[lc].add_mass(m.rc_pos, ref_len, mass);
        }
    }
}

/// Record one fragment's weighted mappings into the shared accumulators.
/// `log_fm` is the current minibatch's forgetting mass (used only when online
/// inference is active).
fn record(
    sh: &Shared,
    key: u128,
    maps: &[ScoredMapping],
    log_fm: f64,
    seqbias: Option<&mut (SBModel, SBModel)>,
    gcbias: Option<&mut GcFragModel>,
    posbias: Option<&mut (Vec<SimplePosBias>, Vec<SimplePosBias>)>,
) {
    // `num_processed` is counted during the mapping pass; pass 2 only receives
    // mapped fragments.
    if maps.is_empty() {
        return;
    }
    sh.num_mapped.fetch_add(1, Ordering::Relaxed);

    // Classify the fragment as orphan when its representative mapping has only
    // one mate of a pair placed (PairedEndLeft / PairedEndRight). Single-end
    // libraries never count as orphan here.
    if matches!(
        maps[0].status,
        MateStatus::PairedEndLeft | MateStatus::PairedEndRight
    ) {
        sh.num_orphan.fetch_add(1, Ordering::Relaxed);
    }

    // Sample the observed library format for auto-detection (auto mode only).
    if let Some(det) = sh.detector {
        if det.is_active() {
            if let Some(m) = maps
                .iter()
                .filter(|m| m.format.is_some())
                .max_by(|a, b| a.weight.total_cmp(&b.weight))
            {
                det.add_sample(m.format.unwrap());
            }
        }
    }

    // Strand-compatibility filtering against the expected format: the explicit
    // `-l` type, or — in auto (`-l A`) mode — the format the detector locks in
    // once it has seen enough of the read prefix (`resolved_format`). Until the
    // prefix is consumed this is `None` and nothing is filtered; afterwards the
    // detected type is applied for the rest of the sample (salmon's behavior).
    let expected = sh
        .expected_format
        .or_else(|| sh.detector.and_then(|d| d.resolved_format()));
    let compat: Vec<(&ScoredMapping, f64)> = maps
        .iter()
        .filter_map(|m| match expected {
            Some(exp) => {
                if is_compatible(exp, m.format, m.is_fw, m.status) {
                    Some((m, m.weight))
                } else if sh.ignore_incompat {
                    None
                } else {
                    Some((m, m.weight * sh.incompat_prior))
                }
            }
            None => Some((m, m.weight)),
        })
        .collect();
    if compat.is_empty() {
        return; // no compatible mapping -> fragment is unassigned
    }

    // Capture ONE immutable FLD snapshot pair (PMF + CMF) for this fragment — a
    // pair of cheap `Arc` clones, no per-fragment allocation — and share it
    // between the online posterior and the persistent eq-class weight. Indexing
    // the same snapshot for every mapping guarantees that transcripts sharing a
    // fragment length (in particular exact-duplicate transcripts) receive an
    // identical fragment-length probability even if another thread refreshes the
    // shared snapshot meanwhile. Reading the live FLD per mapping (racing atomic
    // loads on a concurrently-updated distribution) returned slightly different
    // values for the same length and let the VBEM α<1 prior break duplicate
    // symmetry. Mirrors C++ salmon's per-worker `LogCMFCache`. Empty until the
    // first online refresh (the pre-burn-in window, where the FLD term is not
    // folded in anyway).
    let fld_pmf = sh.fld.online_snapshot();
    let fld_cmf = sh.fld.online_cmf_snapshot();

    // Per-fragment bias-collection weights. With online (dual-phase) inference
    // these are abundance-aware posteriors `softmax(mass_t + log w_t)`, which
    // also advances the online masses by `logForgettingMass + log(posterior)`;
    // without it they fall back to the normalized aux (score) weights.
    let collecting = seqbias.is_some() || gcbias.is_some() || posbias.is_some();
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
    let online_post: Option<Vec<f64>> = sh.online.filter(|o| o.collecting()).map(|online| {
        let use_aux = online.num_assigned() >= sh.pre_burnin;
        let mm: Vec<(u32, f64)> = compat
            .iter()
            .map(|(m, w)| {
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
                    &fld_pmf,
                    &fld_cmf,
                );
                let log_cov = if *w > 0.0 { w.ln() } else { f64::NEG_INFINITY };
                (m.tid, log_cov + start_pos_prob + log_frag_prob)
            })
            .collect();
        online.assign_fragment(&mm, log_fm)
    });
    let bias_w: Vec<f64> = if collecting {
        if let Some(post) = &online_post {
            post.clone()
        } else {
            let wsum: f64 = compat.iter().map(|(_, w)| *w).sum();
            compat
                .iter()
                .map(|(_, w)| if wsum > 0.0 { *w / wsum } else { 0.0 })
                .collect()
        }
    } else {
        Vec::new()
    };
    // After burn-in salmon freezes model collection (still advances masses).
    let collect_now = collecting && sh.online.is_none_or(|o| o.collecting());

    if collect_now {
        if let Some(obs) = seqbias {
            for (i, (m, _)) in compat.iter().enumerate() {
                collect_context(sh.salmon, m, bias_w[i], obs);
            }
        }
        if let Some(gc) = gcbias {
            collect_gc(sh, &compat, &bias_w, gc);
        }
        if let Some(pos) = posbias {
            collect_pos(sh, &compat, &bias_w, pos);
        }
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
    // Burn-in gate. The inline (default) path reads the shared `num_processed`
    // counter as it did before, matching the non-deterministic behaviour. The
    // deterministic pass instead gates on the online inference's own progressive
    // count (the fragments replayed so far in key-sorted order), so the burn-in
    // window is honored in that fixed order rather than against the pre-counted
    // total.
    let use_aux = if sh.deterministic {
        sh.online.map_or(0, |o| o.num_assigned())
    } else {
        sh.num_processed.load(Ordering::Relaxed)
    } >= sh.pre_burnin;
    // Per-mapping FLD log-probability via the shared `frag_log_prob` term — the
    // same one the online posterior uses, reading the same per-fragment snapshot
    // pair captured above (so the two paths share one consistent FLD view, as
    // C++ does via a single `aln.logProb`). Proper pairs get the length-
    // conditioned PMF; orphans / single-end reads get the bounded-CMF ambiguous
    // weight. `LOG_1` (= 0) before the auxiliary model is trained. The snapshot
    // is primed from the prior before the run starts and refreshed at every
    // mini-batch boundary, so it is never empty here.
    let log_fps: Vec<f64> = compat
        .iter()
        .map(|(m, _)| {
            let ref_len = sh.salmon.ref_len(m.tid as usize) as i32;
            frag_log_prob(
                m,
                ref_len,
                use_aux,
                sh.model_single_frag_prob,
                sh.paired_lib,
                sh.no_frag_length_dist,
                &fld_pmf,
                &fld_cmf,
            )
        })
        .collect();
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
        .zip(&log_fps)
        .fold(LOG_0, |acc, ((_, w), &lfp)| log_add(acc, w.ln() + lfp));
    let mut pairs: Vec<(u32, f64)> = compat
        .iter()
        .zip(&log_fps)
        .map(|((m, w), &lfp)| (m.tid, (w.ln() + lfp - log_denom).exp()))
        .collect();

    // Abundance-aware FLD training: accept each concordant compatible pair's
    // fragment length with probability = its abundance-aware online posterior
    // (salmon's `if (r < exp(aln.logProb)) fragLengthDist.addVal(...)`, where
    // aln.logProb includes transcriptLogCount). For reads shared between
    // near-duplicates this preferentially samples the dominant transcript's
    // implied length, concentrating the FLD as salmon does (vs adding every best
    // pair at full weight, which overdisperses it). Frozen after the training
    // window (`online_post` is `None`).
    if let Some(post) = &online_post {
        // Seed the acceptance sampler from the fragment's own identity, so the
        // set of fragments that train the FLD is independent of which worker
        // thread processed this fragment (the documented source of multi-threaded
        // run-to-run wobble). Computed once per fragment, only inside the
        // training window.
        let sampler = FldSampler::from_key(key);
        for (i, (m, _)) in compat.iter().enumerate() {
            let conc = m.status == MateStatus::PairedEndPaired
                && m.fragment_len > 0
                && expected.is_none_or(|exp| is_compatible(exp, m.format, m.is_fw, m.status));
            if conc && sampler.u01(i) < post[i] {
                sh.fld.add_val(m.fragment_len as usize, 0.0);
            }
        }
    }

    // Build the equivalence class: sorted, de-duplicated transcript ids + weights.
    pairs.sort_by_key(|p| p.0);
    pairs.dedup_by_key(|p| p.0);
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
    // The fragment's key fixes the per-class weight-summation order (determinism).
    sh.eq.add_group_keyed(key, group, weights, 1);
}

/// Fold the most recent fragment's selective-alignment [`MapStats`] into the
/// shared meta counters. Call once per mapped/attempted fragment on the SA path.
fn accumulate_vm_stats(sh: &Shared, maps_empty: bool) {
    let s = salmon_map::take_last_map_stats();
    if maps_empty {
        if s.decoy_dominated {
            sh.num_decoy.fetch_add(1, Ordering::Relaxed);
        } else if s.had_candidates {
            sh.num_frags_filtered_vm.fetch_add(1, Ordering::Relaxed);
        }
    } else if s.alns_below_threshold > 0 {
        sh.num_below_threshold_vm
            .fetch_add(s.alns_below_threshold as u64, Ordering::Relaxed);
    }
    if s.dovetail {
        sh.num_dovetail.fetch_add(1, Ordering::Relaxed);
    }
}

/// Merge this thread's observed bias models (sequence and GC) into the shared
/// accumulators.
/// The read name written to unmapped_names.txt: the id up to the first
/// whitespace (the conventional fragment name), lossily decoded.
fn read_name(id: &[u8]) -> String {
    let end = id
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(id.len());
    String::from_utf8_lossy(&id[..end]).into_owned()
}

/// Merge this thread's collected unmapped-fragment names into the shared list.
fn merge_unmapped(proc: &mut QuantProcessor) {
    if let Some(shared) = proc.shared.unmapped_names {
        if !proc.unmapped.is_empty() {
            shared.lock().unwrap().append(&mut proc.unmapped);
        }
    }
}

/// Merge this thread's observed bias models (sequence, GC, positional) into the
/// shared accumulators. Inline (non-deterministic) path only; the deterministic
/// pass collects bias directly into shared models in [`run_inference_serial`].
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

/// Run the per-fragment inference over the mapped fragments in a deterministic,
/// key-sorted, single-threaded order. The parallel mapping pass only maps and
/// stores fragments; this pass replays the online posterior, FLD training, bias
/// collection, and eq-class accumulation in a fixed order, so every
/// order-dependent accumulation — the online masses, the FLD, and the observed
/// bias models — is byte-identical regardless of how many threads mapped the
/// reads.
///
/// `frags` must already yield fragments in ascending key order (the RAD store is
/// externally sorted by [`crate::rad::ExternalSort`], so the whole input never
/// has to reside in memory at once). Each item is a `Result` so a read error in
/// the streamed/merged store propagates out. `batch_size` sets the mini-batch
/// granularity for the forgetting-mass schedule and the FLD-snapshot refresh.
pub(crate) fn run_inference_serial<I>(
    sh: &Shared,
    frags: I,
    batch_size: usize,
) -> anyhow::Result<()>
where
    I: IntoIterator<Item = anyhow::Result<(u128, Vec<ScoredMapping>)>>,
{
    let mut seqbias = sh.collect_seqbias.then(|| (SBModel::new(), SBModel::new()));
    let mut gcbias = sh
        .collect_gcbias
        .then(|| GcFragModel::new(sh.cond_gc_bins, sh.gc_bins));
    let mut posbias = sh.collect_posbias.then(|| {
        (
            (0..NUM_LENGTH_CLASSES)
                .map(|_| SimplePosBias::default())
                .collect::<Vec<_>>(),
            (0..NUM_LENGTH_CLASSES)
                .map(|_| SimplePosBias::default())
                .collect::<Vec<_>>(),
        )
    });

    let bs = batch_size.max(1);
    let mut iter = frags.into_iter();
    // Pull the already-sorted stream one mini-batch at a time: holding only one
    // batch keeps pass-2 memory bounded even for very large stores.
    loop {
        let mut batch: Vec<(u128, Vec<ScoredMapping>)> = Vec::with_capacity(bs);
        for _ in 0..bs {
            match iter.next() {
                Some(item) => batch.push(item?),
                None => break,
            }
        }
        if batch.is_empty() {
            break;
        }
        let last = batch.len() < bs;
        // one forgetting-mass timestep per minibatch (online inference)
        let log_fm = sh.online.map_or(0.0, |o| o.next_log_fm());
        for (key, maps) in &batch {
            record(
                sh,
                *key,
                maps,
                log_fm,
                seqbias.as_mut(),
                gcbias.as_mut(),
                posbias.as_mut(),
            );
        }
        // refresh the FLD snapshot at the minibatch boundary (mirrors salmon)
        sh.fld.refresh_online();
        if last {
            break;
        }
    }

    // Merge the collected bias models into the shared accumulators that the
    // effective-length step reads.
    if let (Some(local), Some(shared_mtx)) = (seqbias.as_ref(), sh.seqbias_obs) {
        let mut g = shared_mtx.lock().unwrap();
        g.0.combine_counts(&local.0);
        g.1.combine_counts(&local.1);
    }
    if let (Some(local), Some(shared_mtx)) = (gcbias.as_ref(), sh.gcbias_obs) {
        shared_mtx.lock().unwrap().combine_counts(local);
    }
    if let (Some(local), Some(shared_mtx)) = (posbias.as_ref(), sh.posbias_obs) {
        let mut g = shared_mtx.lock().unwrap();
        for (a, b) in g.0.iter_mut().zip(&local.0) {
            a.combine(b);
        }
        for (a, b) in g.1.iter_mut().zip(&local.1) {
            a.combine(b);
        }
    }
    Ok(())
}

impl<'a, 'r> PairedParallelProcessor<RefRecord<'r>> for QuantProcessor<'a> {
    fn process_record_pair_batch(
        &mut self,
        pairs: impl Iterator<Item = (RefRecord<'r>, RefRecord<'r>)>,
    ) -> paraseq::Result<()> {
        let QuantProcessor {
            shared,
            hs,
            sketch_scratch,
            seqbias,
            gcbias,
            posbias,
            unmapped,
            sam_buf,
            frag_buf,
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
        // One forgetting-mass timestep per mapped batch for the inline online
        // posterior; the deterministic pass advances it per mini-batch in pass 2
        // instead, so leave the schedule untouched here when buffering.
        let log_fm = if sh.deterministic {
            0.0
        } else {
            sh.online.map_or(0.0, |o| o.next_log_fm())
        };
        for (r1, r2) in pairs {
            let s1 = r1.seq();
            let s2 = r2.seq();
            let mut maps = if sh.sketch {
                map_read_pair_sketch(
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
                )
            } else {
                map_read_pair(idx, hs, sh.salmon, s1.as_ref(), s2.as_ref(), sh.map_cfg)
            };
            // Sketch mappings carry no per-hit decoy flag and bypass the
            // selective-alignment finalize, so decoys would otherwise leak into
            // the eq-classes. Apply the same decoy policy here (drop decoy tids,
            // decoy-domination, --allowDecoyOrphans) and account decoy-dominated
            // fragments. SA mode already handled decoys inside finalize.
            if sh.sketch && sh.salmon.info().num_decoys > 0 {
                let decoy_dominated =
                    salmon_map::filter_sketch_decoys(&mut maps, sh.salmon, &sh.map_cfg.score);
                if decoy_dominated {
                    sh.num_decoy.fetch_add(1, Ordering::Relaxed);
                }
            }
            // A fragment mapping to too many places is discarded (salmon's
            // tooManyHits / maxReadOccs): treat as unmapped everywhere below.
            if maps.len() > sh.max_read_occ {
                maps.clear();
            }
            if maps.is_empty() && sh.unmapped_names.is_some() {
                unmapped.push(format!("{} u", read_name(r1.id())));
            }
            if !sh.sketch {
                accumulate_vm_stats(&sh, maps.is_empty());
            }
            if sh.sam.is_some() && !maps.is_empty() {
                crate::sam::write_fragment(
                    sam_buf,
                    sh.salmon,
                    r1.id(),
                    s1.as_ref(),
                    Some((r2.id(), s2.as_ref())),
                    &maps,
                );
            }
            sh.num_processed.fetch_add(1, Ordering::Relaxed);
            if maps.is_empty() {
                continue;
            }
            let key = xxhash_rust::xxh3::xxh3_128(r1.id());
            if sh.deterministic {
                // Pass 1 only maps + buffers; the inference runs deterministically
                // in pass 2 (see `run_inference_serial`).
                frag_buf.push((key, maps));
            } else {
                // Inline (default): run the per-fragment inference now, in parallel.
                record(
                    &sh,
                    key,
                    &maps,
                    log_fm,
                    seqbias.as_mut(),
                    gcbias.as_mut(),
                    posbias.as_mut(),
                );
            }
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        if self.shared.deterministic {
            // Stream this batch's mapped fragments to the RAD store.
            flush_rad(self);
        } else {
            // Inline path refreshes the FLD online snapshot at the mini-batch
            // boundary (see `record`); the deterministic pass trains it in pass 2.
            self.shared.fld.refresh_online();
        }
        flush_sam(self);
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        if self.shared.deterministic {
            // Flush any fragments not yet written by a batch boundary.
            flush_rad(self);
        } else {
            merge_bias(self);
        }
        merge_unmapped(self);
        flush_sam(self);
        Ok(())
    }
}

impl<'a, 'r> ParallelProcessor<RefRecord<'r>> for QuantProcessor<'a> {
    fn process_record_batch(
        &mut self,
        records: impl Iterator<Item = RefRecord<'r>>,
    ) -> paraseq::Result<()> {
        let QuantProcessor {
            shared,
            hs,
            sketch_scratch,
            seqbias,
            gcbias,
            posbias,
            unmapped,
            sam_buf,
            frag_buf,
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
        let log_fm = if sh.deterministic {
            0.0
        } else {
            sh.online.map_or(0.0, |o| o.next_log_fm())
        };
        for rec in records {
            let s = rec.seq();
            let mut maps = if sh.sketch {
                map_single_read_sketch(
                    idx,
                    hs,
                    sketch_scratch.as_mut().unwrap(),
                    s.as_ref(),
                    sh.skip,
                    sh.map_cfg.collect.max_hit_occ,
                    sh.max_read_occ,
                )
            } else {
                map_single_read(idx, hs, sh.salmon, s.as_ref(), sh.map_cfg)
            };
            // Sketch decoy policy (see the paired-end branch): drop decoy tids /
            // decoy-dominated fragments that SA mode handles inside finalize.
            if sh.sketch && sh.salmon.info().num_decoys > 0 {
                let decoy_dominated =
                    salmon_map::filter_sketch_decoys(&mut maps, sh.salmon, &sh.map_cfg.score);
                if decoy_dominated {
                    sh.num_decoy.fetch_add(1, Ordering::Relaxed);
                }
            }
            if maps.len() > sh.max_read_occ {
                maps.clear();
            }
            if maps.is_empty() && sh.unmapped_names.is_some() {
                unmapped.push(format!("{} u", read_name(rec.id())));
            }
            if !sh.sketch {
                accumulate_vm_stats(&sh, maps.is_empty());
            }
            if sh.sam.is_some() && !maps.is_empty() {
                crate::sam::write_fragment(sam_buf, sh.salmon, rec.id(), s.as_ref(), None, &maps);
            }
            sh.num_processed.fetch_add(1, Ordering::Relaxed);
            if maps.is_empty() {
                continue;
            }
            let key = xxhash_rust::xxh3::xxh3_128(rec.id());
            if sh.deterministic {
                frag_buf.push((key, maps));
            } else {
                record(
                    &sh,
                    key,
                    &maps,
                    log_fm,
                    seqbias.as_mut(),
                    gcbias.as_mut(),
                    posbias.as_mut(),
                );
            }
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        if self.shared.deterministic {
            flush_rad(self);
        } else {
            self.shared.fld.refresh_online();
        }
        flush_sam(self);
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        if self.shared.deterministic {
            flush_rad(self);
        } else {
            merge_bias(self);
        }
        merge_unmapped(self);
        flush_sam(self);
        Ok(())
    }
}

/// Write this thread's buffered mapped fragments to the RAD store as one chunk
/// and clear the buffer (deterministic path only). Called at batch boundaries to
/// stream the store to disk rather than holding every fragment in memory.
fn flush_rad(proc: &mut QuantProcessor) {
    if let Some(rad) = proc.shared.rad {
        if !proc.frag_buf.is_empty() {
            rad.write_chunk(&proc.frag_buf);
            proc.frag_buf.clear();
        }
    }
}

/// Flush this thread's accumulated SAM buffer to the shared writer (under lock).
fn flush_sam(proc: &mut QuantProcessor) {
    if let Some(sw) = proc.shared.sam {
        if !proc.sam_buf.is_empty() {
            let _ = sw.write_block(&proc.sam_buf);
            proc.sam_buf.clear();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reproduce the FLD-acceptance rule exactly: accept mapping `i` iff
    /// `u01(i) < post[i]`, evaluating the indices in the given `order`.
    fn accepts(key: u128, posts: &[f64], order: &[usize]) -> Vec<bool> {
        let s = FldSampler::from_key(key);
        let mut out = vec![false; posts.len()];
        for &i in order {
            out[i] = s.u01(i) < posts[i];
        }
        out
    }

    #[test]
    fn fld_acceptance_is_order_and_thread_independent() {
        let key: u128 = 0x1234_5678_9ABC_DEF0_0FED_CBA9_8765_4321;
        // Mix certain-accept (1.0), certain-reject (0.0), and middling probs so
        // the decision actually depends on the draws.
        let posts = [1.0, 0.0, 0.5, 0.5, 0.0, 1.0, 0.3, 0.7, 0.9, 0.1];

        let forward: Vec<usize> = (0..posts.len()).collect();
        let mut reversed = forward.clone();
        reversed.reverse();
        let shuffled = [3usize, 7, 0, 9, 1, 5, 8, 2, 6, 4];

        let a = accepts(key, &posts, &forward);
        let b = accepts(key, &posts, &reversed);
        let c = accepts(key, &posts, &shuffled);

        // The accept decision for index i is a pure function of (key, i, post[i]);
        // permuting evaluation order (i.e. thread scheduling) must not change any
        // decision.
        assert_eq!(a, b, "acceptance must not depend on evaluation order");
        assert_eq!(a, c, "acceptance must not depend on evaluation order");

        // post == 1.0 always accepts, post == 0.0 never accepts (u01 in [0,1)).
        assert!(a[0] && a[5]);
        assert!(!a[1] && !a[4]);
    }

    #[test]
    fn fld_sampler_is_deterministic_and_seed_robust() {
        // Same key -> identical draws (run-to-run determinism).
        let s1 = FldSampler::from_key(0xAAAA_BBBB_CCCC_DDDD);
        let s2 = FldSampler::from_key(0xAAAA_BBBB_CCCC_DDDD);
        for i in 0..16 {
            assert_eq!(s1.u01(i), s2.u01(i));
        }
        // Different keys -> different seeds (no collapse onto a global constant).
        assert_ne!(FldSampler::from_key(1).seed, FldSampler::from_key(2).seed);
        // A zero key still yields a usable, nonzero-seeded, in-range stream.
        let se = FldSampler::from_key(0);
        assert_ne!(se.seed, 0);
        for i in 0..16 {
            let u = se.u01(i);
            assert!((0.0..1.0).contains(&u));
        }
    }

    #[test]
    fn fld_u01_in_unit_interval() {
        let s = FldSampler::from_key(0xDEAD_BEEF_F00D);
        for i in 0..10_000 {
            let u = s.u01(i);
            assert!((0.0..1.0).contains(&u), "u01({i}) = {u} out of [0,1)");
        }
    }
}
