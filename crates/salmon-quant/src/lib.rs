//! `salmon-quant`: the reads-mode quantification driver.
//!
//! Ties the pipeline together for FASTQ input: parse reads in parallel
//! (paraseq) → map each fragment ([`salmon_map`]) → accumulate weighted
//! equivalence classes ([`salmon_eqclass`]) and fragment lengths
//! ([`salmon_model`]) → compute effective lengths → estimate abundances by EM
//! ([`salmon_infer`]) → write `quant.sf` and the JSON metadata
//! ([`output`]).
//!
//! Selective alignment is the default; set [`QuantOptions::sketch`] for the
//! alignment-free pseudoalignment path.

mod output;
mod processor;
mod sam;

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;

use piscem_rs::io::fastx::{reader_with_batch_size, Collection, CollectionType};
use piscem_rs::mapping::hit_searcher::SkippingStrategy;

use salmon_core::{LibraryFormat, ReadType};
use salmon_eqclass::EquivalenceClassBuilder;
use salmon_index::SalmonIndex;
use salmon_infer::{optimize, optimize_with_init, EmOptions};
use salmon_map::MapConfig;
use salmon_model::dumps::BiasDump;
use salmon_model::{FragmentLengthDistribution, LibraryTypeDetector, SBModel};

use processor::{QuantProcessor, Shared};

pub use output::write_outputs;

/// EM burn-in length before the in-loop bias correction, matching salmon's
/// `targetIt` (the abundance estimate after this many iterations is what weights
/// the expected bias models; see the bias-correction block in [`run`]).
const BIAS_PRELIM_ITERS: u32 = 11;

pub use salmon_core::ProgressCounters;

/// Options for a reads-mode quantification run.
#[derive(Debug, Clone)]
pub struct QuantOptions {
    /// directory of a salmon index (produced by `salmon-index`)
    pub index_dir: PathBuf,
    /// paired-end mate-1 files (parallel to `mates2`)
    pub mates1: Vec<PathBuf>,
    /// paired-end mate-2 files
    pub mates2: Vec<PathBuf>,
    /// single-end read files
    pub unmated: Vec<PathBuf>,
    /// output directory
    pub output_dir: PathBuf,
    /// library type string (e.g. `"IU"`, `"ISR"`, `"A"` for auto)
    pub lib_type: String,
    /// worker threads (`0` = all cores)
    pub num_threads: usize,
    /// use the alignment-free pseudoalignment path
    pub sketch: bool,
    /// sketch mode: restrict orphan emission to pairs whose other mate had no
    /// matching k-mers at all (`--sketchStrictOrphans`); default false uses the
    /// relaxed empty-accepted-target rule
    pub sketch_strict_orphan: bool,
    /// selective-alignment mapping configuration
    pub map_config: MapConfig,
    /// EM/VBEM inference options
    pub em: EmOptions,
    /// range-factorization bins for equivalence classes (0 = disabled, salmon default 4)
    pub range_factorization_bins: u32,
    /// prior weight for strand-incompatible mappings; `0` drops them (salmon default)
    pub incompat_prior: f64,
    /// write `aux_info/eq_classes.txt.gz` (salmon format: transcript list + index
    /// classes). `dump_eq_weights` additionally writes per-transcript weights.
    pub dump_eq: bool,
    pub dump_eq_weights: bool,
    /// write `aux_info/unmapped_names.txt` (names of unmapped fragments + status)
    pub write_unmapped_names: bool,
    /// write per-mapping SAM records to this path (`--writeMappings`); spoofed
    /// CIGAR (`<readLen>M` + end soft-clips), matching salmon's standard output
    pub write_mappings: Option<PathBuf>,
    /// enable sequence-specific bias correction (`--seqBias`)
    pub seq_bias: bool,
    /// enable fragment-GC bias correction (`--gcBias`)
    pub gc_bias: bool,
    /// enable positional bias correction (`--posBias`)
    pub pos_bias: bool,
    /// number of bootstrap replicates (`--numBootstraps`); 0 = off
    pub num_bootstraps: u32,
    /// number of Gibbs posterior samples (`--numGibbsSamples`); 0 = off
    pub num_gibbs_samples: u32,
    /// Gibbs thinning factor (`--thinningFactor`, salmon default 16)
    pub thinning_factor: u32,
    /// disable effective-length correction; use the raw reference length
    /// (`--noLengthCorrection`)
    pub no_length_correction: bool,
    /// fragment-length distribution prior mean, SD, and max tracked length
    /// (`--fldMean` / `--fldSD` / `--fldMax`)
    pub fld_mean: f64,
    pub fld_sd: f64,
    pub fld_max: usize,
    /// online-phase forgetting factor (`--forgettingFactor`, salmon default 0.65)
    pub forgetting_factor: f64,
    /// significant digits for the EffectiveLength and NumReads columns of
    /// `quant.sf` (`--sigDigits`, salmon default 3)
    pub sig_digits: u32,
    /// discard a fragment that maps to more than this many places (salmon's
    /// `--maxReadOcc`, default 200)
    pub max_read_occ: usize,
    /// fragment-length sampling stride for the GC bias convolution (salmon's
    /// `--biasSpeedSamp`, default 5)
    pub bias_speed_samp: usize,
    /// online-phase model-training window: the first N fragments train the
    /// auxiliary models (FLD/bias/error), after which they are fixed (salmon's
    /// `--numAuxModelSamples`, default 5,000,000)
    pub num_aux_model_samples: u64,
    /// disable the lower barrier on how short bias correction may make an
    /// effective length (salmon's `--noBiasLengthThreshold`)
    pub no_bias_length_threshold: bool,
    /// number of fragment-GC bins in the GC bias model (salmon's `--numGCBins`,
    /// default 25)
    pub gc_bins: usize,
    /// number of conditioning (context) bins in the GC bias model (salmon's
    /// `--conditionalGCBins`, default 3)
    pub cond_gc_bins: usize,
    /// skip the abundance estimation (EM + Gibbs/bootstrap) and `quant.sf`,
    /// emitting only equivalence classes, library type, and metadata
    /// (salmon's `--skipQuant`)
    pub skip_quant: bool,
    /// fragments processed before the auxiliary (FLD) model is folded into the
    /// online posterior (salmon's `--numPreAuxModelSamples`; salmon's default is
    /// 1,000,000, this port's prior hardcoded value is 5,000)
    pub num_pre_aux_model_samples: u64,
    /// Optional shared progress counters. When `Some`, [`quantify`] reports
    /// processed/mapped fragment counts here as it runs so the caller can drive
    /// a live progress display. `None` (the default) disables sharing.
    pub progress: Option<std::sync::Arc<ProgressCounters>>,
}

impl QuantOptions {
    /// Options with salmon-like defaults for the given index and output dir.
    pub fn new(index_dir: PathBuf, output_dir: PathBuf) -> Self {
        Self {
            index_dir,
            mates1: Vec::new(),
            mates2: Vec::new(),
            unmated: Vec::new(),
            output_dir,
            lib_type: "A".to_string(),
            num_threads: 0,
            sketch: false,
            sketch_strict_orphan: false,
            map_config: MapConfig::default(),
            em: EmOptions::default(),
            range_factorization_bins: 4,
            incompat_prior: 0.0,
            dump_eq: false,
            dump_eq_weights: false,
            write_unmapped_names: false,
            write_mappings: None,
            seq_bias: false,
            gc_bias: false,
            pos_bias: false,
            num_bootstraps: 0,
            num_gibbs_samples: 0,
            thinning_factor: 16,
            no_length_correction: false,
            fld_mean: 250.0,
            fld_sd: 25.0,
            fld_max: 1000,
            forgetting_factor: 0.65,
            sig_digits: 3,
            max_read_occ: 200,
            bias_speed_samp: 5,
            num_aux_model_samples: 5_000_000,
            no_bias_length_threshold: false,
            gc_bins: salmon_model::gcbias::DEFAULT_GC_BINS,
            cond_gc_bins: salmon_model::gcbias::DEFAULT_COND_BINS,
            skip_quant: false,
            num_pre_aux_model_samples: processor::NUM_PRE_BURNIN,
            progress: None,
        }
    }

    fn is_paired(&self) -> bool {
        !self.mates1.is_empty()
    }
}

/// Quantification results (also written to disk by [`write_outputs`]).
#[derive(Debug, Clone)]
pub struct QuantResult {
    pub names: Vec<String>,
    pub lengths: Vec<u32>,
    pub eff_lengths: Vec<f64>,
    pub tpm: Vec<f64>,
    /// estimated read counts per transcript
    pub counts: Vec<f64>,
    pub num_processed: u64,
    pub num_mapped: u64,
    /// mapped fragments placed as orphans (only one mate of a pair mapped)
    pub num_orphan: u64,
    pub num_eq_classes: usize,
    /// index of the first decoy reference (`None` if the index has no decoys);
    /// references at/after this are excluded from `quant.sf` and counted as decoys
    pub first_decoy_index: Option<usize>,
    /// whether the index collapsed duplicate sequences (for meta_info)
    pub keep_duplicates: bool,
    /// fragments dropped because their best alignment was to a decoy
    pub num_decoy_fragments: u64,
    /// fragments whose mates dovetail (overlap past each other)
    pub num_dovetail_fragments: u64,
    /// fragments that had candidate mappings but none survived validation/filtering
    pub num_fragments_filtered_vm: u64,
    /// per-fragment candidate alignments dropped for scoring below threshold,
    /// summed over fragments that nonetheless mapped
    pub num_alignments_below_threshold_for_mapped_fragments_vm: u64,
    pub frag_len_mean: f64,
    /// standard deviation of the observed fragment-length distribution
    pub frag_len_sd: f64,
    /// transcript length-class boundaries (salmon's `length_classes`)
    pub length_classes: Vec<u32>,
    /// normalized fragment-length PMF over lengths `0..=max` (for `flenDist.txt`)
    pub frag_len_dist: Vec<f64>,
    /// quant start timestamp (asctime-style), captured when the run began
    pub start_time: String,
    /// salmon-compatible index seq/name SHA-256/512 hashes (from the index's
    /// info.json), echoed into meta_info.json for downstream provenance tools
    pub index_seq_hash: String,
    pub index_name_hash: String,
    pub index_seq_hash512: String,
    pub index_name_hash512: String,
    pub index_decoy_seq_hash: String,
    pub index_decoy_name_hash: String,
    /// the library type used: the detected format when `-l A`, else the
    /// user-specified one
    pub library_type: String,
    /// posterior samples (bootstrap or Gibbs), one abundance vector each; empty
    /// when neither was requested
    pub bootstraps: Vec<Vec<f64>>,
    /// per-transcript (unique, ambiguous) fragment counts for `ambig_info.tsv`
    pub ambig: (Vec<u32>, Vec<u32>),
    /// observed/expected bias-model tables for the aux dumps; each component is
    /// empty unless the corresponding `--seqBias`/`--gcBias`/`--posBias` ran
    pub bias_dump: BiasDump,
}

/// Run quantification end-to-end, writing outputs and returning the results.
pub fn quantify(opts: &QuantOptions) -> Result<QuantResult> {
    let start_time = asctime_now();
    let salmon = SalmonIndex::load(&opts.index_dir)
        .with_context(|| format!("loading index {}", opts.index_dir.display()))?;
    let num_refs = salmon.num_refs();

    let eq_builder = EquivalenceClassBuilder::new();
    // Gaussian-prior fragment length distribution (mean 250, sd 25); empirical
    // observations from concordant pairs refine it.
    let mut fld =
        FragmentLengthDistribution::new(1.0, opts.fld_max, opts.fld_mean, opts.fld_sd, 4, 0.5, 1);
    // `processed`/`mapped` live in a (possibly caller-shared) `ProgressCounters`
    // so a CLI progress bar can poll them live; the rest are local.
    let progress = opts.progress.clone().unwrap_or_default();
    let num_processed = &progress.processed;
    let num_mapped = &progress.mapped;
    let num_orphan = AtomicU64::new(0);
    let num_decoy = AtomicU64::new(0);
    let num_dovetail = AtomicU64::new(0);
    let num_frags_filtered_vm = AtomicU64::new(0);
    let num_below_threshold_vm = AtomicU64::new(0);
    // collector for `--writeUnmappedNames` (names of unmapped fragments)
    let unmapped_collector: Option<std::sync::Mutex<Vec<String>>> = opts
        .write_unmapped_names
        .then(|| std::sync::Mutex::new(Vec::new()));
    // SAM sink for `--writeMappings` (header written here).
    let sam_writer: Option<sam::SamWriter> = match &opts.write_mappings {
        Some(path) => {
            let cmd = format!("salmon quant -i {}", opts.index_dir.display());
            Some(sam::SamWriter::create(path, &salmon, &cmd).context("opening SAM output")?)
        }
        None => None,
    };
    let nthreads = if opts.num_threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        opts.num_threads
    };

    // Auto-detect the library type when requested (`-l A`).
    let auto_detect = LibraryFormat::is_auto(&opts.lib_type);
    let read_type = if opts.is_paired() {
        ReadType::PairedEnd
    } else {
        ReadType::SingleEnd
    };
    let detector = auto_detect.then(|| LibraryTypeDetector::new(read_type));

    // Strand-compatibility filtering applies only with an explicit library
    // type; in auto mode the type is unknown during the pass.
    let expected_format = if auto_detect {
        None
    } else {
        Some(
            LibraryFormat::parse(&opts.lib_type)
                .with_context(|| format!("invalid library type '{}'", opts.lib_type))?,
        )
    };
    let ignore_incompat = opts.incompat_prior <= 0.0;

    // Shared observed sequence-bias models (merged per-thread) for `--seqBias`.
    let seqbias_obs = opts
        .seq_bias
        .then(|| std::sync::Mutex::new((SBModel::new(), SBModel::new())));

    // For `--gcBias`: per-transcript cumulative G+C counts via a single rank
    // bitvector over the concatenated references (~1 bit/base; salmon's
    // `--reduceGCMemory`). Benchmarked faster and ~2x leaner than the old dense
    // `Vec<Vec<u32>>` (4 bytes/base) with effectively identical results, so it
    // is now the default; `--reduceGCMemory` is accepted as a no-op. `gc_store`
    // presents per-transcript [`GcView`]s.
    let gc_rank: Option<salmon_model::GcRank> = opts
        .gc_bias
        .then(|| salmon_model::GcRank::new(salmon.refseq_concat()));
    let gc_store: Option<salmon_model::GcStore> =
        gc_rank.as_ref().map(|r| salmon_model::GcStore::Rank {
            rank: r,
            offsets: salmon.ref_offsets(),
        });
    let gcbias_obs = opts.gc_bias.then(|| {
        std::sync::Mutex::new(salmon_model::GcFragModel::new(
            opts.cond_gc_bins,
            opts.gc_bins,
        ))
    });

    // For `--posBias`: transcript length-class quantiles + per-transcript class,
    // and the shared observed (5', 3') positional-bias models per length class.
    let length_quantiles: Option<Vec<u32>> = opts.pos_bias.then(|| {
        let lens: Vec<u32> = (0..num_refs).map(|t| salmon.ref_len(t) as u32).collect();
        salmon_model::compute_length_quantiles(&lens, salmon_model::NUM_LENGTH_CLASSES)
    });
    let length_class: Option<Vec<usize>> = length_quantiles.as_ref().map(|q| {
        (0..num_refs)
            .map(|t| salmon_model::length_class_index(q, salmon.ref_len(t) as u32))
            .collect()
    });
    let posbias_obs = opts.pos_bias.then(|| {
        let mk = || {
            (0..salmon_model::NUM_LENGTH_CLASSES)
                .map(|_| salmon_model::SimplePosBias::default())
                .collect::<Vec<_>>()
        };
        std::sync::Mutex::new((mk(), mk()))
    });

    // Online (dual-phase) inference: when any bias model is collected, develop
    // running per-transcript abundances during the mapping pass so the observed
    // bias models are weighted by abundance-aware posteriors (salmon's online
    // phase), not score-only weights. The offline EM still gives the final
    // point estimate.
    // The online estimate runs unconditionally: besides weighting the observed
    // bias models (when bias correction is on), it provides the abundance-aware
    // posterior used to train the fragment-length distribution (salmon's
    // `r < exp(aln.logProb)` acceptance). The offline EM does not read it.
    let online = {
        let ref_lens: Vec<u64> = (0..num_refs).map(|t| salmon.ref_len(t)).collect();
        Some(salmon_infer::OnlineInference::new(
            &ref_lens,
            0.05,
            opts.forgetting_factor,
            opts.num_aux_model_samples,
        ))
    };

    // ---- parallel mapping pass (borrows the accumulators) -------------------
    {
        let shared = Shared {
            salmon: &salmon,
            eq: &eq_builder,
            fld: &fld,
            detector: detector.as_ref(),
            map_cfg: &opts.map_config,
            sketch: opts.sketch,
            sketch_strict_orphan: opts.sketch_strict_orphan,
            max_read_occ: opts.max_read_occ,
            pre_burnin: opts.num_pre_aux_model_samples,
            skip: SkippingStrategy::Permissive,
            range_factorization_bins: opts.range_factorization_bins,
            expected_format,
            incompat_prior: opts.incompat_prior,
            ignore_incompat,
            collect_seqbias: opts.seq_bias,
            seqbias_obs: seqbias_obs.as_ref(),
            collect_gcbias: opts.gc_bias,
            cond_gc_bins: opts.cond_gc_bins,
            gc_bins: opts.gc_bins,
            gc_store,
            gcbias_obs: gcbias_obs.as_ref(),
            collect_posbias: opts.pos_bias,
            length_class: length_class.as_deref(),
            posbias_obs: posbias_obs.as_ref(),
            online: online.as_ref(),
            paired_lib: opts.is_paired(),
            num_processed,
            num_mapped,
            num_orphan: &num_orphan,
            num_decoy: &num_decoy,
            num_dovetail: &num_dovetail,
            num_frags_filtered_vm: &num_frags_filtered_vm,
            num_below_threshold_vm: &num_below_threshold_vm,
            unmapped_names: unmapped_collector.as_ref(),
            sam: sam_writer.as_ref(),
        };
        let mut proc = QuantProcessor::new(shared);
        tracing::info!(
            "mapping reads ({} mode, {nthreads} threads)",
            if opts.sketch {
                "sketch"
            } else {
                "selective-alignment"
            }
        );
        if opts.is_paired() {
            run_paired(&opts.mates1, &opts.mates2, &mut proc, nthreads)?;
        } else {
            run_single(&opts.unmated, &mut proc, nthreads)?;
        }
    }
    {
        let p = num_processed.load(Ordering::Relaxed);
        let m = num_mapped.load(Ordering::Relaxed);
        let pct = if p > 0 {
            100.0 * m as f64 / p as f64
        } else {
            0.0
        };
        tracing::info!("mapped {m} / {p} fragments ({pct:.2}%)");
    }
    if let Some(sw) = &sam_writer {
        sw.flush().context("flushing SAM output")?;
    }

    // Write aux_info/unmapped_names.txt ("<name> <status>" per line; the port maps
    // unmapped fragments as "u" — orphan/decoy sub-codes await mapper-reason
    // reporting, tracked with the *_vm counters).
    if let Some(collector) = &unmapped_collector {
        use std::io::Write as _;
        let names = collector.lock().unwrap();
        std::fs::create_dir_all(opts.output_dir.join("aux_info")).context("creating aux_info")?;
        let mut w = std::io::BufWriter::new(
            std::fs::File::create(opts.output_dir.join("aux_info").join("unmapped_names.txt"))
                .context("creating unmapped_names.txt")?,
        );
        for line in names.iter() {
            writeln!(w, "{line}")?;
        }
        w.flush()?;
    }

    // Resolve the library type: the detected format (when auto), else the
    // user-specified string. Fall back to a sensible default if detection saw
    // no usable samples.
    let library_type = if let Some(det) = &detector {
        det.most_likely_type()
            .map(|f| f.canonical().to_string())
            .unwrap_or_else(|| if opts.is_paired() { "IU" } else { "U" }.to_string())
    } else {
        opts.lib_type.clone()
    };

    // ---- effective lengths from the (now frozen) FLD ------------------------
    // salmon's base effective length is `refLen − E[L | L ≤ refLen]`
    // (`computeSmoothedEffectiveLengths` via `correctionFactorsFromMass`), NOT the
    // truncated-PMF `Σ pmf(l)·(refLen−l+1)` estimate, which falls back to the raw
    // reference length for any transcript shorter than the FLD mean.
    fld.cache();
    let log_pmf = fld.log_pmf();
    let cond_means = fld.conditional_means();
    let mut eff_lengths = vec![0f64; num_refs];
    for tid in 0..num_refs {
        eff_lengths[tid] = if opts.no_length_correction {
            salmon.ref_len(tid) as f64
        } else {
            salmon_model::smoothed_effective_length(&cond_means, salmon.ref_len(tid) as usize)
        };
    }

    // ---- inference ----------------------------------------------------------
    let mut collapsed = eq_builder.finish();
    collapsed.update_eff_lengths(&eff_lengths);
    let num_eq_classes = collapsed.len();

    if opts.dump_eq || opts.dump_eq_weights {
        dump_eq_classes(&opts.output_dir, &salmon, &collapsed, opts.dump_eq_weights)
            .context("writing eq_classes.txt.gz")?;
    }
    let bias_on = (opts.seq_bias || opts.gc_bias || opts.pos_bias) && !opts.no_length_correction;
    // salmon runs a *single* offline EM: after a short burn-in (`targetIt = 10`
    // iterations) it corrects effective lengths in place using that early
    // abundance estimate to weight the expected bias models, then continues the
    // same alpha vector to convergence. We mirror that: when bias-correcting,
    // run EM only `BIAS_PRELIM_ITERS` iterations here to weight the expected
    // models, then warm-start the single full convergence after correction
    // (below). Without bias, this is the final EM and runs to convergence.
    // Avoids a wasteful second full convergence (~20s on the 36M PE set).
    if !opts.skip_quant {
        tracing::info!(
            "estimating abundances ({})",
            if opts.em.use_vbem { "VBEM" } else { "EM" }
        );
    }
    let mut em = if opts.skip_quant {
        // `--skipQuant`: emit equivalence classes + library type + metadata but
        // skip the optimizer (and Gibbs/bootstrap below, and quant.sf). Leave
        // abundances at zero.
        salmon_infer::EmResult {
            alphas: vec![0.0; num_refs],
            iters: 0,
            converged: true,
        }
    } else if bias_on {
        let mut pre = opts.em.clone();
        pre.min_iter = BIAS_PRELIM_ITERS;
        pre.max_iter = BIAS_PRELIM_ITERS;
        // No min-alpha truncation: these alphas continue into the warm-started
        // EM, and a zeroed alpha can never recover under the multiplicative
        // update (salmon keeps one continuous, untruncated vector).
        pre.min_alpha = 0.0;
        optimize(&collapsed, num_refs, &pre, Some(&eff_lengths))
    } else {
        optimize(&collapsed, num_refs, &opts.em, Some(&eff_lengths))
    };

    // ---- bias correction: re-estimate effective lengths --------------------
    // Build the expected (background) bias models from the initial abundance
    // estimate, recompute bias-corrected effective lengths, then re-run
    // inference. Composes `--seqBias`, `--gcBias`, and `--posBias` in any
    // combination via the unified convolution (salmon's updateEffectiveLengths).
    let mut bias_dump = BiasDump::default();
    if bias_on && !opts.skip_quant {
        use rayon::prelude::*;
        let pmf_lin: Vec<f64> = log_pmf.iter().map(|lp| lp.exp()).collect();
        let (fld_cdf, fld_low, fld_high) = salmon_model::seqbias::fld_cdf_and_bounds(&pmf_lin);
        // K excludes the leading sequence context only when seq-correcting.
        let k = if opts.seq_bias {
            salmon_model::seqbias::CONTEXT_LENGTH
        } else {
            1
        };

        // Sequence-bias observed + expected models.
        let seq = seqbias_obs.map(|m| {
            let (mut obs_fw, mut obs_rc) = m.into_inner().unwrap();
            obs_fw.normalize();
            obs_rc.normalize();
            let (exp_fw, exp_rc) = salmon_model::build_expected(
                num_refs,
                |t| salmon.ref_seq(t as u32),
                &em.alphas,
                &eff_lengths,
                &fld_cdf,
            );
            (obs_fw, obs_rc, exp_fw, exp_rc)
        });
        if let Some((ofw, orc, efw, erc)) = seq.as_ref() {
            bias_dump.obs5_seq = ofw.dump().to_vec();
            bias_dump.obs3_seq = orc.dump().to_vec();
            bias_dump.exp5_seq = efw.dump().to_vec();
            bias_dump.exp3_seq = erc.dump().to_vec();
        }

        // Fragment-GC observed + expected models -> clamped ratio model.
        let store = gc_store;
        let gc_ratio_model = if let Some(m) = gcbias_obs {
            let mut obs = m.into_inner().unwrap();
            let mut exp = salmon_model::build_expected_gc(
                num_refs,
                |t| salmon.ref_seq(t as u32),
                |t| store.unwrap().view(t),
                &em.alphas,
                &eff_lengths,
                &fld_cdf,
                fld_low,
                fld_high,
                opts.cond_gc_bins,
                opts.gc_bins,
                k,
                opts.bias_speed_samp,
            );
            // Capture the (normalized) observed/expected GC tables for the dumps
            // before gc_ratio consumes them (normalize is idempotent).
            obs.normalize();
            exp.normalize();
            bias_dump.obs_gc = obs.dump().to_vec();
            bias_dump.exp_gc = exp.dump().to_vec();
            Some(salmon_model::gc_ratio(
                &mut obs,
                &mut exp,
                salmon_model::gcbias::GC_MAX_RATIO,
            ))
        } else {
            None
        };

        // Positional-bias observed (finalized) + expected (5'/3' per length class).
        let pos_models = posbias_obs.map(|m| {
            let (mut obs_fw, mut obs_rc) = m.into_inner().unwrap();
            for x in obs_fw.iter_mut().chain(obs_rc.iter_mut()) {
                x.finalize();
            }
            let (exp_fw, exp_rc) = salmon_model::build_expected_pos(
                num_refs,
                |t| salmon.ref_len(t) as usize,
                &em.alphas,
                &eff_lengths,
                &fld_cdf,
                length_quantiles.as_ref().unwrap(),
                k,
            );
            // Debug: dump finalized pos models for parity comparison with salmon.
            if std::env::var("SALMON_DEBUG_POS").is_ok() {
                use std::io::Write;
                let mut f =
                    std::fs::File::create(opts.output_dir.join("rust_pos_models.txt")).unwrap();
                for (name, models) in [
                    ("obs5", &obs_fw),
                    ("obs3", &obs_rc),
                    ("exp5", &exp_fw),
                    ("exp3", &exp_rc),
                ] {
                    for (lc, m) in models.iter().enumerate() {
                        write!(f, "{name} {lc}").unwrap();
                        for v in m.masses() {
                            write!(f, " {v:.6}").unwrap();
                        }
                        writeln!(f).unwrap();
                    }
                }
            }
            (obs_fw, obs_rc, exp_fw, exp_rc)
        });
        if let Some((ofw, orc, efw, erc)) = pos_models.as_ref() {
            let masses =
                |v: &[salmon_model::SimplePosBias]| v.iter().map(|m| m.masses().to_vec()).collect();
            bias_dump.obs5_pos = masses(ofw);
            bias_dump.obs3_pos = masses(orc);
            bias_dump.exp5_pos = masses(efw);
            bias_dump.exp3_pos = masses(erc);
        }

        let corrected: Vec<f64> = (0..num_refs)
            .into_par_iter()
            .map(|tid| {
                if em.alphas[tid] < 1e-8 {
                    return eff_lengths[tid];
                }
                let s = salmon.ref_seq(tid as u32);
                let ref_len = s.len();
                // Per-position 5'/3' positional-bias factors (obs/exp projected).
                let pos_vecs: Option<(Vec<f64>, Vec<f64>)> =
                    pos_models.as_ref().map(|(ofw, orc, efw, erc)| {
                        let lc = length_class.as_ref().unwrap()[tid];
                        let (mut o5, mut e5) = (vec![0.0; ref_len], vec![0.0; ref_len]);
                        let (mut o3, mut e3) = (vec![0.0; ref_len], vec![0.0; ref_len]);
                        ofw[lc].project_weights(&mut o5);
                        efw[lc].project_weights(&mut e5);
                        orc[lc].project_weights(&mut o3);
                        erc[lc].project_weights(&mut e3);
                        // additively-smoothed obs/exp (no hard floor; neutral tails)
                        let pf = salmon_model::positional_factor(&o5, &e5);
                        let pr = salmon_model::positional_factor(&o3, &e3);
                        (pf, pr)
                    });
                let bias = salmon_model::BiasInputs {
                    seq: seq.as_ref().map(|(of, or, ef, er)| (of, ef, or, er)),
                    gc: gc_ratio_model
                        .as_ref()
                        .map(|g| (g, store.unwrap().view(tid))),
                    pos: pos_vecs
                        .as_ref()
                        .map(|(pf, pr)| (pf.as_slice(), pr.as_slice())),
                };
                salmon_model::corrected_effective_length_full(
                    s,
                    &fld_cdf,
                    fld_low,
                    fld_high,
                    &bias,
                    eff_lengths[tid],
                    opts.bias_speed_samp,
                    opts.no_bias_length_threshold,
                )
            })
            .collect();
        eff_lengths = corrected;
        collapsed.update_eff_lengths(&eff_lengths);
        // Warm-start the single full convergence from the burn-in alphas, as
        // salmon continues the same vector after its in-loop correction.
        let warm = std::mem::take(&mut em.alphas);
        em = optimize_with_init(
            &collapsed,
            num_refs,
            &opts.em,
            Some(&warm),
            Some(&eff_lengths),
        );
    }
    let counts = em.alphas;

    // ---- posterior uncertainty (bootstrap / Gibbs) + ambiguity --------------
    // The packed CSR layout (piscem-infer style) makes these parallel-friendly.
    let packed = salmon_infer::PackedEqClasses::from_collapsed(&collapsed, num_refs);
    let ambig = salmon_infer::ambiguity_counts(&packed);
    let num_mapped_frags = num_mapped.load(Ordering::Relaxed);
    let bootstraps: Vec<Vec<f64>> = if opts.skip_quant {
        Vec::new()
    } else if opts.num_bootstraps > 0 {
        salmon_infer::bootstrap(
            &packed,
            &opts.em,
            opts.num_bootstraps,
            num_mapped_frags,
            true, // useScaledCounts (selective-alignment, no orphan-only quasi)
            0x5A13_0000,
        )
    } else if opts.num_gibbs_samples > 0 {
        // Gibbs prior follows the main optimizer (salmon): with VBEM and a
        // per-transcript prior it is `max(1.0, vbPrior)`; with plain EM it is
        // 1e-3 per transcript. The rust VBEM uses a constant per-transcript prior.
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
        salmon_infer::gibbs_sample(
            &packed,
            &eff_lengths,
            &counts,
            &gopts,
            num_mapped_frags,
            0x6217_0000,
        )
    } else {
        Vec::new()
    };

    // ---- TPM ----------------------------------------------------------------
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

    let result = QuantResult {
        names: (0..num_refs)
            .map(|i| salmon.ref_name(i).to_string())
            .collect(),
        lengths: (0..num_refs).map(|i| salmon.ref_len(i) as u32).collect(),
        eff_lengths,
        tpm,
        counts,
        num_processed: num_processed.load(Ordering::Relaxed),
        num_mapped: num_mapped.load(Ordering::Relaxed),
        num_orphan: num_orphan.load(Ordering::Relaxed),
        num_eq_classes,
        first_decoy_index: salmon.info().first_decoy_index,
        keep_duplicates: false,
        num_decoy_fragments: num_decoy.load(Ordering::Relaxed),
        num_dovetail_fragments: num_dovetail.load(Ordering::Relaxed),
        num_fragments_filtered_vm: num_frags_filtered_vm.load(Ordering::Relaxed),
        num_alignments_below_threshold_for_mapped_fragments_vm: num_below_threshold_vm
            .load(Ordering::Relaxed),
        frag_len_mean: fld.mean(),
        frag_len_sd: fld.sd(),
        length_classes: salmon_model::compute_length_quantiles(
            &(0..num_refs)
                .map(|i| salmon.ref_len(i) as u32)
                .collect::<Vec<_>>(),
            salmon_model::NUM_LENGTH_CLASSES,
        ),
        frag_len_dist: fld.log_pmf().iter().map(|lp| lp.exp()).collect(),
        start_time,
        index_seq_hash: salmon.info().seq_hash.clone(),
        index_name_hash: salmon.info().name_hash.clone(),
        index_seq_hash512: salmon.info().seq_hash512.clone(),
        index_name_hash512: salmon.info().name_hash512.clone(),
        index_decoy_seq_hash: salmon.info().decoy_seq_hash.clone(),
        index_decoy_name_hash: salmon.info().decoy_name_hash.clone(),
        library_type,
        bootstraps,
        ambig,
        bias_dump,
    };

    tracing::info!("writing results to {}", opts.output_dir.display());
    write_outputs(opts, &result)?;
    Ok(result)
}

/// Current local time as an asctime-style string (`Wed Jun 10 20:34:58 2026`),
/// matching salmon's `start_time`/`end_time` format.
pub(crate) fn asctime_now() -> String {
    jiff::Zoned::now()
        .strftime("%a %b %e %H:%M:%S %Y")
        .to_string()
}

/// Write the naive equivalence classes (collapsing any range-factorized
/// sub-classes back to their transcript set) for comparison/diagnostics.
/// Write `aux_info/eq_classes.txt.gz` in salmon's format: the transcript count,
/// the equivalence-class count, all transcript names (one per line, the index
/// order the classes reference), then one line per class. Without `with_weights`
/// classes are collapsed by transcript set (`groupSize`, tids…, summed count),
/// matching salmon's `--dumpEq`; with `with_weights` each class is emitted as-is
/// with its per-transcript combined weights before the count (`--dumpEqWeights`).
fn dump_eq_classes(
    dir: &Path,
    salmon: &SalmonIndex,
    collapsed: &salmon_eqclass::CollapsedEqClasses,
    with_weights: bool,
) -> Result<()> {
    use std::collections::BTreeMap;
    use std::io::Write as _;
    let num_txps = salmon.num_refs();
    std::fs::create_dir_all(dir.join("aux_info"))?;
    let f = std::fs::File::create(dir.join("aux_info").join("eq_classes.txt.gz"))?;
    let mut w = flate2::write::GzEncoder::new(f, flate2::Compression::new(6));

    if with_weights {
        writeln!(w, "{num_txps}")?;
        writeln!(w, "{}", collapsed.classes.len())?;
        for t in 0..num_txps {
            writeln!(w, "{}", salmon.ref_name(t))?;
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
        for t in 0..num_txps {
            writeln!(w, "{}", salmon.ref_name(t))?;
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

/// Open a (optionally gzipped) read file as a boxed reader.
fn open_reader(path: &Path) -> Result<Box<dyn std::io::Read + Send>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Ok(Box::new(MultiGzDecoder::new(f)))
    } else {
        Ok(Box::new(f))
    }
}

fn run_single(paths: &[PathBuf], proc: &mut QuantProcessor, nthreads: usize) -> Result<()> {
    let mut readers = Vec::with_capacity(paths.len());
    for p in paths {
        readers.push(
            reader_with_batch_size(open_reader(p)?)
                .map_err(|e| anyhow::anyhow!("opening {}: {e}", p.display()))?,
        );
    }
    let collection = Collection::new(readers, CollectionType::Single)
        .map_err(|e| anyhow::anyhow!("building read collection: {e}"))?;
    collection
        .process_parallel(proc, nthreads, None)
        .map_err(|e| anyhow::anyhow!("mapping failed: {e}"))?;
    Ok(())
}

fn run_paired(
    mates1: &[PathBuf],
    mates2: &[PathBuf],
    proc: &mut QuantProcessor,
    nthreads: usize,
) -> Result<()> {
    anyhow::ensure!(
        mates1.len() == mates2.len(),
        "mates1 and mates2 must have the same number of files"
    );
    let mut readers = Vec::with_capacity(mates1.len() * 2);
    for (a, b) in mates1.iter().zip(mates2.iter()) {
        readers.push(
            reader_with_batch_size(open_reader(a)?)
                .map_err(|e| anyhow::anyhow!("opening {}: {e}", a.display()))?,
        );
        readers.push(
            reader_with_batch_size(open_reader(b)?)
                .map_err(|e| anyhow::anyhow!("opening {}: {e}", b.display()))?,
        );
    }
    let collection = Collection::new(readers, CollectionType::Paired)
        .map_err(|e| anyhow::anyhow!("building paired read collection: {e}"))?;
    collection
        .process_parallel_paired(proc, nthreads, None)
        .map_err(|e| anyhow::anyhow!("mapping failed: {e}"))?;
    Ok(())
}
