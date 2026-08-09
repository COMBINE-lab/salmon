//! `salmon-quant`: the reads-mode quantification driver.
//!
//! Ties the pipeline together for FASTQ input: parse reads in parallel
//! (paraseq) → map each fragment ([`salmon_map`]) → accumulate weighted
//! equivalence classes ([`salmon_eqclass`]) and fragment lengths
//! ([`salmon_model`]) → compute effective lengths → estimate abundances by EM
//! ([`salmon_infer`]) → write `quant.sf` and the JSON metadata
//! (the `output` module).
//!
//! Selective alignment is the default; set [`QuantOptions::sketch`] for the
//! alignment-free pseudoalignment path.
//!
//! # The shape of a run
//!
//! [`quantify`] is one long function because a quant run *is* one long sequence
//! of phases, each depending on the last:
//!
//! 1. **Load the index** (optionally without the reference sequence, which is
//!    only needed for alignment and sequence-dependent bias).
//! 2. **Map**, in parallel over the input reads. This is the expensive phase and
//!    the one that populates everything else: equivalence classes, the
//!    fragment-length distribution, the library-type tally, the observed bias
//!    models, and any SAM/BAM/RAD output.
//! 3. **Finalize the models** learned during mapping, and from them compute each
//!    transcript's effective length.
//! 4. **Infer** abundances by EM/VBEM over the equivalence classes, optionally
//!    with a bias-correction round and posterior sampling.
//! 5. **Write** `quant.sf` and the metadata.
//!
//! Bias correction adds a wrinkle to step 4: the bias models must be weighted by
//! abundance, but abundance is what we are estimating. salmon resolves this by
//! running a rough EM first, using its output to weight the bias models,
//! recomputing effective lengths, and then running the real EM.

pub mod decode;
mod bam;
mod mapping_record;
mod output;
mod processor;
mod sam;

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{Context, Result};

use piscem_rs::io::fastx::{Collection, CollectionType};
use piscem_rs::mapping::hit_searcher::SkippingStrategy;

use salmon_core::{LibraryFormat, PhaseTimer, ReadType};
use salmon_eqclass::EquivalenceClassBuilder;
use salmon_index::SalmonIndex;
pub use salmon_infer::EmAccel;
use salmon_infer::{optimize_packed_with_init, EmOptions};
use salmon_map::MapConfig;
use salmon_model::dumps::{dump_bias_models_to_file, BiasDump};
use salmon_model::{FragmentLengthDistribution, LibraryTypeDetector, SBModel};
pub use salmon_rad::ChunkCodec;

use processor::{QuantProcessor, Shared};

pub use output::write_outputs;

pub use salmon_core::ProgressCounters;

/// Options for a reads-mode quantification run.
///
/// The library-level equivalent of salmon's command line: every field here
/// corresponds to a flag, and `salmon-cli` does little more than parse arguments
/// into this struct. Defaults come from [`QuantOptions::new`] and match salmon's.
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
    /// Which gzip decoder to use. `Auto` lets the engagement policy decide
    /// whether the parallel decoder can pay for the mapping threads it costs.
    pub decoder: piscem_rs::io::calibrate::DecoderPreference,
    /// Optional JSON overriding thread/decoder policy. Same format piscem uses.
    pub thread_policy: Option<PathBuf>,
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
    /// write per-mapping BAM records to this path (`--writeBam`); mutually
    /// exclusive with [`Self::write_mappings`]
    pub write_bam: Option<PathBuf>,
    /// override the number of BGZF compression workers used by `--writeBam`
    /// (`--bamCompressThreads`). `None` derives it from `num_threads` so that
    /// compression just keeps up with record production; set it only to
    /// benchmark or to compensate for unusually slow/fast storage.
    pub bam_compress_threads: Option<std::num::NonZeroUsize>,
    /// write per-fragment mappings to this RAD file (`--writeRad`); sketch or
    /// selective-alignment profile is chosen from `sketch`. Quantification still
    /// runs; combine with `skip_quant` to map only.
    pub write_rad: Option<PathBuf>,
    /// `--deterministic` mapping pass: derive the fragment-length distribution and
    /// library format **order-independently** (integer count histograms over
    /// uniquely-mapped proper pairs; see [`salmon_model::DiscreteFld`]) instead of
    /// the online log-space FLD + prefix detector, then bake them. This makes the
    /// baked values — and hence the downstream requant — byte-identical across
    /// thread counts.
    pub deterministic_fld: bool,
    /// enable sequence-specific bias correction (`--seqBias`)
    pub seq_bias: bool,
    /// enable fragment-GC bias correction (`--gcBias`)
    pub gc_bias: bool,
    /// enable positional bias correction (`--posBias`)
    pub pos_bias: bool,
    /// EM iterations for the preliminary abundance estimate that weights bias
    /// collection (`--biasSeedEMIters`, hidden); default 11. Rough on purpose —
    /// fully-converged uncorrected estimates over-fit the bias models.
    pub bias_seed_em_iters: u32,
    /// dump observed+expected seq/GC/positional bias models to `bias_models.txt`
    /// in the output dir (`--dumpBiasModels`, hidden; debugging/parity).
    pub dump_bias_models: bool,
    /// chunk compression codec for `--writeRad` output (and the `--deterministic`
    /// intermediate RAD). [`ChunkCodec::None`] writes uncompressed (the default
    /// for library/test callers; the CLI defaults this to lz4).
    pub rad_codec: ChunkCodec,
    /// number of bootstrap replicates (`--numBootstraps`); 0 = off
    pub num_bootstraps: u32,
    /// number of Gibbs posterior samples (`--numGibbsSamples`); 0 = off
    pub num_gibbs_samples: u32,
    /// Gibbs thinning factor (`--thinningFactor`, salmon default 16)
    pub thinning_factor: u32,
    /// disable effective-length correction; use the raw reference length
    /// (`--noLengthCorrection`)
    pub no_length_correction: bool,
    /// model the fragment-length probability of orphan / single-end mappings via
    /// the bounded-CMF "ambiguous" weight (salmon default `true`); `false` =
    /// `--noSingleFragProb`.
    pub model_single_frag_prob: bool,
    /// disable the fragment-length distribution in the per-fragment assignment
    /// probability (`--noFragLengthDist`); default `false`.
    pub no_frag_length_dist: bool,
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
            decoder: piscem_rs::io::calibrate::DecoderPreference::Auto,
            thread_policy: None,
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
            write_bam: None,
            bam_compress_threads: None,
            write_rad: None,
            deterministic_fld: false,
            seq_bias: false,
            gc_bias: false,
            pos_bias: false,
            bias_seed_em_iters: 11,
            dump_bias_models: false,
            rad_codec: ChunkCodec::None,
            num_bootstraps: 0,
            num_gibbs_samples: 0,
            thinning_factor: 16,
            no_length_correction: false,
            model_single_frag_prob: true,
            no_frag_length_dist: false,
            fld_mean: 250.0,
            fld_sd: 25.0,
            fld_max: 1000,
            forgetting_factor: 0.65,
            sig_digits: 3,
            max_read_occ: 250,
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

pub use salmon_core::{Diagnostic, MissingMetaField};

/// Quantification results (also written to disk by [`write_outputs`]).
///
/// The abundance vectors are indexed by reference id and cover *all* references
/// including decoys; the output writer selects the rows that belong in
/// `quant.sf`. The remaining fields are the counters and provenance that end up
/// in `meta_info.json`.
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
    /// fragment mass (sum of equivalence-class counts) that could not be assigned
    /// to any transcript in the final min-alpha redistribution because every
    /// member of the class was truncated; reported, not rescaled. Normally 0.
    pub inference_truncated_mass: f64,
    pub num_eq_classes: usize,
    /// index of the first decoy reference (`None` if the index has no decoys).
    /// The decoy block is `[first_decoy_index, first_decoy_index + num_decoys)`;
    /// any references beyond it are short (sub-`k`, never-seeded) transcripts that
    /// are still reported in `quant.sf` with 0 reads.
    pub first_decoy_index: Option<usize>,
    /// number of decoy references in the contiguous decoy block (see
    /// [`Self::first_decoy_index`]); 0 when the index has no decoys
    pub num_decoys: usize,
    /// whether the index collapsed duplicate sequences (for meta_info)
    /// whether the index was built with `--keepDuplicates`; `None` when the
    /// index predates recording it, which is reported as a missing field rather
    /// than guessed as `false`
    pub keep_duplicates: Option<bool>,
    /// fragments dropped because their best alignment was to a decoy
    pub num_decoy_fragments: u64,
    /// fragments whose mates dovetail (overlap past each other)
    pub num_dovetail_fragments: u64,
    /// fragments that had candidate mappings but none survived validation/filtering
    pub num_fragments_filtered_vm: u64,
    /// per-fragment candidate alignments dropped for scoring below threshold,
    /// summed over fragments that nonetheless mapped
    pub num_alignments_below_threshold_for_mapped_fragments_vm: u64,
    /// where the fragment-length distribution came from; reported as
    /// `frag_length_source` in `meta_info.json`. Always
    /// `FragLengthSource::Reads` here — reads mode trains the FLD during the
    /// mapping pass — but carried so every mode reports the field.
    pub frag_len_source: salmon_model::FragLengthSource,
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
    /// EM/VBEM convergence: iterations actually run and whether the
    /// relative-difference criterion was met before `max_iter`
    pub em_iters: u32,
    pub em_converged: bool,
    /// total wall-clock seconds for the quantification call
    pub total_seconds: f64,
    /// peak resident set size in KiB (Linux `VmHWM`); 0 if unavailable
    pub peak_rss_kb: u64,
    /// library format observed by the auto-detector (when it ran with enough
    /// samples), for provenance and the strandedness sanity check; `None` when
    /// detection did not run or saw too few samples
    pub detected_library_type: Option<String>,
    /// structured run diagnostics / bad-input warnings (also emitted to the log)
    pub diagnostics: Vec<Diagnostic>,
}

/// Run quantification end-to-end, writing outputs and returning the results.
///
/// See the module documentation for the phases; the `PhaseTimer` marks below
/// delimit them in the log.
pub fn quantify(opts: &QuantOptions) -> Result<QuantResult> {
    anyhow::ensure!(
        opts.write_mappings.is_none() || opts.write_bam.is_none(),
        "--writeMappings and --writeBam are mutually exclusive"
    );
    let start_time = asctime_now();
    let run_timer = std::time::Instant::now();
    let mut timer = PhaseTimer::new();
    // The raw reference nucleotides are only needed for selective-alignment
    // extension (non-sketch mapping) and for sequence-dependent bias correction
    // (seq/GC; positional bias is sequence-free). In sketch mode without seq/GC
    // bias we can skip loading them entirely — with a genome decoy this is a
    // multi-GB resident-memory saving.
    let need_refseq = !opts.sketch || opts.seq_bias || opts.gc_bias;
    let salmon = SalmonIndex::load_with_opts(&opts.index_dir, need_refseq)
        .with_context(|| format!("loading index {}", opts.index_dir.display()))?;
    let num_refs = salmon.num_refs();
    timer.mark("index_load");

    let eq_builder = EquivalenceClassBuilder::new();
    // Gaussian-prior fragment length distribution (mean 250, sd 25); empirical
    // observations from concordant pairs refine it.
    let mut fld =
        FragmentLengthDistribution::new(1.0, opts.fld_max, opts.fld_mean, opts.fld_sd, 4, 0.5, 1);
    // Prime the online PMF snapshot from the prior so it is never empty: once
    // `use_aux` turns on (after the pre-burn-in count) a fragment could otherwise
    // be processed before the first mini-batch boundary refresh and fall back to a
    // racy live read. With it primed, every per-fragment read is from a stable
    // snapshot, keeping exact-duplicate transcripts symmetric.
    fld.refresh_online();
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
    let nthreads = if opts.num_threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        opts.num_threads
    };
    // SAM sink for `--writeMappings` (header written here).
    let sam_writer: Option<sam::SamWriter> = match &opts.write_mappings {
        Some(path) => {
            let cmd = format!("salmon quant -i {}", opts.index_dir.display());
            Some(sam::SamWriter::create(path, &salmon, &cmd).context("opening SAM output")?)
        }
        None => None,
    };
    let bam_output: Option<bam::BamOutput> = match &opts.write_bam {
        Some(path) => {
            let cmd = format!("salmon quant -i {}", opts.index_dir.display());
            Some(
                bam::BamOutput::create(path, &salmon, &cmd, nthreads, opts.bam_compress_threads)
                    .context("opening BAM output")?,
            )
        }
        None => None,
    };
    // RAD sink for `--writeRad` (prelude + file tags written here). Sketch vs
    // selective-alignment profile follows the mapping mode salmon is already in.
    let mut rad_writer: Option<salmon_rad::RadOutputWriter> = match &opts.write_rad {
        Some(path) => {
            let names: Vec<&str> = (0..num_refs).map(|t| salmon.ref_name(t)).collect();
            let lens: Vec<u32> = (0..num_refs).map(|t| salmon.ref_len(t) as u32).collect();
            let profile = if opts.sketch {
                salmon_rad::RadProfile::Sketch
            } else {
                salmon_rad::RadProfile::SelectiveAlignment
            };
            // Reserve an FLD slot matching the cached log-PMF length (max_val+1).
            let fld_len = opts.fld_max + 1;
            Some(
                salmon_rad::RadOutputWriter::create(
                    path,
                    &names,
                    &lens,
                    opts.is_paired(),
                    profile,
                    fld_len,
                    opts.rad_codec,
                    // Identity of the index these mappings were made against, so a
                    // requant reports the provenance of the run that produced the
                    // records rather than of whatever index it was handed.
                    &salmon_rad::WriterProvenance {
                        mapping_type: if opts.sketch {
                            salmon_rad::MappingType::Pseudo
                        } else {
                            salmon_rad::MappingType::Mapping
                        },
                        index: Some(salmon_rad::IndexProvenance {
                            seq_hash: salmon.info().seq_hash.clone(),
                            name_hash: salmon.info().name_hash.clone(),
                            seq_hash512: salmon.info().seq_hash512.clone(),
                            name_hash512: salmon.info().name_hash512.clone(),
                            decoy_seq_hash: salmon.info().decoy_seq_hash.clone(),
                            decoy_name_hash: salmon.info().decoy_name_hash.clone(),
                            // `None` for an index built before this was
                            // recorded; the RAD then omits the tag, so a
                            // requant reports it unknown rather than `false`.
                            keep_duplicates: salmon.info().keep_duplicates,
                        }),
                        // No BAM behind reads; nothing an aligner said to carry.
                        source_programs: Vec::new(),
                    },
                )
                .context("opening RAD output")?,
            )
        }
        None => None,
    };
    // Auto-detect the library type when requested (`-l A`).
    let auto_detect = LibraryFormat::is_auto(&opts.lib_type);
    let read_type = if opts.is_paired() {
        ReadType::PairedEnd
    } else {
        ReadType::SingleEnd
    };
    // Run the library-type detector ALWAYS (not just under `-l A`): under an
    // explicit `-l` it still observes the true orientation distribution (sampled
    // pre-strand-filter in `processor::record`), which powers the
    // `library_type_mismatch` diagnostic. It never overrides an explicit `-l` —
    // the resolution below gates its format on `auto_detect`. The cost is bounded
    // (sampling stops after the detector locks in) → effectively free.
    let detector = Some(LibraryTypeDetector::new(read_type));

    // `--deterministic`: an order-independent FLD / library-format accumulator,
    // populated during the mapping pass instead of the online FLD + prefix
    // detector (see `processor::record_discrete`), then baked at end of pass.
    let discrete_fld = (opts.deterministic_fld && opts.is_paired())
        .then(|| salmon_model::DiscreteFld::new(opts.fld_max));
    // `--deterministic` + bias: also build naive eq-classes during mapping so a
    // rough end-of-mapping EM can bake `initial_abundances` to seed the requant's
    // bias model — letting the requant be a single fused RAD read.
    let want_rough_abund = discrete_fld.is_some()
        && !opts.no_length_correction
        && (opts.seq_bias || opts.gc_bias || opts.pos_bias);
    // Orientation-tagged naive eq-classes for the rough seed EM (filtered by the
    // resolved library type at end of mapping; see `processor::record_discrete`).
    let naive_eq = want_rough_abund.then(salmon_eqclass::NaiveEqBuilder::new);

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
        // Built before `Shared` because `Shared` borrows the plan's busy meter,
        // and before the processor because opening the inputs is what decides
        // the decode/map split the pool is sized to.
        let groups: Vec<Vec<PathBuf>> = if opts.is_paired() {
            opts.mates1
                .iter()
                .zip(opts.mates2.iter())
                .map(|(a, b)| vec![a.clone(), b.clone()])
                .collect()
        } else {
            opts.unmated.iter().map(|p| vec![p.clone()]).collect()
        };
        let mut plan = decode::plan(
            &groups,
            nthreads,
            opts.decoder,
            decode::engagement_policy(opts.sketch, opts.thread_policy.as_deref())?,
        )?;

        let shared = Shared {
            busy: &plan.busy,
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
            model_single_frag_prob: opts.model_single_frag_prob,
            no_frag_length_dist: opts.no_frag_length_dist,
            // Assemble eq-classes only if something consumes them: the EM (unless
            // `--skipQuant`) or an eq-class dump. A map-only `--writeRad` run skips
            // assembly but still trains the FLD (to bake) — see `processor::record`.
            build_eq: !opts.skip_quant || opts.dump_eq || opts.dump_eq_weights,
            // The online estimate is consumed only by the EM-bound FLD / bias
            // weighting; with `--skipQuant` it has no consumer, so train the FLD
            // deterministically (also makes a baked FLD reproduce on requant).
            use_online: !opts.skip_quant,
            num_processed,
            num_mapped,
            num_orphan: &num_orphan,
            num_decoy: &num_decoy,
            num_dovetail: &num_dovetail,
            num_frags_filtered_vm: &num_frags_filtered_vm,
            num_below_threshold_vm: &num_below_threshold_vm,
            unmapped_names: unmapped_collector.as_ref(),
            sam: sam_writer.as_ref(),
            bam: bam_output.as_ref(),
            rad: rad_writer.as_ref(),
            discrete_fld: discrete_fld.as_ref(),
            naive_eq: naive_eq.as_ref(),
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
        let readers = std::mem::take(&mut plan.readers);
        let kind = if opts.is_paired() {
            CollectionType::Paired
        } else {
            CollectionType::Single
        };
        let collection = Collection::new(readers, kind)
            .map_err(|e| anyhow::anyhow!("building read collection: {e}"))?;
        let mapped = if opts.is_paired() {
            collection.process_parallel_paired_pool(&mut proc, &plan.map_pool, None)
        } else {
            collection.process_parallel_pool(&mut proc, &plan.map_pool, None)
        };

        // Settle the broker before propagating a mapping error, and never let a
        // broker failure become the run's error. piscem shipped the opposite:
        // a broker error surfaced after mapping completed, unwound past the
        // output finalisation, and left a file that parsed cleanly but declared
        // no content. Failures are recorded, not raised.
        let diagnostics = plan.broker.finish();
        // The engagement threshold is 1 + C/P, i.e. 1/producer_cost_share, so the
        // broker's own solved share *is* the constant this mode should use. Log it
        // rather than inferring the same number from a wall-time crossover sweep:
        // it is the quantity the rule is derived from, it comes out of a single
        // run, and it does not need runs long enough to beat timing noise.
        if let Some(report) = &diagnostics.report {
            if let Some(model) = &report.final_model {
                let share = model.producer_cost_share;
                if share > 0.0 {
                    tracing::info!(
                        producer_cost_share = share,
                        producer_cost_share_uncertainty = model.producer_cost_share_uncertainty,
                        implied_min_threads_per_stream = (1.0 / share).ceil() as usize,
                        mode = if opts.sketch { "sketch" } else { "selective-alignment" },
                        "thread broker cost model"
                    );
                }
            }
        }
        mapped.map_err(|e| anyhow::anyhow!("mapping failed: {e}"))?;
        if let Some(err) = &diagnostics.failure {
            tracing::warn!(
                "thread broker stopped early ({err:?}); mapping completed with the \
                 split it had last applied"
            );
        }
    }
    // `--deterministic`: replace the (untrained) online FLD with one built ONCE,
    // deterministically, from the order-independent accumulator, and take the
    // library format it inferred (overriding the unused prefix detector).
    let det_fmt = if let Some(acc) = &discrete_fld {
        let (df, fmt) = acc.finish(opts.fld_mean, opts.fld_sd);
        fld = df;
        fmt
    } else {
        None
    };
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
    timer.mark("mapping");
    if let Some(sw) = &sam_writer {
        sw.flush().context("flushing SAM output")?;
    }
    if let Some(output) = bam_output {
        output.finish().context("finalizing BAM output")?;
    }
    // NB: `rad_writer` is finalized *after* EM below, so the cached FLD and the
    // abundance estimates can be baked into the header (single-pass deterministic
    // requant + exact FLD parity). The mapping pass scope above has closed, so all
    // worker references to `rad_writer` are already dropped; it sits idle until
    // then.

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
    let library_type = if let Some(f) = det_fmt {
        f.canonical().to_string()
    } else if auto_detect {
        // Auto mode only: adopt the prefix detector's format. Under an explicit
        // `-l` the (now always-on) detector is used purely for the mismatch
        // diagnostic and must NOT override the user's choice.
        match &detector {
            Some(det) => det.final_format().canonical().to_string(),
            None => opts.lib_type.clone(),
        }
    } else {
        opts.lib_type.clone()
    };

    // ---- effective lengths from the (now frozen) FLD ------------------------
    //
    // Effective length is the divisor that turns fragment counts into
    // abundances: how many positions a fragment could have started at, given how
    // long fragments actually are. It can only be computed now, because the
    // fragment-length distribution was being learned throughout the mapping pass.
    // `cache()` freezes it, after which every lookup is an array read.
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
    //
    // `finish` drains the concurrent builder into a flat, sorted vector, which is
    // what makes the class ordering (and hence the EM) reproducible.
    let mut collapsed = eq_builder.finish();
    collapsed.update_eff_lengths(&eff_lengths);
    let num_eq_classes = collapsed.len();
    timer.mark("eff_length_collapse");

    if opts.dump_eq || opts.dump_eq_weights {
        dump_eq_classes(&opts.output_dir, &salmon, &collapsed, opts.dump_eq_weights)
            .context("writing eq_classes.txt.gz")?;
    }
    let bias_on = (opts.seq_bias || opts.gc_bias || opts.pos_bias) && !opts.no_length_correction;
    // salmon runs a *single* offline EM: after a short burn-in (`targetIt = 10`
    // iterations) it corrects effective lengths in place using that early
    // abundance estimate to weight the expected bias models, then continues the
    // same alpha vector to convergence. We mirror that: when bias-correcting,
    // run EM only `bias_seed_em_iters` iterations here to weight the expected
    // models, then warm-start the single full convergence after correction
    // (below). Without bias, this is the final EM and runs to convergence.
    // Avoids a wasteful second full convergence (~20s on the 36M PE set).
    if !opts.skip_quant {
        tracing::info!(
            "estimating abundances ({})",
            if opts.em.use_vbem { "VBEM" } else { "EM" }
        );
    }
    // The packed CSR layout is built ONCE and reused by every consumer below: the
    // burn-in EM, the warm-started EM after bias correction, and bootstrap /
    // Gibbs. `optimize`/`optimize_with_init` each build their own internally, so
    // going through them meant re-copying labels, counts and weights (hundreds of
    // MB on a human-scale run) two to three times over data that is bit-identical
    // apart from `combined`, which `refresh_combined` rewrites in place.
    //
    // `weights` (the raw conditional weights) is only read by Gibbs, so it is
    // populated only when Gibbs will actually run.
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
    let mut em = if opts.skip_quant && want_rough_abund {
        // `--deterministic` + bias map-only pass: run a SHORT EM on the NAIVE
        // eq-classes (uniform-weight, with library-incompatible placements dropped
        // against the resolved type) to produce rough abundances baked as
        // `initial_abundances` — these seed the requant's bias expected model.
        // Deliberately under-converged (a rough seed avoids overfitting the bias
        // model on uncorrected estimates) and computed from the fixed FLD's
        // effective lengths, so it is deterministic.
        let resolved = expected_format.or(det_fmt);
        let mut naive_collapsed = naive_eq.as_ref().unwrap().finish(resolved);
        naive_collapsed.update_eff_lengths(&eff_lengths);
        let mut ro = opts.em.clone();
        ro.min_iter = opts.bias_seed_em_iters;
        ro.max_iter = opts.bias_seed_em_iters;
        ro.min_alpha = 0.0;
        salmon_infer::optimize(&naive_collapsed, num_refs, &ro, Some(&eff_lengths))
    } else if opts.skip_quant {
        // `--skipQuant`: emit equivalence classes + library type + metadata but
        // skip the optimizer (and Gibbs/bootstrap below, and quant.sf). Leave
        // abundances at zero.
        salmon_infer::EmResult {
            alphas: vec![0.0; num_refs],
            iters: 0,
            converged: true,
            dropped_mass: 0.0,
        }
    } else if bias_on {
        let mut pre = opts.em.clone();
        pre.min_iter = opts.bias_seed_em_iters;
        pre.max_iter = opts.bias_seed_em_iters;
        // No min-alpha truncation: these alphas continue into the warm-started
        // EM, and a zeroed alpha can never recover under the multiplicative
        // update (salmon keeps one continuous, untruncated vector).
        pre.min_alpha = 0.0;
        optimize_packed_with_init(&packed, &pre, true, None, Some(&eff_lengths))
    } else {
        optimize_packed_with_init(&packed, &opts.em, true, None, Some(&eff_lengths))
    };

    // ---- bias correction: re-estimate effective lengths --------------------
    // Build the expected (background) bias models from the initial abundance
    // estimate, recompute bias-corrected effective lengths, then re-run
    // inference. Composes `--seqBias`, `--gcBias`, and `--posBias` in any
    // combination via the unified convolution (salmon's updateEffectiveLengths).
    let mut bias_dump = BiasDump::default();
    if bias_on && !opts.skip_quant {
        use rayon::prelude::*;
        // Expected-bias models are built only over quantification targets: decoys
        // (the contiguous tail at `first_decoy_index`) are never expressed and so
        // contribute nothing, but a decoy chromosome would otherwise hand one
        // rayon worker an O(ref_len) sweep over hundreds of Mbp. See issue #1019.
        let num_targets = salmon.info().first_decoy_index.unwrap_or(num_refs);
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
                num_targets,
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
        // Precompute the 5'/3' (obs − exp) log-bias tables once; the per-position
        // factor build in the correction sweep evaluates these instead of
        // re-evaluating both models per context (~1.4× on that build).
        let seq_tab = seq.as_ref().map(|(of, or, ef, er)| {
            (
                salmon_model::LogBiasTable::new(of, ef),
                salmon_model::LogBiasTable::new(or, er),
            )
        });

        // Fragment-GC observed + expected models -> clamped ratio model.
        let store = gc_store;
        let gc_ratio_model = if let Some(m) = gcbias_obs {
            let mut obs = m.into_inner().unwrap();
            let mut exp = salmon_model::build_expected_gc(
                num_targets,
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
                num_targets,
                |t| salmon.ref_len(t) as usize,
                &em.alphas,
                &eff_lengths,
                &fld_cdf,
                length_quantiles.as_ref().unwrap(),
                k,
            );
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

        // `--dumpBiasModels` (hidden): write the observed+expected seq/GC/pos
        // models to a text file for debugging / C++-parity comparison. The output
        // dir is normally created later by `write_outputs`, so ensure it here.
        if opts.dump_bias_models {
            std::fs::create_dir_all(&opts.output_dir)
                .with_context(|| format!("creating {}", opts.output_dir.display()))?;
            dump_bias_models_to_file(&opts.output_dir.join("bias_models.txt"), &bias_dump)
                .context("writing bias_models.txt")?;
        }

        let corrected: Vec<f64> = (0..num_refs)
            .into_par_iter()
            .map(|tid| {
                // Decoys (the contiguous tail at `first_decoy_index`) are never
                // quantified or reported, so their bias-corrected effective length
                // is unused. Skipping them here avoids a second, single-threaded
                // #1019-style stall: a decoy that catches a few non-decoy-dominated
                // fragments gets a tiny non-zero alpha, and the per-reference
                // `corrected_effective_length_full` sweep over a whole genome-decoy
                // chromosome (up to ~250 Mb) then serializes the entire phase while
                // the other worker threads sit idle. Mirrors PR #1020's exclusion of
                // decoys from the expected bias models.
                if tid >= num_targets || em.alphas[tid] < 1e-8 {
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
                    seq: seq_tab.as_ref().map(|(f, r)| (f, r)),
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
        // Only the combined weights changed, so patch them into the existing
        // packed layout rather than rebuilding it.
        packed.refresh_combined(&collapsed);
        // Warm-start the single full convergence from the burn-in alphas, as
        // salmon continues the same vector after its in-loop correction.
        let warm = std::mem::take(&mut em.alphas);
        em = optimize_packed_with_init(&packed, &opts.em, true, Some(&warm), Some(&eff_lengths));
    }
    // The min-alpha truncation is applied inside the EM as a mass-preserving
    // redistribution (see `redistribute_truncated`): the truncated mass flows to
    // eq-class co-members, so `counts` already sum to `total_count` minus only the
    // mass in fully-truncated classes (`dropped_mass`, normally 0). No rescale.
    let inference_truncated_mass = em.dropped_mass;
    let (em_iters, em_converged) = (em.iters, em.converged);
    let counts = em.alphas;

    // Finalize RAD output now that the FLD and abundances are known: bake the
    // cached FLD (always — exact-parity, single-pass requant) and, when this run
    // actually quantified, the abundance estimates (a prior for future bias-aware
    // requant; the precise iteration count to bake is deferred to that work —
    // likely under-converged to avoid overfitting bias models). `--skipQuant`
    // leaves `counts` at zero, so only the FLD is baked.
    if let Some(mut rw) = rad_writer.take() {
        // Bake the FLD (in `--deterministic` mode this is the order-independent
        // `DiscreteFld`-built one; otherwise the online-trained FLD) and, when this
        // run quantified, the abundance estimates (a prior for future bias-aware
        // requant). `--skipQuant` leaves `counts` at zero, so only the FLD is baked.
        rw.set_frag_length_dist(fld.log_pmf());
        // Bake abundances from a full quant, or the rough EM run for `--deterministic`
        // + bias (its `--skipQuant` map pass) — both seed a requant's bias model.
        if !opts.skip_quant || want_rough_abund {
            rw.set_initial_abundances(&counts);
        }
        // Bake the resolved library format so a reader can apply concordance
        // filtering under `-l A` without re-inferring it. `--deterministic`: the
        // order-independent detected format; auto: the prefix detector's; explicit:
        // the parsed `-l` format.
        let resolved_fmt = match (det_fmt, &detector) {
            (Some(f), _) => Some(f),
            (None, Some(d)) => Some(d.final_format()),
            (None, None) => expected_format,
        };
        if let Some(f) = resolved_fmt {
            rw.set_library_format(f.format_id());
        }
        // What the pass *saw*, which the records cannot show: a RAD holds only
        // the fragments that mapped, so without these a requant would report a
        // 100% mapping rate by construction.
        rw.set_map_counters(salmon_rad::MapCounters {
            num_processed: num_processed.load(Ordering::Relaxed),
            num_dovetail: num_dovetail.load(Ordering::Relaxed),
            num_filtered_vm: num_frags_filtered_vm.load(Ordering::Relaxed),
            num_below_threshold_vm: num_below_threshold_vm.load(Ordering::Relaxed),
            num_decoy_fragments: num_decoy.load(Ordering::Relaxed),
        });
        rw.finalize().context("finalizing RAD output")?;
    }

    timer.mark("em_bias");

    // ---- posterior uncertainty (bootstrap / Gibbs) + ambiguity --------------
    // The packed CSR layout (piscem-infer style) makes these parallel-friendly.
    let ambig = salmon_infer::ambiguity_counts(&packed);
    // Same `posterior` the packed layout was built for, so the run cannot decide
    // it needs Gibbs after the weights it reads were skipped.
    let bootstraps: Vec<Vec<f64>> = match posterior {
        salmon_infer::PosteriorMethod::None => Vec::new(),
        salmon_infer::PosteriorMethod::Bootstrap { replicates } => {
            salmon_infer::bootstrap(&packed, &opts.em, replicates, 0x5A13_0000)
        }
        salmon_infer::PosteriorMethod::Gibbs { samples } => {
            // Gibbs prior follows the main optimizer (salmon): with VBEM and a
            // per-transcript prior it is `max(1.0, vbPrior)`; with plain EM it is
            // 1e-3 per transcript. The rust VBEM uses a constant per-transcript prior.
            let prior = if opts.em.use_vbem {
                opts.em.vb_prior.max(1.0)
            } else {
                1e-3
            };
            let gopts = salmon_infer::GibbsOptions {
                num_samples: samples,
                thinning: opts.thinning_factor,
                prior,
                per_transcript_prior: true,
            };
            salmon_infer::gibbs_sample(&packed, &eff_lengths, &counts, &gopts, 0x6217_0000)
        }
    };
    timer.mark("posterior");

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

    // ---- run diagnostics: bad-input red flags from already-computed aggregates
    // (evaluated once here → no per-fragment cost). Emitted to the log AND to
    // meta_info.json's `diagnostics` array for machine-readable consumption.
    let np = num_processed.load(Ordering::Relaxed);
    let nm = num_mapped.load(Ordering::Relaxed);
    let detected_library_type = det_fmt
        .as_ref()
        .map(|f| f.canonical().to_string())
        .or_else(|| {
            detector
                .as_ref()
                .filter(|d| d.can_guess())
                .map(|d| d.final_format().canonical().to_string())
        });
    let diagnostics = salmon_core::input_diagnostics(
        np,
        nm,
        auto_detect,
        &library_type,
        detected_library_type.as_deref(),
    );
    for d in &diagnostics {
        if d.severity == "error" {
            tracing::error!("{}", d.message);
        } else {
            tracing::warn!("{}", d.message);
        }
    }

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
        inference_truncated_mass,
        num_eq_classes,
        first_decoy_index: salmon.info().first_decoy_index,
        num_decoys: salmon.info().num_decoys,
        keep_duplicates: salmon.info().keep_duplicates,
        num_decoy_fragments: num_decoy.load(Ordering::Relaxed),
        num_dovetail_fragments: num_dovetail.load(Ordering::Relaxed),
        num_fragments_filtered_vm: num_frags_filtered_vm.load(Ordering::Relaxed),
        num_alignments_below_threshold_for_mapped_fragments_vm: num_below_threshold_vm
            .load(Ordering::Relaxed),
        frag_len_source: salmon_model::FragLengthSource::Reads,
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
        em_iters,
        em_converged,
        total_seconds: run_timer.elapsed().as_secs_f64(),
        peak_rss_kb: salmon_core::peak_rss_kb(),
        detected_library_type,
        diagnostics,
    };

    tracing::info!("writing results to {}", opts.output_dir.display());
    write_outputs(opts, &result)?;
    timer.mark("output");
    Ok(result)
}

/// Current local time as an asctime-style string (`Wed Jun 10 20:34:58 2026`),
/// matching salmon's `start_time`/`end_time` format.
///
/// Matched exactly because downstream tools parse these fields.
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
/// Serialize the equivalence classes for inspection or for a downstream tool.
///
/// This is the intermediate the EM actually consumes, so dumping it is the way to
/// see what salmon inferred *before* the statistics: which transcript sets the
/// reads were compatible with, and how many reads fell into each.
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
