//! `salmon-quant`: the reads-mode mapping pass.
//!
//! Turns FASTQ input into an intermediate RAD file: parse reads in parallel
//! (paraseq) → map each fragment ([`salmon_map`]) → write its placements to
//! the RAD, while deriving — order-independently — the fragment-length
//! distribution and library format ([`salmon_model`]), which are baked into
//! the RAD header together with the run's mapping counters.
//!
//! [`map_to_rad`] estimates **no abundances**. Quantification is the second
//! phase: `salmon_align::quantify_rad` reads the RAD this pass wrote. That
//! split is what makes `salmon quant` byte-identical at any thread count; the
//! `salmon` binary's reads-mode driver chains the two phases. (Before 2.12 this
//! crate also had a one-pass "online" path that quantified during mapping; it
//! was deprecated in 2.6.0 and removed in 2.12.0.)
//!
//! Selective alignment is the default; set [`QuantOptions::sketch`] for the
//! alignment-free pseudoalignment path.
//!
//! # The shape of a run
//!
//! 1. **Load the index** (optionally without the reference sequence, which is
//!    only needed for alignment and sequence-dependent bias).
//! 2. **Map**, in parallel over the input reads, writing the RAD records and
//!    any SAM/BAM output, and tallying the order-independent fragment-length /
//!    library-format accumulator (paired) or the whole-input library-type
//!    detector (single-end).
//! 3. **Finalize** the fragment-length distribution and, when bias correction
//!    was requested, run a short EM over naive equivalence classes for rough
//!    abundances that seed the requant's bias model.
//! 4. **Bake** the FLD, library format, counters and seed abundances into the
//!    RAD header, and write the mapping-pass metadata (what `--skipQuant`
//!    ships).

mod bam;
pub mod decode;
mod gpu_map;
mod mapping_record;
mod output;
mod processor;
mod sam;

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};

use anyhow::{Context, Result};

use piscem_rs::io::fastx::{Collection, CollectionType};
use piscem_rs::mapping::hit_searcher::SkippingStrategy;

use salmon_core::{LibraryFormat, PhaseTimer, ReadType};
use salmon_index::{RefseqLoad, SalmonIndex};
pub use salmon_infer::{BootstrapMode, EmAccel, EmOptions, EmStop};
use salmon_map::MapConfig;
use salmon_model::{FragmentLengthDistribution, LibraryTypeDetector};
pub use salmon_rad::ChunkCodec;

use processor::{QuantProcessor, Shared, UnmappedWriter};

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
    /// The classes are built by the quantification pass, so these are for the
    /// two-phase driver to forward there; [`map_to_rad`] ignores them.
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
    /// an explicit location for the run's RAD (`--writeRad`), which the
    /// two-phase driver then keeps. [`map_to_rad`] ignores it and writes to the
    /// path it is given; sketch or selective-alignment profile is chosen from
    /// `sketch`.
    pub write_rad: Option<PathBuf>,
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
    /// disable length correction: every transcript gets the same constant
    /// effective length, [`salmon_model::NO_LENGTH_CORRECTION_EFF_LEN`]
    /// (`--noLengthCorrection`)
    pub no_length_correction: bool,
    /// divide by the raw reference length instead of the fragment-length-derived
    /// effective length (`--noEffectiveLengthCorrection`); `no_length_correction`
    /// wins if both are set
    pub no_effective_length_correction: bool,
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
    /// whether the user supplied `--fldMean`/`--fldSD` (rather than taking the
    /// defaults); forwarded to the deterministic requant, which widens a merely
    /// default single-end prior for the single-end length term.
    pub fld_prior_user_supplied: bool,
    /// significant digits for the EffectiveLength and NumReads columns of
    /// `quant.sf` (`--sigDigits`, salmon default 3)
    pub sig_digits: u32,
    /// discard a fragment that maps to more than this many places (salmon's
    /// `--maxReadOcc`, default 250)
    pub max_read_occ: usize,
    /// fragment-length sampling stride for the GC bias convolution (salmon's
    /// `--biasSpeedSamp`, default 5)
    pub bias_speed_samp: usize,
    /// disable the lower barrier on how short bias correction may make an
    /// effective length (salmon's `--noBiasLengthThreshold`)
    pub no_bias_length_threshold: bool,
    /// number of fragment-GC bins in the GC bias model (salmon's `--numGCBins`,
    /// default 25)
    pub gc_bins: usize,
    /// number of conditioning (context) bins in the GC bias model (salmon's
    /// `--conditionalGCBins`, default 3)
    pub cond_gc_bins: usize,
    /// stop after the mapping pass, keeping the RAD as the run's output
    /// (salmon's `--skipQuant`). For the two-phase driver; [`map_to_rad`]
    /// never estimates abundances, so it ignores this.
    pub skip_quant: bool,
    /// Optional shared progress counters. When `Some`, [`map_to_rad`] reports
    /// processed/mapped fragment counts here as it runs so the caller can drive
    /// a live progress display. `None` (the default) disables sharing.
    pub progress: Option<std::sync::Arc<ProgressCounters>>,
}

impl QuantOptions {
    /// The effective-length mode `--noLengthCorrection` and
    /// `--noEffectiveLengthCorrection` select.
    pub fn length_correction(&self) -> salmon_model::LengthCorrection {
        salmon_model::LengthCorrection::from_flags(
            self.no_length_correction,
            self.no_effective_length_correction,
        )
    }

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
            no_effective_length_correction: false,
            model_single_frag_prob: true,
            no_frag_length_dist: false,
            fld_mean: 250.0,
            fld_sd: 25.0,
            fld_max: 1000,
            fld_prior_user_supplied: false,
            sig_digits: 3,
            max_read_occ: 250,
            bias_speed_samp: 5,
            no_bias_length_threshold: false,
            gc_bins: salmon_model::gcbias::DEFAULT_GC_BINS,
            cond_gc_bins: salmon_model::gcbias::DEFAULT_COND_BINS,
            skip_quant: false,
            progress: None,
        }
    }

    fn is_paired(&self) -> bool {
        !self.mates1.is_empty()
    }
}

pub use salmon_core::{Diagnostic, MissingMetaField};

/// Mapping-pass results (also written to disk by [`write_outputs`]): the
/// counters and provenance that end up in `meta_info.json` and
/// `lib_format_counts.json`. There are no abundances here — those come from
/// quantifying the RAD the pass wrote.
///
/// Per-reference vectors are indexed by reference id and cover *all*
/// references including decoys.
#[derive(Debug, Clone)]
pub struct MapResult {
    pub names: Vec<String>,
    pub lengths: Vec<u32>,
    pub num_processed: u64,
    pub num_mapped: u64,
    /// mapped fragments placed as orphans (only one mate of a pair mapped)
    pub num_orphan: u64,
    /// fragments with at least one mapping strand-compatible with the expected
    /// library format (`lib_format_counts.json`'s `num_compatible_fragments`)
    pub num_compatible_fragments: u64,
    /// fragments that mapped but had *no* strand-compatible mapping
    /// (`lib_format_counts.json`'s `num_incompatible_fragments`). On a stranded
    /// library this is the wrong-strand fraction, the cheapest available signal
    /// of double-stranded (genomic DNA) carry-over; see #1130.
    ///
    /// Counted whatever `--incompatPrior` is: under the default 0 these
    /// fragments are dropped and excluded from `num_mapped`, above 0 they are
    /// down-weighted and *included* in it.
    ///
    /// Both tallies cover only fragments filtered against a known expected
    /// format, so under `-l A` fragments consumed before the detector locks in
    /// are in neither, and `num_compatible_fragments + num_incompatible_fragments`
    /// can fall short of the mapped total.
    pub num_incompatible_fragments: u64,
    /// per-observed-format fragment counts (`counts[format_id]`, one count per
    /// distinct format among a mapped fragment's placements), tallied for every
    /// mapped fragment whatever the expected format. The raw histogram behind
    /// `lib_format_counts.json`'s per-format keys and its
    /// concordant/inconsistent/strand-bias derivations.
    pub lib_format_counts: salmon_core::LibFormatCountsArray,
    /// the GPU backend that ran the `--sketch` k-mer lookup, or `None` when
    /// the mapping ran on the CPU only
    pub gpu_backend: Option<&'static str>,
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
    /// iterations of the rough bias-seed EM (0 when it did not run)
    pub em_iters: u32,
    /// whether that seed EM met its convergence criterion (it is cut short on
    /// purpose, so usually not), and its per-component summary
    pub em_converged: bool,
    pub em_components: Option<salmon_infer::ComponentSummary>,
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

/// Map the reads to a RAD file at `rad_path`, baking the fragment-length
/// distribution, library format, mapping counters and (with bias correction
/// requested) rough seed abundances into its header, and write the
/// mapping-pass metadata into `opts.output_dir`.
///
/// This estimates **no abundances**: quantify the RAD with
/// `salmon_align::quantify_rad` for those. [`QuantOptions::write_rad`],
/// [`QuantOptions::skip_quant`] and the eq-class dump flags describe the whole
/// run and are for the caller's two-phase driver; this pass ignores them.
///
/// See the module documentation for the phases; the `PhaseTimer` marks below
/// delimit them in the log.
pub fn map_to_rad(opts: &QuantOptions, rad_path: &Path) -> Result<MapResult> {
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
    // The bias models read whole transcripts as bytes; mapping only reads short
    // windows, so it holds them at two bits per base.
    let refseq_load = if opts.seq_bias || opts.gc_bias {
        RefseqLoad::Bytes
    } else if !opts.sketch {
        RefseqLoad::Packed
    } else {
        RefseqLoad::None
    };
    // shared with the thread that builds the GPU k-mer lookup
    let salmon = std::sync::Arc::new(
        SalmonIndex::load_with_opts(&opts.index_dir, refseq_load)
            .with_context(|| format!("loading index {}", opts.index_dir.display()))?,
    );
    let num_refs = salmon.num_refs();
    timer.mark("index_load");

    // piscem's default; the GPU lookup reproduces it exactly
    let skip = SkippingStrategy::Permissive;
    // `--sketch` k-mer lookup on the GPU (same mappings as the CPU path).
    // It is built on a thread while the CPU maps the first batches, so its
    // build (`gpu_kmer_table` in the timing log) does not delay the run.
    // shortcut: rebuilt every run; persist it beside the index if it ever
    // outlasts the mapping.
    let kmer_gpu = std::sync::Arc::new(std::sync::OnceLock::new());
    let gpu_lookup = opts.sketch
        && opts.em.gpu
        && salmon_infer::kmers::KmerGpu::compiled()
        && salmon.inner().poison_table().is_none();
    if gpu_lookup {
        let (salmon, slot) = (salmon.clone(), kmer_gpu.clone());
        let max_hit_occ = opts.map_config.collect.max_hit_occ;
        let max_read_occ = opts.max_read_occ.min(u32::MAX as usize) as u32;
        std::thread::spawn(move || {
            let t0 = std::time::Instant::now();
            gpu_map::urgent();
            let idx = salmon.inner();
            use salmon_infer::kmers::{KmerGpu, KmerTable};
            use salmon_map::sketch_hits as sh;
            // piscem's sshash dictionary by default; the flat k-mer table is
            // kept behind a switch while the two are compared
            // shortcut: env switch, drop the loser once measured on real data
            let build = || {
                if std::env::var("SALMON_GPU_LOOKUP").as_deref() == Ok("table") {
                    let (kstart, keep, kmer) = sh::unitig_kmers(idx);
                    KmerTable::build(idx.k(), kstart, keep, kmer).map(|t| (t, Vec::new()))
                } else {
                    sh::sshash_parts(idx).and_then(|p| {
                        // `build_sshash` takes only offsets below 2^32
                        let unitig_base = p.soff.iter().map(|&o| o as u32).collect();
                        KmerTable::build_sshash(
                            idx.k(),
                            p.m,
                            p.magic,
                            p.spss,
                            p.soff,
                            sh::kept(idx),
                            p.bstart,
                            p.offs,
                            // its minimizer table is built on the device where it can be
                            !KmerGpu::builds_table(),
                            salmon_infer::kmers::HEAVY_BUCKET,
                        )
                        .map(|t| (t, unitig_base))
                    })
                }
            };
            // its own pool, so its threads can take the build's priority; as
            // many threads as the global pool, which `-p` sized
            // opening the device takes ~50 ms: beside the CPU half of the build
            let both = || rayon::join(build, KmerGpu::available);
            let (table, available) = match rayon::ThreadPoolBuilder::new()
                .num_threads(rayon::current_num_threads())
                .start_handler(|_| gpu_map::urgent())
                .build()
            {
                Ok(pool) => pool.install(both),
                Err(_) => both(),
            };
            if !available {
                return;
            }
            match table
                .ok_or_else(|| "the index does not fit the GPU lookup".to_string())
                .and_then(|(table, unitig_base)| {
                    let mut gpu = KmerGpu::new(&table, skip == SkippingStrategy::Permissive)?;
                    // here, not on the first mapping batches
                    gpu.warm_up(&table);
                    // the target selection too, on the backends that have it
                    let select = sh::occurrences(idx).is_some_and(|(off, ent, ref_len_bits)| {
                        let occ = salmon_infer::kmers::Occurrences {
                            off,
                            ent,
                            ref_len_bits,
                        };
                        let max_occ = max_hit_occ.saturating_sub(1).min(u32::MAX as usize);
                        gpu.set_occurrences(&occ, max_occ as u32, max_read_occ)
                            .map_err(|e| tracing::debug!("GPU target selection off ({e})"))
                            .is_ok()
                    });
                    Ok(gpu_map::GpuLookup {
                        gpu,
                        unitig_base,
                        select,
                    })
                }) {
                Ok(gpu) => {
                    let device_mb = gpu.gpu.device_bytes().map_or(-1.0, |b| b as f64 / 1e6);
                    let select = gpu.select;
                    let _ = slot.set(gpu);
                    let elapsed_s = t0.elapsed().as_secs_f64();
                    tracing::info!(target: "salmon::timing", phase = "gpu_kmer_table", elapsed_s, device_mb, select, "phase complete");
                }
                Err(e) => tracing::warn!("GPU k-mer lookup off ({e}); mapping on the CPU"),
            }
        });
    }
    let (kmer_gpu_off, kmer_gpu_used) = (AtomicBool::new(false), AtomicBool::new(false));

    // The fragment-length prior (mean 250, sd 25 by default). A paired run
    // replaces it at end of pass with the order-independent `DiscreteFld`
    // estimate; a single-end run observes no lengths and bakes the prior.
    let mut fld =
        FragmentLengthDistribution::new(1.0, opts.fld_max, opts.fld_mean, opts.fld_sd, 4, 0.5, 1);
    // `processed`/`mapped` live in a (possibly caller-shared) `ProgressCounters`
    // so a CLI progress bar can poll them live; the rest are local.
    let progress = opts.progress.clone().unwrap_or_default();
    let num_processed = &progress.processed;
    let num_mapped = &progress.mapped;
    let num_orphan = AtomicU64::new(0);
    let num_frags_compat = AtomicU64::new(0);
    let num_frags_incompat = AtomicU64::new(0);
    let lib_format_counts: [AtomicU64; salmon_core::NUM_LIB_FORMATS] =
        std::array::from_fn(|_| AtomicU64::new(0));
    let num_decoy = AtomicU64::new(0);
    let num_dovetail = AtomicU64::new(0);
    let num_frags_filtered_vm = AtomicU64::new(0);
    let num_below_threshold_vm = AtomicU64::new(0);
    // `--writeUnmappedNames`: open the file up front and let workers stream into
    // it per batch. Opening here also surfaces an unwritable output directory
    // before the mapping pass rather than after it.
    let unmapped_collector: Option<std::sync::Mutex<UnmappedWriter>> = if opts.write_unmapped_names
    {
        std::fs::create_dir_all(opts.output_dir.join("aux_info")).context("creating aux_info")?;
        Some(std::sync::Mutex::new(
            UnmappedWriter::create(&opts.output_dir.join("aux_info").join("unmapped_names.txt"))
                .context("creating unmapped_names.txt")?,
        ))
    } else {
        None
    };
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
    // The RAD sink (prelude + file tags written here). Sketch vs
    // selective-alignment profile follows the mapping mode salmon is already in.
    // An `Option` only so it can be `take`n for finalization after the pass.
    let mut rad_writer: Option<salmon_rad::RadOutputWriter> = {
        let path = rad_path;
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
                        // The decoy boundary, so a requant can exclude the
                        // decoy block from quant.sf as this mode does
                        // (#1140). The header must carry every reference
                        // (records index into it), so the boundary is the
                        // only way a reader can tell decoys apart.
                        first_decoy_index: salmon.info().first_decoy_index.map(|i| i as u64),
                        num_decoys: (salmon.info().num_decoys > 0)
                            .then_some(salmon.info().num_decoys as u64),
                    }),
                    // No BAM behind reads; nothing an aligner said to carry.
                    source_programs: Vec::new(),
                },
            )
            .context("opening RAD output")?,
        )
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
    //
    // `--deterministic` single-end: the prefix would depend on the thread order,
    // and the paired accumulator below has no single-end form, so tally everything.
    let whole_input_detector = !opts.is_paired();
    let detector = Some(if whole_input_detector {
        LibraryTypeDetector::over_whole_input(read_type)
    } else {
        LibraryTypeDetector::new(read_type)
    });

    // Paired: an order-independent FLD / library-format accumulator, populated
    // during the mapping pass (see `processor::record_discrete`), then baked at
    // end of pass.
    let discrete_fld = opts
        .is_paired()
        .then(|| salmon_model::DiscreteFld::new(opts.fld_max));
    // `--deterministic` + bias: also build naive eq-classes during mapping so a
    // rough end-of-mapping EM can bake `initial_abundances` to seed the requant's
    // bias model — letting the requant be a single fused RAD read.
    let want_rough_abund = discrete_fld.is_some()
        && opts.length_correction().is_effective()
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
            decode::engagement_policy(
                opts.sketch,
                opts.is_paired(),
                opts.thread_policy.as_deref(),
            )?,
            // shortcut: decided before the device opens, so a build with a GPU
            // backend but no usable device still reads ahead; harmless but
            // two threads over `-p`
            gpu_lookup,
        )?;

        let shared = Shared {
            busy: &plan.busy,
            salmon: &salmon,
            detector: detector.as_ref(),
            map_cfg: &opts.map_config,
            sketch: opts.sketch,
            sketch_strict_orphan: opts.sketch_strict_orphan,
            max_read_occ: opts.max_read_occ,
            skip,
            gpu: Some(gpu_map::GpuMap {
                gpu: &kmer_gpu,
                off: &kmer_gpu_off,
                used: &kmer_gpu_used,
            }),
            expected_format,
            incompat_prior: opts.incompat_prior,
            ignore_incompat,
            num_processed,
            num_mapped,
            num_orphan: &num_orphan,
            num_frags_compat: &num_frags_compat,
            num_frags_incompat: &num_frags_incompat,
            lib_format_counts: &lib_format_counts,
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
        let stages = tracing::enabled!(target: "salmon::timing", tracing::Level::DEBUG);
        processor::stage::ON.store(stages, Ordering::Relaxed);
        let map_start = std::time::Instant::now();
        let mapped = if opts.is_paired() {
            collection.process_parallel_paired_pool(&mut proc, &plan.map_pool, None)
        } else {
            collection.process_parallel_pool(&mut proc, &plan.map_pool, None)
        };
        if stages {
            log_map_stages(map_start.elapsed().as_secs_f64(), nthreads);
        }

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
                        mode = if opts.sketch {
                            "sketch"
                        } else {
                            "selective-alignment"
                        },
                        "thread broker cost model"
                    );
                }
            }
            // The windowed model above is the broker's live estimate, and on a
            // run that never settles (a solved target clamped at the producer
            // floor, for instance) it is not populated at all. The integrated
            // ratio has neither problem: cumulative producer busy time over
            // cumulative total busy time, whole run, no windowing. It is also
            // the same busy-time definition the control law solves with, so
            // this -- not a CPU-seconds ratio, which counts coordinator
            // overhead the broker deliberately excludes -- is the number an
            // engagement threshold should be derived from.
            if let Some(pm) = &report.producer_measurement {
                let p = pm.busy_nanos as f64;
                let c = plan.busy.nanos() as f64;
                if p > 0.0 && c > 0.0 {
                    let share = p / (p + c);
                    tracing::info!(
                        producer_busy_s = p / 1e9,
                        consumer_busy_s = c / 1e9,
                        integrated_producer_cost_share = share,
                        implied_min_threads_per_stream = (1.0 / share).ceil() as usize,
                        mode = if opts.sketch {
                            "sketch"
                        } else {
                            "selective-alignment"
                        },
                        "thread broker integrated cost ratio"
                    );
                }
            }
        }
        // Name the inputs on the way out. A FASTQ parse failure arrives here as
        // paraseq's bare content string ("FASTX error: Invalid FASTQ separator:
        // N, expected '+'") with no idea which file produced it: `ProcessError`
        // carries no path, and paraseq drains every reader concurrently on
        // threads it owns, so the identity is gone before the Result crosses
        // back. `-1/-2/-r` each take many files, so that message on its own
        // names nothing the user can act on — this was the only input-error
        // class in the tree that did not chain a path. This is the last place
        // that still holds the paths, and it runs once per run (not per record),
        // so the context is free. With one input we can name it outright; with
        // several we say plainly that the list is the candidate set and not the
        // culprit, rather than implying we know which.
        if let Err(e) = mapped {
            let files: Vec<String> = groups
                .iter()
                .flatten()
                .map(|p| p.display().to_string())
                .collect();
            let msg = e.to_string();
            let err = anyhow::anyhow!("mapping failed: {msg}");
            // A decoder error (a truncated or corrupt compressed file) already
            // names its file; only a parse error still needs the candidates.
            if files.iter().any(|f| msg.contains(f.as_str())) {
                return Err(err);
            }
            return Err(err.context(match files.as_slice() {
                [one] => format!("while reading {one}"),
                many => format!(
                    "while reading one of the {} read files (the FASTX reader does not \
                     report which): {}",
                    many.len(),
                    many.join(", ")
                ),
            }));
        }
        if let Some(err) = &diagnostics.failure {
            tracing::warn!(
                "thread broker stopped early ({err:?}); mapping completed with the \
                 split it had last applied"
            );
        }
    }
    // Paired: replace the prior with the FLD built ONCE, deterministically, from
    // the order-independent accumulator, and take the library format it inferred
    // (overriding the prefix detector, kept only for the mismatch diagnostic).
    let det_fmt = if let Some(acc) = &discrete_fld {
        let (df, fmt) = acc.finish(opts.fld_mean, opts.fld_sd);
        fld = df;
        fmt
    } else if whole_input_detector {
        detector.as_ref().map(|d| d.final_format())
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
    if let Some(sw) = sam_writer {
        sw.finish().context("finalizing SAM output")?;
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
    // Every worker has joined, so nothing more will be written: flush and close.
    if let Some(collector) = unmapped_collector {
        collector
            .into_inner()
            .unwrap()
            .finish()
            .context("writing unmapped_names.txt")?;
    }

    // Resolve the library type: the detected format (when auto), else the
    // user-specified string. Fall back to a sensible default if detection saw
    // no usable samples.
    // Auto mode only: adopt the resolved format — the deterministic
    // accumulator's answer when it ran, else the prefix detector's. Under an
    // explicit `-l` both are computed purely for the mismatch diagnostic and
    // must NOT override the user's choice (the strand filter used the explicit
    // type, and this string labels `lib_format_counts.json`'s
    // `expected_format`).
    let library_type = if auto_detect {
        match det_fmt {
            Some(f) => f.canonical().to_string(),
            None => match &detector {
                Some(det) => det.final_format().canonical().to_string(),
                None => opts.lib_type.clone(),
            },
        }
    } else {
        LibraryFormat::parse(&opts.lib_type)
            .map(|f| f.canonical().to_string())
            .unwrap_or_else(|_| opts.lib_type.clone())
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
    let cond_means = fld.conditional_means();
    let length_correction = opts.length_correction();
    let eff_lengths: Vec<f64> = (0..num_refs)
        .map(|tid| length_correction.effective_length(&cond_means, salmon.ref_len(tid) as usize))
        .collect();
    timer.mark("eff_length");

    // ---- rough abundance seed (bias only) -----------------------------------
    // With bias correction requested, run a SHORT EM on the NAIVE eq-classes
    // (uniform-weight, with library-incompatible placements dropped against the
    // resolved type) to produce rough abundances baked as `initial_abundances`
    // — these seed the requant's bias expected model. Deliberately
    // under-converged (a rough seed avoids overfitting the bias model on
    // uncorrected estimates) and computed from the fixed FLD's effective
    // lengths, so it is deterministic. Otherwise no abundances are estimated
    // here: quantification happens in the requant pass.
    let (counts, em_iters, em_converged, em_components) = if want_rough_abund {
        let resolved = expected_format.or(det_fmt);
        let mut naive_collapsed = naive_eq.as_ref().unwrap().finish(resolved);
        naive_collapsed.update_eff_lengths(&eff_lengths);
        let mut ro = opts.em.clone();
        ro.min_iter = opts.bias_seed_em_iters;
        ro.max_iter = opts.bias_seed_em_iters;
        ro.min_alpha = 0.0;
        let em = salmon_infer::optimize(
            &naive_collapsed,
            num_refs,
            &ro,
            salmon_infer::EffLens::new(&eff_lengths),
        );
        // Cut short on purpose: its stopping is not a convergence claim.
        tracing::info!("EM: {} iterations to seed the bias models", em.iters);
        (em.alphas, em.iters, em.converged, em.components)
    } else {
        (vec![0.0; num_refs], 0, true, None)
    };

    // Finalize RAD output now that the FLD is known: bake the cached FLD
    // (always — exact-parity, single-pass requant) and, with bias requested, the
    // rough seed abundances.
    if let Some(mut rw) = rad_writer.take() {
        // The FLD is the order-independent `DiscreteFld`-built one (paired) or
        // the prior (single-end).
        rw.set_frag_length_dist(fld.log_pmf());
        // The rough seed EM's abundances seed the requant's bias model.
        if want_rough_abund {
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

    timer.mark("rad_finalize");

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

    let result = MapResult {
        names: (0..num_refs)
            .map(|i| salmon.ref_name(i).to_string())
            .collect(),
        lengths: (0..num_refs).map(|i| salmon.ref_len(i) as u32).collect(),
        num_processed: num_processed.load(Ordering::Relaxed),
        num_mapped: num_mapped.load(Ordering::Relaxed),
        num_orphan: num_orphan.load(Ordering::Relaxed),
        num_compatible_fragments: num_frags_compat.load(Ordering::Relaxed),
        num_incompatible_fragments: num_frags_incompat.load(Ordering::Relaxed),
        lib_format_counts: std::array::from_fn(|i| lib_format_counts[i].load(Ordering::Relaxed)),
        gpu_backend: kmer_gpu_used
            .load(Ordering::Relaxed)
            .then_some(salmon_infer::GPU_BACKEND)
            .flatten(),
        first_decoy_index: salmon.info().first_decoy_index,
        num_decoys: salmon.info().num_decoys,
        keep_duplicates: salmon.info().keep_duplicates,
        num_decoy_fragments: num_decoy.load(Ordering::Relaxed),
        num_dovetail_fragments: num_dovetail.load(Ordering::Relaxed),
        num_fragments_filtered_vm: num_frags_filtered_vm.load(Ordering::Relaxed),
        num_alignments_below_threshold_for_mapped_fragments_vm: num_below_threshold_vm
            .load(Ordering::Relaxed),
        // A single-end library has no fragment lengths to observe: the FLD is
        // the --fldMean/--fldSD prior, and saying "reads" would be a
        // provenance lie in the field that exists for provenance (#1140).
        frag_len_source: if opts.is_paired() {
            salmon_model::FragLengthSource::Reads
        } else {
            salmon_model::FragLengthSource::Prior
        },
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
        em_iters,
        em_converged,
        em_components,
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

pub(crate) use salmon_core::asctime_now;

/// Log the mapping loop's per-stage thread time (paired path) and what is left of `threads × wall`: reading, decompression and
/// waiting in paraseq's reader.
fn log_map_stages(wall_s: f64, threads: usize) {
    use processor::stage;
    let mut sum = 0.0;
    for (name, ns) in stage::NAMES.iter().zip(&stage::NS) {
        let s = ns.swap(0, Ordering::Relaxed) as f64 / 1e9;
        sum += s;
        tracing::debug!(target: "salmon::timing", stage = name, thread_s = s, "map stage");
    }
    tracing::debug!(
        target: "salmon::timing",
        wall_s,
        threads,
        thread_s = wall_s * threads as f64,
        outside_s = wall_s * threads as f64 - sum,
        frags_gpu = stage::FRAGS_GPU.swap(0, Ordering::Relaxed),
        frags_cpu_in_gpu_batch = stage::FRAGS_CPU_IN_GPU_BATCH.swap(0, Ordering::Relaxed),
        frags_selected = stage::FRAGS_SELECTED.swap(0, Ordering::Relaxed),
        "map stages (outside = reading, decompression, waiting)"
    );
}
