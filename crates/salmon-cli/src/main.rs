//! `salmon` command-line interface (Rust port).
//!
//! Provides the two most-used subcommands so far: `index` (build a salmon
//! index over a transcriptome) and `quant` (quantify from FASTQ reads, via
//! selective alignment or the pseudoalignment-only `--sketch` path). Flag
//! names mirror C++ salmon where they overlap. Alignment-based `quant -a`,
//! `quantmerge`, and `alevin` are stubbed for later phases.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};

// mimalloc as the global allocator: the quant hot path is highly multithreaded
// and allocation-heavy, where mimalloc markedly outperforms the system allocator.
#[cfg(not(feature = "sysalloc"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use salmon_align::{quantify_alignments, AlignQuantOptions};
use salmon_index::{build as build_index, IndexBuildOptions};
use salmon_quant::{quantify_with_aligner, ProgressCounters, QuantOptions};

/// Build the optional selective-alignment backend for `--gpu`. Returns `None`
/// for the default CPU path. With the `gpu` feature, `--gpu` acquires a GPU
/// (falling back to a CPU full-length reference backend if none is present);
/// without it, `--gpu` is a hard error.
#[cfg(feature = "gpu")]
fn build_alignment_backend(
    gpu: bool,
) -> anyhow::Result<Option<Box<dyn salmon_map::Aligner + Sync>>> {
    if !gpu {
        return Ok(None);
    }
    match salmon_gpu::GpuAligner::new() {
        Some(g) => {
            tracing::info!("GPU alignment backend ready (full-length mode)");
            Ok(Some(Box::new(g)))
        }
        None => {
            tracing::warn!("--gpu: no GPU adapter found; using the CPU full-length backend");
            Ok(Some(Box::new(salmon_gpu::RefAligner)))
        }
    }
}

#[cfg(not(feature = "gpu"))]
fn build_alignment_backend(
    gpu: bool,
) -> anyhow::Result<Option<Box<dyn salmon_map::Aligner + Sync>>> {
    anyhow::ensure!(
        !gpu,
        "--gpu requested but salmon was built without GPU support; rebuild with `--features gpu`"
    );
    Ok(None)
}

use std::io::IsTerminal;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

/// The active mapping progress bar, if one is being shown. Tracing log writes
/// route through its `suspend` so status lines print cleanly above the live
/// spinner instead of garbling it.
static ACTIVE_PROGRESS: Mutex<Option<indicatif::ProgressBar>> = Mutex::new(None);

/// `tracing` writer that prints through the active progress bar (if any), so
/// log output and the spinner coexist; otherwise writes straight to stderr.
/// Keeps all status/progress on the same handle (stderr), matching C++ salmon.
#[derive(Clone, Copy)]
struct ProgressAwareWriter;
impl std::io::Write for ProgressAwareWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match ACTIVE_PROGRESS.lock().unwrap().as_ref() {
            Some(pb) => pb.suspend(|| std::io::stderr().write_all(buf))?,
            None => std::io::stderr().write_all(buf)?,
        }
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        std::io::stderr().flush()
    }
}
impl tracing_subscriber::fmt::MakeWriter<'_> for ProgressAwareWriter {
    type Writer = ProgressAwareWriter;
    fn make_writer(&self) -> Self::Writer {
        *self
    }
}

/// RAII handle for the live mapping spinner: spawns a thread that polls the
/// shared [`ProgressCounters`] and updates the spinner; on drop it stops the
/// thread, clears the bar, and unregisters it.
struct ProgressGuard {
    bar: indicatif::ProgressBar,
    stop: Arc<AtomicBool>,
    handle: Option<std::thread::JoinHandle<()>>,
}

impl ProgressGuard {
    /// Start a mapping spinner driven by `counters`.
    fn start(counters: Arc<ProgressCounters>) -> Self {
        let bar = indicatif::ProgressBar::new_spinner();
        // A swimming-salmon spinner, a dimmed elapsed clock, and humanized live
        // counts: "<fish>  [elapsed]  N fragments · R/s · M mapped (P%)".
        bar.set_style(
            indicatif::ProgressStyle::with_template(
                "{spinner:.green.bold} {elapsed_precise:.dim}  {msg}",
            )
            .unwrap_or_else(|_| indicatif::ProgressStyle::default_spinner())
            .tick_strings(&[
                "><(((°>      ",
                " ><(((°>     ",
                "  ><(((°>    ",
                "   ><(((°>   ",
                "    ><(((°>  ",
                "     ><(((°> ",
                "      ><(((°>",
            ]),
        );
        // ~4 Hz repaint — lively but light over SSH/tmux.
        bar.enable_steady_tick(std::time::Duration::from_millis(250));
        *ACTIVE_PROGRESS.lock().unwrap() = Some(bar.clone());
        let stop = Arc::new(AtomicBool::new(false));
        let (bc, sc) = (bar.clone(), stop.clone());
        let handle = std::thread::spawn(move || {
            let start = std::time::Instant::now();
            while !sc.load(Ordering::Relaxed) {
                let p = counters.processed.load(Ordering::Relaxed);
                let m = counters.mapped.load(Ordering::Relaxed);
                let secs = start.elapsed().as_secs_f64().max(1e-3);
                let rate = (p as f64 / secs) as u64;
                let pct = if p > 0 {
                    100.0 * m as f64 / p as f64
                } else {
                    0.0
                };
                bc.set_message(format!(
                    "{} fragments · {}/s · {} mapped ({pct:.1}%)",
                    indicatif::HumanCount(p),
                    indicatif::HumanCount(rate),
                    indicatif::HumanCount(m),
                ));
                std::thread::sleep(std::time::Duration::from_millis(500));
            }
        });
        ProgressGuard {
            bar,
            stop,
            handle: Some(handle),
        }
    }
}

impl Drop for ProgressGuard {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Relaxed);
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
        self.bar.finish_and_clear();
        *ACTIVE_PROGRESS.lock().unwrap() = None;
    }
}

#[derive(Parser)]
#[command(
    name = "salmon",
    version,
    about = "RNA-seq transcript quantification (Rust port)"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
    /// Reduce logging to warnings and errors only.
    #[arg(short = 'q', long = "quiet", global = true)]
    quiet: bool,
    /// Accepted for compatibility with C++ salmon. salmon 2.0 never contacts the
    /// network to check for a newer release, so this flag (and the
    /// `SALMON_NO_VERSION_CHECK` environment variable) is a no-op.
    #[arg(long = "no-version-check", global = true)]
    no_version_check: bool,
}

// clap subcommand args structs differ in size (IndexArgs vs the larger QuantArgs);
// this enum is constructed exactly once at startup, so the size delta is irrelevant.
#[allow(clippy::large_enum_variant)]
#[derive(Subcommand)]
enum Command {
    /// Build a salmon index from a transcriptome FASTA.
    Index(IndexArgs),
    /// Quantify transcript abundances from FASTQ reads.
    Quant(QuantArgs),
    /// Merge a column across multiple samples' quant files into a matrix.
    Quantmerge(QuantMergeArgs),
    /// Diagnostic: per-read best-mapping detail (placement, seed coverage, score).
    DebugMap(DebugMapArgs),
    /// Single-cell quantification (removed; redirects to alevin-fry).
    #[command(disable_help_flag = true)]
    Alevin(AlevinArgs),
    /// Perform super-secret operation.
    #[command(hide = true)]
    Swim,
}

/// `salmon alevin` is not part of the Rust release — single-cell quantification
/// lives in the alevin-fry ecosystem. We accept the subcommand (and swallow any
/// arguments) only so that legacy `salmon alevin ...` invocations get a clear,
/// actionable redirect instead of a clap "invalid subcommand" error.
#[derive(Args)]
struct AlevinArgs {
    #[arg(trailing_var_arg = true, allow_hyphen_values = true, num_args = 0.., hide = true)]
    _ignored: Vec<String>,
}

#[derive(Args)]
struct QuantMergeArgs {
    /// Quantification directories (one per sample).
    #[arg(long = "quants", num_args = 1.., required = true)]
    quants: Vec<PathBuf>,
    /// Optional sample names (default: the directory basenames).
    #[arg(long = "names", num_args = 1..)]
    names: Vec<String>,
    /// Column to merge: {len, elen, tpm, numreads}.
    #[arg(short = 'c', long = "column", default_value = "TPM")]
    column: String,
    /// Merge gene-level quant.genes.sf instead of transcript-level quant.sf.
    #[arg(long = "genes")]
    genes: bool,
    /// Value written for targets missing from a sample.
    #[arg(long = "missing", default_value = "NA")]
    missing: String,
    /// Output matrix file.
    #[arg(short = 'o', long = "output", required = true)]
    output: PathBuf,
}

#[derive(Args)]
struct DebugMapArgs {
    #[arg(short = 'i', long = "index", required = true)]
    index: PathBuf,
    /// Single-end FASTQ (plain text).
    #[arg(short = 'r', long = "reads", required = true)]
    reads: PathBuf,
    /// Seed with reference-extended MEMs (cross unitig boundaries).
    #[arg(long = "refMEMs", conflicts_with = "unimems")]
    refmems: bool,
    /// Seed with true unitig-constrained uni-MEMs (pufferfish-style).
    #[arg(long = "uniMEMs", conflicts_with = "refmems")]
    unimems: bool,
}

#[derive(Args)]
struct IndexArgs {
    /// Transcript FASTA file(s).
    #[arg(short = 't', long = "transcripts", num_args = 1.., required = true)]
    transcripts: Vec<PathBuf>,
    /// Output index directory.
    #[arg(short = 'i', long = "index", required = true)]
    index: PathBuf,
    /// K-mer length (odd, ≤ 31 recommended).
    #[arg(short = 'k', long = "kmerLen", default_value_t = 31)]
    kmer_len: usize,
    /// Minimizer length (0 = auto).
    #[arg(short = 'm', long = "minimizerLen", default_value_t = 0)]
    minimizer_len: usize,
    /// Worker threads (0 = all cores).
    #[arg(short = 'p', long = "threads", default_value_t = 0)]
    threads: usize,
    /// Keep the intermediate compacted-dBG files.
    #[arg(long = "keepIntermediate")]
    keep_intermediate: bool,
    /// Keep the post-cleaning ("fixed") reference FASTA (non-ACGT bases replaced)
    /// instead of deleting it after the index is built.
    #[arg(long = "keepFixedFasta")]
    keep_fixed_fasta: bool,
    /// [accepted for salmon compatibility; no effect] TwoPaCo bloom-filter size —
    /// the Rust index builder (cf1-rs) does not expose this knob.
    #[arg(long = "filterSize")]
    filter_size: Option<usize>,
    /// Directory for build intermediates (cleaned FASTA + cDBG files); the final
    /// index still lands in --index. Defaults to the index dir. (salmon's --tmpdir)
    #[arg(long = "tmpdir")]
    tmpdir: Option<PathBuf>,
    /// Retain exact-duplicate transcript sequences instead of collapsing them
    /// (salmon collapses duplicates by default).
    #[arg(long = "keepDuplicates")]
    keep_duplicates: bool,
    /// GENCODE-format references: truncate each transcript name at the first `|`
    /// (e.g. `ENST...|ENSG...|...` -> `ENST...`).
    #[arg(long = "gencode")]
    gencode: bool,
    /// File of decoy sequence names (one per line). The named records (which must
    /// appear last in the FASTA) are indexed but flagged as decoys.
    #[arg(short = 'd', long = "decoys")]
    decoys: Option<PathBuf>,
    /// Don't clip poly-A tails from the ends of target sequences. By default
    /// (matching salmon/pufferfish `FixFasta`), a reference ending in a run of
    /// at least 10 `A`s has its trailing `A`s trimmed before indexing (an
    /// all-`A` reference is dropped).
    #[arg(short = 'n', long = "no-clip")]
    no_clip: bool,
}

#[derive(Args)]
struct QuantArgs {
    /// Salmon index directory (reads mode).
    #[arg(short = 'i', long = "index", required_unless_present = "alignments")]
    index: Option<PathBuf>,
    /// Alignment-based mode: a BAM of reads aligned to the transcriptome.
    #[arg(short = 'a', long = "alignments")]
    alignments: Option<PathBuf>,
    /// Transcriptome FASTA (alignment mode `-a`): enables the alignment error model.
    #[arg(short = 't', long = "targets")]
    targets: Option<PathBuf>,
    /// Disable the alignment error model (alignment mode).
    #[arg(long = "noErrorModel")]
    no_error_model: bool,
    /// Library type (e.g. IU, ISR, A for auto).
    #[arg(short = 'l', long = "libType", default_value = "A")]
    lib_type: String,
    /// Mate-1 FASTQ file(s) for paired-end reads.
    #[arg(short = '1', long = "mates1", num_args = 1..)]
    mates1: Vec<PathBuf>,
    /// Mate-2 FASTQ file(s) for paired-end reads.
    #[arg(short = '2', long = "mates2", num_args = 1..)]
    mates2: Vec<PathBuf>,
    /// Single-end FASTQ file(s).
    #[arg(short = 'r', long = "unmatedReads", num_args = 1..)]
    unmated: Vec<PathBuf>,
    /// Output directory.
    #[arg(short = 'o', long = "output", required = true)]
    output: PathBuf,
    /// Transcript-to-gene map (GTF/GFF, or a 2-column TSV); also writes
    /// gene-level estimates to `quant.genes.sf`.
    #[arg(short = 'g', long = "geneMap", value_name = "FILE")]
    gene_map: Option<PathBuf>,
    /// Worker threads (0 = all cores).
    #[arg(short = 'p', long = "threads", default_value_t = 0)]
    threads: usize,
    /// Use the alignment-free pseudoalignment path.
    #[arg(long = "sketch")]
    sketch: bool,
    /// Sketch mode: only emit an orphan when the mate had no matching k-mers at
    /// all (strict). Default orphans whenever the mate has no consistent target,
    /// which tracks selective alignment more closely.
    #[arg(long = "sketchStrictOrphans")]
    sketch_strict_orphans: bool,
    /// Use the standard EM optimizer instead of VBEM.
    #[arg(long = "useEM")]
    use_em: bool,
    /// Disable the live mapping progress spinner (it is shown only on an
    /// interactive terminal and never when `--quiet`).
    #[arg(long = "no-progress")]
    no_progress: bool,
    /// Metagenomic preset (salmon's `--meta`): use plain EM (not VBEM), disable
    /// range-factorized/rich equivalence classes, and initialize abundances
    /// uniformly — settings better suited to metagenomic references. Overrides
    /// `--useEM`/`--rangeFactorizationBins`.
    #[arg(long = "meta")]
    meta: bool,
    /// Range-factorization bins for equivalence classes (0 disables).
    #[arg(long = "rangeFactorizationBins", default_value_t = 4)]
    range_factorization_bins: u32,
    /// Prior weight for strand-incompatible mappings (0 drops them).
    #[arg(long = "incompatPrior", default_value_t = 0.0)]
    incompat_prior: f64,
    /// Dump equivalence classes to aux_info/eq_classes.txt.gz (salmon format).
    #[arg(long = "dumpEq")]
    dump_eq: bool,
    /// Like --dumpEq, but also write per-transcript weights in each class.
    #[arg(long = "dumpEqWeights")]
    dump_eq_weights: bool,
    /// Write the names of unmapped fragments to aux_info/unmapped_names.txt.
    #[arg(long = "writeUnmappedNames")]
    write_unmapped_names: bool,
    /// Write per-mapping SAM records to this file (spoofed CIGAR, like salmon's
    /// standard `--writeMappings`).
    #[arg(short = 'z', long = "writeMappings")]
    write_mappings: Option<PathBuf>,
    /// Enable sequence-specific bias correction.
    #[arg(long = "seqBias")]
    seq_bias: bool,
    /// Enable fragment-GC bias correction.
    #[arg(long = "gcBias")]
    gc_bias: bool,
    /// Enable positional bias correction.
    #[arg(long = "posBias")]
    pos_bias: bool,
    /// Number of bootstrap replicates for posterior uncertainty.
    #[arg(long = "numBootstraps", default_value_t = 0)]
    num_bootstraps: u32,
    /// Number of Gibbs posterior samples (mutually exclusive with bootstraps).
    #[arg(long = "numGibbsSamples", default_value_t = 0)]
    num_gibbs_samples: u32,
    /// Gibbs thinning factor.
    #[arg(long = "thinningFactor", default_value_t = 16)]
    thinning_factor: u32,
    /// (Long reads) Oxford Nanopore model — not supported; use oarfish instead.
    #[arg(long = "ont")]
    ont: bool,
    /// Disable effective-length correction (use raw reference length).
    #[arg(long = "noLengthCorrection")]
    no_length_correction: bool,
    /// Minimum alignment score as a fraction of the perfect score.
    #[arg(long = "minScoreFraction", default_value_t = 0.65)]
    min_score_fraction: f32,
    /// Orphan chain sub-optimality threshold (salmon's `orphanChainSubThresh`).
    /// `0.0` (default) aligns every orphan candidate (more sensitive than salmon's
    /// 0.95, which prunes low-chain-coverage orphans before alignment).
    #[arg(long = "orphanChainSubThresh", default_value_t = 0.0)]
    orphan_chain_sub_thresh: f32,
    /// Score the full read with one DP instead of PuffAligner-style inter-MEM-gap scoring.
    #[arg(long = "fullLengthAlignment")]
    full_length_alignment: bool,
    /// Score selective alignment on the GPU (requires a build with
    /// `--features gpu`; runs on Metal or Vulkan via wgpu). Implies
    /// `--fullLengthAlignment` and produces output identical to the CPU
    /// full-length path. Falls back to a CPU full-length backend if no GPU
    /// adapter is found.
    #[arg(long = "gpu")]
    gpu: bool,
    /// Seed with reference-extended MEMs (cross unitig boundaries).
    #[arg(long = "refMEMs", conflicts_with = "unimems")]
    refmems: bool,
    /// Seed with true unitig-constrained uni-MEMs (pufferfish-style).
    #[arg(long = "uniMEMs", conflicts_with = "refmems")]
    unimems: bool,
    /// Match score for selective alignment (reads mode).
    #[arg(long = "ma", default_value_t = 2)]
    ma: i32,
    /// Mismatch penalty for selective alignment (reads mode).
    #[arg(long = "mp", default_value_t = 4)]
    mp: i32,
    /// Gap-open penalty for selective alignment (reads mode).
    #[arg(long = "go", default_value_t = 6)]
    go: i32,
    /// Gap-extend penalty for selective alignment (reads mode).
    #[arg(long = "ge", default_value_t = 2)]
    ge: i32,
    /// Consensus slack: a target is kept only if its best chain score is at least
    /// `(1 - slack)` of the max chain score for that mate (salmon default 0.35;
    /// `1.0` keeps every target).
    #[arg(long = "consensusSlack", default_value_t = 0.35)]
    consensus_slack: f32,
    /// Skip k-mers whose unitig occurs in more than this many references
    /// (repetitive-hit guard; salmon's `maxOccsPerHit`).
    #[arg(long = "maxOccsPerHit", default_value_t = 1000)]
    max_occs_per_hit: usize,
    /// VBEM per-feature Dirichlet prior weight.
    #[arg(long = "vbPrior", default_value_t = 1e-2)]
    vb_prior: f64,
    /// Mean of the fragment-length distribution prior.
    #[arg(long = "fldMean", default_value_t = 250.0)]
    fld_mean: f64,
    /// Standard deviation of the fragment-length distribution prior.
    #[arg(long = "fldSD", default_value_t = 25.0)]
    fld_sd: f64,
    /// Maximum fragment length tracked by the fragment-length distribution.
    #[arg(long = "fldMax", default_value_t = 1000)]
    fld_max: usize,
    /// Online-phase forgetting factor in (0.5, 1.0].
    #[arg(short = 'f', long = "forgettingFactor", default_value_t = 0.65)]
    forgetting_factor: f64,
    /// Initialize the optimizer uniformly instead of from the online estimates.
    #[arg(long = "initUniform")]
    init_uniform: bool,
    /// Interpret --vbPrior as a per-transcript prior (salmon's default).
    #[arg(long = "perTranscriptPrior")]
    per_transcript_prior: bool,
    /// Interpret --vbPrior as a per-nucleotide prior (scaled by effective length).
    #[arg(long = "perNucleotidePrior", conflicts_with = "per_transcript_prior")]
    per_nucleotide_prior: bool,
    /// Significant digits for the EffectiveLength and NumReads columns of quant.sf.
    #[arg(long = "sigDigits", default_value_t = 3)]
    sig_digits: u32,
    /// Discard a read/fragment that maps to more than this many places (reads mode).
    #[arg(short = 'w', long = "maxReadOcc", default_value_t = 200)]
    max_read_occ: usize,
    /// Allow dovetailed mappings (mates extending past each other) as concordant.
    #[arg(long = "allowDovetail")]
    allow_dovetail: bool,
    /// ksw2 DP bandwidth for selective alignment (reads mode).
    #[arg(long = "bandwidth", default_value_t = 15)]
    bandwidth: i32,
    /// Attempt to recover the mate of an orphaned read near its mapped partner
    /// (selective-alignment mode). Off by default.
    #[arg(long = "recoverOrphans")]
    recover_orphans: bool,
    /// Discard orphan mappings in selective-alignment (reads) mode: only
    /// concordantly-paired fragments are kept.
    #[arg(long = "discardOrphansQuasi")]
    discard_orphans_quasi: bool,
    /// An alignment to an annotated transcript is invalid if its score is
    /// `< decoyThreshold * bestDecoyScore` (selective-alignment mode).
    #[arg(long = "decoyThreshold", default_value_t = 1.0)]
    decoy_threshold: f64,
    /// Drop any mapping whose alignment probability
    /// `exp(-scoreExp * (best - score))` is below this (selective-alignment mode).
    #[arg(long = "minAlnProb", default_value_t = 1e-5)]
    min_aln_prob: f64,
    /// Soft-weight decay: a mapping's probability is proportional to
    /// `exp(-scoreExp * (bestScore - score))`. Larger downweights sub-optimal
    /// mappings more steeply. (salmon's --scoreExp; default 1.0)
    #[arg(long = "scoreExp", default_value_t = 1.0)]
    score_exp: f64,
    /// Instead of soft-weighting multimapping locations by alignment score, keep
    /// only the best-scoring mapping(s), each with equal weight.
    #[arg(long = "hardFilter")]
    hard_filter: bool,
    /// Allow soft-clipping of read ends during selective alignment: unaligned
    /// read-end bases are clipped rather than penalized. (salmon's --softclip)
    #[arg(long = "softclip")]
    softclip: bool,
    /// Allow soft-clipping only of read bases that overhang a transcript end
    /// (a restricted form of --softclip). (salmon's --softclipOverhangs)
    #[arg(long = "softclipOverhangs")]
    softclip_overhangs: bool,
    /// Per-target pre-merge chain sub-optimality threshold: keep chains scoring
    /// `>= best_chain_score * thresh` (selective-alignment mode). Range [0,1].
    /// Rust default 0.8 (salmon's default is 0.75).
    #[arg(long = "preMergeChainSubThresh", default_value_t = 0.8)]
    pre_merge_chain_sub_thresh: f32,
    /// Post-merge concordant chain-pair sub-optimality threshold: after pairing,
    /// keep pairs whose read-coverage is `>= best_pair_coverage * thresh`
    /// (paired-end, selective-alignment mode). Range [0,1]. Default 0.0 (off);
    /// salmon's default is 0.9, thresholded on chain score rather than coverage.
    #[arg(long = "postMergeChainSubThresh", default_value_t = 0.0)]
    post_merge_chain_sub_thresh: f32,
    /// Fragment-length sampling stride for the GC bias eff-length convolution.
    /// Larger = faster bias correction, coarser. (salmon's --biasSpeedSamp)
    #[arg(long = "biasSpeedSamp", default_value_t = 5)]
    bias_speed_samp: usize,
    /// Number of leading fragments used to train the auxiliary models
    /// (fragment-length / bias / error); they are fixed afterward.
    #[arg(long = "numAuxModelSamples", default_value_t = 5_000_000)]
    num_aux_model_samples: u64,
    /// Disable the lower threshold on how short bias correction may make an
    /// effective length (increases precision, reduces robustness). Experimental.
    #[arg(long = "noBiasLengthThreshold")]
    no_bias_length_threshold: bool,
    /// Number of read-position bins in the alignment error model
    /// (alignment mode only). (salmon's --numErrorBins)
    #[arg(long = "numErrorBins", default_value_t = 4)]
    num_error_bins: usize,
    /// Number of fragment-GC bins for the GC bias model. (salmon's --numGCBins)
    #[arg(long = "numGCBins", default_value_t = 25)]
    num_gc_bins: usize,
    /// Number of conditioning (context) bins for the GC bias model.
    /// (salmon's --conditionalGCBins)
    #[arg(long = "conditionalGCBins", default_value_t = 3)]
    conditional_gc_bins: usize,
    /// Discard orphan (single-mate) alignments in a paired library
    /// (alignment mode). Reads mode uses --discardOrphansQuasi.
    #[arg(long = "discardOrphans")]
    discard_orphans: bool,
    /// [accepted for salmon compatibility; no effect] alignment-mode mapping
    /// cache size — the Rust align path streams with bounded buffers.
    #[arg(long = "mappingCacheMemoryLimit")]
    mapping_cache_memory_limit: Option<usize>,
    /// [accepted for salmon compatibility; no effect] the Rust chainer uses a
    /// loss-less ref-distance early-break, not salmon's numRounds heuristic.
    #[arg(long = "disableChainingHeuristic")]
    disable_chaining_heuristic: bool,
    /// [accepted for salmon compatibility; no effect yet] bases to skip after a
    /// k-mer miss during seeding — lives in the piscem-rs seed-walk; not yet tunable.
    #[arg(long = "mismatchSeedSkip")]
    mismatch_seed_skip: Option<u32>,
    /// [accepted for salmon compatibility; now the default] the rank-bitvector
    /// GC representation it selects is used unconditionally — it is faster and
    /// ~2x leaner than the dense representation with identical results.
    #[arg(long = "reduceGCMemory")]
    reduce_gc_memory: bool,
    /// Skip transcript quantification (EM + Gibbs/bootstrap) and quant.sf, while
    /// still mapping, building equivalence classes (use with --dumpEq), detecting
    /// library type, and writing metadata. (salmon's --skipQuant)
    #[arg(long = "skipQuant")]
    skip_quant: bool,
    /// [not yet implemented] quantify from a previously-written equivalence-class
    /// file instead of mapping/alignments (an input mode, not the --dumpEq output).
    #[arg(long = "eqclasses")]
    eqclasses: Option<PathBuf>,
    /// Fragments processed before the auxiliary (fragment-length) model is
    /// applied during online inference. (salmon's --numPreAuxModelSamples;
    /// salmon's default is 1000000, this port's is 5000.)
    #[arg(long = "numPreAuxModelSamples", default_value_t = 5000)]
    num_pre_aux_model_samples: u64,
    /// [accepted; not yet implemented] disable fragment-length-distribution
    /// concordance in the per-fragment probability.
    #[arg(long = "noFragLengthDist")]
    no_frag_length_dist: bool,
    /// [accepted; not yet implemented] disable the single-end/orphan
    /// fragment-length probability estimate.
    #[arg(long = "noSingleFragProb")]
    no_single_frag_prob: bool,
    /// [accepted; not yet implemented] write a BAM of posterior-sampled
    /// alignments (alignment mode).
    #[arg(short = 's', long = "sampleOut")]
    sample_out: bool,
    /// [accepted; not yet implemented] also include unaligned reads in the
    /// sampled BAM (requires --sampleOut).
    #[arg(short = 'u', long = "sampleUnaligned")]
    sample_unaligned: bool,
    /// [accepted; not yet implemented] write qualities into the sampled BAM
    /// (only meaningful with --sampleOut / --writeMappings).
    #[arg(long = "writeQualities")]
    write_qualities: bool,
    /// [accepted; not yet implemented] hit-filtering policy (BEFORE/AFTER/BOTH/
    /// NONE) for selective alignment; the Rust port filters after chaining.
    #[arg(long = "hitFilterPolicy")]
    hit_filter_policy: Option<String>,
    /// [accepted; not yet implemented] cap on candidate mappings during orphan
    /// recovery (salmon's --maxRecoverReadOcc).
    #[arg(long = "maxRecoverReadOcc")]
    max_recover_read_occ: Option<u32>,
    /// [accepted; no effect] deprecated in salmon too — selective alignment is
    /// the default mapping mode (use --sketch for pseudoalignment).
    #[arg(long = "validateMappings")]
    validate_mappings: bool,
}

/// Resolve the `--uniMEMs` / `--refMEMs` flags into a seeding mode (clap
/// enforces they are mutually exclusive).
fn seed_mode(unimems: bool, refmems: bool) -> salmon_map::SeedMode {
    if unimems {
        salmon_map::SeedMode::UniMem
    } else if refmems {
        salmon_map::SeedMode::RefMem
    } else {
        salmon_map::SeedMode::Sparse
    }
}

fn run_index(args: IndexArgs) -> Result<()> {
    if args.filter_size.is_some() {
        tracing::warn!("--filterSize is accepted for salmon compatibility but has no effect: the Rust index builder (cf1-rs) does not expose the cDBG bloom-filter size.");
    }
    let mut opts = IndexBuildOptions::new(args.transcripts, args.index);
    opts.k = args.kmer_len;
    opts.m = args.minimizer_len;
    opts.threads = args.threads;
    opts.keep_intermediate = args.keep_intermediate;
    opts.keep_fixed_fasta = args.keep_fixed_fasta;
    opts.tmpdir = args.tmpdir;
    opts.keep_duplicates = args.keep_duplicates;
    opts.decoys = args.decoys;
    opts.gencode = args.gencode;
    opts.clip_polya = !args.no_clip;
    let info = build_index(&opts).context("index build failed")?;
    println!(
        "indexed {} references (k={}, m={})",
        info.num_refs, info.k, info.m
    );
    Ok(())
}

/// Long-read (`--ont`) quantification is intentionally not supported: for long
/// reads, oarfish is the recommended tool. Print a pointer and exit rather than
/// silently producing estimates that are inappropriate for long reads.
fn long_read_redirect() -> ! {
    eprintln!(
        "salmon (Rust): long-read quantification (--ont) is not supported.\n\
         For long-read (Oxford Nanopore / PacBio) transcript quantification, please use oarfish,\n\
         our dedicated long-read quantification tool:\n\
         \n    https://github.com/COMBINE-lab/oarfish\n\n\
         oarfish models long-read alignments directly and is the recommended choice for this data."
    );
    std::process::exit(2);
}

/// Aggregate transcript estimates to gene level and write `quant.genes.sf`.
fn write_gene_level(
    out_dir: &std::path::Path,
    gene_map: &std::path::Path,
    names: &[String],
    lengths: &[u32],
    eff_lengths: &[f64],
    tpm: &[f64],
    counts: &[f64],
) -> Result<()> {
    let map = salmon_core::genemap::read_transcript_gene_map(gene_map)
        .with_context(|| format!("reading gene map {}", gene_map.display()))?;
    let unmapped = salmon_core::genemap::write_gene_quant(
        &out_dir.join("quant.genes.sf"),
        names,
        lengths,
        eff_lengths,
        tpm,
        counts,
        &map,
    )
    .context("writing quant.genes.sf")?;
    let n_genes = map.values().collect::<std::collections::HashSet<_>>().len();
    if unmapped > 0 {
        eprintln!(
            "warning: {unmapped} transcript(s) had no entry in the gene map and were omitted from quant.genes.sf"
        );
    }
    println!("wrote gene-level estimates for {n_genes} genes to quant.genes.sf");
    Ok(())
}

fn run_quant(args: QuantArgs, quiet: bool) -> Result<()> {
    if args.ont {
        long_read_redirect();
    }
    // Accepted-for-compatibility flags that have no effect in the Rust port.
    if args.mapping_cache_memory_limit.is_some() {
        tracing::warn!("--mappingCacheMemoryLimit is accepted for salmon compatibility but has no effect: the Rust alignment path streams with bounded buffers (no mass-banking cache).");
    }
    if args.disable_chaining_heuristic {
        tracing::warn!("--disableChainingHeuristic is accepted for salmon compatibility but has no effect: the Rust chainer uses a loss-less ref-distance early-break, not salmon's numRounds heuristic.");
    }
    if args.mismatch_seed_skip.is_some() {
        tracing::warn!("--mismatchSeedSkip is accepted for salmon compatibility but is not yet implemented: the post-miss seed skip lives in the piscem-rs k-mer seed walk, which does not yet expose it.");
    }
    if args.reduce_gc_memory {
        tracing::warn!("--reduceGCMemory is a no-op: the rank-bitvector GC representation it selects is now used by default (faster and ~2x leaner, identical results).");
    }
    if args.eqclasses.is_some() {
        tracing::warn!("--eqclasses (quantify from a precomputed equivalence-class file) is not yet implemented and is ignored; mapping/alignment input is used instead.");
    }
    if args.no_frag_length_dist {
        tracing::warn!("--noFragLengthDist is accepted but not yet implemented and has no effect: the fragment-length distribution is still used in the per-fragment probability.");
    }
    if args.no_single_frag_prob {
        tracing::warn!("--noSingleFragProb is accepted but not yet implemented and has no effect.");
    }
    if args.sample_out || args.sample_unaligned || args.write_qualities {
        tracing::warn!("--sampleOut/--sampleUnaligned/--writeQualities (posterior-sampled BAM output) are accepted but not yet implemented and have no effect.");
    }
    if args.hit_filter_policy.is_some() {
        tracing::warn!("--hitFilterPolicy is accepted but not yet implemented and has no effect: the Rust port filters hits after chaining (salmon's AFTER default).");
    }
    if args.max_recover_read_occ.is_some() {
        tracing::warn!(
            "--maxRecoverReadOcc is accepted but not yet implemented and has no effect."
        );
    }
    if args.validate_mappings {
        tracing::warn!("--validateMappings has no effect (deprecated in salmon too): selective alignment is the default mapping mode; pass --sketch for pseudoalignment.");
    }
    // Bound the global rayon pool (used by the EM and bias passes) to the
    // requested thread count; otherwise it spans every core regardless of -p.
    let nthreads = if args.threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        args.threads
    };
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build_global();
    let out_dir = args.output.clone();
    let gene_map = args.gene_map.clone();
    // `--meta` (metagenomic preset) forces plain EM (no VBEM) and disables rich /
    // range-factorized equivalence classes, matching salmon's --meta. salmon's
    // preset also sets initUniform; the Rust offline EM already initializes
    // uniformly, so that part is inherent. --meta overrides --useEM/--rangeFactorizationBins.
    let use_vbem = !args.use_em && !args.meta;
    let range_factorization_bins = if args.meta {
        0
    } else {
        args.range_factorization_bins
    };
    if args.meta {
        tracing::info!(
            "--meta: metagenomic preset (plain EM, no range-factorized equivalence classes, uniform init)"
        );
    }
    // Alignment-based mode: quantify directly from a BAM.
    if let Some(bam) = args.alignments {
        let mut opts = AlignQuantOptions::new(bam, args.output);
        opts.lib_type = args.lib_type;
        opts.em.use_vbem = use_vbem;
        opts.range_factorization_bins = range_factorization_bins;
        opts.transcripts = args.targets;
        opts.no_error_model = args.no_error_model;
        opts.seq_bias = args.seq_bias;
        opts.gc_bias = args.gc_bias;
        opts.pos_bias = args.pos_bias;
        opts.incompat_prior = args.incompat_prior;
        opts.em.vb_prior = args.vb_prior;
        opts.em.per_nucleotide_prior = args.per_nucleotide_prior;
        opts.sig_digits = args.sig_digits;
        opts.fld_mean = args.fld_mean;
        opts.fld_sd = args.fld_sd;
        opts.fld_max = args.fld_max;
        opts.forgetting_factor = args.forgetting_factor;
        opts.init_uniform = args.init_uniform;
        opts.bias_speed_samp = args.bias_speed_samp;
        opts.num_aux_model_samples = args.num_aux_model_samples;
        opts.no_bias_length_threshold = args.no_bias_length_threshold;
        opts.num_error_bins = args.num_error_bins;
        opts.discard_orphans = args.discard_orphans;
        opts.gc_bins = args.num_gc_bins;
        opts.cond_gc_bins = args.conditional_gc_bins;
        opts.skip_quant = args.skip_quant;
        opts.num_pre_aux_model_samples = args.num_pre_aux_model_samples;
        // --scoreExp is selective-alignment-mode only (it scales the
        // best-minus-score soft weight); alignment mode has no such term.
        // Live progress spinner on an interactive terminal (unless --quiet/--no-progress).
        let progress = Arc::new(ProgressCounters::default());
        let guard = if !quiet && !args.no_progress && std::io::stderr().is_terminal() {
            opts.progress = Some(progress.clone());
            Some(ProgressGuard::start(progress.clone()))
        } else {
            None
        };
        let res = quantify_alignments(&opts).context("alignment-based quantification failed")?;
        drop(guard); // stop + clear the spinner before the summary
        let pct = if res.num_processed > 0 {
            100.0 * res.num_mapped as f64 / res.num_processed as f64
        } else {
            0.0
        };
        tracing::info!(
            "processed {} fragments, mapped {} ({:.2}%), {} equivalence classes",
            res.num_processed,
            res.num_mapped,
            pct,
            res.num_eq_classes
        );
        if let Some(gm) = &gene_map {
            write_gene_level(
                &out_dir,
                gm,
                &res.names,
                &res.lengths,
                &res.eff_lengths,
                &res.tpm,
                &res.counts,
            )?;
        }
        return Ok(());
    }

    // Reads (selective-alignment / pseudoalignment) mode.
    anyhow::ensure!(
        !args.mates1.is_empty() || !args.unmated.is_empty(),
        "no reads provided: pass -1/-2 (paired), -r (single-end), or -a (BAM)"
    );
    anyhow::ensure!(
        args.mates1.len() == args.mates2.len(),
        "the number of -1 and -2 files must match"
    );
    let index = args.index.context("reads mode requires -i/--index")?;

    let mut opts = QuantOptions::new(index, args.output);
    opts.mates1 = args.mates1;
    opts.mates2 = args.mates2;
    opts.unmated = args.unmated;
    opts.lib_type = args.lib_type;
    opts.num_threads = args.threads;
    opts.sketch = args.sketch;
    opts.sketch_strict_orphan = args.sketch_strict_orphans;
    opts.em.use_vbem = use_vbem;
    opts.range_factorization_bins = range_factorization_bins;
    opts.incompat_prior = args.incompat_prior;
    opts.dump_eq = args.dump_eq;
    opts.dump_eq_weights = args.dump_eq_weights;
    opts.write_unmapped_names = args.write_unmapped_names;
    opts.write_mappings = args.write_mappings;
    opts.seq_bias = args.seq_bias;
    opts.gc_bias = args.gc_bias;
    opts.pos_bias = args.pos_bias;
    opts.num_bootstraps = args.num_bootstraps;
    opts.num_gibbs_samples = args.num_gibbs_samples;
    opts.thinning_factor = args.thinning_factor;
    opts.no_length_correction = args.no_length_correction;
    opts.map_config.align.min_score_fraction = args.min_score_fraction;
    opts.map_config.pair.orphan_chain_sub_thresh = args.orphan_chain_sub_thresh;
    opts.map_config.align.full_length_alignment = args.full_length_alignment;
    opts.map_config.align.bandwidth = args.bandwidth;
    opts.map_config.pair.allow_dovetail = args.allow_dovetail;
    opts.map_config.align.softclip = args.softclip;
    opts.map_config.align.softclip_overhangs = args.softclip_overhangs;
    opts.map_config.score.score_exp = args.score_exp;
    // mapping policy (Tier 2): all default to the prior hardcoded behavior
    opts.map_config.recover_orphans = args.recover_orphans;
    if args.discard_orphans_quasi {
        opts.map_config.pair.allow_orphans = false;
    }
    opts.map_config.score.decoy_thresh = args.decoy_threshold;
    opts.map_config.score.min_aln_prob = args.min_aln_prob;
    opts.map_config.score.hard_filter = args.hard_filter;
    // chaining sub-optimality thresholds (Tier 2)
    opts.map_config.collect.chain.chain_subopt_thresh = args.pre_merge_chain_sub_thresh;
    opts.map_config.pair.post_merge_chain_sub_thresh = args.post_merge_chain_sub_thresh;
    // model toggles (Tier 2)
    opts.bias_speed_samp = args.bias_speed_samp;
    opts.num_aux_model_samples = args.num_aux_model_samples;
    opts.no_bias_length_threshold = args.no_bias_length_threshold;
    opts.gc_bins = args.num_gc_bins;
    opts.cond_gc_bins = args.conditional_gc_bins;
    opts.skip_quant = args.skip_quant;
    opts.num_pre_aux_model_samples = args.num_pre_aux_model_samples;
    opts.map_config.seed_mode = seed_mode(args.unimems, args.refmems);
    // alignment scoring (selective alignment)
    opts.map_config.align.match_score = args.ma as i8;
    opts.map_config.align.mismatch_pen = args.mp as i8;
    opts.map_config.align.gap_open_pen = args.go as i8;
    opts.map_config.align.gap_extend_pen = args.ge as i8;
    // chaining consensus + repetitive-hit guard
    opts.map_config.collect.consensus_fraction = 1.0 - args.consensus_slack;
    opts.map_config.collect.max_hit_occ = args.max_occs_per_hit;
    // inference + fragment-length-distribution knobs
    opts.em.vb_prior = args.vb_prior;
    opts.em.per_nucleotide_prior = args.per_nucleotide_prior;
    opts.sig_digits = args.sig_digits;
    opts.max_read_occ = args.max_read_occ;
    opts.fld_mean = args.fld_mean;
    opts.fld_sd = args.fld_sd;
    opts.fld_max = args.fld_max;
    opts.forgetting_factor = args.forgetting_factor;

    // Live mapping spinner on an interactive terminal (unless --quiet/--no-progress).
    let progress = Arc::new(ProgressCounters::default());
    let guard = if !quiet && !args.no_progress && std::io::stderr().is_terminal() {
        opts.progress = Some(progress.clone());
        Some(ProgressGuard::start(progress.clone()))
    } else {
        None
    };
    if args.gpu {
        // --gpu only scores in full-length mode (the batchable banded DP).
        opts.map_config.align.full_length_alignment = true;
    }
    let aligner = build_alignment_backend(args.gpu)?;
    let res = quantify_with_aligner(&opts, aligner.as_deref()).context("quantification failed")?;
    drop(guard); // stop + clear the spinner before the summary
    let pct = if res.num_processed > 0 {
        100.0 * res.num_mapped as f64 / res.num_processed as f64
    } else {
        0.0
    };
    // Status to stderr (same handle as the rest of salmon's logging, matching C++).
    tracing::info!(
        "processed {} fragments, mapped {} ({:.2}%), {} equivalence classes",
        res.num_processed,
        res.num_mapped,
        pct,
        res.num_eq_classes
    );
    if let Some(gm) = &gene_map {
        write_gene_level(
            &out_dir,
            gm,
            &res.names,
            &res.lengths,
            &res.eff_lengths,
            &res.tpm,
            &res.counts,
        )?;
    }
    Ok(())
}

/// Shown when a legacy C++ invocation uses a removed/renamed option or
/// subcommand, so users get an actionable pointer instead of a bare clap error.
const MIGRATION_HINT: &str =
    "\nnote: salmon 2.0 is the Rust rewrite and changed or removed some command-line \
options.\n      See the migration guide (MIGRATION.md) for the C++ → 2.0 flag mapping, \
and note that\n      indices must be rebuilt for 2.0 (the index format changed).";

fn main() -> Result<()> {
    // Parse with a friendly hint: on an unrecognized flag/subcommand (the common
    // failure mode for an old C++ command line), append the migration pointer to
    // clap's error. Help/version/other errors keep clap's normal behavior.
    let cli = match Cli::try_parse() {
        Ok(cli) => cli,
        Err(e) => {
            use clap::error::ErrorKind;
            if matches!(
                e.kind(),
                ErrorKind::UnknownArgument | ErrorKind::InvalidSubcommand
            ) {
                eprint!("{e}");
                eprintln!("{MIGRATION_HINT}");
                std::process::exit(2);
            }
            e.exit();
        }
    };
    // Single-cell quant was removed; redirect to alevin-fry and exit.
    if let Command::Alevin(_) = cli.command {
        eprintln!(
            "salmon alevin is not available in salmon 2.0.\n\n\
             Single-cell quantification has moved to the alevin-fry ecosystem:\n  \
             - piscem / salmon (this tool) for mapping → RAD\n  \
             - alevin-fry for cell-barcode/UMI resolution and counting\n\n\
             See https://github.com/COMBINE-lab/alevin-fry and \
             https://alevin-fry.readthedocs.io"
        );
        std::process::exit(1);
    }
    // RUST_LOG (if set) wins; otherwise `info`, or `warn` under --quiet.
    let default_level = if cli.quiet { "warn" } else { "info" };
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new(default_level)),
        )
        .with_writer(ProgressAwareWriter)
        .init();

    if cli.no_version_check {
        tracing::debug!(
            "--no-version-check accepted; salmon 2.0 performs no startup version check."
        );
    }

    match cli.command {
        Command::Index(args) => run_index(args),
        Command::Quant(args) => run_quant(args, cli.quiet),
        Command::Quantmerge(args) => run_quantmerge(args),
        Command::DebugMap(args) => run_debug_map(args),
        Command::Alevin(_) => unreachable!("alevin handled before logging init"),
        Command::Swim => {
            run_swim();
            Ok(())
        }
    }
}

/// `salmon swim` — the super-secret operation. (An easter egg carried over from
/// C++ salmon: it prints the ASCII banner and swims away.)
fn run_swim() {
    print!(
        r#"
    _____       __
   / ___/____ _/ /___ ___  ____  ____
   \__ \/ __ `/ / __ `__ \/ __ \/ __ \
  ___/ / /_/ / / / / / / / /_/ / / / /
 /____/\__,_/_/_/ /_/ /_/\____/_/ /_/

         ><(((°>   ><(((°>   ><(((°>

"#
    );
}

fn run_quantmerge(args: QuantMergeArgs) -> Result<()> {
    use salmon_core::quantmerge::{merge_quants, MergeColumn};
    let column = MergeColumn::parse(&args.column).with_context(|| {
        format!(
            "unrecognized --column '{}'; expected one of {{len, elen, tpm, numreads}}",
            args.column
        )
    })?;
    // sample names default to the directory basenames
    let names: Vec<String> = if args.names.is_empty() {
        args.quants
            .iter()
            .map(|d| {
                d.file_name()
                    .map(|s| s.to_string_lossy().into_owned())
                    .unwrap_or_default()
            })
            .collect()
    } else {
        anyhow::ensure!(
            args.names.len() == args.quants.len(),
            "--names ({}) must match the number of --quants ({})",
            args.names.len(),
            args.quants.len()
        );
        args.names.clone()
    };
    let quant_file = if args.genes {
        "quant.genes.sf"
    } else {
        "quant.sf"
    };
    let quant_paths: Vec<PathBuf> = args.quants.iter().map(|d| d.join(quant_file)).collect();
    for p in &quant_paths {
        anyhow::ensure!(p.is_file(), "sample quant file not found: {}", p.display());
    }
    let n_missing = merge_quants(&quant_paths, &names, column, &args.missing, &args.output)
        .context("merging quant files")?;
    if n_missing > 0 {
        eprintln!(
            "warning: {n_missing} missing entries written as \"{}\"",
            args.missing
        );
    }
    println!(
        "merged {} samples ({}) -> {}",
        quant_paths.len(),
        args.column,
        args.output.display()
    );
    Ok(())
}

fn run_debug_map(args: DebugMapArgs) -> Result<()> {
    use piscem_rs::mapping::hit_searcher::HitSearcher;
    use salmon_index::SalmonIndex;
    use salmon_map::{debug_best_mapping, MapConfig};
    use std::io::BufRead;

    let idx = SalmonIndex::load(&args.index).context("loading index")?;
    let mut hs = HitSearcher::new(idx.inner());
    let cfg = MapConfig {
        seed_mode: seed_mode(args.unimems, args.refmems),
        ..MapConfig::default()
    };

    let f = std::fs::File::open(&args.reads).context("opening reads")?;
    let mut lines = std::io::BufReader::new(f).lines();
    println!("name\ttid\trefname\tref_pos\tis_fw\tchain_cov\tread_len\tcov_frac\tfull_score\tperfect\tscore_frac\tnum_mems");
    while let Some(Ok(header)) = lines.next() {
        let seq = match lines.next() {
            Some(Ok(s)) => s,
            _ => break,
        };
        let _ = lines.next(); // '+'
        let _ = lines.next(); // qual
        let name = header
            .trim_start_matches('@')
            .split_whitespace()
            .next()
            .unwrap_or("");
        if let Some(d) = debug_best_mapping(idx.inner(), &mut hs, &idx, seq.as_bytes(), &cfg) {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{:.3}\t{}",
                name,
                d.tid,
                idx.ref_name(d.tid as usize),
                d.ref_pos,
                d.is_fw,
                d.chain_cov,
                d.read_len,
                d.chain_cov as f64 / d.read_len as f64,
                d.full_score,
                d.perfect,
                d.full_score as f64 / d.perfect as f64,
                d.num_mems,
            );
        }
    }
    Ok(())
}
