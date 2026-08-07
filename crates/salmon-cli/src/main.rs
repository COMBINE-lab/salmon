//! `salmon` command-line interface (Rust port).
//!
//! # What this layer does
//!
//! Almost nothing, deliberately. Its job is to turn a command line into one of
//! the option structs the library crates take (`IndexBuildOptions`,
//! `QuantOptions`, `AlignQuantOptions`, …), call the corresponding entry point,
//! and report the result. All the science lives in those crates, which is why
//! they can be used as libraries and tested without a command line.
//!
//! What does live here is everything that only matters to a person at a
//! terminal: argument validation and conflict detection, the progress bar,
//! logging setup, and choosing which mode a given combination of flags means.
//!
//! **A note on doc comments below.** The `///` comments on the argument structs
//! are not just documentation — `clap` turns them verbatim into the `--help`
//! text a user sees. So they are written for users, and anything addressed to a
//! reader of the source is a plain `//` comment instead.
//!
//! Provides the two most-used subcommands so far: `index` (build a salmon
//! index over a transcriptome) and `quant` (quantify from FASTQ reads, via
//! selective alignment or the pseudoalignment-only `--sketch` path). Flag
//! names mirror C++ salmon where they overlap. Alignment-based `quant -a` and
//! `quantmerge` are stubbed for later phases; `alevin` is a permanent redirect
//! to the alevin-fry ecosystem.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};

// mimalloc as the global allocator: the quant hot path is highly multithreaded
// and allocation-heavy, where mimalloc markedly outperforms the system allocator.
#[cfg(not(feature = "sysalloc"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use salmon_align::{
    project_genome_bam_to_rad, quantify_alignments, quantify_rad, AlignQuantOptions,
    ExplicitFldArgs, FldPolicy, GenomeProjectOptions,
};
use salmon_index::{build as build_index, IndexBuildOptions};
use salmon_quant::{quantify, ChunkCodec, EmAccel, ProgressCounters, QuantOptions};

use std::io::IsTerminal;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};

// ---------------------------------------------------------------------------
// Progress display
//
// A long run needs to show that it is making progress, but log lines and a live
// spinner both write to the terminal and will garble each other. The two types
// below solve that by routing every log write through the progress bar, which
// knows how to clear and redraw itself around foreign output.
// ---------------------------------------------------------------------------

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
#[command(name = "salmon", version, about = "RNA-seq transcript quantification")]
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

/// CLI surface for [`EmAccel`]: the EM/VBEM convergence acceleration scheme.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
enum EmAccelArg {
    /// Plain fixed-point iteration; output unchanged from historical salmon.
    #[default]
    None,
    /// SQUAREM acceleration (same fixpoint, far fewer M-steps; not byte-identical).
    Squarem,
    /// DAAREM: damped Anderson acceleration over a window of past iterates; faster
    /// than SQUAREM on high-dimensional problems. Same fixpoint; not byte-identical.
    Daarem,
}

impl From<EmAccelArg> for EmAccel {
    fn from(a: EmAccelArg) -> Self {
        match a {
            EmAccelArg::None => EmAccel::None,
            EmAccelArg::Squarem => EmAccel::Squarem,
            EmAccelArg::Daarem => EmAccel::Daarem,
        }
    }
}

/// Defaults for the fragment-length prior flags. Held here (rather than in
/// `default_value_t`) so the flags can stay `Option` and salmon can distinguish
/// a supplied value from an inherited one — the distinction issue #1062 needs in
/// order to warn only when a value the user actually chose is being ignored.
const DEFAULT_FLD_MEAN: f64 = 250.0;
const DEFAULT_FLD_SD: f64 = 25.0;
const DEFAULT_FLD_MAX: usize = 1000;

/// CLI surface for [`FldPolicy`]: where a RAD requant's fragment-length
/// distribution comes from.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
enum FldPolicyArg {
    /// Prefer the distribution baked into the RAD header (exact parity with the
    /// run that wrote it); otherwise derive it, or fall back to the prior.
    #[default]
    Baked,
    /// Ignore any baked distribution; re-derive it from this RAD's own
    /// uniquely-mapped proper pairs.
    Derive,
    /// Ignore both the baked distribution and the RAD's fragment lengths; use
    /// `--fldMean`/`--fldSD` alone.
    Prior,
}

impl From<FldPolicyArg> for FldPolicy {
    fn from(p: FldPolicyArg) -> Self {
        match p {
            FldPolicyArg::Baked => FldPolicy::Baked,
            FldPolicyArg::Derive => FldPolicy::Derive,
            FldPolicyArg::Prior => FldPolicy::Prior,
        }
    }
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
    /// Single-cell quantification: removed; use the alevin-fry ecosystem.
    /// Hidden from help; retained only to redirect legacy `salmon alevin` calls.
    #[command(hide = true, disable_help_flag = true)]
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
    /// Seed with sparse fixed-k k-mer anchors instead of the default
    /// unitig-constrained uni-MEMs.
    #[arg(long = "sparseSeeds", conflicts_with_all = ["unimems", "refmems"])]
    sparse_seeds: bool,
    /// Seed with reference-extended MEMs (cross unitig boundaries).
    #[arg(long = "refMEMs", conflicts_with = "unimems")]
    refmems: bool,
    /// Seed with true unitig-constrained uni-MEMs (pufferfish-style; default).
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
    /// Minimizer length (0 = auto; 19 for the default k=31).
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
    /// Directory for sshash's external minimizer-sort scratch. Defaults to a
    /// `sshash_tmp` subdirectory of --tmpdir (or the index dir when --tmpdir is
    /// unset). Override to place the sort scratch on a separate/fast disk; the
    /// directory is created before the build and removed afterwards.
    #[arg(long = "sshashTmpDir")]
    sshash_tmp_dir: Option<PathBuf>,
    /// RAM ceiling, in GiB, for sshash's external minimizer sort (the main
    /// build-time memory/disk trade-off). Default is 8; a smaller value uses
    /// less RAM but spills to disk sooner.
    #[arg(long = "ramLimit")]
    ram_limit_gib: Option<usize>,
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
    #[arg(
        short = 'i',
        long = "index",
        required_unless_present_any = ["alignments", "rad"]
    )]
    index: Option<PathBuf>,
    /// Alignment-based mode: a BAM of reads aligned to the transcriptome.
    #[arg(short = 'a', long = "alignments")]
    alignments: Option<PathBuf>,
    /// RAD-input mode: a RAD file of mappings (from `salmon quant --writeRad` or
    /// `piscem map-bulk`) to quantify directly, in parallel.
    #[arg(long = "rad", conflicts_with_all = ["alignments", "mates1", "mates2", "unmated"])]
    rad: Option<PathBuf>,
    /// Transcriptome FASTA (alignment mode `-a`): enables the alignment error model.
    #[arg(short = 't', long = "targets")]
    targets: Option<PathBuf>,
    /// Disable the alignment error model (alignment mode).
    #[arg(long = "noErrorModel")]
    no_error_model: bool,
    /// Opt in to the order-independent alignment error model in deterministic
    /// mode (`-a --deterministic`, requires `-t`). Off by default: deterministic
    /// mode scores by the BAM `AS` tag, which benchmarks at least as accurately
    /// against ground truth and needs only a single BAM pass; `--errorModel` adds
    /// a second pass to train the model. Use it for parity with salmon's classic
    /// error-modeled quant, or when the aligner's `AS` is absent/unreliable.
    #[arg(long = "errorModel")]
    error_model: bool,
    /// Genome-alignment mode: a GTF/GFF annotation. Given with `-a <genome.bam>`,
    /// salmon projects the spliced genome alignments into transcriptome
    /// coordinates (via bramble) and quantifies the projected placements. The
    /// transcript set is taken from the annotation.
    #[arg(long = "annotation", value_name = "GTF|GFF")]
    annotation: Option<PathBuf>,
    /// Genome FASTA (genome-alignment mode) — needed only for bias correction, to
    /// reconstruct transcript sequences from the annotation's exons.
    #[arg(long = "genome", value_name = "FASTA")]
    genome: Option<PathBuf>,
    /// Genome-alignment projection: penalty multiplier per junction mismatch
    /// (bramble `junc_miss_discount`; default 1.0 = no penalty).
    #[arg(long = "juncMissDiscount", value_name = "F")]
    junc_miss_discount: Option<f64>,
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
    /// Transcript-to-gene map (GTF/GFF, or a 2-column TSV), optionally
    /// compressed (gzip, BGZF, bzip2, xz, zstd); also writes gene-level
    /// estimates to `quant.genes.sf`.
    #[arg(short = 'g', long = "geneMap", value_name = "FILE")]
    gene_map: Option<PathBuf>,
    /// Ignore trailing transcript-version suffixes when joining `--geneMap`
    /// entries to quantified transcripts, so `ENST00000456328.2` matches
    /// `ENST00000456328`. Matches tximport's option of the same name. Needed
    /// when the index was built from an Ensembl cDNA FASTA, which carries the
    /// version in the identifier, and the annotation is an Ensembl GTF, which
    /// puts it in a separate `transcript_version` attribute.
    #[arg(long = "ignoreTxVersion")]
    ignore_tx_version: bool,
    /// Worker threads (0 = all cores). salmon also runs one auxiliary thread for
    /// FASTQ parsing/decompression, so total CPU use can slightly exceed this
    /// value; on schedulers that enforce the limit, request one extra core.
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
    /// standard `--writeMappings`). `--writeSam` is an accepted alias, naming
    /// the format explicitly to pair with `--writeBam`.
    #[arg(
        short = 'z',
        long = "writeMappings",
        visible_alias = "writeSam",
        conflicts_with = "write_bam"
    )]
    write_mappings: Option<PathBuf>,
    /// Write per-mapping BAM records to this file. Mutually exclusive with
    /// --writeMappings/--writeSam.
    #[arg(long = "writeBam", conflicts_with = "write_mappings")]
    write_bam: Option<PathBuf>,
    /// BGZF compression workers for --writeBam [default: about one per 3
    /// mapping threads, at most 8].
    ///
    /// These threads run alongside the -p mapping threads and do nothing but
    /// compress BGZF blocks. The default balances the two stages from measured
    /// throughput: one worker compresses ~165 MiB/s of BAM records, while one
    /// mapping thread produces at most ~52 MiB/s of them, so roughly one worker
    /// per 3 mapping threads keeps compression from becoming the bottleneck.
    /// The result is capped at 8, which is already about 3x the fastest record
    /// production measured at any thread count.
    ///
    /// The cost is very lopsided, so the default errs upward: one worker too
    /// few can halve throughput, whereas surplus workers just block on an empty
    /// queue and cost nothing measurable. It stops at the balance point rather
    /// than going higher, because past that these cores do more good mapping
    /// reads. Raise it if you are writing BAM to unusually fast storage and can
    /// spare the cores; there is no reason to lower it.
    // `requires` alone would let `--writeSam out.sam --bamCompressThreads 8`
    // through and then quietly ignore it, since SAM has no compression stage.
    #[arg(
        long = "bamCompressThreads",
        requires = "write_bam",
        conflicts_with = "write_mappings"
    )]
    bam_compress_threads: Option<std::num::NonZeroUsize>,
    /// Write per-fragment mappings to a RAD file at this path. Sketch or
    /// selective-alignment profile is chosen automatically from the mapping
    /// mode. Quantification still runs; add --skipQuant to map only. The file is
    /// piscem map-bulk compatible and can be re-quantified with `--rad`.
    #[arg(long = "writeRad")]
    write_rad: Option<PathBuf>,
    /// Deterministic quantification (reads/FASTQ input): map once to an
    /// intermediate RAD, then quantify from it with a fixed fragment-length
    /// distribution, so the result is byte-identical across runs and thread
    /// counts. Avoids a second mapping pass. The intermediate RAD is written under
    /// the output directory and deleted on success unless --keepRad (or use
    /// --writeRad PATH to choose its location and keep it).
    #[arg(long = "deterministic")]
    deterministic: bool,
    /// Keep the intermediate RAD produced by --deterministic (by default it is
    /// deleted once quantification finishes). Ignored without --deterministic.
    #[arg(long = "keepRad")]
    keep_rad: bool,
    /// Enable sequence-specific bias correction.
    #[arg(long = "seqBias")]
    seq_bias: bool,
    /// Enable fragment-GC bias correction.
    #[arg(long = "gcBias")]
    gc_bias: bool,
    /// Enable positional bias correction.
    #[arg(long = "posBias")]
    pos_bias: bool,
    /// EM iterations for the preliminary abundance estimate that weights bias
    /// model collection (expert; rough on purpose to avoid overfitting).
    #[arg(long = "biasSeedEMIters", default_value_t = 11, hide = true)]
    bias_seed_em_iters: u32,
    /// Dump observed+expected seq/GC/positional bias models to
    /// `bias_models.txt` (debugging / parity).
    #[arg(long = "dumpBiasModels", hide = true)]
    dump_bias_models: bool,
    /// Compression codec for RAD output (`--writeRad` and the `--deterministic`
    /// intermediate): `lz4` (default), `zstd` (better ratio), or `none`.
    #[arg(long = "radCompress", default_value = "lz4", value_parser = ["lz4", "zstd", "none"])]
    rad_compress: String,
    /// Write uncompressed RAD chunks (overrides --radCompress).
    #[arg(long = "noCompressRad")]
    no_compress_rad: bool,
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
    /// Seed with sparse fixed-k k-mer anchors instead of the default
    /// unitig-constrained uni-MEMs.
    #[arg(long = "sparseSeeds", conflicts_with_all = ["unimems", "refmems"])]
    sparse_seeds: bool,
    /// Seed with reference-extended MEMs (cross unitig boundaries).
    #[arg(long = "refMEMs", conflicts_with = "unimems")]
    refmems: bool,
    /// Seed with true unitig-constrained uni-MEMs (pufferfish-style). This is the
    /// default; the flag is accepted for explicitness/back-compat.
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
    /// EM/VBEM convergence acceleration. `squarem` reaches the same abundances in
    /// far fewer M-steps but is not byte-identical to the default `none`.
    #[arg(long = "emAccel", value_enum, default_value_t = EmAccelArg::None)]
    em_accel: EmAccelArg,
    /// Mean of the fragment-length distribution prior [default: 250].
    ///
    /// A prior, not an override: wherever fragment lengths are observed, this
    /// seeds the distribution with weight 1 and the observations dominate. With
    /// `--rad` it is ignored entirely unless `--fldPolicy` says otherwise, since
    /// a salmon RAD carries a baked distribution.
    // `Option` rather than `default_value_t` so salmon can tell a supplied value
    // from an inherited one, and warn only when a supplied one is ignored.
    #[arg(long = "fldMean")]
    fld_mean: Option<f64>,
    /// Standard deviation of the fragment-length distribution prior [default: 25].
    ///
    /// Same prior semantics as `--fldMean`.
    #[arg(long = "fldSD")]
    fld_sd: Option<f64>,
    /// Maximum fragment length tracked by the fragment-length distribution
    /// [default: 1000].
    #[arg(long = "fldMax")]
    fld_max: Option<usize>,
    /// Where the fragment-length distribution comes from when quantifying a RAD
    /// (`--rad` only) [default: baked].
    ///
    /// salmon's RAD writer always bakes its fragment-length distribution into
    /// the header, so by default a requant reproduces the writing run exactly
    /// and `--fldMean`/`--fldSD` have no effect. Use `derive` to re-derive the
    /// distribution from the RAD's own fragment lengths, or `prior` to use
    /// `--fldMean`/`--fldSD` alone — `prior` is the setting that makes a
    /// fragment-length sensitivity analysis actually vary anything.
    #[arg(long = "fldPolicy", value_enum, default_value_t = FldPolicyArg::Baked)]
    fld_policy: FldPolicyArg,
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
    /// Discard a read/fragment that maps to more than this many places — counting
    /// the total set of distinct mappings (concordant + orphan union) (reads mode).
    #[arg(short = 'w', long = "maxReadOcc", default_value_t = 250)]
    max_read_occ: usize,
    /// Only emit a single-mate (orphan) mapping when the read's mate is entirely
    /// unmapped. By default, when a pair has no concordant mapping, orphans are
    /// reported for both mates (their union); with this flag a read that maps only
    /// because its mate also mapped (to a disjoint set) is not reported as an orphan.
    #[arg(long = "orphansRequireUnmappedMate")]
    orphans_require_unmapped_mate: bool,
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
    /// When a fragment's best transcript alignment is outscored by a decoy
    /// (genome) alignment, keep the transcript placement instead of discarding the
    /// fragment. Off by default (decoy-dominated fragments are dropped, matching
    /// salmon). Increases the reported mapping rate on decoy-aware indices at the
    /// cost of retaining transcript mappings for fragments better explained by the
    /// genome; only affects fragments that retain a transcript candidate.
    #[arg(long = "allowDecoyOrphans")]
    allow_decoy_orphans: bool,
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

/// Resolve the `--sparseSeeds` / `--uniMEMs` / `--refMEMs` flags into a seeding
/// mode (clap enforces they are mutually exclusive). The default is uni-MEM
/// seeding (`--uniMEMs`): it is faster than sparse k-mer seeding and at least as
/// accurate, and faithfully reproduces pufferfish's seeding. `--sparseSeeds`
/// selects the older sparse fixed-k anchors.
// Resolve the three mutually-exclusive seeding flags into one mode. Kept as a
// function so the precedence between them is stated once rather than inline.
fn seed_mode(unimems: bool, refmems: bool, sparse_seeds: bool) -> salmon_map::SeedMode {
    if sparse_seeds {
        salmon_map::SeedMode::Sparse
    } else if refmems {
        salmon_map::SeedMode::RefMem
    } else {
        // default (and explicit `--uniMEMs`)
        let _ = unimems;
        salmon_map::SeedMode::UniMem
    }
}

// `salmon index`: translate the arguments into `IndexBuildOptions` and build.
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
    opts.sshash_tmp = args.sshash_tmp_dir;
    opts.ram_limit_gib = args.ram_limit_gib;
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
// Long-read quantification is out of scope, as it is upstream. Rather than
// failing obscurely on data salmon cannot model, point the user at the right
// tool. `-> !` means this never returns: it exits the process.
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

/// Where the `--geneMap` annotation is and how to join it to quantified
/// transcripts. Carried together so the option travels with the path through
/// every quantification path.
#[derive(Clone, Debug)]
struct GeneMapOpts {
    path: PathBuf,
    /// `--ignoreTxVersion`
    ignore_tx_version: bool,
}

/// Name of the file under `aux_info/` listing transcripts the gene map did not
/// cover. It sits beside `unmapped_names.txt` — which is about *fragments* that
/// did not map — so the name says transcripts explicitly.
const UNMATCHED_TXP_FILE: &str = "genemap_unmatched_txps.json";

/// Aggregate transcript estimates to gene level and write `quant.genes.sf`.
fn write_gene_level(
    out_dir: &std::path::Path,
    gene_map: &GeneMapOpts,
    names: &[String],
    lengths: &[u32],
    eff_lengths: &[f64],
    tpm: &[f64],
    counts: &[f64],
) -> Result<()> {
    let map = salmon_core::genemap::read_transcript_gene_map(&gene_map.path)
        .with_context(|| format!("reading gene map {}", gene_map.path.display()))?;
    // Listed by name rather than only counted, so the failure survives a lost
    // terminal and can be fed straight back into whatever built the annotation.
    let unmatched_list = out_dir.join("aux_info").join(UNMATCHED_TXP_FILE);
    let summary = salmon_core::genemap::write_gene_quant(
        &out_dir.join("quant.genes.sf"),
        names,
        lengths,
        eff_lengths,
        tpm,
        counts,
        &map,
        gene_map.ignore_tx_version,
        Some(&unmatched_list),
    )
    .context("writing quant.genes.sf")?;

    if summary.unmatched > 0 {
        eprintln!(
            "warning: {} of {} transcript(s) had no entry in the gene map. Each was written to \
             quant.genes.sf as its own single-transcript gene, so those rows are named by \
             transcript, not by gene. Their names are listed in {}",
            summary.unmatched,
            summary.total(),
            unmatched_list.display()
        );
    }
    // A near-total mismatch means the identifiers do not line up, not that the
    // annotation is merely incomplete. With the fallback in place the file is
    // full-looking either way, so the row breakdown is the only thing that
    // distinguishes a real gene-level result from a transcript list wearing its
    // name — state it rather than leave it to be inferred.
    if summary.match_rate_is_low() {
        eprintln!(
            "warning: the gene map matched only {} of {} quantified transcripts ({:.1}%), so of \
             the {} rows in quant.genes.sf, {} are unmatched transcripts standing in as their \
             own gene. Treat this output as transcript-level until the join is fixed.",
            summary.matched,
            summary.total(),
            summary.match_rate() * 100.0,
            summary.genes_written,
            summary.unmatched
        );
        let would_match = salmon_core::genemap::matches_ignoring_tx_version(names, &map);
        if !gene_map.ignore_tx_version && would_match > summary.matched {
            eprintln!(
                "warning: {would_match} would match with --ignoreTxVersion, which compares \
                 identifiers without their trailing version suffix (e.g. {} vs {}). The Ensembl \
                 cDNA FASTA carries transcript versions in the identifier while the GTF puts \
                 them in a separate transcript_version attribute, so an index built from the \
                 FASTA does not join to the GTF as-is.",
                names
                    .iter()
                    .find(|n| !map.contains_key(*n))
                    .map(String::as_str)
                    .unwrap_or("ENST00000456328.2"),
                map.keys()
                    .next()
                    .map(String::as_str)
                    .unwrap_or("ENST00000456328"),
            );
        }
    }
    // Count the rows written, not the genes in the gene map: the two differ
    // whenever the join is imperfect, and reporting the latter is what let a
    // header-only quant.genes.sf be announced as a success (#1075).
    //
    // With the fallback in place the same trap reappears one step along, since
    // a stand-in row is a transcript wearing a gene's place in the file. So the
    // count is only called "genes" when every row is one; otherwise the split is
    // spelled out, and the last line of the run agrees with the warnings above
    // it rather than contradicting them.
    if summary.unmatched == 0 {
        println!(
            "wrote gene-level estimates for {} genes to quant.genes.sf",
            summary.genes_written
        );
    } else {
        println!(
            "wrote {} rows to quant.genes.sf: {} gene(s), plus {} unmatched transcript(s) as \
             their own gene",
            summary.genes_written,
            // saturating: a stand-in row whose name collides with a real gene
            // merges into it, so the two counts can overlap by one per collision.
            summary.genes_written.saturating_sub(summary.unmatched),
            summary.unmatched
        );
    }
    Ok(())
}

/// Deterministic reads-mode quantification: map the FASTQ once to an
/// intermediate RAD (baking a fixed, order-independent FLD and the resolved
/// library format), then quantify from that RAD via [`quantify_rad`]. Because the
/// FLD is fixed before equivalence-class assembly, the result is byte-identical
/// across runs and thread counts. Re-reading the RAD is cheap, so this avoids a
/// second *mapping* pass; the requant derives its own abundances from the RAD.
/// Bias correction needs the reference sequence, but since we control both passes
/// and phase 1 loads the index, we hand the index's sequences to phase 2 — so
/// `--deterministic --seqBias/--gcBias/--posBias` needs no separate `-t`.
// `--deterministic` (reads mode): a two-pass run. The first pass maps to a RAD
// with order-independent models baked into the header; the second quantifies
// from it. Splitting the run this way is what makes the result independent of
// the thread count.
fn run_deterministic(
    mut map_opts: QuantOptions,
    rad_out: Option<PathBuf>,
    keep_rad: bool,
    bias_targets: Option<PathBuf>,
    gene_map: Option<&GeneMapOpts>,
    out_dir: &std::path::Path,
) -> Result<()> {
    let bias_on = map_opts.seq_bias || map_opts.gc_bias || map_opts.pos_bias;
    // An explicit `--writeRad PATH` is honoured (and kept — it was a requested
    // output); otherwise a temp under the output directory, removed on success
    // unless `--keepRad`.
    let explicit = rad_out.is_some();
    let rad_path = rad_out.unwrap_or_else(|| out_dir.join("intermediate_mappings.rad"));
    // The RAD writer opens its file before the mapping pass, so the output
    // directory (where the default intermediate lives) must exist first.
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("creating output directory {}", out_dir.display()))?;

    // Phase 1 — map once, write the RAD (bakes the deterministic FLD + resolved
    // library format), skipping the online EM: quantification happens in phase 2.
    map_opts.write_rad = Some(rad_path.clone());
    map_opts.skip_quant = true;
    // Derive the FLD + library format order-independently during this pass and
    // bake them (see `QuantOptions::deterministic_fld`), so phase 2 is a single
    // pass and the whole result is byte-identical across thread counts.
    map_opts.deterministic_fld = true;
    let map_res = quantify(&map_opts).context("deterministic mapping pass failed")?;
    let pct = if map_res.num_processed > 0 {
        100.0 * map_res.num_mapped as f64 / map_res.num_processed as f64
    } else {
        0.0
    };
    tracing::info!(
        "deterministic: mapped {} / {} fragments ({:.2}%); quantifying from intermediate RAD",
        map_res.num_mapped,
        map_res.num_processed,
        pct
    );

    // Phase 2 — deterministic quant from the RAD (fixed baked FLD + library
    // format ⇒ order-independent eq-classes + EM). Mirrors the `--rad` knob wiring.
    let mut q = AlignQuantOptions::new(rad_path.clone(), out_dir.to_path_buf());
    q.lib_type = map_opts.lib_type.clone();
    q.em = map_opts.em.clone();
    q.range_factorization_bins = map_opts.range_factorization_bins;
    // Bias needs the reference sequence. Prefer the index's own sequences (we
    // already built the index for phase 1), so the user need not pass `-t`; an
    // explicit `-t` is still honoured as a fallback.
    q.transcripts = bias_targets;
    if bias_on {
        q.ref_seqs = Some(
            salmon_index::load_ref_seqs(&map_opts.index_dir)
                .context("loading reference sequences from the index for deterministic bias")?,
        );
    }
    q.seq_bias = map_opts.seq_bias;
    q.gc_bias = map_opts.gc_bias;
    q.pos_bias = map_opts.pos_bias;
    q.bias_seed_em_iters = map_opts.bias_seed_em_iters;
    q.score_exp = map_opts.map_config.score.score_exp;
    q.incompat_prior = map_opts.incompat_prior;
    q.sig_digits = map_opts.sig_digits;
    q.fld_mean = map_opts.fld_mean;
    q.fld_sd = map_opts.fld_sd;
    q.fld_max = map_opts.fld_max;
    q.gc_bins = map_opts.gc_bins;
    q.cond_gc_bins = map_opts.cond_gc_bins;
    q.num_bootstraps = map_opts.num_bootstraps;
    q.num_gibbs_samples = map_opts.num_gibbs_samples;
    q.thinning_factor = map_opts.thinning_factor;
    let res = quantify_rad(&q, &rad_path).context("deterministic requant pass failed")?;
    let pct2 = if res.num_processed > 0 {
        100.0 * res.num_mapped as f64 / res.num_processed as f64
    } else {
        0.0
    };
    tracing::info!(
        "{} fragments from intermediate RAD, {} quantified ({:.2}%); {} equivalence classes",
        res.num_processed,
        res.num_mapped,
        pct2,
        res.num_eq_classes
    );
    if let Some(gm) = gene_map {
        write_gene_level(
            out_dir,
            gm,
            &res.names,
            &res.lengths,
            &res.eff_lengths,
            &res.tpm,
            &res.counts,
        )?;
    }
    // Remove the intermediate unless kept (an explicit --writeRad path is always
    // kept).
    if explicit || keep_rad {
        tracing::info!("kept intermediate RAD at {}", rad_path.display());
    } else if let Err(e) = std::fs::remove_file(&rad_path) {
        tracing::warn!(
            "could not remove intermediate RAD {}: {e}",
            rad_path.display()
        );
    }
    Ok(())
}

/// Alignment-mode `--deterministic`: write the transcriptome BAM's placements to
/// an intermediate salmon RAD (one pass — baking an order-independent FLD +
/// library format, and rough abundances when bias is requested), then quantify
/// from it via the same deterministic [`quantify_rad`] the reads path uses, so
/// the result is byte-identical across thread counts. Score-based (carries the
/// BAM `AS` tag); the online alignment error model is not used in this mode.
// The same two-pass split for alignment input: BAM to RAD, then quantify.
fn run_deterministic_align(
    opts: AlignQuantOptions,
    rad_out: Option<PathBuf>,
    keep_rad: bool,
    gene_map: Option<&GeneMapOpts>,
    codec: ChunkCodec,
) -> Result<()> {
    let out_dir = opts.output_dir.clone();
    // An explicit `--writeRad PATH` is honoured and kept; otherwise a temp under
    // the output dir, removed on success unless `--keepRad` (or `--skipQuant`,
    // where the RAD is the deliverable).
    let explicit = rad_out.is_some();
    let rad_path = rad_out.unwrap_or_else(|| out_dir.join("intermediate_alignments.rad"));
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("creating output directory {}", out_dir.display()))?;

    // Phase 1 — one BAM pass: write placements + bake FLD / library format
    // (+ rough abundances when bias is on) so phase 2 is a single baked pass.
    let summary = salmon_align::write_alignment_rad(&opts, &rad_path, codec)
        .context("deterministic alignment RAD-write pass failed")?;
    let pct = if summary.num_processed > 0 {
        100.0 * summary.num_mapped as f64 / summary.num_processed as f64
    } else {
        0.0
    };
    tracing::info!(
        "deterministic: wrote {} / {} aligned fragments ({:.2}%) to the intermediate RAD; quantifying",
        summary.num_mapped,
        summary.num_processed,
        pct
    );

    // Phase 2 — quantify from the fully-baked RAD (single pass, order-independent).
    let res = quantify_rad(&opts, &rad_path).context("deterministic requant pass failed")?;
    let pct2 = if res.num_processed > 0 {
        100.0 * res.num_mapped as f64 / res.num_processed as f64
    } else {
        0.0
    };
    tracing::info!(
        "{} fragments from intermediate RAD, {} quantified ({:.2}%); {} equivalence classes",
        res.num_processed,
        res.num_mapped,
        pct2,
        res.num_eq_classes
    );
    if let Some(gm) = gene_map {
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
    if explicit || keep_rad || opts.skip_quant {
        tracing::info!("kept intermediate RAD at {}", rad_path.display());
    } else if let Err(e) = std::fs::remove_file(&rad_path) {
        tracing::warn!(
            "could not remove intermediate RAD {}: {e}",
            rad_path.display()
        );
    }
    Ok(())
}

/// Resolve the RAD chunk-compression codec from the CLI flags.
fn rad_codec_from_args(no_compress_rad: bool, rad_compress: &str) -> Result<ChunkCodec> {
    if no_compress_rad {
        return Ok(ChunkCodec::None);
    }
    Ok(match rad_compress {
        "none" => ChunkCodec::None,
        "lz4" => ChunkCodec::Lz4,
        "zstd" => ChunkCodec::Zstd,
        other => anyhow::bail!("unknown --radCompress codec '{other}' (expected lz4|zstd|none)"),
    })
}

/// Genome-alignment mode (`-a <genome.bam> --annotation <gtf>`): project the
/// spliced genome alignments into transcriptome coordinates with bramble, write
/// a salmon RAD, then quantify from it with the deterministic `quantify_rad`.
/// `q` supplies the EM / bias / output / bootstrap knobs for that requant.
// `--annotation`: project a genome BAM into transcript coordinates, then
// quantify the resulting RAD through the ordinary path.
fn run_genome_project(
    gp: GenomeProjectOptions,
    q: AlignQuantOptions,
    rad_out: Option<PathBuf>,
    keep_rad: bool,
    gene_map: Option<&GeneMapOpts>,
) -> Result<()> {
    let out_dir = q.output_dir.clone();
    let explicit = rad_out.is_some();
    let rad_path = rad_out.unwrap_or_else(|| out_dir.join("genome_projected.rad"));
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("creating output directory {}", out_dir.display()))?;

    // Phase 1 — project the genome BAM into a fully-baked transcriptome RAD.
    let art = project_genome_bam_to_rad(&gp, &rad_path).context("genome projection pass failed")?;
    let pct = if art.num_processed > 0 {
        100.0 * art.num_mapped as f64 / art.num_processed as f64
    } else {
        0.0
    };
    tracing::info!(
        "projected {} / {} read groups ({:.2}%) onto {} annotated transcripts; quantifying",
        art.num_mapped,
        art.num_processed,
        pct,
        art.names.len()
    );

    // Phase 2 — quantify from the projected RAD (single baked pass). Hand the
    // reconstructed transcript sequences to bias correction (no `-t` needed).
    let mut q = q;
    if q.ref_seqs.is_none() {
        q.ref_seqs = art.ref_seqs;
    }
    let res = quantify_rad(&q, &rad_path).context("genome-projected quantification failed")?;
    tracing::info!(
        "{} fragments from projected RAD, {} quantified; {} equivalence classes",
        res.num_processed,
        res.num_mapped,
        res.num_eq_classes
    );
    if let Some(gm) = gene_map {
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
    if explicit || keep_rad {
        tracing::info!("kept projected RAD at {}", rad_path.display());
    } else if let Err(e) = std::fs::remove_file(&rad_path) {
        tracing::warn!("could not remove projected RAD {}: {e}", rad_path.display());
    }
    Ok(())
}

// `salmon quant`: the main entry point. Most of its length is deciding which of
// the several input modes the flags describe (reads, BAM, genome BAM, RAD; plain
// or deterministic), validating that the combination makes sense, and reporting
// clearly when it does not — a wrong mode chosen silently would be far worse
// than an error.
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
    if args.sketch && args.decoy_threshold != 1.0 {
        tracing::warn!(
            "--decoyThreshold {} has no effect in --sketch mode: sketch (pseudoalignment) \
             returns only equally-best mappings, so the decoy-domination comparison \
             (bestTxpScore < decoyThreshold * bestDecoyScore) never triggers. A fragment is \
             treated as decoy-dominated only when it maps to decoys *and no transcript*; \
             otherwise decoy hits are dropped and transcript hits kept (use --allowDecoyOrphans \
             to keep transcript hits even when a decoy also matches).",
            args.decoy_threshold
        );
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
    let gene_map = args.gene_map.clone().map(|path| GeneMapOpts {
        path,
        ignore_tx_version: args.ignore_tx_version,
    });
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
        warn_unsupported_mapping_output(
            &args.write_mappings,
            &args.write_bam,
            "alignment mode (-a)",
            "the input alignments already are the per-mapping records",
        );
        warn_unsupported_fld_policy(args.fld_policy, "alignment mode (-a)");
        let mut opts = AlignQuantOptions::new(bam, args.output);
        opts.lib_type = args.lib_type;
        opts.em.use_vbem = use_vbem;
        opts.range_factorization_bins = range_factorization_bins;
        opts.transcripts = args.targets;
        opts.no_error_model = args.no_error_model;
        opts.deterministic_error_model = args.error_model;
        opts.seq_bias = args.seq_bias;
        opts.gc_bias = args.gc_bias;
        opts.pos_bias = args.pos_bias;
        opts.bias_seed_em_iters = args.bias_seed_em_iters;
        opts.incompat_prior = args.incompat_prior;
        opts.em.vb_prior = args.vb_prior;
        opts.em.per_nucleotide_prior = args.per_nucleotide_prior;
        opts.em.accel = args.em_accel.into();
        opts.sig_digits = args.sig_digits;
        opts.fld_mean = args.fld_mean.unwrap_or(DEFAULT_FLD_MEAN);
        opts.fld_sd = args.fld_sd.unwrap_or(DEFAULT_FLD_SD);
        opts.fld_max = args.fld_max.unwrap_or(DEFAULT_FLD_MAX);
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
        opts.num_bootstraps = args.num_bootstraps;
        opts.num_gibbs_samples = args.num_gibbs_samples;
        opts.thinning_factor = args.thinning_factor;
        // --scoreExp is selective-alignment-mode only (it scales the
        // best-minus-score soft weight); alignment mode has no such term.

        // Genome-alignment mode: `-a <genome.bam> --annotation <gtf>` projects the
        // spliced genome alignments into transcriptome coordinates (bramble) and
        // quantifies the projection. Inherently RAD-based / deterministic.
        if let Some(annotation) = args.annotation.clone() {
            // Genome projection is inherently deterministic (RAD-based) and has no
            // error model (bramble exposes no projected CIGAR), so flags that only
            // matter for transcriptomic alignment are accepted but no-ops here.
            if args.deterministic {
                tracing::warn!(
                    "--deterministic is redundant in genome-projection mode (--annotation): \
                     projection is already deterministic; ignoring it"
                );
            }
            if args.error_model {
                tracing::warn!(
                    "--errorModel has no effect in genome-projection mode (--annotation): \
                     projection scores by bramble similarity, not an alignment error model"
                );
            }
            let codec = rad_codec_from_args(args.no_compress_rad, &args.rad_compress)?;
            let gp = GenomeProjectOptions {
                bam: opts.bam.clone(),
                annotation,
                genome_fasta: args.genome.clone(),
                lib_type: opts.lib_type.clone(),
                junc_miss_discount: args.junc_miss_discount,
                bias: opts.seq_bias || opts.gc_bias || opts.pos_bias,
                fld_mean: opts.fld_mean,
                fld_sd: opts.fld_sd,
                fld_max: opts.fld_max,
                bias_seed_em_iters: opts.bias_seed_em_iters,
                em: opts.em.clone(),
                rad_codec: codec,
            };
            return run_genome_project(
                gp,
                opts,
                args.write_rad.clone(),
                args.keep_rad,
                gene_map.as_ref(),
            );
        }

        // `--deterministic`: write the BAM's placements to an intermediate RAD
        // (baking the FLD/library-format + rough abundances) and requantify from
        // it, for results byte-identical across thread counts. Score-based; no
        // online error model.
        if args.deterministic {
            let codec = rad_codec_from_args(args.no_compress_rad, &args.rad_compress)?;
            return run_deterministic_align(
                opts,
                args.write_rad.clone(),
                args.keep_rad,
                gene_map.as_ref(),
                codec,
            );
        }

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
        // Alignment mode: the records are already aligned (by the upstream
        // aligner), so `num_processed` is the count of *aligned* fragments and
        // `num_mapped` is those with a strand-compatible placement we could
        // quantify. Report it that way so a (necessarily) <100% rate on a
        // stranded library reads as "strand-incompatible", not "lost alignments".
        if res.num_processed > res.num_mapped {
            tracing::info!(
                "{} aligned fragments; {} strand-compatible and quantified ({:.2}%); \
                 {} dropped as incompatible with library type {}; {} equivalence classes",
                res.num_processed,
                res.num_mapped,
                pct,
                res.num_processed - res.num_mapped,
                opts.lib_type,
                res.num_eq_classes
            );
        } else {
            tracing::info!(
                "{} aligned fragments, all strand-compatible and quantified; {} equivalence classes",
                res.num_processed,
                res.num_eq_classes
            );
        }
        if res.inference_truncated_mass > 0.0 {
            tracing::warn!(
                "{:.3} fragments of equivalence-class mass could not be assigned (every member \
                 transcript was truncated below the min-alpha threshold); reported as \
                 inference_truncated_mass in meta_info.json",
                res.inference_truncated_mass
            );
        }
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

    // RAD-input mode: quantify directly from a RAD file of mappings, in parallel.
    if let Some(rad_path) = args.rad {
        warn_unsupported_mapping_output(
            &args.write_mappings,
            &args.write_bam,
            "RAD-input mode (--rad)",
            "a RAD carries no read names or sequences, so SAM/BAM records cannot be reconstructed",
        );
        // The RAD reuses alignment mode's inference/output; AlignQuantOptions
        // carries the relevant knobs (lib type, EM/VBEM, score-exp, FLD prior,
        // bootstrap/Gibbs). `bam` is unused by the RAD path.
        let mut opts = AlignQuantOptions::new(rad_path.clone(), args.output);
        opts.lib_type = args.lib_type;
        opts.em.use_vbem = use_vbem;
        opts.range_factorization_bins = range_factorization_bins;
        // Bias correction from RAD: seq/GC use the reference at each fragment's
        // recorded position (the transcriptome `-t`), positional uses only the RAD.
        opts.transcripts = args.targets;
        opts.seq_bias = args.seq_bias;
        opts.gc_bias = args.gc_bias;
        opts.pos_bias = args.pos_bias;
        opts.bias_seed_em_iters = args.bias_seed_em_iters;
        opts.score_exp = args.score_exp;
        opts.incompat_prior = args.incompat_prior;
        opts.em.vb_prior = args.vb_prior;
        opts.em.per_nucleotide_prior = args.per_nucleotide_prior;
        opts.em.accel = args.em_accel.into();
        opts.sig_digits = args.sig_digits;
        opts.fld_mean = args.fld_mean.unwrap_or(DEFAULT_FLD_MEAN);
        opts.fld_sd = args.fld_sd.unwrap_or(DEFAULT_FLD_SD);
        opts.fld_max = args.fld_max.unwrap_or(DEFAULT_FLD_MAX);
        opts.fld_policy = args.fld_policy.into();
        // Only what the user actually typed, so a plain requant is not warned
        // about flags it inherited (see `warn_baked_fld_overrides`).
        opts.explicit_fld_args = ExplicitFldArgs {
            mean: args.fld_mean.is_some(),
            sd: args.fld_sd.is_some(),
            max: args.fld_max.is_some(),
        };
        opts.skip_quant = args.skip_quant;
        opts.num_bootstraps = args.num_bootstraps;
        opts.num_gibbs_samples = args.num_gibbs_samples;
        opts.thinning_factor = args.thinning_factor;
        let res = quantify_rad(&opts, &rad_path).context("RAD-input quantification failed")?;
        let pct = if res.num_processed > 0 {
            100.0 * res.num_mapped as f64 / res.num_processed as f64
        } else {
            0.0
        };
        tracing::info!(
            "{} fragments from RAD, {} quantified ({:.2}%); {} equivalence classes",
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
        "no reads provided: pass -1/-2 (paired), -r (single-end), -a (BAM), or --rad (RAD)"
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
    opts.write_bam = args.write_bam;
    opts.bam_compress_threads = args.bam_compress_threads;
    opts.write_rad = args.write_rad;
    opts.seq_bias = args.seq_bias;
    opts.gc_bias = args.gc_bias;
    opts.pos_bias = args.pos_bias;
    opts.bias_seed_em_iters = args.bias_seed_em_iters;
    opts.dump_bias_models = args.dump_bias_models;
    // RAD chunk compression (only consumed when a RAD is actually written:
    // --writeRad output or the --deterministic intermediate). --noCompressRad
    // overrides --radCompress.
    opts.rad_codec = if args.no_compress_rad {
        ChunkCodec::None
    } else {
        match args.rad_compress.as_str() {
            "none" => ChunkCodec::None,
            "lz4" => ChunkCodec::Lz4,
            "zstd" => ChunkCodec::Zstd,
            other => {
                anyhow::bail!("unknown --radCompress codec '{other}' (expected lz4|zstd|none)")
            }
        }
    };
    opts.num_bootstraps = args.num_bootstraps;
    opts.num_gibbs_samples = args.num_gibbs_samples;
    opts.thinning_factor = args.thinning_factor;
    opts.no_length_correction = args.no_length_correction;
    opts.model_single_frag_prob = !args.no_single_frag_prob;
    opts.no_frag_length_dist = args.no_frag_length_dist;
    opts.map_config.align.min_score_fraction = args.min_score_fraction;
    opts.map_config.pair.orphan_chain_sub_thresh = args.orphan_chain_sub_thresh;
    opts.map_config.align.full_length_alignment = args.full_length_alignment;
    opts.map_config.align.bandwidth = args.bandwidth;
    opts.map_config.pair.allow_dovetail = args.allow_dovetail;
    opts.map_config.pair.orphans_require_unmapped_mate = args.orphans_require_unmapped_mate;
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
    opts.map_config.score.allow_decoy_orphans = args.allow_decoy_orphans;
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
    opts.map_config.seed_mode = seed_mode(args.unimems, args.refmems, args.sparse_seeds);
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
    opts.em.accel = args.em_accel.into();
    opts.sig_digits = args.sig_digits;
    opts.max_read_occ = args.max_read_occ;
    opts.fld_mean = args.fld_mean.unwrap_or(DEFAULT_FLD_MEAN);
    opts.fld_sd = args.fld_sd.unwrap_or(DEFAULT_FLD_SD);
    opts.fld_max = args.fld_max.unwrap_or(DEFAULT_FLD_MAX);
    warn_unsupported_fld_policy(args.fld_policy, "reads mode");
    opts.forgetting_factor = args.forgetting_factor;

    // Deterministic mode: map once to an intermediate RAD, then quantify from it
    // with a fixed FLD (byte-identical output, no second mapping pass).
    if args.deterministic {
        let rad_out = opts.write_rad.take(); // honour an explicit --writeRad path
        return run_deterministic(
            opts,
            rad_out,
            args.keep_rad,
            args.targets.clone(),
            gene_map.as_ref(),
            &out_dir,
        );
    }

    // Live mapping spinner on an interactive terminal (unless --quiet/--no-progress).
    let progress = Arc::new(ProgressCounters::default());
    let guard = if !quiet && !args.no_progress && std::io::stderr().is_terminal() {
        opts.progress = Some(progress.clone());
        Some(ProgressGuard::start(progress.clone()))
    } else {
        None
    };
    let res = quantify(&opts).context("quantification failed")?;
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
    if res.inference_truncated_mass > 0.0 {
        tracing::warn!(
            "{:.3} fragments of equivalence-class mass could not be assigned (every member \
             transcript was truncated below the min-alpha threshold); reported as \
             inference_truncated_mass in meta_info.json",
            res.inference_truncated_mass
        );
    }
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
// An easter egg, and a deliberate one: salmon has always had it.
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
        seed_mode: seed_mode(args.unimems, args.refmems, args.sparse_seeds),
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

/// `--fldPolicy` selects between a RAD's baked fragment-length distribution and
/// recomputing one, so it means nothing where there is no RAD header to read.
/// Say so rather than accepting it silently.
// Warn rather than ignore. A flag that has no effect in the chosen mode is a
// misunderstanding on the user's part, and silently discarding it would leave
// them believing they had configured something they had not.
fn warn_unsupported_fld_policy(policy: FldPolicyArg, mode: &str) {
    if policy == FldPolicyArg::default() {
        return;
    }
    tracing::warn!(
        "--fldPolicy has no effect in {mode}: it chooses whether to use the \
         fragment-length distribution baked into a RAD header, which only \
         applies to --rad; ignoring it"
    );
}

/// `--writeMappings`/`--writeSam`/`--writeBam` only apply to the read-mapping
/// path. Rather than ignoring them silently in the alignment and RAD-input
/// modes, say so — the repo's accept-and-warn policy for flags that are no-ops
/// in the selected mode.
fn warn_unsupported_mapping_output(
    write_mappings: &Option<PathBuf>,
    write_bam: &Option<PathBuf>,
    mode: &str,
    why: &str,
) {
    let flag = match (write_mappings, write_bam) {
        (Some(_), _) => "--writeMappings/--writeSam",
        (_, Some(_)) => "--writeBam",
        _ => return,
    };
    tracing::warn!("{flag} has no effect in {mode}: {why}; ignoring it");
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;

    /// Parse a minimal `salmon quant` invocation plus `extra`, returning the
    /// mapping-output flags only (`QuantArgs` itself is not `Debug`).
    /// Parse a `quant` command line and return just the two mapping-output paths.
    fn write_flags(extra: &[&str]) -> Result<(Option<PathBuf>, Option<PathBuf>), clap::Error> {
        quant_args(extra).map(|args| (args.write_mappings, args.write_bam))
    }

    /// Parse a minimal `salmon quant` invocation plus `extra`.
    /// Parse a `quant` command line with the mandatory arguments filled in, so
    /// each test only states the flags it is actually about.
    fn quant_args(extra: &[&str]) -> Result<QuantArgs, clap::Error> {
        let mut argv = vec![
            "salmon", "quant", "-i", "idx", "-l", "A", "-r", "r.fq", "-o", "out",
        ];
        argv.extend_from_slice(extra);
        match Cli::try_parse_from(argv)?.command {
            Command::Quant(args) => Ok(args),
            _ => unreachable!("parsed a non-quant subcommand"),
        }
    }

    /// `--writeSam` is a spelling of `--writeMappings`, giving SAM output a
    /// format-named flag to match `--writeBam`.
    #[test]
    /// `--writeSam` is kept as a visible alias for compatibility, so it must land
    /// in the same field rather than being a second, separate option.
    fn write_sam_is_an_alias_for_write_mappings() {
        let expected = Some(PathBuf::from("m.sam"));
        for sam_flag in ["--writeMappings", "--writeSam", "-z"] {
            let (mappings, bam) = write_flags(&[sam_flag, "m.sam"]).unwrap();
            assert_eq!(mappings, expected, "{sam_flag} should set write_mappings");
            assert_eq!(bam, None, "{sam_flag} should not set write_bam");
        }
    }

    /// The alias shares `write_mappings`' arg id, so it inherits the
    /// `--writeBam` conflict rather than silently allowing both formats.
    #[test]
    /// Asking for both mapping-output formats at once is a usage error, and must
    /// be reported as one instead of silently honouring whichever came last.
    fn write_sam_conflicts_with_write_bam() {
        for sam_flag in ["--writeMappings", "--writeSam", "-z"] {
            let err = write_flags(&[sam_flag, "m.sam", "--writeBam", "m.bam"])
                .expect_err("should conflict with --writeBam");
            assert_eq!(
                err.kind(),
                clap::error::ErrorKind::ArgumentConflict,
                "{sam_flag} + --writeBam"
            );
        }
    }

    /// `--bamCompressThreads` is user-facing but only meaningful with
    /// `--writeBam`: unset means "derive it", and it must reject 0 rather than
    /// silently starting no compression workers.
    #[test]
    fn bam_compress_threads_is_optional_and_positive() {
        let derived = quant_args(&["--writeBam", "m.bam"]).unwrap();
        assert_eq!(
            derived.bam_compress_threads, None,
            "unset must stay None so the derived default applies"
        );

        let explicit = quant_args(&["--writeBam", "m.bam", "--bamCompressThreads", "4"]).unwrap();
        assert_eq!(
            explicit.bam_compress_threads.map(NonZeroUsize::get),
            Some(4)
        );

        let zero = quant_args(&["--writeBam", "m.bam", "--bamCompressThreads", "0"])
            .map(|_| ())
            .expect_err("0 compression workers is not a valid pool size");
        assert_eq!(zero.kind(), clap::error::ErrorKind::ValueValidation);
    }

    /// Without `--writeBam` there is nothing to compress, so the flag is a
    /// usage error rather than a silently ignored value.
    #[test]
    fn bam_compress_threads_requires_write_bam() {
        let err = quant_args(&["--bamCompressThreads", "4"])
            .map(|_| ())
            .expect_err("--bamCompressThreads alone should be rejected");
        assert_eq!(err.kind(), clap::error::ErrorKind::MissingRequiredArgument);

        // SAM output has no compression stage, so pairing the two is a
        // conflict, not a silently ignored value.
        for sam_flag in ["--writeMappings", "--writeSam", "-z"] {
            let err = quant_args(&[sam_flag, "m.sam", "--bamCompressThreads", "4"])
                .map(|_| ())
                .expect_err("--bamCompressThreads needs --writeBam");
            assert_eq!(
                err.kind(),
                clap::error::ErrorKind::ArgumentConflict,
                "{sam_flag} + --bamCompressThreads"
            );
        }
    }

    /// #1062: salmon can only warn "your --fldMean was ignored" if it can tell a
    /// supplied value from an inherited default, which is why these are `Option`
    /// rather than `default_value_t`.
    #[test]
    /// The warning about ignored fragment-length arguments must fire only for
    /// values the user actually typed — warning about a default the user never
    /// set would be noise.
    fn fld_prior_args_distinguish_supplied_from_default() {
        let bare = quant_args(&[]).unwrap();
        assert_eq!(bare.fld_mean, None);
        assert_eq!(bare.fld_sd, None);
        assert_eq!(bare.fld_max, None);

        let given = quant_args(&["--fldMean", "150", "--fldSD", "60"]).unwrap();
        assert_eq!(given.fld_mean, Some(150.0));
        assert_eq!(given.fld_sd, Some(60.0));
        assert_eq!(given.fld_max, None, "--fldMax was not supplied");

        // Supplying the default value explicitly still counts as supplying it:
        // the user asked for it, so an ignored-value warning is warranted.
        let explicit_default = quant_args(&["--fldMean", "250"]).unwrap();
        assert_eq!(explicit_default.fld_mean, Some(DEFAULT_FLD_MEAN));
    }

    /// The policy must default to `baked` so an ordinary requant keeps
    /// reproducing the run that wrote the RAD.
    #[test]
    fn fld_policy_defaults_to_baked_and_parses_all_variants() {
        assert_eq!(quant_args(&[]).unwrap().fld_policy, FldPolicyArg::Baked);
        for (arg, expected) in [
            ("baked", FldPolicyArg::Baked),
            ("derive", FldPolicyArg::Derive),
            ("prior", FldPolicyArg::Prior),
        ] {
            let parsed = quant_args(&["--fldPolicy", arg]).unwrap().fld_policy;
            assert_eq!(parsed, expected, "--fldPolicy {arg}");
            assert_eq!(FldPolicy::from(parsed), FldPolicy::from(expected));
        }
        assert!(
            quant_args(&["--fldPolicy", "bogus"]).map(|_| ()).is_err(),
            "unknown policies must be rejected"
        );
    }
}
