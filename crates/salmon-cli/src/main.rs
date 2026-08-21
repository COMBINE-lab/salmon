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
//! names mirror C++ salmon where they overlap. Alignment-based `quant -a`
//! (online and deterministic, including genome-BAM projection via
//! `--annotation`) and `quantmerge` are fully implemented; `alevin` is a
//! permanent redirect to the alevin-fry ecosystem.

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
    /// directory is created before the build and removed afterwards unless
    /// --keepIntermediate is set.
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
    /// Mutually exclusive with the read inputs (-1/-2/-r): quantifying a BAM
    /// while silently dropping supplied FASTQs was the "wrong mode chosen
    /// silently" failure the mode dispatch exists to prevent (#1140).
    #[arg(
        short = 'a',
        long = "alignments",
        conflicts_with_all = ["mates1", "mates2", "unmated"]
    )]
    alignments: Option<PathBuf>,
    /// RAD-input mode: a RAD file of mappings (from `salmon quant --writeRad` or
    /// `piscem map-bulk`) to quantify directly, in parallel.
    #[arg(long = "rad", conflicts_with_all = ["alignments", "mates1", "mates2", "unmated"])]
    rad: Option<PathBuf>,
    /// Transcriptome FASTA for alignment mode (`-a`) and RAD-input mode
    /// (`--rad`): supplies the reference sequences that seq/GC bias correction
    /// needs, and (under `--online -a`) enables the alignment error model — the
    /// default alignment path takes `--errorModel` instead. Ignored in reads
    /// mode, where the index already carries the sequences.
    #[arg(short = 't', long = "targets")]
    targets: Option<PathBuf>,
    /// Disable the alignment error model (alignment mode).
    /// Direct opposite of `--errorModel`; passing both was resolved by silent
    /// precedence (`--noErrorModel` won) and salmon then advised passing
    /// `--errorModel`, which was already on the command line (#1140).
    #[arg(long = "noErrorModel", conflicts_with = "error_model")]
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
    /// Genome projection is reachable only through `-a <genome.bam>`, so
    /// without it the whole request — annotation, genome and junction discount
    /// — was dropped without a word, and the paths were never even checked for
    /// existence. `--genome`/`--juncMissDiscount` already required this flag;
    /// this closes the chain (#1140).
    ///
    /// The `requires` below is necessary but not sufficient: clap skips a
    /// `requires` when the required argument conflicts with one already
    /// present, and `-1/-2`, `-r` and `--rad` all conflict with `-a` — which is
    /// precisely how a user gets here. `run_quant` therefore checks it
    /// explicitly as well.
    #[arg(long = "annotation", value_name = "GTF|GFF", requires = "alignments")]
    annotation: Option<PathBuf>,
    // `--genome` / `--juncMissDiscount` are read only by the genome-projection
    // branch, which `--annotation` selects; without it they were silently
    // accepted and dead (#1140), so clap now enforces the dependency.
    /// Genome FASTA (genome-alignment mode) — needed only for bias correction, to
    /// reconstruct transcript sequences from the annotation's exons.
    #[arg(long = "genome", value_name = "FASTA", requires = "annotation")]
    genome: Option<PathBuf>,
    /// Genome-alignment projection: penalty multiplier per junction mismatch
    /// (bramble `junc_miss_discount`; default 1.0 = no penalty).
    #[arg(long = "juncMissDiscount", value_name = "F", requires = "annotation")]
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
    /// Single-end FASTQ file(s). Mutually exclusive with the paired inputs:
    /// supplying both silently quantified only the pairs while `cmd_info.json`
    /// still recorded the single-end files, which is the same "wrong mode
    /// chosen silently" failure `-a`'s conflict exists to prevent (#1140).
    #[arg(short = 'r', long = "unmatedReads", num_args = 1.., conflicts_with_all = ["mates1", "mates2"])]
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
    /// Which gzip decoder to use: `auto` (default), `serial`, `parallel`, or
    /// `parallel=N` to pin N slots per decodable input.
    ///
    /// `-p` is one budget shared between inflating the input and mapping it, so
    /// engaging the parallel decoder costs mapping threads. `auto` engages it
    /// only when the budget is large enough for it to add decode concurrency
    /// rather than merely move threads around -- a threshold that is higher in
    /// selective-alignment mode than in sketch, because more work per fragment
    /// means a smaller share of the budget should be decoding.
    #[arg(long = "decoder", value_parser = ["auto", "serial", "parallel", "libz", "libdeflate", "rapidgzip"])]
    decoder: Option<String>,
    /// JSON file overriding thread and decoder policy.
    ///
    /// Every field is optional and an unknown field is an error rather than a
    /// silent no-op, so a typo refuses to load instead of leaving the run
    /// looking configured while behaving as if it were not. Example:
    /// `{"parallel_decode": {"min_threads_per_stream": 16}}`.
    #[arg(long = "threadPolicy", value_name = "FILE")]
    thread_policy: Option<std::path::PathBuf>,
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
    /// Write per-mapping SAM records to this file. Records carry base qualities,
    /// a base-level CIGAR with NM and MD, and the mate and placement tags.
    /// `--writeSam` is an accepted alias, naming the format explicitly to pair
    /// with `--writeBam`.
    #[arg(
        short = 'z',
        long = "writeMappings",
        visible_alias = "writeSam",
        conflicts_with = "write_bam"
    )]
    write_mappings: Option<PathBuf>,
    /// Write per-mapping BAM records to this file: the same records as
    /// --writeMappings, BGZF-compressed. Mutually exclusive with
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
    /// Deterministic quantification: map once to an intermediate RAD, then
    /// quantify from it with a fixed fragment-length distribution, so the
    /// result is byte-identical across runs and thread counts. **The default
    /// since 2.6.0** — the flag is accepted as a no-op; use --online for the
    /// deprecated one-pass path. The intermediate RAD is written under
    /// the output directory and deleted on success unless --keepRad (or use
    /// --writeRad PATH to choose its location and keep it). With --skipQuant the
    /// run stops after mapping and the RAD is kept as its output, in the output
    /// directory even when --radScratchDir is given.
    #[arg(long = "deterministic")]
    deterministic: bool,
    /// Quantify with the pre-2.6 online one-pass path instead of the default
    /// deterministic two-phase flow. Deprecated: results depend on the thread
    /// count (the default is byte-identical at any -p); scheduled for removal
    /// in 2.7.0.
    #[arg(long = "online", conflicts_with = "deterministic")]
    online: bool,
    /// Keep the intermediate RAD produced by --deterministic (by default it is
    /// deleted once quantification finishes). Ignored without --deterministic.
    #[arg(long = "keepRad")]
    keep_rad: bool,
    /// Directory for the intermediate RAD written by --deterministic (default:
    /// the output directory). Point this at node-local or memory-backed
    /// storage (e.g. /dev/shm) when the output volume is slow, networked, or
    /// space-constrained. Ignored when --writeRad names an explicit path.
    #[arg(long = "radScratchDir")]
    rad_scratch_dir: Option<PathBuf>,
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
    /// Option-typed (default applied at use: `lz4`) so passing it on a path
    /// that writes no RAD can warn instead of being silently inert (#1140).
    #[arg(long = "radCompress", value_parser = ["lz4", "zstd", "none"])]
    rad_compress: Option<String>,
    /// Write uncompressed RAD chunks (overrides --radCompress).
    #[arg(long = "noCompressRad")]
    no_compress_rad: bool,
    /// Number of bootstrap replicates for posterior uncertainty.
    #[arg(long = "numBootstraps", default_value_t = 0)]
    num_bootstraps: u32,
    /// Number of Gibbs posterior samples (mutually exclusive with bootstraps).
    ///
    /// `conflicts_with` makes that exclusivity real. It was previously only
    /// stated in this doc comment and resolved silently by dispatch order, so
    /// passing both ran the bootstrap and ignored the Gibbs request without
    /// saying so.
    #[arg(
        long = "numGibbsSamples",
        default_value_t = 0,
        conflicts_with = "num_bootstraps"
    )]
    num_gibbs_samples: u32,
    /// Gibbs thinning factor.
    /// Option-typed (default 16 at use) so passing it without Gibbs sampling
    /// can be reported rather than silently doing nothing — the sibling
    /// `--bamCompressThreads` without `--writeBam` already hard-errors.
    #[arg(long = "thinningFactor")]
    thinning_factor: Option<u32>,
    /// (Long reads) Oxford Nanopore model — not supported; use oarfish instead.
    #[arg(long = "ont")]
    ont: bool,
    /// Disable effective-length correction (use raw reference length).
    #[arg(long = "noLengthCorrection")]
    no_length_correction: bool,
    /// Minimum alignment score as a fraction of the perfect score (default 0.65).
    /// Option-typed so passing it under `--sketch` — which computes no alignment
    /// scores — can warn instead of being silently inert (#1140).
    #[arg(long = "minScoreFraction")]
    min_score_fraction: Option<f32>,
    /// Orphan chain sub-optimality threshold (salmon's `orphanChainSubThresh`).
    /// `0.0` (default) aligns every orphan candidate (more sensitive than salmon's
    /// 0.95, which prunes low-chain-coverage orphans before alignment).
    #[arg(long = "orphanChainSubThresh")]
    orphan_chain_sub_thresh: Option<f32>,
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
    /// Match score for selective alignment (reads mode; default 2).
    #[arg(long = "ma")]
    ma: Option<i32>,
    /// Mismatch penalty for selective alignment (reads mode; default 4).
    #[arg(long = "mp")]
    mp: Option<i32>,
    /// Gap-open penalty for selective alignment (reads mode; default 6).
    #[arg(long = "go")]
    go: Option<i32>,
    /// Gap-extend penalty for selective alignment (reads mode; default 2).
    #[arg(long = "ge")]
    ge: Option<i32>,
    /// Consensus slack: a target is kept only if its best chain score is at least
    /// `(1 - slack)` of the max chain score for that mate (salmon default 0.35;
    /// `1.0` keeps every target).
    #[arg(long = "consensusSlack")]
    consensus_slack: Option<f32>,
    /// Skip k-mers whose unitig occurs in more than this many references
    /// (repetitive-hit guard; salmon's `maxOccsPerHit`).
    #[arg(long = "maxOccsPerHit")]
    max_occs_per_hit: Option<usize>,
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
    /// Online-phase forgetting factor in (0.5, 1.0]; default 0.65. Online-path
    /// only (--online, deprecated): ignored with a warning under the default
    /// deterministic mode.
    #[arg(short = 'f', long = "forgettingFactor")]
    forgetting_factor: Option<f64>,
    /// Initialize the optimizer uniformly. In reads mode (one-pass and
    /// --deterministic) the optimizer already starts uniform, so this is the
    /// existing behavior; in online alignment mode (-a) it replaces the
    /// warm start from the online estimates.
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
    #[arg(short = 'w', long = "maxReadOcc")]
    max_read_occ: Option<usize>,
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
    #[arg(long = "bandwidth")]
    bandwidth: Option<i32>,
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
    #[arg(long = "decoyThreshold")]
    decoy_threshold: Option<f64>,
    /// Drop any mapping whose alignment probability
    /// `exp(-scoreExp * (best - score))` is below this (selective-alignment mode).
    #[arg(long = "minAlnProb")]
    min_aln_prob: Option<f64>,
    /// Soft-weight decay: a mapping's probability is proportional to
    /// `exp(-scoreExp * (bestScore - score))`. Larger downweights sub-optimal
    /// mappings more steeply. (salmon's --scoreExp; default 1.0)
    #[arg(long = "scoreExp")]
    score_exp: Option<f64>,
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
    #[arg(long = "preMergeChainSubThresh")]
    pre_merge_chain_sub_thresh: Option<f32>,
    /// Post-merge concordant chain-pair sub-optimality threshold: after pairing,
    /// keep pairs whose read-coverage is `>= best_pair_coverage * thresh`
    /// (paired-end, selective-alignment mode). Range [0,1]. Default 0.0 (off);
    /// salmon's default is 0.9, thresholded on chain score rather than coverage.
    #[arg(long = "postMergeChainSubThresh")]
    post_merge_chain_sub_thresh: Option<f32>,
    /// Fragment-length sampling stride for the GC bias eff-length convolution.
    /// Larger = faster bias correction, coarser. (salmon's --biasSpeedSamp)
    #[arg(long = "biasSpeedSamp")]
    bias_speed_samp: Option<usize>,
    /// Number of leading fragments used to train the auxiliary models
    /// (fragment-length / bias / error); they are fixed afterward; default
    /// 5000000. Online-path only (--online, deprecated): ignored with a
    /// warning under the default deterministic mode.
    #[arg(long = "numAuxModelSamples")]
    num_aux_model_samples: Option<u64>,
    /// Disable the lower threshold on how short bias correction may make an
    /// effective length (increases precision, reduces robustness). Experimental.
    #[arg(long = "noBiasLengthThreshold")]
    no_bias_length_threshold: bool,
    /// Number of read-position bins in the alignment error model
    /// (alignment mode only). (salmon's --numErrorBins)
    #[arg(long = "numErrorBins")]
    num_error_bins: Option<usize>,
    /// Number of fragment-GC bins for the GC bias model. (salmon's --numGCBins)
    #[arg(long = "numGCBins")]
    num_gc_bins: Option<usize>,
    /// Number of conditioning (context) bins for the GC bias model.
    /// (salmon's --conditionalGCBins)
    #[arg(long = "conditionalGCBins")]
    conditional_gc_bins: Option<usize>,
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
    /// Skip transcript quantification (EM + Gibbs/bootstrap) and quant.sf: map
    /// (or read the alignments), detect the library type, and write metadata.
    /// The run's deliverable is the RAD of the mappings, kept in the output
    /// directory; requantify it later with `--rad <file>`. On the default reads
    /// path the run stops after mapping, so equivalence classes come from that
    /// requant (`--rad <file> --dumpEq`), not from this run. (salmon's --skipQuant)
    #[arg(long = "skipQuant")]
    skip_quant: bool,
    /// [not yet implemented] quantify from a previously-written equivalence-class
    /// file instead of mapping/alignments (an input mode, not the --dumpEq output).
    #[arg(long = "eqclasses")]
    eqclasses: Option<PathBuf>,
    /// Fragments processed before the auxiliary (fragment-length) model is
    /// applied during online inference. (salmon's --numPreAuxModelSamples;
    /// salmon's default is 1000000, this port's is 5000.) Online-path only
    /// (--online, deprecated): ignored with a warning under the default
    /// deterministic mode.
    #[arg(long = "numPreAuxModelSamples")]
    num_pre_aux_model_samples: Option<u64>,
    /// Disable fragment-length-distribution concordance in the per-fragment
    /// probability. Honoured on every quantification path except the deprecated
    /// online alignment path (--online -a), which warns that it lacks it.
    #[arg(long = "noFragLengthDist")]
    no_frag_length_dist: bool,
    /// Disable the single-end/orphan fragment-length probability estimate:
    /// orphans take a flat penalty instead of the bounded-CMF ambiguous weight,
    /// which can change which transcript wins a contested orphan. Honoured on
    /// the reads paths and on every RAD-based path (deterministic reads phase 2,
    /// -a, --rad). Redundant on the deprecated online alignment path
    /// (--online -a), which weights every orphan flat by construction.
    #[arg(long = "noSingleFragProb")]
    no_single_frag_prob: bool,
    /// [accepted; not yet implemented] write a BAM of posterior-sampled
    /// alignments (alignment mode).
    #[arg(short = 's', long = "sampleOut")]
    sample_out: bool,
    /// Also write a record for every fragment that did not map, flagged 0x4, so
    /// the mapping output covers the whole library rather than only the part
    /// that mapped. Applies to --writeMappings/--writeSam and --writeBam.
    #[arg(short = 'u', long = "sampleUnaligned")]
    sample_unaligned: bool,
    /// [accepted for salmon compatibility; has no effect] base qualities are not
    /// optional output. They are written into every record whenever the input
    /// carries them, with or without this flag.
    #[arg(long = "writeQualities")]
    write_qualities: bool,
    /// Read group to declare in the mapping output's header, as tab-separated
    /// KEY:value fields with a mandatory ID (for example
    /// `ID:sample1\tSM:sample1\tPL:ILLUMINA`). A literal \t is accepted in place
    /// of a tab. Every record is then tagged RG:Z:<ID>.
    #[arg(long = "rgLine", value_name = "LINE")]
    rg_line: Option<String>,
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
    /// The parsed transcript → gene map, read before quantification starts.
    ///
    /// Held rather than re-read at the end because gene-level aggregation is the
    /// last thing a run does: reading it there meant an unreadable annotation
    /// failed the run *after* mapping, inference and `quant.sf` had all completed
    /// (#1074). Loading it up front costs one parse of a file the run needs
    /// anyway, and turns minutes of wasted work into an immediate error.
    map: std::collections::HashMap<String, String>,
}

/// Read and validate `--geneMap` before any quantification work begins.
///
/// Returns `None` when the option was not given. Any failure — missing file,
/// unreadable, malformed — surfaces here rather than after the run (#1074).
fn load_gene_map(path: Option<PathBuf>, ignore_tx_version: bool) -> Result<Option<GeneMapOpts>> {
    let Some(path) = path else {
        return Ok(None);
    };
    tracing::warn!(
        "tximport (R) and pytximport (Python) are the preferred way to obtain \
         gene-level estimates from salmon output: they offer several ways to \
         aggregate transcript abundances to the gene level, including one that \
         computes an offset accounting for changes in average transcript length \
         between samples. --geneMap may be removed in a future release."
    );
    let map = salmon_core::genemap::read_transcript_gene_map(&path)
        .with_context(|| format!("reading gene map {}", path.display()))?;
    tracing::info!(
        "read {} transcript to gene mappings from {}",
        map.len(),
        path.display()
    );
    Ok(Some(GeneMapOpts {
        path,
        ignore_tx_version,
        map,
    }))
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
    // Parsed up front by `load_gene_map`, so nothing here can fail on a bad path.
    let map = &gene_map.map;
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
        map,
        gene_map.ignore_tx_version,
        Some(&unmatched_list),
    )
    .context("writing quant.genes.sf")?;

    if summary.unmatched > 0 {
        eprintln!(
            "warning: {} of {} transcript(s) had no entry in the gene map {}. Each was \
             written to quant.genes.sf as its own single-transcript gene, so those rows are \
             named by transcript, not by gene. Their names are listed in {}",
            summary.unmatched,
            summary.total(),
            gene_map.path.display(),
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
        let would_match = salmon_core::genemap::matches_ignoring_tx_version(names, map);
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

/// The phase-2 options for a deterministic reads-mode run: every knob from the
/// reads-mode request that the RAD quantifier can honour, carried across.
///
/// This is a separate function so the forwarding can be asserted field by field
/// in a test. A knob that is simply forgotten here does not fail loudly: the
/// second pass runs with the default, the run succeeds, and the flag the user
/// passed silently did nothing. That is how `--biasSpeedSamp` and
/// `--noBiasLengthThreshold` came to be ignored under `--deterministic`.
///
/// # Deliberately not forwarded
///
/// The point of the function is to make "audited and correctly omitted"
/// distinguishable from "forgotten", which a test cannot do for an absence. Every
/// field of `AlignQuantOptions` this leaves alone is here, with why:
///
/// * `ref_seqs`: needs to read the index from disk, and only the caller knows
///   whether bias correction asked for it. Set immediately after this returns.
/// * `explicit_fld_args`, `fld_policy`: phase 1 already folded the user's FLD
///   prior into the bake. Forwarding them would have phase 2 warn that a baked
///   FLD overrides the run's own `--fldMean`, against a bake that came from that
///   very argument.
/// * `forgetting_factor`, `num_aux_model_samples`, `num_pre_aux_model_samples`:
///   online-inference knobs. Phase 2 has no online phase to apply them to.
/// * `num_error_bins`, `discard_orphans`, `no_error_model`,
///   `deterministic_error_model`: BAM-side knobs. Phase 2 reads a RAD whose
///   placements were already scored and filtered during mapping.
/// * `bam`, `output_dir`, `progress`: set by the caller, which owns the paths and
///   the run's lifecycle.
/// * `prior_seconds`, `external_start_time`: computed from the caller's clock
///   after phase 1 finishes, so they cannot be part of a mapping done before it
///   starts; the caller sets them just before phase 2.
///
/// One absence here was a bug rather than a decision, found *from* this list:
/// `skip_quant` is set on phase 1 unconditionally (its job is to map and write
/// the RAD) and never reached phase 2, so `--skipQuant --deterministic`
/// quantified anyway — fixed in #1149 by stopping after phase 1. Being able to
/// see that from this list is the argument for the list.
fn requant_options(
    map_opts: &QuantOptions,
    rad_path: &std::path::Path,
    out_dir: &std::path::Path,
    bias_targets: Option<PathBuf>,
) -> AlignQuantOptions {
    let mut q = AlignQuantOptions::new(rad_path.to_path_buf(), out_dir.to_path_buf());
    q.lib_type = map_opts.lib_type.clone();
    q.em = map_opts.em.clone();
    q.range_factorization_bins = map_opts.range_factorization_bins;
    q.transcripts = bias_targets;
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
    // The three below were missing until #1140's audit. They are not
    // online-only knobs: the RAD quantifier honours all three, so dropping them
    // meant the flag was accepted and then ignored.
    q.bias_speed_samp = map_opts.bias_speed_samp;
    q.no_bias_length_threshold = map_opts.no_bias_length_threshold;
    q.init_uniform = map_opts.init_uniform;
    // Phase 2 is where the equivalence classes exist, so it is where they are
    // dumped. Phase 1 has none to dump (COMBINE-lab/salmon#1140), which is why
    // this mapping has to be built from the request as the user made it, before
    // the driver clears the flags for phase 1.
    q.dump_eq = map_opts.dump_eq;
    q.dump_eq_weights = map_opts.dump_eq_weights;
    // Honoured by the RAD quantifier's orphan weighting, so the flag survives
    // the phase-1/phase-2 split instead of dying with the mapping pass
    // (#1140, audit D14).
    q.model_single_frag_prob = map_opts.model_single_frag_prob;
    // Phase 1's mapper already applied the hard filter when choosing what to
    // write; phase 2 applies it again over the RAD's per-hit scores, which is
    // where the eq-class weights are actually formed.
    q.hard_filter = map_opts.map_config.score.hard_filter;
    // Same shape as the eq dumps: phase 2 owns the run's final bias models
    // (phase 1's dump site is gated off by skip_quant), so `--dumpBiasModels`
    // is honoured by the requant's shared writer (#1140: the flag silently
    // vanished when deterministic became the default).
    q.dump_bias_models = map_opts.dump_bias_models;
    // The RAD path gained these two with #1148; before that they existed only in
    // reads mode, so `--deterministic` accepted them and quantified as though
    // they had not been passed.
    q.no_length_correction = map_opts.no_length_correction;
    q.no_frag_length_dist = map_opts.no_frag_length_dist;
    q
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
/// Where a two-phase run's salmon-owned intermediate RAD lives: in
/// `--radScratchDir` when given (created if needed; pid-suffixed so concurrent
/// runs sharing one scratch volume cannot collide at the finalize rename), else
/// under the output directory as before.
fn intermediate_rad_path(
    out_dir: &std::path::Path,
    scratch_dir: Option<&std::path::Path>,
    stem: &str,
) -> Result<PathBuf> {
    match scratch_dir {
        Some(d) => {
            std::fs::create_dir_all(d)
                .with_context(|| format!("creating --radScratchDir {}", d.display()))?;
            Ok(d.join(format!("{stem}-{}.rad", std::process::id())))
        }
        None => Ok(out_dir.join(format!("{stem}.rad"))),
    }
}

/// Best-effort removal of a salmon-owned intermediate RAD after a failed phase:
/// the finalized file if it exists, and any `.partial-*` sibling a failed
/// mapping pass left behind. The failure that fills a disk is precisely the one
/// where the user needs the space back (#1140); an explicit `--writeRad` path
/// is user-owned evidence and is never cleaned here.
fn cleanup_intermediate_rad(path: &std::path::Path) {
    let _ = std::fs::remove_file(path);
    if let (Some(dir), Some(name)) = (path.parent(), path.file_name().and_then(|n| n.to_str())) {
        let prefix = format!("{name}.partial-");
        if let Ok(rd) = std::fs::read_dir(dir) {
            for e in rd.flatten() {
                if e.file_name()
                    .to_str()
                    .is_some_and(|n| n.starts_with(&prefix))
                {
                    let _ = std::fs::remove_file(e.path());
                }
            }
        }
    }
}

/// Total size of the given input files, for the intermediate-RAD space
/// preflight; unreadable paths contribute 0 (the mapping pass will report them
/// properly).
fn total_file_bytes<'a>(paths: impl IntoIterator<Item = &'a PathBuf>) -> u64 {
    paths
        .into_iter()
        .filter_map(|p| std::fs::metadata(p).ok())
        .map(|m| m.len())
        .sum()
}

/// Warn up front when the volume that will hold the intermediate RAD has less
/// free space than the total compressed input. The lz4 intermediate is
/// typically of the same order as the compressed input (and grows with the
/// multimapping rate, so no fixed per-fragment constant is safe — #1140); this
/// is a warning rather than an error because the estimate is a proxy, and the
/// RAD writer's periodic low-space check turns a genuinely filling disk into a
/// prompt, actionable failure.
fn preflight_rad_space(rad_path: &std::path::Path, input_bytes: u64) {
    let dir = rad_path.parent().unwrap_or(std::path::Path::new("."));
    if let Some(free) = salmon_core::free_disk_bytes(dir) {
        if free < input_bytes {
            let gib = |b: u64| b as f64 / (1024.0 * 1024.0 * 1024.0);
            tracing::warn!(
                "{} has {:.1} GiB free, less than the {:.1} GiB of input; the \
                 deterministic intermediate RAD is typically of the same order as \
                 the compressed input and may not fit. Consider --radScratchDir \
                 <dir> (e.g. /dev/shm or node-local scratch) or --radCompress zstd.",
                dir.display(),
                gib(free),
                gib(input_bytes)
            );
        }
    }
}

// `--deterministic` (reads mode): a two-pass run. The first pass maps to a RAD
// with order-independent models baked into the header; the second quantifies
// from it. Splitting the run this way is what makes the result independent of
// the thread count.
fn run_deterministic(
    mut map_opts: QuantOptions,
    rad_out: Option<PathBuf>,
    scratch_dir: Option<&std::path::Path>,
    keep_rad: bool,
    gene_map: Option<&GeneMapOpts>,
    out_dir: &std::path::Path,
    // The live mapping spinner, owned here so it can be stopped the moment
    // phase 1 finishes (its counters do not advance during the requant).
    spinner: Option<ProgressGuard>,
) -> Result<()> {
    // Wall-clock accounting spans both phases: phase 2 rewrites
    // `meta_info.json` into the final output directory, and its
    // `total_time_seconds` must cover the whole run, not just the requant
    // (#1140: a 10M-pair run otherwise reports ~0.4s against minutes of wall).
    let run_start = std::time::Instant::now();
    let bias_on = map_opts.seq_bias || map_opts.gc_bias || map_opts.pos_bias;
    // An explicit `--writeRad PATH` is honoured (and kept — it was a requested
    // output); otherwise a temp under the output directory, removed on success
    // unless `--keepRad`.
    let explicit = rad_out.is_some();
    if explicit && scratch_dir.is_some() {
        tracing::warn!("--radScratchDir is ignored when --writeRad names an explicit path");
    }
    // Scratch placement is for intermediates. Under `--skipQuant` there is no
    // phase 2 and the RAD is the run's output, so it belongs in the output
    // directory under its stable name rather than pid-suffixed on a scratch
    // volume nothing was meant to be kept in. An explicit `--writeRad` still
    // wins over both.
    let deliverable = map_opts.skip_quant;
    if deliverable && !explicit && scratch_dir.is_some() {
        tracing::info!(
            "--skipQuant makes the RAD this run's output, so it is written to the output \
             directory rather than to --radScratchDir"
        );
    }
    let rad_path = match rad_out {
        Some(p) => p,
        None => intermediate_rad_path(
            out_dir,
            if deliverable { None } else { scratch_dir },
            "intermediate_mappings",
        )?,
    };
    // The RAD writer opens its file before the mapping pass, so the output
    // directory (where the default intermediate lives) must exist first.
    std::fs::create_dir_all(out_dir)
        .with_context(|| format!("creating output directory {}", out_dir.display()))?;
    preflight_rad_space(
        &rad_path,
        total_file_bytes(
            map_opts
                .mates1
                .iter()
                .chain(map_opts.mates2.iter())
                .chain(map_opts.unmated.iter()),
        ),
    );

    // Phase 2's options, built from the request exactly as the user made it and
    // before phase 1 rewrites parts of it for its own purposes. Everything it
    // needs (the RAD path, the bias targets) is known by now, and building it
    // here is what keeps "the phase-2 options" meaning the user's request rather
    // than whatever phase 1 left behind.
    // No `-t` fallback here: reads mode warns and ignores the flag (the index
    // carries the reference sequences; see the bias load below). The
    // `bias_targets` seam on `requant_options` exists for the alignment/RAD
    // drivers, which have no index to fall back on.
    let mut q = requant_options(&map_opts, &rad_path, out_dir, None);

    // Phase 1 — map once, write the RAD (bakes the deterministic FLD + resolved
    // library format), skipping the online EM: quantification happens in phase 2.
    //
    // `skip_quant` is about to be set for phase 1 whatever the user asked, since
    // phase 1's job is to map and write the RAD. Take the request first: with
    // `--skipQuant` there is no phase 2, and the RAD is the deliverable rather
    // than an intermediate. Alignment mode already worked this way.
    let stop_after_mapping = map_opts.skip_quant;
    map_opts.write_rad = Some(rad_path.clone());
    map_opts.skip_quant = true;
    // The classes live in phase 2, which already has these flags. Left on here,
    // phase 1 pays for `build_eq` and then writes a well-formed
    // `eq_classes.txt.gz` declaring zero classes, which phase 2 never overwrites
    // (COMBINE-lab/salmon#1140).
    map_opts.dump_eq = false;
    map_opts.dump_eq_weights = false;
    // Derive the FLD + library format order-independently during this pass and
    // bake them (see `QuantOptions::deterministic_fld`), so phase 2 is a single
    // pass and the whole result is byte-identical across thread counts.
    map_opts.deterministic_fld = true;
    let map_res = match quantify(&map_opts) {
        Ok(r) => r,
        Err(e) => {
            if !explicit {
                cleanup_intermediate_rad(&rad_path);
            }
            return Err(e).context("deterministic mapping pass failed");
        }
    };
    drop(spinner); // stop + clear the mapping spinner before phase-2 logging
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

    if stop_after_mapping {
        // The classes are built by phase 2, which is exactly the pass
        // `--skipQuant` asks not to run, so there is nothing to dump. Saying so
        // beats writing no file and no explanation: one-pass `--skipQuant
        // --dumpEq` does produce a dump, so the difference would otherwise look
        // like the flag being dropped.
        if q.dump_eq || q.dump_eq_weights {
            tracing::warn!(
                "--dumpEq/--dumpEqWeights has nothing to write under --skipQuant --deterministic: \
                 equivalence classes are built by the quantification pass, which --skipQuant \
                 skips. Quantify the RAD that was kept with `--rad {} --dumpEq` to produce them.",
                rad_path.display()
            );
        }
        tracing::info!(
            "--skipQuant: mapping only; the RAD at {} is the output",
            rad_path.display()
        );
        return Ok(());
    }

    // Phase 2 — deterministic quant from the RAD (fixed baked FLD + library
    // format ⇒ order-independent eq-classes + EM). Mirrors the `--rad` knob wiring.
    // Bias needs the reference sequence, and in reads mode it always comes
    // from the index (which carries the sequences by construction — we just
    // mapped against it); `-t` is warned-and-ignored at dispatch. A previous
    // comment here claimed `-t` was "honoured as a fallback", but the fallback
    // was unreachable: with bias on this load always wins, with bias off the
    // sequences are never read (#1140).
    if bias_on {
        q.ref_seqs = Some(
            salmon_index::load_ref_seqs(&map_opts.index_dir)
                .context("loading reference sequences from the index for deterministic bias")?,
        );
    }
    q.prior_seconds = run_start.elapsed().as_secs_f64();
    q.external_start_time = Some(map_res.start_time.clone());
    // Phase 1 wrote the true invocation record and named the read files; keep
    // both instead of letting the requant's RAD-centric view replace them
    // (#1140).
    q.preserve_cmd_info = true;
    q.read_files = map_opts
        .mates1
        .iter()
        .chain(map_opts.mates2.iter())
        .chain(map_opts.unmated.iter())
        .map(|p| p.display().to_string())
        .collect();
    let res = match quantify_rad(&q, &rad_path) {
        Ok(r) => r,
        Err(e) => {
            if !explicit && !keep_rad {
                cleanup_intermediate_rad(&rad_path);
            }
            return Err(e).context("deterministic requant pass failed");
        }
    };
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
    // No `--skipQuant` guard here, deliberately: this driver returns early for
    // that case (see `stop_after_mapping` above), so reaching this point means
    // phase 2 ran and produced real abundances. Guarding on
    // `map_opts.skip_quant` instead — which phase 1 forces true regardless of
    // what the user asked — made `-g` a silent no-op in the DEFAULT reads mode
    // and printed a "--skipQuant" explanation on runs that never passed it.
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
    mut opts: AlignQuantOptions,
    rad_out: Option<PathBuf>,
    scratch_dir: Option<&std::path::Path>,
    keep_rad: bool,
    gene_map: Option<&GeneMapOpts>,
    codec: ChunkCodec,
) -> Result<()> {
    // Both phases count toward the reported run time (see run_deterministic).
    let run_start = std::time::Instant::now();
    let start_time = salmon_align::asctime_now();
    // This driver converts the user's BAM into the intermediate RAD itself, so
    // the fragment-length distribution really was derived from alignments —
    // phase 2 must not report the internal handoff (`rad_baked`), which is the
    // right answer only for a RAD the user supplied (#1140).
    opts.input_is_alignments = true;
    let out_dir = opts.output_dir.clone();
    // An explicit `--writeRad PATH` is honoured and kept; otherwise a temp under
    // the output dir, removed on success unless `--keepRad` (or `--skipQuant`,
    // where the RAD is the deliverable).
    let explicit = rad_out.is_some();
    if explicit && scratch_dir.is_some() {
        tracing::warn!("--radScratchDir is ignored when --writeRad names an explicit path");
    }
    // Same rule as the reads path: `--skipQuant` makes the RAD this run's
    // output, and an output belongs in the output directory under a stable name
    // rather than pid-suffixed on a scratch volume.
    let deliverable = opts.skip_quant;
    if deliverable && !explicit && scratch_dir.is_some() {
        tracing::info!(
            "--skipQuant makes the RAD this run's output, so it is written to the output \
             directory rather than to --radScratchDir"
        );
    }
    let rad_path = match rad_out {
        Some(p) => p,
        None => intermediate_rad_path(
            &out_dir,
            if deliverable { None } else { scratch_dir },
            "intermediate_alignments",
        )?,
    };
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("creating output directory {}", out_dir.display()))?;
    preflight_rad_space(&rad_path, total_file_bytes(std::iter::once(&opts.bam)));

    // Phase 1 — one BAM pass: write placements + bake FLD / library format
    // (+ rough abundances when bias is on) so phase 2 is a single baked pass.
    let summary = match salmon_align::write_alignment_rad(&opts, &rad_path, codec) {
        Ok(r) => r,
        Err(e) => {
            if !explicit {
                cleanup_intermediate_rad(&rad_path);
            }
            return Err(e).context("deterministic alignment RAD-write pass failed");
        }
    };
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
    // The unusable-`AS` verdict is born here in phase 1, but `meta_info.json`
    // is written by phase 2 — hand it across so pipelines that only keep the
    // output directory see it in `diagnostics`, not just the log (#1140).
    if let Some(w) = &summary.score_tag_warning {
        opts.extra_diagnostics.push(salmon_core::Diagnostic::new(
            "bam_score_tags_unusable",
            "warning",
            w.clone(),
        ));
    }

    // Phase 2 — quantify from the fully-baked RAD (single pass, order-independent).
    opts.prior_seconds = run_start.elapsed().as_secs_f64();
    opts.external_start_time = Some(start_time);
    let res = match quantify_rad(&opts, &rad_path) {
        Ok(r) => r,
        Err(e) => {
            if !explicit && !keep_rad && !opts.skip_quant {
                cleanup_intermediate_rad(&rad_path);
            }
            return Err(e).context("deterministic requant pass failed");
        }
    };
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
    // `--skipQuant` estimates no abundances; aggregating the zero rows to
    // gene level would announce success over an empty result (#1140 — the
    // deterministic reads driver already skipped it, the others wrote
    // zeros silently).
    if opts.skip_quant && gene_map.is_some() {
        tracing::info!(
            "--skipQuant: gene-level aggregation (-g) skipped — no abundances were estimated"
        );
    } else if let Some(gm) = gene_map {
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
    scratch_dir: Option<&std::path::Path>,
    keep_rad: bool,
    gene_map: Option<&GeneMapOpts>,
) -> Result<()> {
    // Both phases count toward the reported run time (see run_deterministic).
    let run_start = std::time::Instant::now();
    let start_time = salmon_align::asctime_now();
    let out_dir = q.output_dir.clone();
    let explicit = rad_out.is_some();
    if explicit && scratch_dir.is_some() {
        tracing::warn!("--radScratchDir is ignored when --writeRad names an explicit path");
    }
    // Deliverable-vs-intermediate rule, same as the reads/alignment drivers
    // (#1149/#1140): under `--skipQuant` the projected RAD is the run's only
    // output, so it belongs in the output directory under its stable name —
    // not pid-suffixed in scratch, where a "successful" run would leave the
    // output directory empty.
    let deliverable = q.skip_quant;
    if deliverable && scratch_dir.is_some() {
        tracing::info!(
            "--skipQuant makes the projected RAD this run's output, so it is written to \
             the output directory rather than to --radScratchDir"
        );
    }
    let rad_path = match rad_out {
        Some(p) => p,
        None => intermediate_rad_path(
            &out_dir,
            if deliverable { None } else { scratch_dir },
            "genome_projected",
        )?,
    };
    std::fs::create_dir_all(&out_dir)
        .with_context(|| format!("creating output directory {}", out_dir.display()))?;
    preflight_rad_space(&rad_path, total_file_bytes(std::iter::once(&gp.bam)));

    // Phase 1 — project the genome BAM into a fully-baked transcriptome RAD.
    let art = match project_genome_bam_to_rad(&gp, &rad_path) {
        Ok(r) => r,
        Err(e) => {
            if !explicit {
                cleanup_intermediate_rad(&rad_path);
            }
            return Err(e).context("genome projection pass failed");
        }
    };
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
    // Genome projection is alignment input too: the FLD is observed while
    // projecting the genome BAM, not carried in from a user-supplied RAD, so
    // `frag_length_source` must say `alignments` rather than `rad_baked`.
    q.input_is_alignments = true;
    if q.ref_seqs.is_none() {
        q.ref_seqs = art.ref_seqs;
    }
    q.prior_seconds = run_start.elapsed().as_secs_f64();
    q.external_start_time = Some(start_time);
    let res = match quantify_rad(&q, &rad_path) {
        Ok(r) => r,
        Err(e) => {
            if !explicit && !keep_rad {
                cleanup_intermediate_rad(&rad_path);
            }
            return Err(e).context("genome-projected quantification failed");
        }
    };
    tracing::info!(
        "{} fragments from projected RAD, {} quantified; {} equivalence classes",
        res.num_processed,
        res.num_mapped,
        res.num_eq_classes
    );
    // `--skipQuant` estimates no abundances; aggregating the zero rows to
    // gene level would announce success over an empty result (#1140 — the
    // deterministic reads driver already skipped it, the others wrote
    // zeros silently).
    if q.skip_quant && gene_map.is_some() {
        tracing::info!(
            "--skipQuant: gene-level aggregation (-g) skipped — no abundances were estimated"
        );
    } else if let Some(gm) = gene_map {
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
    if explicit || keep_rad || deliverable {
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
/// Reject an obviously incomplete uncompressed FASTQ before mapping it.
///
/// A truncated FASTQ — interrupted download, aborted rsync, full disk — was
/// read up to its last complete record and quantified silently: 8 KB of an
/// 84 KB file produced a clean `quant.sf` over 38 of 400 reads, exit 0, with
/// nothing said. That is a wrong answer presented as a right one, and salmon
/// already rejects the equivalent for its other inputs (a truncated RAD says
/// "the RAD file may be truncated" and exits 1).
///
/// Only uncompressed inputs need this. A truncated gzip/zstd/bzip2 stream
/// already fails in the decompressor ("unexpected end of file"), which covers
/// the form real sequencing data almost always arrives in. A truncated *plain*
/// file has no such envelope, so the completeness has to be checked here.
///
/// The check is exact rather than heuristic: a valid FASTQ is a whole number of
/// 4-line records, so a line count that is not a multiple of 4 is incomplete,
/// full stop. It cannot reject a valid file — in particular a final record with
/// no trailing newline still counts as a line. The cost is one sequential pass
/// over uncompressed input; large FASTQ is essentially always compressed, and
/// mapping will read the file again anyway.
fn preflight_fastq_complete(paths: &[PathBuf]) -> Result<()> {
    use std::io::Read;
    for path in paths {
        let mut f = match std::fs::File::open(path) {
            Ok(f) => f,
            // Missing/unreadable inputs are reported by the reader, with its
            // own wording; do not pre-empt it with a worse message.
            Err(_) => continue,
        };
        let mut magic = [0u8; 4];
        let n = f.read(&mut magic).unwrap_or(0);
        // Compressed: the decompressor owns this check (see above).
        let compressed = (n >= 2 && magic[0] == 0x1f && magic[1] == 0x8b)          // gzip
            || (n >= 4 && magic[..4] == [0x28, 0xb5, 0x2f, 0xfd])                  // zstd
            || (n >= 3 && &magic[..3] == b"BZh")                                   // bzip2
            || (n >= 4 && magic[..4] == [0xfd, 0x37, 0x7a, 0x58]); // xz
        if compressed || n == 0 {
            continue;
        }
        use std::io::{BufReader, Seek, SeekFrom};
        f.seek(SeekFrom::Start(0))?;
        let mut reader = BufReader::with_capacity(1 << 20, f);
        let mut buf = vec![0u8; 1 << 20];
        let (mut lines, mut last) = (0u64, 0u8);
        loop {
            let got = reader.read(&mut buf)?;
            if got == 0 {
                break;
            }
            lines += bytecount_newlines(&buf[..got]);
            last = buf[got - 1];
        }
        // A final record with no trailing newline is still a line.
        if last != b'\n' {
            lines += 1;
        }
        if !lines.is_multiple_of(4) {
            anyhow::bail!(
                "{} looks truncated: it holds {lines} lines, which is not a whole number of \
                 4-line FASTQ records. Quantifying it would silently report results for only \
                 the part that survived. Re-transfer or re-generate the file.",
                path.display()
            );
        }

        // The line count alone is not enough: a byte-level cut can land such
        // that the surviving lines still number a multiple of 4, with the last
        // one a partial quality string. So check the final record's shape,
        // which is exact and does not depend on a trailing newline (some valid
        // files have none): the separator line starts with `+`, and the
        // sequence and quality lines are the same length.
        let tail_len = 1 << 20;
        let size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        let mut f = std::fs::File::open(path)?;
        f.seek(SeekFrom::Start(size.saturating_sub(tail_len)))?;
        let mut tail = Vec::new();
        f.take(tail_len).read_to_end(&mut tail)?;
        let tail_lines: Vec<&[u8]> = tail
            .split(|&c| c == b'\n')
            .filter(|l| !l.is_empty())
            .collect();
        if tail_lines.len() >= 4 {
            let n = tail_lines.len();
            let (seq, sep, qual) = (tail_lines[n - 3], tail_lines[n - 2], tail_lines[n - 1]);
            if !sep.starts_with(b"+") || seq.len() != qual.len() {
                anyhow::bail!(
                    "{}: its final FASTQ record is incomplete — the `+` separator is \
                     missing, or the sequence and quality lines differ in length. That is \
                     what a truncated transfer looks like, but a corrupted or mis-generated \
                     file looks the same from here, so check the file rather than assuming \
                     either. Quantifying it would silently report results for only the part \
                     that parsed.",
                    path.display()
                );
            }
        }
    }
    Ok(())
}

/// Count `b'\n'` in a slice. Split out so the hot loop above stays readable.
#[inline]
fn bytecount_newlines(b: &[u8]) -> u64 {
    b.iter().filter(|&&c| c == b'\n').count() as u64
}

/// Warn, by name, about flags the selected mode cannot act on.
///
/// The rule this serves: a flag that does nothing must say so. Salmon already
/// applied it for sketch-vs-selective-alignment and for the online-only knobs;
/// the readiness sweep found it had never been written for reads-vs-`-a`/`--rad`,
/// where ~30 mapping flags were accepted in total silence (#1140).
///
/// `mode` names the selected mode as the user would recognise it, and the `why_*`
/// pair explains in one clause why the flags cannot apply — a warning that only
/// says "ignored" leaves the user unsure whether they typed the wrong flag or
/// chose the wrong mode.
///
/// Agreement is this function's job, not the callers'. It owns *both* halves:
/// the subject verb ("`--writeRad` has" / "`--x`, `--y` have") and the pronoun
/// and verb opening the reason ("it shapes" / "they shape"). Callers used to
/// pass the reason with a hard-coded "they …", so a lone flag read as
/// "`--decoyThreshold` has no effect …: they configure salmon's read mapper" —
/// the exact defect the doc above this function already claimed was handled,
/// because only the subject half had been (#1162).
///
/// `why_verb` is `Some((third-person singular, plural))` when the reason's
/// subject is the flags themselves, so its pronoun and verb have to agree with
/// how many were passed; the pair is given explicitly rather than by appending
/// an "s" to a stem, so an irregular verb cannot be silently mangled later.
/// It is `None` when the reason supplies its own subject — the mode ("it writes
/// no intermediate RAD") or a noun phrase ("posterior sampling happens …") — in
/// which case `why_rest` is used verbatim and no agreement applies.
fn warn_inert_in_mode(
    mode: &str,
    why_verb: Option<(&str, &str)>,
    why_rest: &str,
    flags: &[(&str, bool)],
) {
    if let Some(msg) = inert_in_mode_message(mode, why_verb, why_rest, flags) {
        tracing::warn!("{msg}");
    }
}

/// The message [`warn_inert_in_mode`] logs, split out so the agreement rules can
/// be asserted directly. Returns `None` when no listed flag was passed.
///
/// Kept separate because the defect this guards was invisible at the call sites:
/// the warning read correctly in every test that passed two flags and wrongly in
/// every case that passed one, and nothing asserted the singular form (#1162).
fn inert_in_mode_message(
    mode: &str,
    why_verb: Option<(&str, &str)>,
    why_rest: &str,
    flags: &[(&str, bool)],
) -> Option<String> {
    let named: Vec<&str> = flags
        .iter()
        .filter_map(|&(name, passed)| passed.then_some(name))
        .collect();
    if named.is_empty() {
        return None;
    }
    let single = named.len() == 1;
    let (subject, verb) = if single {
        (named[0].to_string(), "has")
    } else {
        (named.join(", "), "have")
    };
    let why = match why_verb {
        Some((singular_verb, plural_verb)) => {
            let (pronoun, agreed) = if single {
                ("it", singular_verb)
            } else {
                ("they", plural_verb)
            };
            format!("{pronoun} {agreed} {why_rest}")
        }
        None => why_rest.to_string(),
    };
    Some(format!(
        "{subject} {verb} no effect in {mode}: {why}. Ignored."
    ))
}

/// The reads-mode mapping/scoring/seeding knobs, paired with whether the user
/// passed each. Every one of these configures salmon's own mapper, so none can
/// apply to input that arrives already aligned (`-a`) or already mapped
/// (`--rad`). Kept in one place so the `-a` and `--rad` branches cannot drift.
fn reads_only_mapping_flags(args: &QuantArgs) -> Vec<(&'static str, bool)> {
    vec![
        ("--sketch", args.sketch),
        ("--sketchStrictOrphans", args.sketch_strict_orphans),
        ("--ma", args.ma.is_some()),
        ("--mp", args.mp.is_some()),
        ("--go", args.go.is_some()),
        ("--ge", args.ge.is_some()),
        ("--minScoreFraction", args.min_score_fraction.is_some()),
        ("--fullLengthAlignment", args.full_length_alignment),
        ("--bandwidth", args.bandwidth.is_some()),
        ("--consensusSlack", args.consensus_slack.is_some()),
        ("--maxOccsPerHit", args.max_occs_per_hit.is_some()),
        ("--maxReadOcc", args.max_read_occ.is_some()),
        ("--minAlnProb", args.min_aln_prob.is_some()),
        ("--softclip", args.softclip),
        ("--softclipOverhangs", args.softclip_overhangs),
        ("--sparseSeeds", args.sparse_seeds),
        ("--refMEMs", args.refmems),
        ("--uniMEMs", args.unimems),
        ("--mismatchSeedSkip", args.mismatch_seed_skip.is_some()),
        ("--recoverOrphans", args.recover_orphans),
        ("--allowDovetail", args.allow_dovetail),
        ("--discardOrphansQuasi", args.discard_orphans_quasi),
        (
            "--orphansRequireUnmappedMate",
            args.orphans_require_unmapped_mate,
        ),
        (
            "--preMergeChainSubThresh",
            args.pre_merge_chain_sub_thresh.is_some(),
        ),
        (
            "--postMergeChainSubThresh",
            args.post_merge_chain_sub_thresh.is_some(),
        ),
        (
            "--orphanChainSubThresh",
            args.orphan_chain_sub_thresh.is_some(),
        ),
        ("--decoyThreshold", args.decoy_threshold.is_some()),
        // NB --hardFilter and --discardOrphans are deliberately absent: both
        // are honoured on these paths (see the -a/--rad option blocks).
        ("--allowDecoyOrphans", args.allow_decoy_orphans),
        ("--decoder", args.decoder.is_some()),
        ("--threadPolicy", args.thread_policy.is_some()),
        ("--writeUnmappedNames", args.write_unmapped_names),
    ]
}

fn run_quant(args: QuantArgs, quiet: bool) -> Result<()> {
    if args.ont {
        long_read_redirect();
    }
    // Validate `-l` here, before any mode dispatch, so every path rejects a
    // bad library type identically. Reads mode parsed it with `?` deep inside
    // `quantify`, while every alignment/RAD call site discarded the error with
    // `.ok()` and quantified as though unstranded — so `-l XYZ` was a hard
    // error on one path and a silent wrong answer on the others (#1140
    // readiness sweep). One check up front is also the only way genome
    // projection gets the same treatment.
    if !salmon_core::LibraryFormat::is_auto(&args.lib_type)
        && salmon_core::LibraryFormat::parse(&args.lib_type).is_err()
    {
        anyhow::bail!(
            "invalid library type '{}'. Expected `A` (auto-detect) or a format string: \
             an optional relative orientation `I`/`O`/`M` (inward/outward/matching), then \
             strandedness `S`/`U` (stranded/unstranded), then for stranded libraries the \
             read-1 strand `F`/`R` (forward/reverse) — e.g. `IU`, `ISF`, `ISR`, `U`, `SF`, `SR`.",
            args.lib_type
        );
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
    if args.sample_out {
        tracing::warn!("--sampleOut (a BAM of posterior-sampled alignments) is accepted but not yet implemented and has no effect. --sampleUnaligned now works on its own, against --writeMappings/--writeBam.");
    }
    if args.sample_unaligned && args.write_mappings.is_none() && args.write_bam.is_none() {
        tracing::warn!("--sampleUnaligned has no effect without --writeMappings/--writeSam or --writeBam: there is no mapping output to add the unaligned records to.");
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
    // Genome projection lives inside the `-a` branch, so `--annotation` without
    // it silently discards the entire request. clap's `requires` does not fire
    // for the cases that matter (see the flag's doc comment), so check here.
    if args.annotation.is_some() && args.alignments.is_none() {
        anyhow::bail!(
            "--annotation selects genome-projection mode, which reads a genome BAM: pass it \
             together with `-a <genome.bam>`. With FASTQ or --rad input there is nothing to \
             project, and the annotation would be ignored."
        );
    }

    // `--noLengthCorrection` and bias correction cannot both apply: every
    // driver computes `bias_on = (seq||gc||pos) && !no_length_correction`, so
    // the bias request was silently discarded — and `meta_info.json` went on
    // reporting `gc_bias_correct: true` for a run that did no such thing. The
    // conflict is fully known at invocation, so say so instead (#1140).
    if args.no_length_correction && (args.seq_bias || args.gc_bias || args.pos_bias) {
        let asked: Vec<&str> = [
            ("--seqBias", args.seq_bias),
            ("--gcBias", args.gc_bias),
            ("--posBias", args.pos_bias),
        ]
        .iter()
        .filter_map(|&(n, on)| on.then_some(n))
        .collect();
        anyhow::bail!(
            "--noLengthCorrection cannot be combined with {}: bias correction works by \
             adjusting effective lengths, which is exactly what --noLengthCorrection turns \
             off. Drop whichever you did not mean.",
            asked.join(" / ")
        );
    }

    // Bias correction reads the reference sequence at each alignment's position,
    // and a BAM carries only reference names and lengths — so alignment input
    // needs `-t`. The default (deterministic) `-a` path only discovered this
    // inside `quantify_rad`, i.e. after the whole BAM had been read and written
    // to an intermediate RAD (hours on real data), and raised `--rad`'s wording
    // for an invocation containing no `--rad`. The precondition is fully known
    // at invocation, so check it before any work (#1140).
    // Genome projection (`-a --annotation`) is exempt: it reconstructs the
    // transcript sequences from `--genome` plus the annotation and never reads `-t`.
    if args.alignments.is_some()
        && args.annotation.is_none()
        && args.targets.is_none()
        && (args.seq_bias || args.gc_bias || args.pos_bias)
    {
        let asked: Vec<&str> = [
            ("--seqBias", args.seq_bias),
            ("--gcBias", args.gc_bias),
            ("--posBias", args.pos_bias),
        ]
        .iter()
        .filter_map(|&(n, on)| on.then_some(n))
        .collect();
        anyhow::bail!(
            "{} with alignment input (-a) requires -t/--targets: bias correction is computed \
             from the reference sequence at each alignment's position, and a BAM carries only \
             reference names and lengths. Pass the transcriptome FASTA the alignments were \
             produced against.",
            asked.join(" / ")
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
    // Clamp to what the machine can actually give. An unsatisfiable `-p` is an
    // ordinary user error — an unset `$SLURM_CPUS_PER_TASK`, a typo'd digit —
    // and rayon's failure for it is not a diagnosis: the thread pool aborts the
    // process with a raw panic (exit 134) whose message asserts the global pool
    // both was and was not initialized, and names neither `-p` nor the cause.
    // The threshold scales with the machine, so a small node trips far below a
    // big one. Say what happened and continue with a workable value.
    let available = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let nthreads = if nthreads > available {
        tracing::warn!(
            "-p {nthreads} exceeds the {available} hardware thread(s) this machine reports; \
             using {available}. (If this came from a scheduler variable, it was probably unset.)"
        );
        available
    } else {
        nthreads
    };
    // Propagate the failure rather than discarding it: a swallowed error here
    // resurfaces much later as SIGABRT from an unrelated-looking place.
    if let Err(e) = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build_global()
    {
        anyhow::bail!("could not create a thread pool of {nthreads} thread(s): {e}");
    }
    let out_dir = args.output.clone();
    let gene_map = load_gene_map(args.gene_map.clone(), args.ignore_tx_version)?;

    // 2.6.0 flips the default to the deterministic two-phase flow; `--online`
    // selects the pre-2.6 one-pass path for one deprecation cycle. The knobs
    // that only exist on the online path warn when passed without it, rather
    // than being silently inert under the new default — silent inertness is
    // the exact failure mode #1140's audit spent a release hunting down.
    // In RAD-input mode there is no online path at all — a RAD is always
    // quantified with the offline estimator — so `--online` there must not
    // claim the one-pass path will run, and must not suppress the inert-knob
    // warnings below (`--rad X --online -f 0.9` used to drop `-f` silently
    // behind a misleading deprecation notice, #1140).
    // `--deterministic` is the default since 2.6.0, so passing it explicitly is
    // a no-op. Say so at INFO rather than staying silent: the flag lives in
    // scripts and CI from the opt-in period, and a one-line note tells those
    // users they can drop it without implying anything is wrong. (`--online`
    // conflicts with `--deterministic` at the clap layer, so the two branches
    // are mutually exclusive.)
    if args.deterministic {
        tracing::info!(
            "--deterministic is the default since 2.6.0 and no longer needs to be passed; \
             accepting it as a no-op. Pass --online for the deprecated one-pass path \
             (removal planned for 2.7.0)."
        );
    }
    let online_path_exists = args.rad.is_none();
    if args.online && !online_path_exists {
        tracing::warn!(
            "--online has no effect in RAD-input mode (--rad): a RAD is always quantified \
             with the offline (deterministic) estimator"
        );
    }
    if args.online && online_path_exists {
        tracing::warn!(
            "--online selects the deprecated one-pass quantification path (removal planned \
             for 2.7.0). Its results depend on the thread count; the default deterministic \
             mode is byte-identical at any -p and was measured as fast or faster."
        );
    } else {
        let mut inert: Vec<&str> = Vec::new();
        if args.forgetting_factor.is_some() {
            inert.push("--forgettingFactor");
        }
        if args.num_aux_model_samples.is_some() {
            inert.push("--numAuxModelSamples");
        }
        if args.num_pre_aux_model_samples.is_some() {
            inert.push("--numPreAuxModelSamples");
        }
        if !inert.is_empty() {
            // Agree with how many were passed: one flag reads "--forgettingFactor
            // configures … to use it", several read "…, … configure … to use
            // them". Same singular/plural care as `warn_inert_in_mode`.
            let (subject, verb, pronoun) = if inert.len() == 1 {
                (inert[0].to_string(), "configures", "it")
            } else {
                (inert.join(", "), "configure", "them")
            };
            tracing::warn!(
                "{subject} {verb} the online inference path, which the default deterministic \
                 mode does not run: ignored. Pass --online (deprecated, removal planned \
                 for 2.7.0) to use {pronoun}."
            );
        }
    }
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
    // Snapshot before the mode branches move fields out of `args`; both `-a`
    // and `--rad` report the same set.
    let reads_only_flags = reads_only_mapping_flags(&args);
    if let Some(bam) = args.alignments {
        warn_inert_in_mode(
            "alignment mode (-a)",
            Some(("configures", "configure")),
            "salmon's read mapper, and these alignments were produced by another aligner",
            &reads_only_flags,
        );
        if args.index.is_some() {
            // Accepted-but-inert flags must say so (#1140): the BAM is
            // quantified as given; no index is consulted in this mode.
            tracing::warn!(
                "-i/--index is ignored in alignment mode (-a): the input alignments are \
                 quantified directly and the index is never consulted"
            );
        }
        if args.write_unmapped_names {
            // Only the reads-mode mapper tracks per-read mapping failures; an
            // alignment/RAD input carries only the records that aligned.
            tracing::warn!(
                "--writeUnmappedNames is ignored in alignment mode (-a): unmapped-read \
                 tracking exists only in reads mode"
            );
        }
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
        opts.forgetting_factor = args.forgetting_factor.unwrap_or(0.65);
        // `--meta` promises a uniform start, and this is the one path where
        // that is a real choice: the online alignment optimizer warm-starts
        // from its online estimates unless told otherwise. Every other path
        // starts uniform by construction, which is why the preset's own log
        // line was true everywhere except here (#1140, audit D15).
        opts.init_uniform = args.init_uniform || args.meta;
        opts.model_single_frag_prob = !args.no_single_frag_prob;
        // The other half of the --scoreExp policy, which this path already
        // honours; accepting one and ignoring the other left it half-applied.
        opts.hard_filter = args.hard_filter;
        opts.dump_eq = args.dump_eq;
        opts.dump_eq_weights = args.dump_eq_weights;
        opts.dump_bias_models = args.dump_bias_models;
        opts.no_length_correction = args.no_length_correction;
        opts.no_frag_length_dist = args.no_frag_length_dist;
        opts.bias_speed_samp = args.bias_speed_samp.unwrap_or(5);
        opts.num_aux_model_samples = args.num_aux_model_samples.unwrap_or(5_000_000);
        opts.no_bias_length_threshold = args.no_bias_length_threshold;
        opts.num_error_bins = args.num_error_bins.unwrap_or(4);
        opts.discard_orphans = args.discard_orphans;
        opts.gc_bins = args.num_gc_bins.unwrap_or(25);
        opts.cond_gc_bins = args.conditional_gc_bins.unwrap_or(3);
        opts.skip_quant = args.skip_quant;
        opts.num_pre_aux_model_samples = args.num_pre_aux_model_samples.unwrap_or(5_000);
        opts.num_bootstraps = args.num_bootstraps;
        opts.num_gibbs_samples = args.num_gibbs_samples;
        opts.thinning_factor = args.thinning_factor.unwrap_or(16);
        // --scoreExp scales the best-minus-score soft weight, and the RAD
        // quantifier applies it to AS-scored placements exactly as it does to
        // selective-alignment ones — so the deterministic `-a` path honours it.
        // A previous comment here claimed alignment mode "has no such term",
        // which was true only of the online path (#1140, audit C11).
        opts.score_exp = args.score_exp.unwrap_or(1.0);

        // Genome-alignment mode: `-a <genome.bam> --annotation <gtf>` projects the
        // spliced genome alignments into transcriptome coordinates (bramble) and
        // quantifies the projection. Inherently RAD-based / deterministic.
        if let Some(annotation) = args.annotation.clone() {
            // Genome projection is inherently deterministic (RAD-based) and has no
            // error model (bramble exposes no projected CIGAR), so flags that only
            // matter for transcriptomic alignment are accepted but no-ops here.
            if args.online {
                tracing::warn!(
                    "--online has no effect in genome-projection mode (--annotation): \
                     projection is inherently RAD-based and deterministic; ignoring it"
                );
            }
            if args.error_model {
                tracing::warn!(
                    "--errorModel has no effect in genome-projection mode (--annotation): \
                     projection scores by bramble similarity, not an alignment error model"
                );
            }
            let codec = rad_codec_from_args(
                args.no_compress_rad,
                args.rad_compress.as_deref().unwrap_or("lz4"),
            )?;
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
                args.rad_scratch_dir.as_deref(),
                args.keep_rad,
                gene_map.as_ref(),
            );
        }

        // Default (2.6.0): deterministic alignment quant — write the BAM's
        // placements to an intermediate RAD (baking the FLD/library-format +
        // rough abundances) and requantify from it, for results byte-identical
        // across thread counts. Score-based (`AS`); `--errorModel` opts into
        // the order-independent error model. `--online` selects the deprecated
        // one-pass path below, which is the only place the online error model
        // still lives.
        if !args.online {
            let codec = rad_codec_from_args(
                args.no_compress_rad,
                args.rad_compress.as_deref().unwrap_or("lz4"),
            )?;
            return run_deterministic_align(
                opts,
                args.write_rad.clone(),
                args.rad_scratch_dir.as_deref(),
                args.keep_rad,
                gene_map.as_ref(),
                codec,
            );
        }

        // The online alignment path never supported these two (C++ did); the
        // deterministic default honours them, so under `--online` they warn
        // rather than silently doing nothing (#1140: no parity work for a path
        // scheduled for removal).
        if args.no_length_correction || args.no_frag_length_dist {
            tracing::warn!(
                "--noLengthCorrection/--noFragLengthDist are not supported on the deprecated \
                 online alignment path: ignored. The default (deterministic) mode honours them."
            );
        }
        // The online alignment path writes no RAD, so the flags that place or
        // shape one are inert here; the deterministic default honours them
        // (#1140: a flag that does nothing must say so).
        warn_inert_in_mode(
            "the deprecated online alignment path (--online -a)",
            Some(("shapes", "shape")),
            "the intermediate RAD, which that path never writes (the default deterministic \
             mode honours them)",
            &[
                ("--writeRad", args.write_rad.is_some()),
                ("--keepRad", args.keep_rad),
                ("--radScratchDir", args.rad_scratch_dir.is_some()),
                ("--radCompress", args.rad_compress.is_some()),
                ("--noCompressRad", args.no_compress_rad),
            ],
        );
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
        // The truncated-mass warning now fires in the shared output writer,
        // covering this path and every quantify_rad driver alike (#1140).
        // `--skipQuant` estimates no abundances; aggregating the zero rows to
        // gene level would announce success over an empty result (#1140 — the
        // deterministic reads driver already skipped it, the others wrote
        // zeros silently).
        if args.skip_quant && gene_map.is_some() {
            tracing::info!(
                "--skipQuant: gene-level aggregation (-g) skipped — no abundances were estimated"
            );
        } else if let Some(gm) = &gene_map {
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
        warn_inert_in_mode(
            "RAD-input mode (--rad)",
            Some(("configures", "configure")),
            "salmon's read mapper, and a RAD already holds the mappings",
            &reads_only_flags,
        );
        if args.index.is_some() {
            // Accepted-but-inert flags must say so (#1140), same as `-a` above.
            tracing::warn!(
                "-i/--index is ignored in RAD-input mode (--rad): the RAD's own reference \
                 table is used and the index is never consulted"
            );
        }
        if args.write_unmapped_names {
            tracing::warn!(
                "--writeUnmappedNames is ignored in RAD-input mode (--rad): a RAD carries \
                 only the fragments that mapped, and no read names"
            );
        }
        // The input *is* the RAD; the flags that would place or shape a new
        // one have nothing to act on (#1140).
        warn_inert_in_mode(
            "RAD-input mode (--rad)",
            Some(("shapes", "shape")),
            "an intermediate RAD, and the --rad file is the input",
            &[
                ("--writeRad", args.write_rad.is_some()),
                ("--keepRad", args.keep_rad),
                ("--radScratchDir", args.rad_scratch_dir.is_some()),
                ("--radCompress", args.rad_compress.is_some()),
                ("--noCompressRad", args.no_compress_rad),
            ],
        );
        warn_unsupported_mapping_output(
            &args.write_mappings,
            &args.write_bam,
            "RAD-input mode (--rad)",
            "a RAD carries no read names or sequences, so SAM/BAM records cannot be reconstructed",
        );
        // Fail on an unwritable output directory *before* reading the RAD and
        // running the EM, as the reads, alignment and genome-projection drivers
        // already do. Without this the run did the whole job and then discarded
        // it: on the 26M-fragment benchmark the failure arrived at 12.74 s
        // against a 12.81 s successful run, i.e. after the entire inference,
        // and the message named no path (#1162).
        //
        // `aux_info` rather than the output directory alone, because that is the
        // one the write actually failed on. Creating only the parent would still
        // succeed for an existing-but-read-only `-o`, which is the more likely
        // way to hit this than a missing parent.
        let aux_dir = args.output.join("aux_info");
        std::fs::create_dir_all(&aux_dir).with_context(|| {
            format!(
                "creating output directory {} (needed before quantification; \
                 check that it is writable)",
                args.output.display()
            )
        })?;
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
        opts.score_exp = args.score_exp.unwrap_or(1.0);
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
        opts.model_single_frag_prob = !args.no_single_frag_prob;
        opts.hard_filter = args.hard_filter;
        // A RAD records orphan placements, so the flag that governs them means
        // something here — without this, requantifying a RAD contradicted the
        // `-a` run that wrote it (#1140).
        opts.discard_orphans = args.discard_orphans;
        opts.dump_eq = args.dump_eq;
        opts.dump_eq_weights = args.dump_eq_weights;
        opts.dump_bias_models = args.dump_bias_models;
        // The GC/bias shape knobs quantify_rad honours everywhere else; these
        // four were the exact requant_options omission recurring in this block
        // (#1140): --gcBias was wired while its tuning flags silently stayed
        // at defaults.
        opts.gc_bins = args.num_gc_bins.unwrap_or(25);
        opts.cond_gc_bins = args.conditional_gc_bins.unwrap_or(3);
        opts.bias_speed_samp = args.bias_speed_samp.unwrap_or(5);
        opts.no_bias_length_threshold = args.no_bias_length_threshold;
        opts.no_length_correction = args.no_length_correction;
        opts.no_frag_length_dist = args.no_frag_length_dist;
        opts.num_bootstraps = args.num_bootstraps;
        opts.num_gibbs_samples = args.num_gibbs_samples;
        opts.thinning_factor = args.thinning_factor.unwrap_or(16);
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
        // `--skipQuant` estimates no abundances; aggregating the zero rows to
        // gene level would announce success over an empty result (#1140 — the
        // deterministic reads driver already skipped it, the others wrote
        // zeros silently).
        if args.skip_quant && gene_map.is_some() {
            tracing::info!(
                "--skipQuant: gene-level aggregation (-g) skipped — no abundances were estimated"
            );
        } else if let Some(gm) = &gene_map {
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
    // `-2` alone is reads input with a missing half, not "no reads": the generic
    // message below sent users hunting for an input they had already typed. The
    // count-mismatch check further down only fires once `-1` is present, so this
    // case had no message of its own.
    anyhow::ensure!(
        args.mates2.is_empty() || !args.mates1.is_empty(),
        "-2 given without -1: mates must be supplied as a pair (-1 <mates_1> -2 <mates_2>). \
         For single-end reads use -r."
    );

    // Parse the read group up front: a malformed --rgLine should fail before any
    // mapping work, not after it, and the header is written the moment the
    // output file is opened.
    let read_group = args
        .rg_line
        .as_deref()
        .map(salmon_quant::mapping_record::ReadGroup::parse)
        .transpose()
        .context("parsing --rgLine")?;
    anyhow::ensure!(
        !args.mates1.is_empty() || !args.unmated.is_empty(),
        "no reads provided: pass -1/-2 (paired), -r (single-end), -a (BAM), or --rad (RAD)"
    );
    // Alignment-input knobs have nothing to act on when salmon does its own
    // mapping. Their help text says "alignment mode", but saying it only in the
    // help means a user who passes one here learns nothing.
    warn_inert_in_mode(
        "reads mode",
        Some(("configures", "configure")),
        "how salmon reads *existing* alignments, and reads mode produces its \
         own (reads mode uses --discardOrphansQuasi for the orphan policy)",
        &[
            ("--errorModel", args.error_model),
            ("--noErrorModel", args.no_error_model),
            ("--numErrorBins", args.num_error_bins.is_some()),
            ("--discardOrphans", args.discard_orphans),
        ],
    );
    // Orphans exist only relative to a mate, and dovetailing and post-merge
    // pairing describe two mates' geometry, so none of these mean anything for
    // single-end input.
    if !args.unmated.is_empty() {
        warn_inert_in_mode(
            "single-end mode (-r)",
            Some(("describes", "describe")),
            "how two mates relate, and single-end reads have no mate",
            &[
                ("--recoverOrphans", args.recover_orphans),
                ("--allowDovetail", args.allow_dovetail),
                ("--discardOrphansQuasi", args.discard_orphans_quasi),
                (
                    "--orphansRequireUnmappedMate",
                    args.orphans_require_unmapped_mate,
                ),
                (
                    "--postMergeChainSubThresh",
                    args.post_merge_chain_sub_thresh.is_some(),
                ),
            ],
        );
    }
    // Sketch mode returns only equally-best mappings, so the decoy-domination
    // comparison never runs and --decoyThreshold is inert. This lives here, not
    // before mode dispatch: with `-a`/`--rad` input `--sketch` is itself inert
    // (both are reported by `reads_only_mapping_flags`), so warning up front
    // told a user that sketch was running on a run that never mapped anything.
    if args.sketch && args.decoy_threshold.unwrap_or(1.0) != 1.0 {
        tracing::warn!(
            "--decoyThreshold {} has no effect in --sketch mode: sketch (pseudoalignment) \
             returns only equally-best mappings, so the decoy-domination comparison \
             (bestTxpScore < decoyThreshold * bestDecoyScore) never triggers. A fragment is \
             treated as decoy-dominated only when it maps to decoys *and no transcript*; \
             otherwise decoy hits are dropped and transcript hits kept (use --allowDecoyOrphans \
             to keep transcript hits even when a decoy also matches).",
            args.decoy_threshold.unwrap_or(1.0)
        );
    }
    // The one-pass reads path writes no intermediate RAD, so the flags that
    // place or shape one are inert — `-a --online` and `--rad` already say so.
    if args.online {
        warn_inert_in_mode(
            "the deprecated one-pass reads path (--online)",
            None,
            "it writes no intermediate RAD",
            &[
                ("--writeRad", args.write_rad.is_some()),
                ("--keepRad", args.keep_rad),
                ("--radScratchDir", args.rad_scratch_dir.is_some()),
                ("--radCompress", args.rad_compress.is_some()),
                ("--noCompressRad", args.no_compress_rad),
            ],
        );
    }
    // Posterior sampling happens inside the quantification pass, so
    // `--skipQuant` and it cannot both be honoured. The sibling cases
    // (`--skipQuant --dumpEq`, `--skipQuant -g`) both explain themselves;
    // these two were discarded in silence (#1140).
    if args.skip_quant {
        warn_inert_in_mode(
            "a --skipQuant run",
            None,
            "posterior sampling happens inside the quantification pass that --skipQuant \
             skips (quantify the kept RAD with --rad to draw samples from it)",
            &[
                ("--numBootstraps", args.num_bootstraps > 0),
                ("--numGibbsSamples", args.num_gibbs_samples > 0),
            ],
        );
    }
    // `--thinningFactor` thins a Gibbs chain; with no chain there is nothing to
    // thin. The analogous `--bamCompressThreads` without `--writeBam` errors,
    // and this said nothing at all.
    if args.thinning_factor.is_some() && args.num_gibbs_samples == 0 {
        warn_inert_in_mode(
            "a run without Gibbs sampling",
            Some(("thins", "thin")),
            "a Gibbs chain, which needs --numGibbsSamples",
            &[("--thinningFactor", true)],
        );
    }
    // Bias sub-knobs tune a model that only exists when bias correction is on.
    if !(args.seq_bias || args.gc_bias || args.pos_bias) {
        warn_inert_in_mode(
            "a run without bias correction",
            Some(("tunes", "tune")),
            "the bias model, which needs --seqBias, --gcBias or --posBias",
            &[
                ("--numGCBins", args.num_gc_bins.is_some()),
                ("--conditionalGCBins", args.conditional_gc_bins.is_some()),
                ("--biasSpeedSamp", args.bias_speed_samp.is_some()),
                ("--noBiasLengthThreshold", args.no_bias_length_threshold),
            ],
        );
    }
    preflight_fastq_complete(&args.mates1)?;
    preflight_fastq_complete(&args.mates2)?;
    preflight_fastq_complete(&args.unmated)?;
    anyhow::ensure!(
        args.mates1.len() == args.mates2.len(),
        "the number of -1 and -2 files must match"
    );
    let index = args.index.context("reads mode requires -i/--index")?;
    if args.targets.is_some() {
        // Accepted-but-inert flags must say so (#1140). Both reads paths ignore
        // `-t`: bias correction loads the reference sequences from the index
        // (which carries them by construction), and the alignment error model
        // it enables belongs to alignment input (`-a`).
        tracing::warn!(
            "-t/--targets is ignored in reads mode: bias correction uses the index's own \
             reference sequences, and the alignment error model applies only to -a input"
        );
    }

    let mut opts = QuantOptions::new(index, args.output);
    opts.mates1 = args.mates1;
    opts.mates2 = args.mates2;
    opts.unmated = args.unmated;
    opts.lib_type = args.lib_type;
    opts.num_threads = args.threads;
    opts.decoder = piscem_rs::io::calibrate::DecoderPreference::parse(
        args.decoder.as_deref().unwrap_or("auto"),
    )
    .map_err(|e| {
        anyhow::anyhow!(
            "--decoder {}: {e}",
            args.decoder.as_deref().unwrap_or("auto")
        )
    })?;
    opts.thread_policy = args.thread_policy.clone();
    opts.sketch = args.sketch;
    opts.sketch_strict_orphan = args.sketch_strict_orphans;
    opts.em.use_vbem = use_vbem;
    opts.range_factorization_bins = range_factorization_bins;
    opts.incompat_prior = args.incompat_prior;
    opts.dump_eq = args.dump_eq;
    opts.dump_eq_weights = args.dump_eq_weights;
    opts.write_unmapped_names = args.write_unmapped_names;
    opts.write_unaligned = args.sample_unaligned;
    opts.read_group = read_group;
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
        match args.rad_compress.as_deref().unwrap_or("lz4") {
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
    opts.thinning_factor = args.thinning_factor.unwrap_or(16);
    opts.no_length_correction = args.no_length_correction;
    opts.model_single_frag_prob = !args.no_single_frag_prob;
    opts.no_frag_length_dist = args.no_frag_length_dist;
    opts.map_config.align.min_score_fraction = args.min_score_fraction.unwrap_or(0.65);
    opts.map_config.pair.orphan_chain_sub_thresh = args.orphan_chain_sub_thresh.unwrap_or(0.0);
    opts.map_config.align.full_length_alignment = args.full_length_alignment;
    opts.map_config.align.bandwidth = args.bandwidth.unwrap_or(15);
    opts.map_config.pair.allow_dovetail = args.allow_dovetail;
    opts.map_config.pair.orphans_require_unmapped_mate = args.orphans_require_unmapped_mate;
    opts.map_config.align.softclip = args.softclip;
    opts.map_config.align.softclip_overhangs = args.softclip_overhangs;
    opts.map_config.score.score_exp = args.score_exp.unwrap_or(1.0);
    // mapping policy (Tier 2): all default to the prior hardcoded behavior
    opts.map_config.recover_orphans = args.recover_orphans;
    if args.discard_orphans_quasi {
        opts.map_config.pair.allow_orphans = false;
    }
    opts.map_config.score.decoy_thresh = args.decoy_threshold.unwrap_or(1.0);
    opts.map_config.score.min_aln_prob = args.min_aln_prob.unwrap_or(1e-5);
    opts.map_config.score.hard_filter = args.hard_filter;
    opts.map_config.score.allow_decoy_orphans = args.allow_decoy_orphans;
    // chaining sub-optimality thresholds (Tier 2)
    opts.map_config.collect.chain.chain_subopt_thresh =
        args.pre_merge_chain_sub_thresh.unwrap_or(0.8);
    opts.map_config.pair.post_merge_chain_sub_thresh =
        args.post_merge_chain_sub_thresh.unwrap_or(0.0);
    // model toggles (Tier 2)
    opts.bias_speed_samp = args.bias_speed_samp.unwrap_or(5);
    opts.num_aux_model_samples = args.num_aux_model_samples.unwrap_or(5_000_000);
    opts.no_bias_length_threshold = args.no_bias_length_threshold;
    opts.gc_bins = args.num_gc_bins.unwrap_or(25);
    opts.cond_gc_bins = args.conditional_gc_bins.unwrap_or(3);
    opts.skip_quant = args.skip_quant;
    opts.num_pre_aux_model_samples = args.num_pre_aux_model_samples.unwrap_or(5_000);
    opts.map_config.seed_mode = seed_mode(args.unimems, args.refmems, args.sparse_seeds);
    // alignment scoring (selective alignment)
    opts.map_config.align.match_score = args.ma.unwrap_or(2) as i8;
    opts.map_config.align.mismatch_pen = args.mp.unwrap_or(4) as i8;
    opts.map_config.align.gap_open_pen = args.go.unwrap_or(6) as i8;
    opts.map_config.align.gap_extend_pen = args.ge.unwrap_or(2) as i8;

    // The sketch (pseudoalignment) path computes no alignment scores, uses its
    // own hardcoded seeding, and governs orphans solely through
    // `--sketchStrictOrphans` — so every selective-alignment knob below is
    // inert under `--sketch`. A flag that does nothing must say so (#1140);
    // this is the mode×flag applicability check for the sketch/SA split.
    if args.sketch {
        let inert: Vec<&str> = [
            ("--ma", args.ma.is_some()),
            ("--mp", args.mp.is_some()),
            ("--go", args.go.is_some()),
            ("--ge", args.ge.is_some()),
            ("--minScoreFraction", args.min_score_fraction.is_some()),
            ("--fullLengthAlignment", args.full_length_alignment),
            (
                "--orphanChainSubThresh",
                args.orphan_chain_sub_thresh.is_some(),
            ),
            ("--discardOrphansQuasi", args.discard_orphans_quasi),
            (
                "--orphansRequireUnmappedMate",
                args.orphans_require_unmapped_mate,
            ),
            ("--sparseSeeds", args.sparse_seeds),
            ("--refMEMs", args.refmems),
            ("--uniMEMs", args.unimems),
            // The ten below are equally inert under --sketch — the sketch entry
            // point takes only strict_orphan/allow_dovetail/strat — but were
            // missing from this list, so they were silently accepted while
            // their neighbours warned (#1140 readiness sweep).
            ("--recoverOrphans", args.recover_orphans),
            ("--softclip", args.softclip),
            ("--softclipOverhangs", args.softclip_overhangs),
            ("--hardFilter", args.hard_filter),
            ("--minAlnProb", args.min_aln_prob.is_some()),
            ("--scoreExp", args.score_exp.is_some()),
            ("--bandwidth", args.bandwidth.is_some()),
            ("--consensusSlack", args.consensus_slack.is_some()),
            (
                "--preMergeChainSubThresh",
                args.pre_merge_chain_sub_thresh.is_some(),
            ),
            (
                "--postMergeChainSubThresh",
                args.post_merge_chain_sub_thresh.is_some(),
            ),
        ]
        .iter()
        .filter_map(|&(name, passed)| passed.then_some(name))
        .collect();
        if !inert.is_empty() {
            tracing::warn!(
                "{} {} no effect under --sketch: pseudoalignment computes no alignment \
                 scores, seeds with its own fixed sketch parameters, and controls orphans \
                 only via --sketchStrictOrphans",
                inert.join(", "),
                if inert.len() == 1 { "has" } else { "have" }
            );
        }
    } else if args.sketch_strict_orphans {
        tracing::warn!(
            "--sketchStrictOrphans has no effect without --sketch: selective alignment \
             governs orphans via --discardOrphansQuasi / --orphansRequireUnmappedMate"
        );
    }
    // chaining consensus + repetitive-hit guard
    opts.map_config.collect.consensus_fraction = 1.0 - args.consensus_slack.unwrap_or(0.35);
    opts.map_config.collect.max_hit_occ = args.max_occs_per_hit.unwrap_or(1000);
    // inference + fragment-length-distribution knobs
    opts.em.vb_prior = args.vb_prior;
    opts.em.per_nucleotide_prior = args.per_nucleotide_prior;
    opts.em.accel = args.em_accel.into();
    opts.sig_digits = args.sig_digits;
    opts.max_read_occ = args.max_read_occ.unwrap_or(250);
    opts.fld_mean = args.fld_mean.unwrap_or(DEFAULT_FLD_MEAN);
    opts.fld_sd = args.fld_sd.unwrap_or(DEFAULT_FLD_SD);
    opts.fld_max = args.fld_max.unwrap_or(DEFAULT_FLD_MAX);
    warn_unsupported_fld_policy(args.fld_policy, "reads mode");
    opts.forgetting_factor = args.forgetting_factor.unwrap_or(0.65);
    opts.init_uniform = args.init_uniform;

    // Default (2.6.0): deterministic mode — map once to an intermediate RAD,
    // then quantify from it with a fixed FLD (byte-identical output at any
    // thread count, no second mapping pass). `--online` selects the deprecated
    // one-pass path below.
    if !args.online {
        let rad_out = opts.write_rad.take(); // honour an explicit --writeRad path
                                             // The phase-1 mapping spinner, wired exactly as the one-pass path wires
                                             // its own: a default mode cannot be a silent black box for the minutes
                                             // the mapping pass takes. Handed into the driver so it stops when
                                             // phase 1 does, rather than spinning frozen through phase 2.
        let progress = Arc::new(ProgressCounters::default());
        let spinner = if !quiet && !args.no_progress && std::io::stderr().is_terminal() {
            opts.progress = Some(progress.clone());
            Some(ProgressGuard::start(progress))
        } else {
            None
        };
        return run_deterministic(
            opts,
            rad_out,
            args.rad_scratch_dir.as_deref(),
            args.keep_rad,
            gene_map.as_ref(),
            &out_dir,
            spinner,
        );
    }
    // `--initUniform` asks for a uniform optimizer start, and in one-pass reads
    // mode that is already the only initialization there is: the offline EM
    // always starts uniform (the C++ online-blended warm start was deliberately
    // not carried over — measured not to change the converged estimate, and
    // dropping it removed an order-dependent step; see the `--meta` note above,
    // which records the same fact). The same holds under `--deterministic`,
    // whose RAD quantifier likewise always starts uniform. So the flag is
    // honoured by construction in both modes and is accepted silently; the one
    // path where it toggles anything is online alignment mode (`-a`), whose
    // offline EM otherwise warm-starts from the online estimates.

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
    // `--skipQuant` estimates no abundances; aggregating the zero rows to
    // gene level would announce success over an empty result (#1140 — the
    // deterministic reads driver already skipped it, the others wrote
    // zeros silently).
    if args.skip_quant && gene_map.is_some() {
        tracing::info!(
            "--skipQuant: gene-level aggregation (-g) skipped — no abundances were estimated"
        );
    } else if let Some(gm) = &gene_map {
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

    /// One inert flag must not be described in the plural.
    ///
    /// The warning read "`--decoyThreshold` has no effect in RAD-input mode
    /// (--rad): **they configure** salmon's read mapper" — the subject agreed
    /// while the reason did not. It survived because the helper's own doc
    /// comment claimed singular and plural were handled, and every test that
    /// exercised it happened to pass more than one flag (#1162).
    #[test]
    fn inert_flag_reason_agrees_with_how_many_were_passed() {
        let verb = Some(("configures", "configure"));
        let rest = "salmon's read mapper, and a RAD already holds the mappings";

        let one = inert_in_mode_message(
            "RAD-input mode (--rad)",
            verb,
            rest,
            &[("--decoyThreshold", true), ("--sketch", false)],
        )
        .expect("one flag was passed");
        assert!(
            one.starts_with("--decoyThreshold has no effect"),
            "subject should be singular: {one}"
        );
        assert!(
            one.contains("it configures salmon's read mapper"),
            "reason should be singular too: {one}"
        );
        assert!(!one.contains("they"), "no plural anywhere: {one}");

        let many = inert_in_mode_message(
            "RAD-input mode (--rad)",
            verb,
            rest,
            &[("--sketch", true), ("--decoyThreshold", true)],
        )
        .expect("two flags were passed");
        assert!(
            many.starts_with("--sketch, --decoyThreshold have no effect"),
            "subject should be plural: {many}"
        );
        assert!(
            many.contains("they configure salmon's read mapper"),
            "reason should be plural: {many}"
        );
    }

    /// A reason that supplies its own subject is used verbatim: forcing a
    /// pronoun onto it would produce "it it writes no intermediate RAD".
    #[test]
    fn inert_reason_with_its_own_subject_is_left_alone() {
        let msg = inert_in_mode_message(
            "the deprecated one-pass reads path (--online)",
            None,
            "it writes no intermediate RAD",
            &[("--keepRad", true)],
        )
        .expect("one flag was passed");
        assert!(
            msg.ends_with("(--online): it writes no intermediate RAD. Ignored."),
            "reason should be verbatim: {msg}"
        );
    }

    /// Nothing passed, nothing said.
    #[test]
    fn inert_message_is_silent_when_no_flag_was_passed() {
        assert!(inert_in_mode_message(
            "RAD-input mode (--rad)",
            Some(("configures", "configure")),
            "salmon's read mapper",
            &[("--sketch", false)],
        )
        .is_none());
    }

    /// `-a` with read inputs must be a parse error, not a silent choice of the
    /// BAM: the pre-#1140 behavior quantified the alignments and dropped the
    /// FASTQs without a word.
    #[test]
    fn alignments_conflict_with_read_inputs() {
        for reads in [vec!["-1", "r1.fq", "-2", "r2.fq"], vec!["-r", "r.fq"]] {
            let mut argv = vec!["salmon", "quant", "-l", "A", "-o", "out", "-a", "aln.bam"];
            argv.extend(reads.iter());
            let err = match Cli::try_parse_from(&argv) {
                Err(e) => e,
                Ok(_) => panic!("-a plus reads must not parse: {argv:?}"),
            };
            assert_eq!(err.kind(), clap::error::ErrorKind::ArgumentConflict);
        }
        // ...and -a alone still parses.
        assert!(
            Cli::try_parse_from(["salmon", "quant", "-l", "A", "-o", "out", "-a", "aln.bam"])
                .is_ok()
        );
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

    use std::path::Path;

    /// Every knob phase 2 can honour has to survive the hand-off, and the only
    /// way to be sure is to check them one by one.
    ///
    /// A forgotten field is invisible at runtime: the requant runs with the
    /// default, the run succeeds, and the flag the user passed did nothing.
    /// `--biasSpeedSamp` and `--noBiasLengthThreshold` were dropped exactly that
    /// way (COMBINE-lab/salmon#1140), and `--initUniform` never reached reads
    /// mode at all. Each value below is deliberately different from the
    /// `AlignQuantOptions` default, so a field that is not forwarded fails here
    /// rather than merely looking plausible.
    #[test]
    fn deterministic_requant_forwards_every_supported_knob() {
        let mut map_opts = QuantOptions::new(PathBuf::from("idx"), PathBuf::from("out"));
        map_opts.lib_type = "ISR".to_string();
        map_opts.em.vb_prior = 0.017;
        map_opts.em.use_vbem = false;
        map_opts.range_factorization_bins = 7;
        map_opts.seq_bias = true;
        map_opts.gc_bias = true;
        map_opts.pos_bias = true;
        map_opts.bias_seed_em_iters = 13;
        map_opts.map_config.score.score_exp = 1.75;
        map_opts.incompat_prior = 1e-7;
        map_opts.sig_digits = 5;
        map_opts.fld_mean = 271.0;
        map_opts.fld_sd = 33.0;
        map_opts.fld_max = 777;
        map_opts.gc_bins = 17;
        map_opts.cond_gc_bins = 5;
        map_opts.num_bootstraps = 11;
        map_opts.num_gibbs_samples = 23;
        map_opts.thinning_factor = 9;
        map_opts.bias_speed_samp = 25;
        map_opts.no_bias_length_threshold = true;
        map_opts.init_uniform = true;
        map_opts.dump_eq = true;
        map_opts.dump_eq_weights = true;
        map_opts.dump_bias_models = true;
        map_opts.model_single_frag_prob = false;
        map_opts.no_length_correction = true;
        map_opts.no_frag_length_dist = true;

        let q = requant_options(
            &map_opts,
            Path::new("inter.rad"),
            Path::new("out"),
            Some(PathBuf::from("targets.fa")),
        );

        assert_eq!(q.bam, PathBuf::from("inter.rad"));
        assert_eq!(q.output_dir, PathBuf::from("out"));
        assert_eq!(q.lib_type, "ISR");
        assert_eq!(q.em.vb_prior, 0.017);
        assert!(!q.em.use_vbem);
        assert_eq!(q.range_factorization_bins, 7);
        assert_eq!(q.transcripts, Some(PathBuf::from("targets.fa")));
        assert!(q.seq_bias && q.gc_bias && q.pos_bias);
        assert_eq!(q.bias_seed_em_iters, 13);
        assert_eq!(q.score_exp, 1.75);
        assert_eq!(q.incompat_prior, 1e-7);
        assert_eq!(q.sig_digits, 5);
        assert_eq!(q.fld_mean, 271.0);
        assert_eq!(q.fld_sd, 33.0);
        assert_eq!(q.fld_max, 777);
        assert_eq!(q.gc_bins, 17);
        assert_eq!(q.cond_gc_bins, 5);
        assert_eq!(q.num_bootstraps, 11);
        assert_eq!(q.num_gibbs_samples, 23);
        assert_eq!(q.thinning_factor, 9);
        assert_eq!(q.bias_speed_samp, 25);
        assert!(q.no_bias_length_threshold);
        assert!(q.init_uniform);
        assert!(q.dump_eq);
        assert!(q.dump_eq_weights);
        assert!(q.dump_bias_models);
        assert!(!q.model_single_frag_prob);
        assert!(q.no_length_correction);
        assert!(q.no_frag_length_dist);
    }

    /// The three knobs #1140 found dropped, checked against a default request
    /// rather than against each other: with none of them set, phase 2 must see
    /// the `AlignQuantOptions` defaults, so the test above is measuring
    /// forwarding and not a coincidence.
    #[test]
    fn deterministic_requant_leaves_unset_knobs_at_their_defaults() {
        let map_opts = QuantOptions::new(PathBuf::from("idx"), PathBuf::from("out"));
        let defaults = AlignQuantOptions::new(PathBuf::from("inter.rad"), PathBuf::from("out"));
        let q = requant_options(&map_opts, Path::new("inter.rad"), Path::new("out"), None);
        assert_eq!(q.bias_speed_samp, defaults.bias_speed_samp);
        assert_eq!(
            q.no_bias_length_threshold,
            defaults.no_bias_length_threshold
        );
        assert_eq!(q.init_uniform, defaults.init_uniform);
        assert_eq!(q.dump_eq, defaults.dump_eq);
        assert_eq!(q.dump_eq_weights, defaults.dump_eq_weights);
        assert_eq!(q.no_length_correction, defaults.no_length_correction);
        assert_eq!(q.no_frag_length_dist, defaults.no_frag_length_dist);
    }

    /// `--rgLine` reaches the parser as a plain string, so the shell-friendly
    /// escaped-tab spelling has to survive to where it is parsed rather than
    /// being mangled by clap.
    #[test]
    fn rg_line_is_taken_verbatim() {
        let args = quant_args(&["--rgLine", "ID:s1\\tSM:sample"]).unwrap();
        assert_eq!(args.rg_line.as_deref(), Some("ID:s1\\tSM:sample"));
    }

    /// `--sampleUnaligned` keeps its salmon short form, since it is the flag
    /// users already have in their scripts.
    #[test]
    fn sample_unaligned_has_both_spellings() {
        for flag in ["--sampleUnaligned", "-u"] {
            assert!(quant_args(&[flag]).unwrap().sample_unaligned, "{flag}");
        }
        assert!(!quant_args(&[]).unwrap().sample_unaligned);
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

    /// `--geneMap` is read before any quantification work, so an unusable
    /// annotation fails immediately rather than after mapping, inference and
    /// `quant.sf` have all completed (#1074).
    #[test]
    fn gene_map_is_validated_before_the_run() {
        let dir = tempfile::tempdir().unwrap();

        // Absent option: nothing to load, and nothing to complain about.
        assert!(load_gene_map(None, false).unwrap().is_none());

        // A path that does not exist fails here, where no work has been done yet.
        let missing = dir.path().join("nope.gtf");
        let err = load_gene_map(Some(missing.clone()), false).unwrap_err();
        assert!(
            format!("{err:#}").contains("nope.gtf"),
            "the error should name the annotation, got: {err:#}"
        );

        // A usable annotation is parsed once, up front, and carried on the opts.
        let gtf = dir.path().join("genes.gtf");
        std::fs::write(
            &gtf,
            "chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
             chr1\tsrc\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";\n",
        )
        .unwrap();
        let opts = load_gene_map(Some(gtf.clone()), true).unwrap().unwrap();
        assert_eq!(opts.path, gtf);
        assert!(opts.ignore_tx_version);
        assert_eq!(opts.map.get("T1").map(String::as_str), Some("G1"));
        assert_eq!(opts.map.get("T2").map(String::as_str), Some("G1"));
    }

    /// The salmon-owned intermediate goes to the scratch dir (pid-suffixed so
    /// concurrent runs sharing one scratch volume cannot collide) when given,
    /// else under the output directory with the historical fixed name.
    #[test]
    fn intermediate_rad_path_placement() {
        let tmp = tempfile::tempdir().unwrap();
        let out = tmp.path().join("out");
        let scratch = tmp.path().join("scratch");

        let p = intermediate_rad_path(&out, None, "intermediate_mappings").unwrap();
        assert_eq!(p, out.join("intermediate_mappings.rad"));

        let p = intermediate_rad_path(&out, Some(&scratch), "intermediate_mappings").unwrap();
        assert!(scratch.is_dir(), "scratch dir must be created");
        assert_eq!(p.parent().unwrap(), scratch);
        let name = p.file_name().unwrap().to_str().unwrap();
        assert!(
            name.starts_with("intermediate_mappings-") && name.ends_with(".rad"),
            "pid-suffixed name, got {name}"
        );
    }

    /// After a failed pass, cleanup removes the finalized intermediate and any
    /// `.partial-*` sibling, and touches nothing else in the directory.
    #[test]
    fn cleanup_intermediate_rad_removes_file_and_partials() {
        let tmp = tempfile::tempdir().unwrap();
        let rad = tmp.path().join("intermediate_mappings.rad");
        let partial = tmp.path().join("intermediate_mappings.rad.partial-12345");
        let bystander = tmp.path().join("quant.sf");
        std::fs::write(&rad, b"x").unwrap();
        std::fs::write(&partial, b"x").unwrap();
        std::fs::write(&bystander, b"x").unwrap();

        cleanup_intermediate_rad(&rad);
        assert!(!rad.exists());
        assert!(!partial.exists());
        assert!(bystander.exists());

        // Idempotent on an already-clean directory.
        cleanup_intermediate_rad(&rad);
    }
}
