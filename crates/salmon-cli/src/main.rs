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

use salmon_align::{quantify_alignments, AlignQuantOptions};
use salmon_index::{build as build_index, IndexBuildOptions};
use salmon_quant::{quantify, QuantOptions};

#[derive(Parser)]
#[command(name = "salmon", version, about = "RNA-seq transcript quantification (Rust port)")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Build a salmon index from a transcriptome FASTA.
    Index(IndexArgs),
    /// Quantify transcript abundances from FASTQ reads.
    Quant(QuantArgs),
    /// Diagnostic: per-read best-mapping detail (placement, seed coverage, score).
    DebugMap(DebugMapArgs),
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
}

#[derive(Args)]
struct QuantArgs {
    /// Salmon index directory (reads mode).
    #[arg(short = 'i', long = "index", required_unless_present = "alignments")]
    index: Option<PathBuf>,
    /// Alignment-based mode: a BAM of reads aligned to the transcriptome.
    #[arg(short = 'a', long = "alignments")]
    alignments: Option<PathBuf>,
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
    /// Worker threads (0 = all cores).
    #[arg(short = 'p', long = "threads", default_value_t = 0)]
    threads: usize,
    /// Use the alignment-free pseudoalignment path.
    #[arg(long = "sketch")]
    sketch: bool,
    /// Use the standard EM optimizer instead of VBEM.
    #[arg(long = "useEM")]
    use_em: bool,
    /// Range-factorization bins for equivalence classes (0 disables).
    #[arg(long = "rangeFactorizationBins", default_value_t = 4)]
    range_factorization_bins: u32,
    /// Prior weight for strand-incompatible mappings (0 drops them).
    #[arg(long = "incompatPrior", default_value_t = 0.0)]
    incompat_prior: f64,
    /// Dump naive equivalence classes to aux_info/eq_classes.txt.
    #[arg(long = "dumpEq")]
    dump_eq: bool,
    /// Enable sequence-specific bias correction.
    #[arg(long = "seqBias")]
    seq_bias: bool,
    /// Enable fragment-GC bias correction.
    #[arg(long = "gcBias")]
    gc_bias: bool,
    /// Enable positional bias correction.
    #[arg(long = "posBias")]
    pos_bias: bool,
    /// Minimum alignment score as a fraction of the perfect score.
    #[arg(long = "minScoreFraction", default_value_t = 0.65)]
    min_score_fraction: f32,
    /// Score the full read with one DP instead of PuffAligner-style inter-MEM-gap scoring.
    #[arg(long = "fullLengthAlignment")]
    full_length_alignment: bool,
    /// Seed with reference-extended MEMs (cross unitig boundaries).
    #[arg(long = "refMEMs", conflicts_with = "unimems")]
    refmems: bool,
    /// Seed with true unitig-constrained uni-MEMs (pufferfish-style).
    #[arg(long = "uniMEMs", conflicts_with = "refmems")]
    unimems: bool,
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
    let mut opts = IndexBuildOptions::new(args.transcripts, args.index);
    opts.k = args.kmer_len;
    opts.m = args.minimizer_len;
    opts.threads = args.threads;
    opts.keep_intermediate = args.keep_intermediate;
    let info = build_index(&opts).context("index build failed")?;
    println!("indexed {} references (k={}, m={})", info.num_refs, info.k, info.m);
    Ok(())
}

fn run_quant(args: QuantArgs) -> Result<()> {
    // Alignment-based mode: quantify directly from a BAM.
    if let Some(bam) = args.alignments {
        let mut opts = AlignQuantOptions::new(bam, args.output);
        opts.lib_type = args.lib_type;
        opts.em.use_vbem = !args.use_em;
        opts.range_factorization_bins = args.range_factorization_bins;
        let res = quantify_alignments(&opts).context("alignment-based quantification failed")?;
        let pct = if res.num_processed > 0 {
            100.0 * res.num_mapped as f64 / res.num_processed as f64
        } else {
            0.0
        };
        println!(
            "processed {} fragments, mapped {} ({:.2}%), {} equivalence classes",
            res.num_processed, res.num_mapped, pct, res.num_eq_classes
        );
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
    opts.em.use_vbem = !args.use_em;
    opts.range_factorization_bins = args.range_factorization_bins;
    opts.incompat_prior = args.incompat_prior;
    opts.dump_eq = args.dump_eq;
    opts.seq_bias = args.seq_bias;
    opts.gc_bias = args.gc_bias;
    opts.pos_bias = args.pos_bias;
    opts.map_config.align.min_score_fraction = args.min_score_fraction;
    opts.map_config.align.full_length_alignment = args.full_length_alignment;
    opts.map_config.seed_mode = seed_mode(args.unimems, args.refmems);

    let res = quantify(&opts).context("quantification failed")?;
    let pct = if res.num_processed > 0 {
        100.0 * res.num_mapped as f64 / res.num_processed as f64
    } else {
        0.0
    };
    println!(
        "processed {} fragments, mapped {} ({:.2}%), {} equivalence classes",
        res.num_processed, res.num_mapped, pct, res.num_eq_classes
    );
    Ok(())
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .with_writer(std::io::stderr)
        .init();

    let cli = Cli::parse();
    match cli.command {
        Command::Index(args) => run_index(args),
        Command::Quant(args) => run_quant(args),
        Command::DebugMap(args) => run_debug_map(args),
    }
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
        let name = header.trim_start_matches('@').split_whitespace().next().unwrap_or("");
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
