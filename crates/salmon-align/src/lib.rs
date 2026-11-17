//! `salmon-align`: alignment-based (BAM) quantification.
//!
//! # When you would use this
//!
//! salmon's usual input is raw reads, which it maps itself. But sometimes the
//! alignments already exist: another tool produced them, a pipeline requires a
//! specific aligner, or the reads were aligned to the genome and projected. In
//! those cases the expensive half of the work is already done, and salmon only
//! needs to do the statistics.
//!
//! So this crate is the mapping stage replaced by a reader. Everything after it —
//! equivalence classes, the fragment-length distribution, bias models, the EM,
//! the output files — is the same code the reads path uses. The three entry
//! points differ only in where the placements come from:
//!
//! * a transcriptome BAM (this module),
//! * a genome BAM plus an annotation, projected into transcript coordinates
//!   (the `genome_project` module),
//! * a RAD file of mappings (the `rad` module).
//!
//! The BAM must be *name-grouped*: all records for one read adjacent, so a
//! fragment's placements can be gathered in one pass without holding the whole
//! file in memory.
//!
//! The alternative to mapping FASTQ reads: take a BAM of reads already aligned
//! to the *transcriptome* (each `@SQ` is a transcript) and quantify directly
//! from those alignments. References and their lengths come from the BAM
//! header, so no index or FASTA is required for basic quantification.
//!
//! For each fragment (the contiguous run of records sharing a read name) the
//! set of transcripts it aligns to becomes an equivalence class, weighted by
//! alignment score (`AS`) via the same `exp(-scoreExp·(best−score))` rule used
//! in the mapping path. Fragment lengths (`TLEN`) feed the fragment-length
//! distribution. The class set then flows through the shared EM
//! ([`salmon_infer`]) to `quant.sf`. Mirrors salmon's `quant -a` mode (the
//! position-binned alignment error model is a later refinement).

mod bam_rad;
mod error_model;
mod genome_project;
mod rad;

pub use bam_rad::{write_alignment_rad, AlignRadSummary};
pub use genome_project::{project_genome_bam_to_rad, GenomeProjectOptions, ProjectionArtifacts};
pub use rad::quantify_rad;

use std::collections::HashMap;
use std::io::{BufRead, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_sam::Header;
use rayon::prelude::*;
use serde::Serialize;

use error_model::AlnOp;
use salmon_core::{LibraryFormat, MateStatus, ReadOrientation, ReadStrandedness, ReadType};
use salmon_infer::EmOptions;
use salmon_model::FragmentLengthDistribution;

const SALMON_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Where the fragment-length distribution comes from when quantifying a RAD
/// (`--fldPolicy`).
///
/// salmon's own RAD writer always bakes its FLD into the header so a requant
/// reproduces the writing run exactly. That default is right for reproduction
/// but leaves no way to ask a different question of the same file, which is
/// what the other two variants restore.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum FldPolicy {
    /// Use the FLD baked into the RAD header when present (exact parity with
    /// the run that wrote the file); otherwise derive it, or fall back to the
    /// prior for single-end libraries. The default.
    #[default]
    Baked,
    /// Ignore any baked FLD and re-derive it from this RAD's uniquely-mapped
    /// proper pairs, seeded by `--fldMean`/`--fldSD`.
    Derive,
    /// Ignore both the baked FLD and the observations; use `--fldMean`/`--fldSD`
    /// alone. This is the only setting under which those flags fully determine
    /// the distribution, so it is the one to use for a fragment-length
    /// sensitivity analysis.
    Prior,
}

/// Which fragment-length prior flags the user supplied explicitly, as opposed
/// to inheriting their defaults.
///
/// salmon cannot warn "your `--fldMean` was ignored" without knowing the user
/// actually typed it — every run has *some* value for these.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ExplicitFldArgs {
    pub mean: bool,
    pub sd: bool,
    pub max: bool,
}

impl ExplicitFldArgs {
    /// True when at least one fragment-length prior flag was supplied.
    pub fn any(self) -> bool {
        self.mean || self.sd || self.max
    }

    /// The supplied flags, formatted for a message: `--fldMean/--fldSD`.
    pub fn names(self) -> String {
        let mut names = Vec::new();
        if self.mean {
            names.push("--fldMean");
        }
        if self.sd {
            names.push("--fldSD");
        }
        if self.max {
            names.push("--fldMax");
        }
        names.join("/")
    }
}

/// Options for alignment-based quantification.
#[derive(Debug, Clone)]
pub struct AlignQuantOptions {
    /// BAM of alignments to the transcriptome (records grouped by read name)
    pub bam: PathBuf,
    /// output directory
    pub output_dir: PathBuf,
    /// library type string (recorded in output)
    pub lib_type: String,
    /// EM/VBEM options
    pub em: EmOptions,
    /// range-factorization bins (0 disables)
    pub range_factorization_bins: u32,
    /// soft-weight decay applied to alignment-score differences
    pub score_exp: f64,
    /// transcriptome FASTA (`-t`); required to train the alignment error model
    pub transcripts: Option<PathBuf>,
    /// reference sequences (transcript-id order, decoys included) supplied
    /// directly instead of loading `transcripts` — used by `--deterministic`,
    /// which already has the index and can hand its sequences to the requant
    /// pass so bias correction needs no separate `-t`. Takes precedence over
    /// [`Self::transcripts`] when set.
    pub ref_seqs: Option<salmon_core::RefSeqs>,
    /// disable the alignment error model (salmon's `--noErrorModel`)
    pub no_error_model: bool,
    /// opt in to the order-independent error model in `--deterministic` alignment
    /// mode (`--errorModel`). Off by default: deterministic mode scores by the BAM
    /// `AS` tag, which benchmarks at least as accurately against truth and runs a
    /// single BAM pass; enabling this trains the model in a second BAM pass. Only
    /// consulted by the `--deterministic` BAM→RAD producer, and needs `-t`.
    pub deterministic_error_model: bool,
    /// enable sequence-specific bias correction (`--seqBias`)
    pub seq_bias: bool,
    /// enable fragment-GC bias correction (`--gcBias`)
    pub gc_bias: bool,
    /// enable positional bias correction (`--posBias`)
    pub pos_bias: bool,
    /// EM iterations for the rough abundance seed that weights bias collection
    /// when no baked prior is present (`--biasSeedEMIters`, hidden); default 11.
    pub bias_seed_em_iters: u32,
    /// weight multiplier for orientation-incompatible alignments; `0` drops them
    /// (salmon's default `ignoreIncompat` behavior)
    pub incompat_prior: f64,
    /// fragment-length distribution prior mean, SD, and max tracked length
    /// (`--fldMean` / `--fldSD` / `--fldMax`)
    pub fld_mean: f64,
    pub fld_sd: f64,
    pub fld_max: usize,
    /// RAD-input mode: where the fragment-length distribution comes from
    /// (`--fldPolicy`). Ignored outside `--rad`, which has no baked FLD to
    /// prefer.
    pub fld_policy: FldPolicy,
    /// which of `--fldMean`/`--fldSD`/`--fldMax` the user actually supplied.
    /// Used only to warn when a baked FLD makes them ineffective, so that a
    /// default-valued run stays quiet.
    pub explicit_fld_args: ExplicitFldArgs,
    /// `--fldMean`/`--fldSD` were supplied by the user -- to this run or, for a
    /// deterministic run's phase 2, to the mapping pass that baked the RAD.
    /// When `false` and the library is single-end, the FLD is only salmon's
    /// default guess, and the RAD quantifier widens it for the single-end
    /// length term.
    pub fld_prior_user_supplied: bool,
    /// wall-clock seconds a driver already spent before this call (phase 1 of
    /// `--deterministic`, or genome projection); added to the reported
    /// `total_time_seconds` so `meta_info.json` covers the whole run rather
    /// than only the quantification pass
    pub prior_seconds: f64,
    /// start time (asctime) of that earlier phase; when set, reported as the
    /// run's `start_time` instead of this pass's own
    pub external_start_time: Option<String>,
    /// keep the `cmd_info.json` an earlier phase already wrote instead of
    /// overwriting it: reads-mode `--deterministic` phase 1 records the real
    /// invocation (index, read files, threads); phase 2's own record would
    /// replace it with a RAD-centric one and lose all of that (#1140)
    pub preserve_cmd_info: bool,
    /// this run converted alignments into the RAD it is quantifying (`-a`
    /// `--deterministic`, genome projection) rather than being handed one.
    /// Only metadata depends on it: `frag_length_source` must say `alignments`
    /// — the distribution came from the input alignments — instead of leaking
    /// `rad_baked`, which names an intermediate the user never asked for. It
    /// cannot be inferred: `preserve_cmd_info` marks the *reads* driver, and a
    /// standalone `--rad` (where `rad_baked` is correct) looks identical here
    /// otherwise. Left false by anything reading a user-supplied RAD.
    pub input_is_alignments: bool,
    /// the read files behind this RAD, when the driver knows them (reads-mode
    /// `--deterministic`), for `lib_format_counts.json`'s `read_files` — which
    /// otherwise reports `[]` because a RAD does not name its reads (#1140)
    pub read_files: Vec<String>,
    /// `--dumpEq`: write `aux_info/eq_classes.txt.gz`, the equivalence classes
    /// the EM consumed, collapsed by transcript set.
    pub dump_eq: bool,
    /// `--dumpEqWeights`: the same file with each class's per-transcript
    /// combined weights, and without collapsing range-factorized sub-classes.
    pub dump_eq_weights: bool,
    /// `--noLengthCorrection`: give every transcript the same constant effective
    /// length ([`salmon_model::NO_LENGTH_CORRECTION_EFF_LEN`], as salmon 1.x
    /// did), instead of the fragment-length-derived one.
    ///
    /// For protocols whose read count does not grow with transcript length (3'
    /// tagged-end libraries above all), any length-dependent divisor, the raw
    /// length included, puts back a dependence the data does not have (#1243).
    /// Bias correction is disabled with it, for the same reason: its whole job
    /// is to adjust an effective length this mode is not computing.
    pub no_length_correction: bool,
    /// `--noEffectiveLengthCorrection`: divide by the raw reference length
    /// instead of the fragment-length-derived effective length (C++ salmon's
    /// flag of the same name). `no_length_correction` wins if both are set.
    /// Bias correction is disabled with it too.
    pub no_effective_length_correction: bool,
    /// `--noFragLengthDist`: drop the fragment-length term from a placement's
    /// weight, leaving score and compatibility to decide it.
    pub no_frag_length_dist: bool,
    /// `--hardFilter`: keep only a fragment's best-scoring placement(s) at
    /// weight 1 instead of soft-weighting the rest by
    /// `exp(-scoreExp * (best - score))`. The two are halves of one policy, and
    /// `--scoreExp` was already honoured here, so accepting `--hardFilter` and
    /// ignoring it left that policy half-implemented (#1140 readiness sweep).
    pub hard_filter: bool,
    /// `--noSingleFragProb` (inverted): model an orphan's fragment length with
    /// the bounded-CMF ambiguous weight. When false, orphans in a paired
    /// library take a flat `LOG_EPSILON` penalty instead, matching the reads
    /// path (#1140, D14).
    pub model_single_frag_prob: bool,
    /// diagnostics born in an earlier phase the driver ran (e.g. the
    /// deterministic alignment pass's unusable-`AS` verdict), appended to this
    /// pass's own in `meta_info.json` — so a pipeline that only keeps the
    /// output directory sees them, not just whoever read the log (#1140)
    pub extra_diagnostics: Vec<salmon_core::Diagnostic>,
    /// `--dumpBiasModels` (hidden): write the observed+expected seq/GC/pos
    /// bias models to `bias_models.txt` for debugging / C++-parity comparison.
    /// Written by the shared output writer from the pass's [`BiasDump`](salmon_model::dumps::BiasDump), so it
    /// works on every path that runs bias correction (#1140: it used to exist
    /// only on the deprecated one-pass reads path).
    pub dump_bias_models: bool,
    /// significant digits for the EffectiveLength and NumReads columns of
    /// `quant.sf` (`--sigDigits`, salmon default 3)
    pub sig_digits: u32,
    /// number of read-position bins in the alignment error model
    /// (`--numErrorBins`, salmon default 4)
    pub num_error_bins: usize,
    /// drop orphan (single-mate) placements in a paired library instead of
    /// fragment-length-penalizing them (`--discardOrphans`, alignment mode)
    pub discard_orphans: bool,
    /// fragment-length sampling stride for the GC bias convolution
    /// (`--biasSpeedSamp`, default 5)
    pub bias_speed_samp: usize,
    /// disable the lower barrier on bias-corrected effective lengths
    /// (`--noBiasLengthThreshold`)
    pub no_bias_length_threshold: bool,
    /// GC bias model bin counts (`--numGCBins` × `--conditionalGCBins`, 25×3)
    pub gc_bins: usize,
    pub cond_gc_bins: usize,
    /// skip abundance estimation + `quant.sf`, emitting only eq-classes,
    /// library type, and metadata (salmon's `--skipQuant`)
    pub skip_quant: bool,
    /// number of posterior bootstrap replicates to draw (`--numBootstraps`, 0 = off)
    pub num_bootstraps: u32,
    /// number of posterior Gibbs samples to draw (`--numGibbsSamples`, 0 = off);
    /// ignored when `num_bootstraps > 0` (bootstrap takes precedence, like reads mode)
    pub num_gibbs_samples: u32,
    /// Gibbs thinning factor (`--thinningFactor`, salmon default 16)
    pub thinning_factor: u32,
    /// Optional shared progress counters. When `Some`, the BAM pass reports
    /// processed/mapped fragment counts here as it runs so the caller can drive
    /// a live progress display. `None` (the default) disables sharing.
    pub progress: Option<std::sync::Arc<salmon_core::ProgressCounters>>,
}

impl AlignQuantOptions {
    /// The effective-length mode `--noLengthCorrection` and
    /// `--noEffectiveLengthCorrection` select.
    pub fn length_correction(&self) -> salmon_model::LengthCorrection {
        salmon_model::LengthCorrection::from_flags(
            self.no_length_correction,
            self.no_effective_length_correction,
        )
    }

    pub fn new(bam: PathBuf, output_dir: PathBuf) -> Self {
        Self {
            bam,
            output_dir,
            lib_type: "IU".to_string(),
            em: EmOptions::default(),
            range_factorization_bins: 4,
            score_exp: 1.0,
            transcripts: None,
            ref_seqs: None,
            no_error_model: false,
            deterministic_error_model: false,
            seq_bias: false,
            gc_bias: false,
            pos_bias: false,
            bias_seed_em_iters: 11,
            incompat_prior: 0.0,
            fld_mean: 250.0,
            fld_sd: 25.0,
            fld_max: 1000,
            fld_policy: FldPolicy::default(),
            explicit_fld_args: ExplicitFldArgs::default(),
            fld_prior_user_supplied: false,
            prior_seconds: 0.0,
            external_start_time: None,
            preserve_cmd_info: false,
            input_is_alignments: false,
            read_files: Vec::new(),
            dump_eq: false,
            dump_eq_weights: false,
            no_length_correction: false,
            no_effective_length_correction: false,
            no_frag_length_dist: false,
            hard_filter: false,
            model_single_frag_prob: true,
            extra_diagnostics: Vec::new(),
            dump_bias_models: false,
            sig_digits: 3,
            num_error_bins: 4,
            discard_orphans: false,
            bias_speed_samp: 5,
            no_bias_length_threshold: false,
            gc_bins: salmon_model::gcbias::DEFAULT_GC_BINS,
            cond_gc_bins: salmon_model::gcbias::DEFAULT_COND_BINS,
            skip_quant: false,
            num_bootstraps: 0,
            num_gibbs_samples: 0,
            thinning_factor: 16,
            progress: None,
        }
    }
}

/// Quantification results.
/// What a quantification run read its fragments from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FragmentSource {
    /// an alignment file given with `-a`
    Bam,
    /// a RAD file given with `--rad`, or written as a `--deterministic` intermediate
    Rad,
}

#[derive(Debug, Clone)]
pub struct AlignQuantResult {
    pub names: Vec<String>,
    pub lengths: Vec<u32>,
    pub eff_lengths: Vec<f64>,
    pub tpm: Vec<f64>,
    pub counts: Vec<f64>,
    pub num_processed: u64,
    pub num_mapped: u64,
    /// mapped fragments placed as orphans (only one mate mapped).
    ///
    /// `Some` whenever the records were tallied, which is every path here:
    /// any RAD input — orphan status is in the records themselves, so this
    /// works for a piscem RAD too, independently of whether the file carries
    /// salmon's provenance — and BAM input, which tallies orphans while
    /// streaming. The `Option` remains because a future producer could leave
    /// the count unknown, and `meta_info.json` must be able to say so rather
    /// than report a fabricated 0.
    pub num_orphan: Option<u64>,
    /// fragments judged against a known expected library format that had at
    /// least one strand-compatible placement / had none (#1130). An unstranded
    /// format is judged like any other and simply finds every fragment
    /// compatible. Both zero only when no placement was ever compared to an
    /// expected format, where the writer falls back to the historical
    /// `num_mapped` / ratio-1.0 form. Feed `lib_format_counts.json`, mirroring
    /// the reads-mode fields in `salmon-quant`'s `QuantResult`.
    pub num_compatible_fragments: u64,
    pub num_incompatible_fragments: u64,
    /// per-observed-format fragment counts (`counts[format_id]`, one count per
    /// distinct format among a fragment's placements), tallied for every
    /// fragment with placements whatever the expected format. The raw histogram
    /// behind `lib_format_counts.json`'s per-format keys and its
    /// concordant/inconsistent/strand-bias derivations.
    pub lib_format_counts: salmon_core::LibFormatCountsArray,
    /// what the RAD says about the mapping pass that produced it; empty for a
    /// BAM input, which has no RAD provenance to read
    pub provenance: crate::rad::RadProvenance,
    /// verbatim `@PG` lines describing how the fragments were produced: read
    /// from a BAM's header directly, or carried through a BAM-derived RAD
    pub source_programs: Vec<String>,
    /// what this run actually read, which decides which metadata is even
    /// applicable — an index hash is missing from a RAD requant but simply does
    /// not apply to a BAM, and blaming a RAD that was never read would be wrong
    pub source: FragmentSource,
    /// fragment mass dropped in the final min-alpha redistribution (every member
    /// of an equivalence class truncated); normally 0.
    pub inference_truncated_mass: f64,
    /// the GPU backend that warm-started the EM/VBEM point estimate, or `None`
    /// when it ran on the CPU only
    pub gpu_backend: Option<&'static str>,
    pub num_eq_classes: usize,
    pub frag_len_mean: f64,
    pub frag_len_sd: f64,
    /// where the fragment-length distribution came from; reported as
    /// `frag_length_source` in `meta_info.json` so a result's FLD provenance
    /// stays auditable after the run
    pub frag_len_source: salmon_model::FragLengthSource,
    pub length_classes: Vec<u32>,
    pub frag_len_dist: Vec<f64>,
    pub start_time: String,
    pub bias_dump: salmon_model::dumps::BiasDump,
    /// per-transcript (unique, ambiguous) fragment counts for `ambig_info.tsv`
    pub ambig: (Vec<u64>, Vec<u64>),
    /// posterior samples (bootstrap or Gibbs), one abundance vector each; empty
    /// when neither was requested. Length is `num_refs`, matching `quant.sf` rows.
    pub bootstraps: Vec<Vec<f64>>,
    /// EM/VBEM convergence: iterations run and whether the relative-difference
    /// criterion was met before `max_iter`
    pub em_iters: u32,
    pub em_converged: bool,
    /// how each connected component's iteration ended (see
    /// [`salmon_infer::ComponentSummary`]); `None` off the component engine
    pub em_components: Option<salmon_infer::ComponentSummary>,
    /// library format observed for these alignments (baked/auto-detected), for
    /// provenance and the strandedness sanity check; `None` if not observed
    pub detected_library_type: Option<String>,
    /// total wall-clock seconds for the quantification call
    pub total_seconds: f64,
    /// peak resident set size in KiB (Linux `VmHWM`); 0 if unavailable
    pub peak_rss_kb: u64,
    /// structured run diagnostics / bad-input warnings (also emitted to the log)
    pub diagnostics: Vec<salmon_core::Diagnostic>,
}

pub use salmon_core::asctime_now;

/// Extract an integer tag value (e.g. `AS`) as `i32`.
fn value_as_i32(v: &Value) -> Option<i32> {
    match v {
        Value::Int8(x) => Some(*x as i32),
        Value::UInt8(x) => Some(*x as i32),
        Value::Int16(x) => Some(*x as i32),
        Value::UInt16(x) => Some(*x as i32),
        Value::Int32(x) => Some(*x),
        Value::UInt32(x) => Some(*x as i32),
        _ => None,
    }
}

/// One placement of the fragment on a single transcript: the record indices the
/// *aligner* reported together — a proper pair (two mates that point at each
/// other) or a single orphan record. Indices are sorted by reference position
/// (so index 0 is the left mate).
struct Placement {
    tid: u32,
    idxs: Vec<usize>,
}

/// Resolve a fragment's records into the placements the aligner actually
/// intended, rather than cross-producting every read1 with every read2 on a
/// transcript.
///
/// A permissive aligner (e.g. `bowtie2 -k`) reports many alignment records per
/// fragment; the two mates of one reported pair are linked by their mate fields
/// — each record's `RNEXT`/`PNEXT` points at the other's `(tid, pos)` — and, when
/// present, share a hit index (`HI`). Pairing by transcript co-occurrence instead
/// fabricates concordant pairs the aligner never reported, which keeps a fragment
/// artificially orientation-compatible and defeats salmon's dropping of
/// protocol-inconsistent fragments. Here we pair a read1 record with the read2
/// record that reciprocally references it (and matches its `HI` when both carry
/// one); records left unpaired — including all records of a single-end library —
/// become orphan placements.
fn pair_records(recs: &[FragRecord]) -> Vec<Placement> {
    let n = recs.len();
    let mut used = vec![false; n];
    let mut placements: Vec<Placement> = Vec::with_capacity(n);

    for i in 0..n {
        if used[i] || !recs[i].is_read1 {
            continue;
        }
        let r1 = &recs[i];
        // Only mates aligned to the *same* transcript form a single-transcript
        // pair placement; a mate on another transcript leaves r1 an orphan.
        let (Some(mtid), Some(mpos)) = (r1.mate_tid, r1.mate_pos) else {
            continue;
        };
        if mtid != r1.tid {
            continue;
        }
        let mut mate: Option<usize> = None;
        for j in 0..n {
            if used[j] || recs[j].is_read1 {
                continue;
            }
            let r2 = &recs[j];
            if r2.tid != mtid
                || r2.pos != mpos
                || r2.mate_tid != Some(r1.tid)
                || r2.mate_pos != Some(r1.pos)
            {
                continue;
            }
            // A reciprocal coordinate match; prefer one whose HI agrees.
            let hi_ok = matches!((r1.hi, r2.hi), (Some(a), Some(b)) if a == b)
                || r1.hi.is_none()
                || r2.hi.is_none();
            if hi_ok {
                mate = Some(j);
                break;
            }
            if mate.is_none() {
                mate = Some(j);
            }
        }
        if let Some(j) = mate {
            used[i] = true;
            used[j] = true;
            let mut idxs = vec![i, j];
            idxs.sort_by_key(|&k| recs[k].pos);
            placements.push(Placement { tid: r1.tid, idxs });
        }
    }
    for i in 0..n {
        if !used[i] {
            placements.push(Placement {
                tid: recs[i].tid,
                idxs: vec![i],
            });
        }
    }
    placements
}

/// `log(Σ exp(xs))` over an iterator, numerically stable. `xs` is non-empty.
///
/// Same two passes and the same summation order as [`logsumexp`], so the two
/// agree bit for bit; this form exists so a caller holding a run of records can
/// avoid copying their weights into a scratch slice first.
pub(crate) fn logsumexp_iter<I>(xs: I) -> f64
where
    I: Iterator<Item = f64> + Clone,
{
    let m = xs.clone().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    m + xs.map(|x| (x - m).exp()).sum::<f64>().ln()
}

/// [`logsumexp_iter`] of a single value, without its `exp` and `ln`.
///
/// For one finite `x` the general form is `x + ln(exp(x − x)) = x + ln(1)`,
/// and `exp(0) = 1` and `ln(1) = +0` exactly (IEEE 754 / C Annex F special
/// cases, which every libm honours), so it is `x + 0.0` — written as such,
/// since that also turns `-0.0` into `+0.0` the way the general form does. The
/// non-finite cases follow the general form's own: no finite maximum gives
/// `-inf` (`f64::max` skips a NaN), and `+inf` gives `inf − inf`, a NaN.
///
/// A RAD fragment almost always has one placement per transcript, so this is
/// the common case of the per-transcript grouping, where it saved an `exp` and
/// an `ln` per placement.
pub(crate) fn logsumexp_one(x: f64) -> f64 {
    if x.is_nan() || x == f64::NEG_INFINITY {
        f64::NEG_INFINITY
    } else if x == f64::INFINITY {
        f64::NAN
    } else {
        x + 0.0
    }
}

/// `log(Σ exp(xs))`, numerically stable. `xs` is non-empty.
fn logsumexp(xs: &[f64]) -> f64 {
    let m = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    m + xs.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
}

/// Derive the observed library format, the single-read forward flag, and the
/// mate status for one reported pair / orphan (`recs[idxs]`), mirroring salmon's
/// `hitType` (`src/util/SalmonUtils.cpp`).
///
/// salmon classifies orientation from the **leftmost** reference coordinate
/// (`bam_pos`) of each mate — *not* their 5' ends: a forward/reverse pair is
/// inward (TOWARD) iff the forward mate starts at or before the reverse mate, and
/// outward (AWAY) otherwise. For a dovetailed/overlapping short fragment where the
/// reverse mate's leftmost precedes the forward mate's, salmon therefore reports
/// the pair as outward and (under a strict library type) drops it as a
/// zero-probability fragment. We match that convention exactly. Strandedness is
/// keyed on which mate is forward (read 1 forward → SA, read 2 forward → AS).
fn frag_format(recs: &[FragRecord], idxs: &[usize]) -> (Option<LibraryFormat>, bool, MateStatus) {
    if idxs.len() >= 2 {
        let r1 = idxs
            .iter()
            .map(|&i| &recs[i])
            .find(|r| r.is_read1)
            .unwrap_or(&recs[idxs[0]]);
        let r2 = idxs
            .iter()
            .map(|&i| &recs[i])
            .find(|r| !r.is_read1)
            .unwrap_or(&recs[idxs[1]]);
        let (r1_fw, r2_fw) = (!r1.is_reverse, !r2.is_reverse);
        let (orientation, strandedness) = if r1_fw != r2_fw {
            let (fw, rc) = if r1_fw { (r1, r2) } else { (r2, r1) };
            let orientation = if fw.pos <= rc.pos {
                ReadOrientation::Toward
            } else {
                ReadOrientation::Away
            };
            let strandedness = if r1_fw {
                ReadStrandedness::SA
            } else {
                ReadStrandedness::AS
            };
            (orientation, strandedness)
        } else {
            (
                ReadOrientation::Same,
                if r1_fw {
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
            r1_fw,
            MateStatus::PairedEndPaired,
        )
    } else {
        let r = &recs[idxs[0]];
        // A record that is not part of a multi-segment template (`0x1` unset) is a
        // genuine single-end read → SingleEnd. A paired read reported alone here is
        // an orphan (its mate unmapped / not grouped), classified left/right by the
        // `0x40` first-segment flag. Mirrors salmon's `hitType` and the genome
        // projection path (`genome_project.rs`). Conflating single-end reads with
        // PairedEndRight made strand filters (SF/SR) drop every single-end
        // alignment (issue #1057).
        let status = if !r.is_paired {
            MateStatus::SingleEnd
        } else if r.is_read1 {
            MateStatus::PairedEndLeft
        } else {
            MateStatus::PairedEndRight
        };
        (None, !r.is_reverse, status)
    }
}

/// 2-bit encode a base (`A=0, C=1, G=2, T=3`; anything else → `0`).
#[inline]
fn base_2bit(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// Load a (optionally compressed) transcriptome FASTA and return the (ASCII)
/// base sequences aligned to the BAM's reference order (`names`); a name absent
/// from the FASTA yields an empty sequence (its model contributions are
/// skipped). The same bytes feed both the error model (2-bit on the fly) and the
/// bias models.
///
/// Compression is detected from content, so gzip/BGZF, bzip2, xz and zstd all
/// work regardless of how the file is named.
///
/// `lengths` are the input's own reference lengths: a same-named sequence of
/// another length is another sequence, so positions would read the wrong bases.
fn load_ref_bytes(path: &Path, names: &[String], lengths: &[u32]) -> Result<Vec<Vec<u8>>> {
    let reader: Box<dyn BufRead> = salmon_core::compress::open_maybe_compressed(path)
        .with_context(|| format!("opening {}", path.display()))?;

    let mut by_name: HashMap<String, Vec<u8>> = HashMap::new();
    let mut cur_name: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        // Named here as well as at open: decompression failures surface while
        // reading, not while opening, so without this a truncated file reports a
        // bare "unexpected end of file" with nothing to say which file.
        let line = line.with_context(|| format!("reading {}", path.display()))?;
        if let Some(stripped) = line.strip_prefix('>') {
            if let Some(n) = cur_name.take() {
                by_name.insert(n, std::mem::take(&mut cur_seq));
            }
            cur_name = Some(stripped.split_whitespace().next().unwrap_or("").to_string());
        } else if cur_name.is_some() {
            cur_seq.extend(line.trim_end().bytes());
        }
    }
    if let Some(n) = cur_name.take() {
        by_name.insert(n, cur_seq);
    }

    // This is only reached when the user asked for the error or bias models and
    // supplied `-t`, so a file with no records in it is a mistake, not a choice.
    // Left unchecked it is silent: every reference resolves to an empty sequence,
    // the models quietly contribute nothing, and the run reports success with
    // numbers that differ from what was asked for.
    anyhow::ensure!(
        !by_name.is_empty(),
        "no sequences found in {}: expected FASTA records beginning with '>'. \
         The error and sequence-bias models are built from this file, so an empty \
         or non-FASTA one would silently disable them rather than fail.",
        path.display()
    );

    let mut seqs: Vec<Vec<u8>> = names
        .iter()
        .map(|n| by_name.remove(n).unwrap_or_default())
        .collect();
    // A poly-A tail past the input's length is one the salmon index clipped:
    // drop it, or bias correction would count fragment positions on it (#1314).
    for (seq, &len) in seqs.iter_mut().zip(lengths) {
        let l = len as usize;
        if seq.len() > l && seq[l..].iter().all(|b| b.eq_ignore_ascii_case(&b'A')) {
            seq.truncate(l);
        }
    }
    // Not an error: a transcriptome legitimately need not cover every reference
    // in the alignment header. But *no* overlap means the wrong file, which
    // otherwise presents as models that do nothing rather than as a mistake.
    if seqs.iter().all(|s| s.is_empty()) {
        tracing::warn!(
            "none of the {} reference(s) in the alignment input were found in {}. \
             The error and sequence-bias models will contribute nothing. This \
             usually means the FASTA does not match the alignments — check that it \
             is the transcriptome they were aligned against, and that identifiers \
             match (versioned vs unversioned, or a GENCODE-style `|` header).",
            names.len(),
            path.display()
        );
    }
    let (wrong, example) = length_mismatches(names, &seqs, lengths);
    if let Some((name, fasta_len, input_len)) = example {
        tracing::warn!(
            "{wrong} reference(s) in {} have a different length from the alignment input \
             (e.g. {name}: {fasta_len} in the FASTA, {input_len} in the input). The error \
             and bias models would read the wrong bases: check that this is the \
             transcriptome the input was aligned against.",
            path.display()
        );
    }
    Ok(seqs)
}

/// How many same-named references differ in length, and the first of them as
/// `(name, FASTA length, input length)`. An absent reference is not a mismatch
/// (see above), nor is a FASTA sequence longer only by a poly-A tail: the salmon
/// index clips those, so a RAD's lengths legitimately lack them.
fn length_mismatches<'a>(
    names: &'a [String],
    seqs: &[Vec<u8>],
    lengths: &[u32],
) -> (usize, Option<(&'a str, usize, u32)>) {
    let mut wrong = 0;
    let mut first = None;
    for ((name, seq), &len) in names.iter().zip(seqs).zip(lengths) {
        let l = len as usize;
        let clipped_tail = seq.len() > l && seq[l..].iter().all(|b| b.eq_ignore_ascii_case(&b'A'));
        if !seq.is_empty() && seq.len() != l && !clipped_tail {
            wrong += 1;
            first.get_or_insert((name.as_str(), seq.len(), len));
        }
    }
    (wrong, first)
}

/// One alignment record needed by the error model / weighting.
struct FragRecord {
    tid: u32,
    pos: usize,
    read_2bit: Vec<u8>,
    ops: Vec<(AlnOp, usize)>,
    /// aligner `AS` tag; parsed for completeness but not used for weighting —
    /// the error model's `errLike` drives the conditional weight, matching
    /// salmon's behavior for CIGAR-bearing aligners (see [`salmon-align-as-weighting`]).
    #[allow(dead_code)]
    score: i32,
    /// whether the record actually carried an `AS` tag. A missing tag reads as
    /// `score: 0`, which is indistinguishable from a genuine zero, so presence
    /// has to be recorded where the tag is read rather than inferred later.
    has_score: bool,
    frag_len: i32,
    /// reverse-strand alignment (BAM `0x10` flag)
    is_reverse: bool,
    /// first mate of the pair (BAM `0x40` flag)
    is_read1: bool,
    /// the read is part of a multi-segment template (BAM `0x1` "paired" flag).
    /// `false` marks a genuine single-end read, distinguishing it from a
    /// paired-end orphan (mate unmapped / absent) — the two are classified
    /// differently by [`frag_format`] (SingleEnd vs PairedEndLeft/Right).
    is_paired: bool,
    /// reference span (Σ ref-consuming CIGAR op lengths); the read's 3' end on
    /// the reference is `pos + ref_span − 1`
    ref_span: usize,
    /// mate's transcript id (`RNEXT`), if the mate is mapped
    mate_tid: Option<u32>,
    /// mate's 0-based alignment start (`PNEXT`), if the mate is mapped
    mate_pos: Option<usize>,
    /// hit index (`HI` tag): the aligner's pairing id, used to disambiguate which
    /// mate records form one reported alignment when several share coordinates
    hi: Option<i32>,
}

impl FragRecord {}

/// The canonical read name: a trailing `/1` or `/2` mate suffix stripped (as
/// salmon's `getPairedNameLen` does), so the two mates of a fragment group
/// together even when the aligner kept the suffix in the QNAME.
#[inline]
fn canonical_name(name: &[u8]) -> &[u8] {
    let l = name.len();
    if l > 2 && name[l - 2] == b'/' && matches!(name[l - 1], b'1' | b'2') {
        &name[..l - 2]
    } else {
        name
    }
}

/// Groups an alignment stream's consecutive records by read, for the reader
/// threads. Only mapped records are kept: a read with none is counted here and
/// never reaches the workers, so it costs no more than its decode, and the
/// minibatches hold the same reads they would without it.
pub(crate) struct ReadGrouper<R> {
    name: Vec<u8>,
    group: Vec<R>,
    started: bool,
    /// Reads seen with no mapped record: processed, never mapped.
    pub(crate) unmapped_reads: u64,
    /// Whether the mates of each read are next to each other.
    pub(crate) grouping: GroupingCheck,
}

impl<R: sam::alignment::Record> ReadGrouper<R> {
    pub(crate) fn new() -> Self {
        Self {
            name: Vec::new(),
            group: Vec::new(),
            started: false,
            unmapped_reads: 0,
            grouping: GroupingCheck::default(),
        }
    }

    /// Add the next record. Returns the previous read's mapped records when this
    /// record starts a new read.
    pub(crate) fn push(&mut self, rec: R) -> Option<Vec<R>> {
        let name = rec.name()?;
        let name = canonical_name(name.as_ref());
        let mut done = None;
        if !self.started || name != self.name.as_slice() {
            if self.started {
                done = self.close();
            }
            self.name.clear();
            self.name.extend_from_slice(name);
            self.started = true;
        }
        match rec.flags() {
            Ok(flags) if !flags.is_unmapped() => {
                self.grouping.record(flags);
                self.group.push(rec);
            }
            _ => {}
        }
        done
    }

    /// The last read's mapped records, at the end of the stream.
    pub(crate) fn finish(&mut self) -> Option<Vec<R>> {
        if self.started {
            self.started = false;
            self.close()
        } else {
            None
        }
    }

    fn close(&mut self) -> Option<Vec<R>> {
        if self.group.is_empty() {
            self.unmapped_reads += 1;
            None
        } else {
            self.grouping.end_read();
            Some(std::mem::take(&mut self.group))
        }
    }
}

/// Give records without bases (SEQ `*`, usual on secondary alignments) those of
/// another record of the same read. Scored as nothing, such a record gets
/// log-likelihood 0, which beats any real alignment. Needs the same query length,
/// so a hard-clipped record whose bases cannot be lined up stays empty.
fn fill_missing_bases(frags: &mut [FragRecord]) {
    for i in 0..frags.len() {
        if !frags[i].read_2bit.is_empty() {
            continue;
        }
        let query_len: usize = frags[i]
            .ops
            .iter()
            .filter(|(op, _)| op.consume_seq())
            .map(|(_, len)| len)
            .sum();
        let donor = frags.iter().find(|d| {
            d.is_read1 == frags[i].is_read1 && d.read_2bit.len() == query_len && query_len > 0
        });
        let Some(donor) = donor else { continue };
        // Bases are stored in reference orientation, so the other strand needs
        // the reverse complement (`3 - b` in this encoding).
        let bases = if donor.is_reverse == frags[i].is_reverse {
            donor.read_2bit.clone()
        } else {
            donor.read_2bit.iter().rev().map(|b| 3 - b).collect()
        };
        frags[i].read_2bit = bases;
    }
}

/// Map a noodles CIGAR op kind to our `AlnOp`.
fn kind_to_op(k: Kind) -> AlnOp {
    match k {
        Kind::Match => AlnOp::Match,
        Kind::SequenceMatch => AlnOp::SeqMatch,
        Kind::SequenceMismatch => AlnOp::SeqMismatch,
        Kind::Insertion => AlnOp::Ins,
        Kind::Deletion => AlnOp::Del,
        Kind::Skip => AlnOp::RefSkip,
        Kind::SoftClip => AlnOp::SoftClip,
        Kind::HardClip => AlnOp::HardClip,
        Kind::Pad => AlnOp::Pad,
    }
}

/// Bytes read from the front of an alignment input to tell BAM from SAM: at
/// least one whole BGZF block (64 KiB at most), so a BAM's first block always
/// decompresses from them.
const SNIFF_BYTES: u64 = 1 << 16;

type SamReader = sam::io::Reader<Box<dyn BufRead + Send>>;
type BamReader = bam::io::Reader<noodles_bgzf::io::MultithreadedReader<Box<dyn Read + Send>>>;

/// The records of an alignment input whose header has been read, plus any
/// records [`AlignmentRecords::peek_is_paired`] read ahead (replayed first).
///
/// Every alignment path makes one pass over the records, and the input is
/// opened exactly once, by [`open_alignment_input`]. So the input can be a
/// pipe: a FIFO, `-a <(samtools view -b ...)`, or `-a /dev/stdin`.
pub(crate) enum AlignmentRecords {
    Sam(SamReader, Vec<std::io::Result<sam::Record>>),
    Bam(BamReader, Vec<std::io::Result<bam::Record>>),
}

/// Bind `$recs` to the records of an [`AlignmentRecords`] (read-ahead records
/// first) as an `Iterator<Item = io::Result<R>>` and evaluate `$body`, once per
/// record type: the passes are generic over SAM and BAM records.
macro_rules! with_alignment_records {
    ($input:expr, $recs:ident => $body:expr) => {
        match $input {
            $crate::AlignmentRecords::Sam(mut reader, ahead) => {
                let $recs = ahead.into_iter().chain(reader.records());
                $body
            }
            $crate::AlignmentRecords::Bam(mut reader, ahead) => {
                let $recs = ahead.into_iter().chain(reader.records());
                $body
            }
        }
    };
}
pub(crate) use with_alignment_records;

/// Open an alignment input once and read its header.
///
/// SAM or BAM is decided by content, not by name: a BGZF stream whose data
/// starts with `BAM\1` is BAM; anything else is SAM, decompressed by content
/// (plain, gzip/BGZF, bzip2, xz, zstd). A name says nothing for `/dev/stdin` or
/// `/dev/fd/63`, and a pipe cannot be reopened after a wrong guess.
fn open_alignment_input(path: &Path) -> Result<(Header, AlignmentRecords)> {
    let mut file = std::fs::File::open(path)
        .with_context(|| format!("opening alignment input {}", path.display()))?;
    let mut head = Vec::new();
    (&mut file)
        .take(SNIFF_BYTES)
        .read_to_end(&mut head)
        .with_context(|| format!("reading alignment input {}", path.display()))?;
    anyhow::ensure!(
        !head.is_empty(),
        "alignment input {} is empty",
        path.display()
    );
    let mut magic = [0u8; 4];
    let is_bam = noodles_bgzf::io::Reader::new(&head[..])
        .read_exact(&mut magic)
        .is_ok()
        && &magic == b"BAM\x01";
    if is_bam {
        let stream: Box<dyn Read + Send> = Box::new(std::io::Cursor::new(head).chain(file));
        // BGZF decompression on a small worker pool; framing/parsing stays serial
        // on the reader thread, but the decompression overlaps it.
        let workers = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4)
            .clamp(1, 4);
        let workers = std::num::NonZeroUsize::new(workers).unwrap();
        let decoder = noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, stream);
        let mut reader = bam::io::Reader::from(decoder);
        let header = reader.read_header().context("reading BAM header")?;
        Ok((header, AlignmentRecords::Bam(reader, Vec::new())))
    } else {
        let inner = salmon_core::compress::maybe_decompress_stream(head, Box::new(file))
            .with_context(|| format!("reading alignment input {}", path.display()))?;
        let inner: Box<dyn BufRead + Send> = Box::new(std::io::BufReader::new(inner));
        let mut reader = sam::io::Reader::new(inner);
        let header = reader.read_header().context("reading SAM header")?;
        Ok((header, AlignmentRecords::Sam(reader, Vec::new())))
    }
}

impl AlignmentRecords {
    /// Whether the first record belongs to a paired-end template (C++ salmon's
    /// `peekBAMIsPaired`): under `-l A` the user has not said, and the
    /// library-type detector needs to know which formats are candidates. `None`
    /// for an input with no records. The records read to answer are kept and
    /// replayed by the pass, so this costs no second open.
    fn peek_is_paired(&mut self) -> Option<bool> {
        fn first_flag<R: sam::alignment::Record>(
            records: impl Iterator<Item = std::io::Result<R>>,
            ahead: &mut Vec<std::io::Result<R>>,
        ) -> Option<bool> {
            for rec in records {
                // A read error ends the peek; the pass reports it.
                let paired = match &rec {
                    Ok(r) => r.flags().ok().map(|f| f.is_segmented()),
                    Err(_) => {
                        ahead.push(rec);
                        return None;
                    }
                };
                ahead.push(rec);
                if paired.is_some() {
                    return paired;
                }
            }
            None
        }
        match self {
            Self::Sam(reader, ahead) => first_flag(reader.records(), ahead),
            Self::Bam(reader, ahead) => first_flag(reader.records(), ahead),
        }
    }
}

/// Convert one mapped alignment record into a `FragRecord`. `None` skips the
/// record (unmapped / unnamed / no reference). Works over any noodles
/// `alignment::Record`, so the same logic serves BAM and SAM input.
fn record_to_frag<R: sam::alignment::Record>(
    record: &R,
    header: &Header,
    need_seq: bool,
) -> Option<FragRecord> {
    let flags = record.flags().ok()?;
    if flags.is_unmapped() {
        return None;
    }
    // Grouping by name happened on the reader thread; only presence matters here.
    record.name()?;

    let tid = record.reference_sequence_id(header)?.ok()?;
    let pos = record
        .alignment_start()
        .and_then(|r| r.ok())
        .map(|p| p.get() - 1)
        .unwrap_or(0);
    let read_2bit = if need_seq {
        record.sequence().iter().map(base_2bit).collect()
    } else {
        Vec::new()
    };
    let cigar = record.cigar();
    let all_ops = || {
        cigar
            .iter()
            .filter_map(|r| r.ok())
            .map(|op| (kind_to_op(op.kind()), op.len()))
    };
    let ref_span: usize = all_ops()
        .filter(|(o, _)| {
            matches!(
                o,
                AlnOp::Match | AlnOp::SeqMatch | AlnOp::SeqMismatch | AlnOp::Del | AlnOp::RefSkip
            )
        })
        .map(|(_, l)| l)
        .sum();
    // Only the error model reads the op list.
    let ops: Vec<(AlnOp, usize)> = if need_seq {
        all_ops().collect()
    } else {
        Vec::new()
    };
    let parsed_score = record
        .data()
        .get(&Tag::ALIGNMENT_SCORE)
        .and_then(|r| r.ok())
        .and_then(|v| value_as_i32(&v));
    let score = parsed_score.unwrap_or(0);
    let has_score = parsed_score.is_some();
    let frag_len = record.template_length().map(|t| t.abs()).unwrap_or(0);
    let is_reverse = flags.is_reverse_complemented();
    let is_read1 = flags.is_first_segment();
    let is_paired = flags.is_segmented();
    // Mate linkage as the aligner recorded it (RNEXT/PNEXT); a mate that is
    // unmapped or absent leaves these `None`, making the record an orphan.
    let mate_tid = (!flags.is_mate_unmapped())
        .then(|| {
            record
                .mate_reference_sequence_id(header)
                .and_then(|r| r.ok())
        })
        .flatten()
        .map(|t| t as u32);
    let mate_pos = record
        .mate_alignment_start()
        .and_then(|r| r.ok())
        .map(|p| p.get() - 1);
    let hi = record
        .data()
        .get(&Tag::HIT_INDEX)
        .and_then(|r| r.ok())
        .and_then(|v| value_as_i32(&v));
    Some(FragRecord {
        tid: tid as u32,
        pos,
        read_2bit,
        ops,
        score,
        has_score,
        frag_len,
        is_reverse,
        is_read1,
        is_paired,
        ref_span,
        mate_tid,
        mate_pos,
        hi,
    })
}

/// Whether a RAD writer's library is paired: under `-l A` the input decides,
/// otherwise a single-end `-l` says no and anything else says yes.
pub(crate) fn rad_library_is_paired(records: &mut AlignmentRecords, lib_type: &str) -> bool {
    let peeked = if LibraryFormat::is_auto(lib_type) {
        records.peek_is_paired()
    } else {
        None
    };
    peeked.unwrap_or(!matches!(lib_type, "U" | "SF" | "SR" | "S"))
}

/// Is the input coordinate-sorted and *not* grouped by read name?
///
/// Read once from the `@HD` header (no per-record cost). `SO:coordinate` orders
/// records by position, scattering a read's alignments — unusable here. But a file
/// can be coordinate-sorted and then re-grouped by name (`samtools collate`, which
/// sets `GO:query`), or carry a stale `SO:coordinate` after several samtools steps;
/// the reliable signal that records are usably grouped is `GO:query`. So we only
/// reject when coordinate-sorted AND not query-grouped. (Query-name *sorted* files
/// report `SO:queryname`, not coordinate, and pass.)
/// One `@PG` header record: what a program said about how it touched the file.
///
/// `command_line` is the field that actually answers "how was this BAM made";
/// the rest identify the program well enough to make sense of it.
#[derive(Clone, Debug, Default, Serialize, PartialEq, Eq)]
pub struct SourceProgram {
    /// `ID` — unique within the file, and the link target of `PP`
    pub id: String,
    /// `PN` — program name
    #[serde(skip_serializing_if = "Option::is_none")]
    pub program_name: Option<String>,
    /// `VN` — program version
    #[serde(skip_serializing_if = "Option::is_none")]
    pub version: Option<String>,
    /// `CL` — the command line the program was invoked with
    #[serde(skip_serializing_if = "Option::is_none")]
    pub command_line: Option<String>,
    /// `PP` — the previous program in the chain, so the order is recoverable
    #[serde(skip_serializing_if = "Option::is_none")]
    pub previous_id: Option<String>,
    /// `DS` — free-text description
    #[serde(skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
}

/// Render a BAM header's `@PG` chain as SAM lines.
///
/// Kept as text on the way into a RAD so nothing the aligner recorded is dropped
/// in transit; [`parse_source_programs`] reads them back. Empty when the header
/// has no `@PG` records at all, which is common enough — plenty of tools write
/// none — that callers must treat it as ordinary rather than exceptional.
///
/// Values are escaped rather than passed through raw. RAD imposes nothing here:
/// it stores a string as an explicit length followed by bytes, so it would carry
/// any content. The constraint is entirely from *this* representation — records
/// joined by newlines, fields by tabs — so those characters are escaped to keep
/// the framing unambiguous, and [`parse_source_programs`] restores them. SAM
/// forbids them in header values anyway, so it only matters for a malformed
/// input, where preserving the value beats silently rewriting it.
pub fn source_program_lines(header: &noodles_sam::Header) -> Vec<String> {
    let mut lines: Vec<String> = header
        .programs()
        .as_ref()
        .iter()
        .map(|(id, program)| {
            let mut line = format!("@PG\tID:{}", escape_header_value(id.as_ref()));
            // Every tag the aligner wrote is carried through, standard or not:
            // this is provenance, not a schema to validate against.
            for (tag, value) in program.other_fields() {
                let t: &[u8; 2] = tag.as_ref();
                line.push_str(&format!(
                    "\t{}:{}",
                    String::from_utf8_lossy(t),
                    escape_header_value(value.as_ref())
                ));
            }
            line
        })
        .collect();

    // libradicl bounds an over-long tag value rather than letting it wrap the
    // length field, but it bounds by *bytes*, which here would cut a record in
    // half and leave a `@PG` line whose command is silently wrong. Dropping whole
    // trailing records instead keeps every record that survives intact and
    // truthful, so ask whether the chain fits and shorten it ourselves.
    while !lines.is_empty() && !chain_fits(&lines) {
        let dropped = lines.pop().unwrap_or_default();
        tracing::warn!(
            "the @PG chain is longer than a RAD string tag can hold; dropping its \
             trailing record ({}...) so no record is left half-written",
            dropped.chars().take(60).collect::<String>()
        );
    }
    lines
}

/// Whether the joined `@PG` chain fits a RAD string tag.
///
/// Asks the format rather than restating its width: the limit belongs to
/// libradicl, and hard-coding it here is how the two would drift.
fn chain_fits(lines: &[String]) -> bool {
    salmon_rad::value_fits(
        &libradicl::rad_types::TagValue::String(lines.join("\n")),
        &libradicl::rad_types::RadType::String,
    )
}

/// Escape the characters that frame the joined `@PG` representation.
///
/// Backslash is escaped too, so the mapping is reversible: without it a value
/// containing a literal backslash-n (two characters) would be indistinguishable
/// from an escaped newline.
fn escape_header_value(value: &[u8]) -> String {
    let mut out = String::new();
    for c in String::from_utf8_lossy(value).chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            _ => out.push(c),
        }
    }
    out
}

/// Reverse [`escape_header_value`]. An unrecognized escape is left as written,
/// so text that was never escaped survives unchanged.
fn unescape_header_value(value: &str) -> String {
    let mut out = String::new();
    let mut chars = value.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('\\') => out.push('\\'),
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('r') => out.push('\r'),
            Some(other) => {
                out.push('\\');
                out.push(other);
            }
            None => out.push('\\'),
        }
    }
    out
}

/// Parse `@PG` lines produced by [`source_program_lines`] back into records.
///
/// Unknown tags are ignored rather than rejected: a `@PG` record may carry
/// anything, and provenance that fails to parse is worse than provenance that
/// reports only the fields it understands.
pub fn parse_source_programs(lines: &[String]) -> Vec<SourceProgram> {
    lines
        .iter()
        .filter(|l| l.starts_with("@PG\t"))
        .map(|line| {
            let mut prog = SourceProgram::default();
            for field in line.split('\t').skip(1) {
                let Some((tag, value)) = field.split_once(':') else {
                    continue;
                };
                let value = unescape_header_value(value);
                match tag {
                    "ID" => prog.id = value,
                    "PN" => prog.program_name = Some(value),
                    "VN" => prog.version = Some(value),
                    "CL" => prog.command_line = Some(value),
                    "PP" => prog.previous_id = Some(value),
                    "DS" => prog.description = Some(value),
                    _ => {}
                }
            }
            prog
        })
        .collect()
}

/// Spots input whose mates are not next to each other, which the header check
/// below cannot do when the header says nothing (`SO:unknown`, or no `@HD`). Each
/// mate then forms a fragment of its own, so every count doubles while each record
/// looks like an orphan. Tallied on the reader thread, from the flags
/// [`ReadGrouper`] has already read.
#[derive(Default)]
pub(crate) struct GroupingCheck {
    fragments: u64,
    /// fragments with a record whose mate is marked mapped, but no such mate
    split: u64,
    first: bool,
    last: bool,
    expects_mate: bool,
}

impl GroupingCheck {
    /// One mapped record of the current read.
    fn record(&mut self, flags: sam::alignment::record::Flags) {
        self.first |= flags.is_first_segment();
        self.last |= flags.is_last_segment();
        self.expects_mate |= flags.is_segmented() && !flags.is_mate_unmapped();
    }

    /// The current read, which had at least one mapped record, is complete.
    fn end_read(&mut self) {
        self.fragments += 1;
        self.split += u64::from(self.expects_mate && !(self.first && self.last));
        (self.first, self.last, self.expects_mate) = (false, false, false);
    }

    /// Refuses input where at least half the fragments are split: that is not a
    /// filtered BAM but one in the wrong order, and its numbers would be wrong.
    /// A smaller share can be real (a filter that kept one mate), so it only warns.
    pub(crate) fn finish(&self) -> Result<()> {
        let Self {
            fragments, split, ..
        } = *self;
        let hint = "Records of one read must be consecutive: group them by name with \
                    `samtools collate` or `samtools sort -n`";
        anyhow::ensure!(
            split * 2 < fragments || fragments == 0,
            "{split} of {fragments} fragments have a mate marked as mapped that is not next \
             to them, so the input is not grouped by read name. {hint}."
        );
        if split * 20 >= fragments && split > 0 {
            tracing::warn!(
                "{split} of {fragments} fragments have a mate marked as mapped that is not next \
                 to them. They are quantified as orphans. If the mates are really there: {hint}."
            );
        }
        Ok(())
    }
}

fn coordinate_sorted_unusable(header: &noodles_sam::Header) -> bool {
    let Some(hd) = header.header() else {
        return false;
    };
    let mut so_coord = false;
    let mut go_query = false;
    // `SO`/`GO` are non-standard tags in this noodles version → read from other_fields.
    for (tag, value) in hd.other_fields() {
        let t: &[u8; 2] = tag.as_ref();
        let v: &[u8] = value.as_ref();
        if t == b"SO" {
            so_coord = v == &b"coordinate"[..];
        } else if t == b"GO" {
            go_query = v == &b"query"[..];
        }
    }
    so_coord && !go_query
}

/// Build expected bias models from the collected seq/GC/positional observations
/// plus abundance estimates, fold them into bias-corrected **effective lengths**
/// (mutating `eff_lengths` in place), and return the [`BiasDump`] for output.
///
/// Called by the RAD quantifier, which collects the per-fragment observations
/// weighted by abundance-aware posteriors over fixed baked/derived abundances.
/// The caller runs the subsequent re-EM with the corrected lengths.
#[allow(clippy::too_many_arguments)]
/// Bias-correct the effective lengths in place.
///
/// `num_targets` is the end of the quantification targets: the decoy block,
/// when there is one, is the contiguous tail `[num_targets, num_refs)`. Every
/// expected-bias model and the per-reference effective-length sweep stop
/// there.
///
/// Decoys are never expressed, so an abundance guard already skips them in
/// practice; bounding the range as well is what makes that a *guarantee*
/// rather than a consequence — [`salmon_model::build_expected`]'s own contract
/// says as much.
/// Without it, a decoy that picked up any abundance would drag a whole
/// genome-scale chromosome (up to ~250 Mb) through an O(refLen) sweep, which
/// single-threads the dominant phase of the run: issue #1019, fixed by PR
/// #1020 (#1140 follow-up).
fn apply_bias_correction(
    num_refs: usize,
    num_targets: usize,
    ref_bytes: &salmon_core::RefSeqs,
    gc_store: &salmon_model::GcStore<'_>,
    lengths: &[u32],
    length_class: Option<&[usize]>,
    length_quantiles: Option<&[u32]>,
    fld: &FragmentLengthDistribution,
    alphas: &[f64],
    eff_lengths: &mut [f64],
    seq_obs: Option<(salmon_model::SBModel, salmon_model::SBModel)>,
    gc_obs: Option<salmon_model::GcFragModel>,
    pos_obs: Option<(
        Vec<salmon_model::SimplePosBias>,
        Vec<salmon_model::SimplePosBias>,
    )>,
    seq_bias: bool,
    cond_gc_bins: usize,
    gc_bins: usize,
    bias_speed_samp: usize,
    no_bias_length_threshold: bool,
) -> salmon_model::dumps::BiasDump {
    use salmon_model::seqbias::CONTEXT_LENGTH;
    let mut bias_dump = salmon_model::dumps::BiasDump::default();
    let log_pmf = fld.log_pmf();
    let pmf_lin: Vec<f64> = log_pmf.iter().map(|lp| lp.exp()).collect();
    let (fld_cdf, fld_low, fld_high) = salmon_model::seqbias::fld_cdf_and_bounds(&pmf_lin);
    let k = if seq_bias { CONTEXT_LENGTH } else { 1 };
    let refseq_of = |t: usize| &ref_bytes[t];

    let seq = seq_obs.map(|(mut of, mut or)| {
        of.normalize();
        or.normalize();
        let (ef, er) =
            salmon_model::build_expected(num_targets, refseq_of, alphas, eff_lengths, &fld_cdf);
        (of, or, ef, er)
    });
    if let Some((of, or, ef, er)) = seq.as_ref() {
        bias_dump.obs5_seq = of.dump().to_vec();
        bias_dump.obs3_seq = or.dump().to_vec();
        bias_dump.exp5_seq = ef.dump().to_vec();
        bias_dump.exp3_seq = er.dump().to_vec();
    }
    // Precompute the 5'/3' (obs − exp) log-bias tables once (see the reads-mode
    // path in salmon-quant): the per-position factor build evaluates these
    // rather than both models per context.
    let seq_tab = seq.as_ref().map(|(of, or, ef, er)| {
        (
            salmon_model::LogBiasTable::new(of, ef),
            salmon_model::LogBiasTable::new(or, er),
        )
    });
    let gc_ratio_model = if let Some(mut obs) = gc_obs.filter(gc_obs_has_data) {
        let mut exp = salmon_model::build_expected_gc(
            num_targets,
            refseq_of,
            |t| gc_store.view(t),
            alphas,
            eff_lengths,
            &fld_cdf,
            fld_low,
            fld_high,
            cond_gc_bins,
            gc_bins,
            k,
            bias_speed_samp,
        );
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
    let pos_models = pos_obs.map(|(mut of, mut or)| {
        for x in of.iter_mut().chain(or.iter_mut()) {
            x.finalize();
        }
        let (ef, er) = salmon_model::build_expected_pos(
            num_targets,
            |t| lengths[t] as usize,
            alphas,
            eff_lengths,
            &fld_cdf,
            length_quantiles.expect("positional bias requires length quantiles"),
            k,
        );
        (of, or, ef, er)
    });
    if let Some((ofw, orc, efw, erc)) = pos_models.as_ref() {
        let masses =
            |v: &[salmon_model::SimplePosBias]| v.iter().map(|m| m.masses().to_vec()).collect();
        bias_dump.obs5_pos = masses(ofw);
        bias_dump.obs3_pos = masses(orc);
        bias_dump.exp5_pos = masses(efw);
        bias_dump.exp3_pos = masses(erc);
    }

    // Each transcript's correction is independent. Reads mode already performs
    // this sweep with rayon; doing the same here prevents RAD/alignment mode from
    // serializing its dominant bias-correction phase. Each worker writes one
    // disjoint effective length, and the arithmetic within a transcript is
    // unchanged, so thread count cannot alter the numerical result.
    eff_lengths[..num_refs]
        .par_iter_mut()
        .enumerate()
        .for_each(|(tid, eff_length)| {
            // Decoys keep their uncorrected effective length: they are never
            // reported, so the corrected value would be unused, and computing
            // it is the #1019 stall.
            if tid >= num_targets || alphas[tid] < 1e-8 {
                return;
            }
            let s = &ref_bytes[tid];
            let pos_vecs: Option<(Vec<f64>, Vec<f64>)> =
                pos_models.as_ref().map(|(ofw, orc, efw, erc)| {
                    let lc = length_class.expect("positional bias requires length classes")[tid];
                    let rl = s.len();
                    let (mut o5, mut e5) = (vec![0.0; rl], vec![0.0; rl]);
                    let (mut o3, mut e3) = (vec![0.0; rl], vec![0.0; rl]);
                    ofw[lc].project_weights(&mut o5);
                    efw[lc].project_weights(&mut e5);
                    orc[lc].project_weights(&mut o3);
                    erc[lc].project_weights(&mut e3);
                    (
                        salmon_model::positional_factor(&o5, &e5),
                        salmon_model::positional_factor(&o3, &e3),
                    )
                });
            let bias = salmon_model::BiasInputs {
                seq: seq_tab.as_ref().map(|(f, r)| (f, r, s)),
                gc: gc_ratio_model.as_ref().map(|g| (g, gc_store.view(tid))),
                pos: pos_vecs
                    .as_ref()
                    .map(|(pf, pr)| (pf.as_slice(), pr.as_slice())),
            };
            *eff_length = salmon_model::corrected_effective_length_full(
                s.len(),
                &fld_cdf,
                fld_low,
                fld_high,
                &bias,
                *eff_length,
                bias_speed_samp,
                no_bias_length_threshold,
            );
        });
    bias_dump
}

/// Whether the observed fragment-GC model saw any fragment; warns when it did not.
/// Only proper pairs feed it, so a single-end run leaves it empty, and a ratio over
/// that flat table would "correct" every transcript by the expected model alone.
fn gc_obs_has_data(obs: &salmon_model::GcFragModel) -> bool {
    let has = !obs.is_empty();
    if !has {
        tracing::warn!(
            "--gcBias: no paired fragment was observed (single-end reads?), so fragment-GC \
             correction is skipped"
        );
    }
    has
}

/// Run alignment-based quantification end-to-end.
/// Warn that `-l A` over alignment input can only sample what the aligner
/// chose to report.
///
/// Detection can only see the reported alignments. In reads mode salmon maps
/// everything itself, so the sampled orientations are the library's; here they
/// are whatever survived the aligner's own settings, and an upstream
/// orientation/strand filter skews the sample in a way detection cannot
/// recover from.
///
/// Every alignment-input producer that detects a library type calls this — the
/// RAD writer and genome projection — because the caveat is a property of the
/// input, not of which estimator reads it (COMBINE-lab/salmon#1140, audit C12).
pub(crate) fn warn_auto_detect_from_alignments() {
    tracing::warn!(
        "`-l A` with alignment input infers the library type from the alignments \
         the aligner reported, treating them as an unfiltered sample. If the \
         aligner was configured to report only one orientation or strand, \
         detection will mirror that filter rather than the library — e.g. it \
         cannot conclude `IU` when wrong-strand alignments were already excluded \
         upstream. If the input BAM/SAM is filtered this way, pass the library \
         type explicitly instead."
    );
}

/// Guard the assumption that the references we omit from `quant.sf` hold no
/// abundance.
///
/// Every reference excluded here is a decoy, and decoys are supposed to be
/// unreachable by the estimator: reads mode drops decoy placements before the
/// RAD is written, and piscem never indexes decoys as references at all (its
/// `--decoy-paths` build a poison-k-mer table used to reject reads, so no RAD
/// record can point at a decoy). Both are *upstream* invariants, enforced by
/// code this function cannot see.
///
/// If either is ever violated, the failure is silent and nasty: TPM is
/// normalized over every reference and the rows are filtered afterward, so mass
/// landing on an excluded reference simply disappears — `quant.sf`'s TPM stops
/// summing to 1e6 and its NumReads stops summing to the mapped count, with
/// nothing in the output saying why. This exists so that shows up as a loud
/// warning naming the references, instead of as a quiet discrepancy someone
/// has to reverse-engineer months later (#1140 follow-up; the reasoning is in
/// `rad.rs`'s `num_targets`).
fn warn_mass_on_unreported_refs(res: &AlignQuantResult, rows: &[usize]) {
    if rows.len() == res.names.len() {
        return; // nothing excluded
    }
    let reported: std::collections::HashSet<usize> = rows.iter().copied().collect();
    let mut held = 0.0f64;
    let mut worst: Vec<(String, f64)> = Vec::new();
    for (i, &c) in res.counts.iter().enumerate() {
        if reported.contains(&i) || c <= 0.0 {
            continue;
        }
        held += c;
        if worst.len() < 5 {
            worst.push((res.names[i].clone(), c));
        }
    }
    // A hair of mass can land on a decoy through floating-point redistribution
    // without meaning the invariant broke; a whole fragment cannot.
    if held < 1.0 {
        return;
    }
    let shown = worst
        .iter()
        .map(|(n, c)| format!("{n} ({c:.2})"))
        .collect::<Vec<_>>()
        .join(", ");
    tracing::warn!(
        "{held:.2} fragments of abundance landed on {} reference(s) that are excluded from \
         quant.sf (decoys), so the reported NumReads and TPM columns do not account for them \
         (e.g. {shown}). This should be impossible — decoy placements are dropped before the \
         RAD is written, and piscem does not index decoys as references — so it means an \
         upstream invariant changed. Please report this, with the command line, at \
         https://github.com/COMBINE-lab/salmon/issues",
        res.names.len() - rows.len()
    );
}

impl AlignQuantResult {
    /// The references that get an output row: everything except the
    /// RAD-recorded decoy block.
    pub fn quant_rows(&self) -> Vec<usize> {
        let (first_decoy, num_decoys) = self
            .provenance
            .index
            .as_ref()
            .map(|ix| {
                (
                    ix.first_decoy_index.map(|v| v as usize),
                    ix.num_decoys.unwrap_or(0) as usize,
                )
            })
            .unwrap_or((None, 0));
        salmon_core::quant_row_indices(self.names.len(), first_decoy, num_decoys).collect()
    }
}

fn write_outputs(opts: &AlignQuantOptions, res: &AlignQuantResult) -> Result<()> {
    let dir = &opts.output_dir;
    std::fs::create_dir_all(dir.join("aux_info")).context("creating output dirs")?;

    // Which references get output rows: everything except the RAD-recorded
    // decoy block. A reads-mode RAD's header must carry every reference —
    // records index into it — so the writer, not the header, is where decoys
    // are excluded, exactly as reads mode excludes them (#1140). A RAD without
    // the recorded boundary (piscem, BAM-derived, or written before the
    // boundary tags existed) has no decoys to exclude and emits every row, as
    // before.
    let rows = res.quant_rows();
    warn_mass_on_unreported_refs(res, &rows);

    // quant.sf (EffectiveLength + NumReads at --sigDigits decimals; TPM at 6).
    // Skipped under --skipQuant (no abundances), like salmon.
    if !opts.skip_quant {
        let sd = opts.sig_digits as usize;
        salmon_core::write_via_rename(&dir.join("quant.sf"), |w| {
            writeln!(w, "Name\tLength\tEffectiveLength\tTPM\tNumReads")?;
            for &i in &rows {
                writeln!(
                    w,
                    "{}\t{}\t{:.*}\t{:.6}\t{:.*}",
                    res.names[i],
                    res.lengths[i],
                    sd,
                    res.eff_lengths[i],
                    res.tpm[i],
                    sd,
                    res.counts[i]
                )?;
            }
            Ok(())
        })
        .context("writing quant.sf")?;
    }

    /// Shape of the positional-bias model: one model per transcript length
    /// class, each binned over the transcript's relative position. Reported as
    /// both numbers because either alone is ambiguous.
    #[derive(Serialize)]
    struct PosBiasShape {
        length_classes: usize,
        bins_per_class: usize,
    }

    // aux_info/meta_info.json — full field set, matching salmon's alignment mode
    // (no index hashes since there is no index; no keep_duplicates).
    #[derive(Serialize)]
    struct MetaInfo {
        salmon_version: String,
        samp_type: String,
        opt_type: String,
        quant_errors: Vec<String>,
        num_libraries: usize,
        library_types: Vec<String>,
        frag_dist_length: usize,
        frag_length_mean: f64,
        frag_length_sd: f64,
        /// provenance of the two fields above: whether they came from observed
        /// data or from the `--fldMean`/`--fldSD` prior (see `FragLengthSource`)
        frag_length_source: &'static str,
        seq_bias_correct: bool,
        gc_bias_correct: bool,
        pos_bias_correct: bool,
        num_bias_bins: usize,
        // salmon-rs extensions: the shapes of the bias models that actually
        // ran. `num_bias_bins` above cannot carry these — it names C++'s legacy
        // 4096-bin k-mer counter, which this port replaced — so a reader had no
        // way to learn what the running models looked like. `None` when that
        // correction was not requested, which is how a consumer tells "off"
        // from "on with N bins" (#1140 follow-up).
        #[serde(skip_serializing_if = "Option::is_none")]
        seq_bias_context_length: Option<usize>,
        #[serde(skip_serializing_if = "Option::is_none")]
        gc_bias_bins: Option<usize>,
        #[serde(skip_serializing_if = "Option::is_none")]
        pos_bias_bins: Option<PosBiasShape>,
        mapping_type: String,
        /// how the index was built; known only from a salmon RAD's provenance,
        /// omitted (rather than guessed) for a BAM input
        #[serde(skip_serializing_if = "Option::is_none")]
        keep_duplicates: Option<bool>,
        // empty in alignment mode (no salmon index); recovered from the RAD's
        // provenance when quantifying one salmon wrote
        index_seq_hash: String,
        index_name_hash: String,
        index_seq_hash512: String,
        index_name_hash512: String,
        index_decoy_seq_hash: String,
        index_decoy_name_hash: String,
        num_valid_targets: usize,
        num_decoy_targets: usize,
        num_eq_classes: usize,
        serialized_eq_classes: bool,
        eq_class_properties: Vec<String>,
        length_classes: Vec<u32>,
        num_processed: u64,
        num_mapped: u64,
        num_dovetail_fragments: u64,
        num_fragments_filtered_vm: u64,
        num_alignments_below_threshold_for_mapped_fragments_vm: u64,
        percent_mapped: f64,
        num_decoy_fragments: u64,
        inference_truncated_mass: f64,
        num_bootstraps: u32,
        /// mapped fragments placed as orphans; counted from RAD records, so
        /// omitted for a BAM input where this path does not tally them
        #[serde(skip_serializing_if = "Option::is_none")]
        num_orphan: Option<u64>,
        range_factorization_bins: u32,
        /// whether every field above is a true observation of the run.
        ///
        /// `false` means the RAD quantified from did not record what its mapping
        /// pass observed, so the fields named in `incomplete_meta_info_fields`
        /// are placeholders — most visibly `percent_mapped`, which then reads
        /// 100% because a RAD holds only the fragments that mapped.
        meta_info_complete: bool,
        #[serde(skip_serializing_if = "Vec::is_empty")]
        incomplete_meta_info_fields: Vec<salmon_core::MissingMetaField>,
        /// what the source BAM's `@PG` chain said about how it was produced;
        /// omitted when the fragments did not come from an alignment file
        #[serde(skip_serializing_if = "Vec::is_empty")]
        alignment_provenance: Vec<SourceProgram>,
        num_em_iterations: u32,
        em_converged: bool,
        #[serde(skip_serializing_if = "Option::is_none")]
        em_components: Option<salmon_core::EmComponentsMeta>,
        /// GPU backend that warm-started the EM/VBEM; null when CPU only
        gpu_backend: Option<&'static str>,
        detected_library_type: Option<String>,
        /// which inference path produced these results: `deterministic` (RAD
        /// input — every quantification since 2.12.0). Pre-2.12 one-pass runs
        /// wrote `online`, so the field keeps old outputs distinguishable
        /// (#1140).
        inference_path: &'static str,
        total_time_seconds: f64,
        peak_rss_kb: u64,
        diagnostics: Vec<salmon_core::Diagnostic>,
        call: String,
        start_time: String,
        end_time: String,
    }
    let pct = if res.num_processed > 0 {
        100.0 * res.num_mapped as f64 / res.num_processed as f64
    } else {
        0.0
    };
    let eq_class_properties = if opts.range_factorization_bins > 0 {
        vec!["range_factorized".to_string(), "gzipped".to_string()]
    } else {
        vec!["gzipped".to_string()]
    };
    // A salmon-written RAD records what its mapping pass saw; anything else
    // leaves these unknown, which the reader has already warned about and which
    // `meta_info_complete` records for whoever reads the results later.
    let counters = res.provenance.counters;
    let index_prov = res.provenance.index.clone().unwrap_or_default();
    // Only a RAD requant can be short of the mapping pass's own record of itself.
    // Alignment mode never had index hashes — that is by design, not a gap — and
    // it counts its own fragments, so nothing there is missing for want of an
    // upstream.
    let mut incomplete_meta_info_fields = match res.source {
        FragmentSource::Rad => res.provenance.missing_meta_info_fields(),
        FragmentSource::Bam => Vec::new(),
    };
    // What the aligner recorded about itself, carried through a BAM-derived RAD.
    // A BAM with no `@PG` chain is itself worth recording: the gap is upstream,
    // and saying so is more useful than an empty field the reader must interpret.
    let alignment_provenance = crate::parse_source_programs(&res.source_programs);
    let from_alignments = res.source == FragmentSource::Bam
        || res.provenance.mapping_type == Some(salmon_rad::MappingType::Alignment);
    if from_alignments {
        // `@PG` is optional in SAM and plenty of tools omit it, so its absence is
        // an ordinary outcome — but it is still the record of how the alignments
        // came about, and a reader deserves to be told it was unavailable rather
        // than left to wonder whether salmon simply did not look.
        if alignment_provenance.is_empty() {
            incomplete_meta_info_fields.push(salmon_core::MissingMetaField::bam(
                "alignment_provenance",
                "the source BAM's header carries no @PG records, so how its \
                 alignments were produced is not recorded anywhere in the file",
            ));
        } else if alignment_provenance
            .iter()
            .all(|p| p.command_line.is_none())
        {
            // The chain identifies the programs but not how they were invoked,
            // which is the part that actually reproduces the alignments.
            incomplete_meta_info_fields.push(salmon_core::MissingMetaField::bam(
                "alignment_provenance[].command_line",
                "the source BAM's @PG records name the programs but none carries \
                 a CL field, so the commands that produced the alignments are not \
                 recorded in the file",
            ));
        }
    }
    // lib_format_counts.json — reads mode writes one, so alignment mode and a
    // RAD requant should too, or a pipeline that reads it breaks on the requant
    // of its own output. The struct and its C++-matching field semantics live
    // in `salmon_core::LibFormatCountsFile`, shared with reads mode; the
    // per-format histogram and the judged compat/incompat tallies are measured
    // in both input paths here. For phase 2 of `--deterministic` this rewrite
    // is what puts *measured* values in the final output directory — the
    // mapping pass's file is overwritten here, so these fields must be real,
    // not placeholders.
    let lib_diags: Vec<salmon_core::Diagnostic> = {
        // The expected format is the one the strand filter actually used: the
        // detected type only under `-l A` — the RAD derive pass also detects
        // under an explicit `-l` (for the mismatch diagnostic), and reporting
        // that instead would mislabel the file.
        let expected = if LibraryFormat::is_auto(&opts.lib_type) {
            res.detected_library_type
                .clone()
                .unwrap_or_else(|| opts.lib_type.clone())
        } else {
            LibraryFormat::parse(&opts.lib_type)
                .map(|f| f.canonical().to_string())
                .unwrap_or_else(|_| opts.lib_type.clone())
        };
        let counts = salmon_core::LibFormatCountsFile::new(
            // The reads behind the RAD when the driver knows them (reads-mode
            // --deterministic); a standalone --rad/-a input has none to name.
            opts.read_files.clone(),
            expected,
            res.num_mapped,
            res.num_compatible_fragments,
            res.num_incompatible_fragments,
            &res.lib_format_counts,
        );
        counts.log_warnings();
        std::fs::write(
            dir.join("lib_format_counts.json"),
            serde_json::to_string_pretty(&counts)?,
        )?;
        counts.diagnostics()
    };

    // Warn wherever the mass was actually dropped — this writer is shared by
    // every quantify_rad driver, so the default surfaces the condition it computes instead of only recording it
    // silently in meta_info.json (#1140: the warning previously lived only in
    // the deprecated drivers).
    if res.inference_truncated_mass > 0.0 {
        tracing::warn!(
            "{:.3} fragments of equivalence-class mass could not be assigned (every member \
             transcript was truncated below the min-alpha threshold); reported as \
             inference_truncated_mass in meta_info.json",
            res.inference_truncated_mass
        );
    }
    let meta = MetaInfo {
        salmon_version: SALMON_VERSION.to_string(),
        // What kind of posterior samples were actually produced, so tximport
        // (which keys off this together with num_bootstraps) loads them
        // (#1140; previously hardcoded "none").
        samp_type: if res.bootstraps.is_empty() {
            "none"
        } else if opts.num_bootstraps > 0 {
            "bootstrap"
        } else {
            "gibbs"
        }
        .to_string(),
        opt_type: if opts.em.use_vbem { "vb" } else { "em" }.to_string(),
        quant_errors: vec![],
        num_libraries: 1,
        // The format the strand filter actually applied — same rule as
        // `lib_format_counts.json`'s `expected` above, and for the same reason.
        // The trap: detection now runs under an explicit `-l` too (it powers the
        // `library_type_mismatch` diagnostic), so `detected_library_type` is
        // always populated; taking it unconditionally made `-l ISR` on ISF data
        // report ISF here while `cmd_info.json`, `lib_format_counts.json`, and
        // the (zero) `num_mapped` all said ISR. Metadata reports *both* formats
        // deliberately: the applied one here, the inferred one in
        // `detected_library_type` below — neither is redundant, do not collapse
        // them into one field.
        library_types: vec![if LibraryFormat::is_auto(&opts.lib_type) {
            res.detected_library_type
                .clone()
                .unwrap_or_else(|| opts.lib_type.clone())
        } else {
            LibraryFormat::parse(&opts.lib_type)
                .map(|f| f.canonical().to_string())
                .unwrap_or_else(|_| opts.lib_type.clone())
        }],
        frag_dist_length: res.frag_len_dist.len(),
        frag_length_mean: res.frag_len_mean,
        frag_length_sd: res.frag_len_sd,
        // In the integrated two-phase flows provenance describes the user's
        // run, not the internal RAD handoff: phase 1 of this same run observed
        // the FLD (paired) or defaulted it (single-end, no lengths to observe);
        // the intermediate RAD merely carried it, and the user never asked for
        // that RAD. Which input phase 1 read decides the spelling —
        // `preserve_cmd_info` marks the reads driver, `input_is_alignments` the
        // `-a`/genome-projection drivers. The rad_baked* spellings stay for
        // standalone `--rad`, where the RAD really is the source and neither
        // flag is set (#1140).
        frag_length_source: match (
            opts.preserve_cmd_info,
            opts.input_is_alignments,
            res.frag_len_source,
        ) {
            (true, _, salmon_model::FragLengthSource::RadBaked) => {
                salmon_model::FragLengthSource::Reads
            }
            (_, true, salmon_model::FragLengthSource::RadBaked) => {
                salmon_model::FragLengthSource::Alignments
            }
            // Single-end either way: no fragment lengths existed to observe, so
            // what the RAD carries is this run's own `--fldMean`/`--fldSD`.
            (true, _, salmon_model::FragLengthSource::RadBakedPrior)
            | (_, true, salmon_model::FragLengthSource::RadBakedPrior) => {
                salmon_model::FragLengthSource::Prior
            }
            (_, _, s) => s,
        }
        .as_str(),
        seq_bias_correct: opts.seq_bias,
        gc_bias_correct: opts.gc_bias && !res.bias_dump.obs_gc.is_empty(),
        pos_bias_correct: opts.pos_bias,
        // Not the running bias models' bin counts: C++ writes
        // `readBias(FORWARD).counts.size()` (GZipWriter.cpp), the legacy fixed
        // 4096-bin k-mer counter behind `observed_bias.gz`, which this port
        // replaced with SBModel/GC/pos models and stubs out in the aux dumps.
        // There is no such bin count to report here, and the models that did
        // run publish their own sizes in their dumps, so 0 stands (matching the
        // reads writer) rather than a number of a different model's bins.
        num_bias_bins: 0,
        // Measured from the models this run produced, not from the requested
        // options: a model's own array length is the truth, where a
        // configured-but-unbuilt model would otherwise be described as though
        // it had run.
        seq_bias_context_length: (!res.bias_dump.obs5_seq.is_empty())
            .then_some(salmon_model::seqbias::CONTEXT_LENGTH),
        gc_bias_bins: (!res.bias_dump.obs_gc.is_empty()).then_some(res.bias_dump.obs_gc.len()),
        pos_bias_bins: res.bias_dump.obs5_pos.first().map(|first| PosBiasShape {
            length_classes: res.bias_dump.obs5_pos.len(),
            bins_per_class: first.len(),
        }),
        // What the RAD says its producer was; a BAM input has no RAD to ask, and
        // an older or foreign RAD does not record it, so both fall back to the
        // value this path has always reported.
        mapping_type: res
            .provenance
            .mapping_type
            .map_or("alignment", |m| m.as_str())
            .to_string(),
        keep_duplicates: res
            .provenance
            .index
            .as_ref()
            .and_then(|p| p.keep_duplicates),
        index_seq_hash: index_prov.seq_hash,
        index_name_hash: index_prov.name_hash,
        index_seq_hash512: index_prov.seq_hash512,
        index_name_hash512: index_prov.name_hash512,
        index_decoy_seq_hash: index_prov.decoy_seq_hash,
        index_decoy_name_hash: index_prov.decoy_name_hash,
        num_valid_targets: rows.len(),
        // The decoy block this writer actually held back, derived from the same
        // `rows` split so the two counts always partition the references
        // instead of contradicting each other (previously hardcoded 0 even for
        // a decoy-aware index, while `num_decoy_fragments` reported the real
        // tally). A RAD without the recorded boundary (BAM input, piscem, or
        // written before the tags existed) excludes nothing, so 0 there is this
        // output's true decoy count, not a placeholder — any decoys it does
        // carry are indistinguishable and are counted in `num_valid_targets`,
        // exactly as they appear in `quant.sf`.
        num_decoy_targets: res.names.len() - rows.len(),
        num_eq_classes: res.num_eq_classes,
        serialized_eq_classes: opts.dump_eq || opts.dump_eq_weights,
        eq_class_properties,
        length_classes: res.length_classes.clone(),
        num_processed: res.num_processed,
        num_mapped: res.num_mapped,
        num_dovetail_fragments: counters.map_or(0, |c| c.num_dovetail),
        num_fragments_filtered_vm: counters.map_or(0, |c| c.num_filtered_vm),
        num_alignments_below_threshold_for_mapped_fragments_vm: counters
            .map_or(0, |c| c.num_below_threshold_vm),
        percent_mapped: pct,
        // From the RAD-baked mapping counters, like its three siblings above;
        // 0 when the producing pass predates the counter (#1140; previously
        // hardcoded 0 while the loaded value sat unused).
        num_decoy_fragments: res.provenance.counters.map_or(0, |c| c.num_decoy_fragments),
        inference_truncated_mass: res.inference_truncated_mass,
        num_bootstraps: res.bootstraps.len() as u32,
        num_orphan: res.num_orphan,
        meta_info_complete: incomplete_meta_info_fields.is_empty(),
        incomplete_meta_info_fields,
        alignment_provenance,
        range_factorization_bins: opts.range_factorization_bins,
        num_em_iterations: res.em_iters,
        em_converged: res.em_converged,
        em_components: res.em_components.as_ref().map(|c| c.to_meta()),
        gpu_backend: res.gpu_backend,
        detected_library_type: res.detected_library_type.clone(),
        inference_path: match res.source {
            FragmentSource::Rad => "deterministic",
            FragmentSource::Bam => "online",
        },
        total_time_seconds: res.total_seconds,
        peak_rss_kb: res.peak_rss_kb,
        // The run's own diagnostics, the library-format warnings, and any
        // driver-supplied diagnostics (`extra_diagnostics`, e.g. phase 1's
        // AS-tag warning), so `meta_info.json` carries every machine-readable
        // concern (#1140). Appended here, in the shared writer, so the field's
        // promise holds on every driver.
        diagnostics: {
            let mut d = res.diagnostics.clone();
            d.extend(lib_diags.iter().cloned());
            d.extend(opts.extra_diagnostics.iter().cloned());
            d
        },
        call: "quant".to_string(),
        start_time: res.start_time.clone(),
        end_time: asctime_now(),
    };
    std::fs::write(
        dir.join("aux_info").join("meta_info.json"),
        serde_json::to_string_pretty(&meta)?,
    )?;

    // aux_info/bootstrap/{names.tsv.gz, bootstraps.gz}: posterior samples in the
    // same layout reads mode uses (gzipped tab-separated names; gzipped raw
    // little-endian f64, each sample's `num_refs` values contiguously). Alignment
    // mode has no decoy/short references, so every `quant.sf` row is emitted,
    // keeping the two files aligned positionally with `quant.sf`.
    if !res.bootstraps.is_empty() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        let bdir = dir.join("aux_info").join("bootstrap");
        std::fs::create_dir_all(&bdir).context("creating bootstrap dir")?;

        let f = std::fs::File::create(bdir.join("names.tsv.gz"))?;
        let mut enc = GzEncoder::new(f, Compression::new(6));
        for (j, &i) in rows.iter().enumerate() {
            if j > 0 {
                enc.write_all(b"\t")?;
            }
            enc.write_all(res.names[i].as_bytes())?;
        }
        enc.write_all(b"\n")?;
        enc.finish()?;

        // Restricted to the same rows as quant.sf so the two stay aligned
        // positionally, which is what tximport assumes.
        salmon_model::dumps::write_gz_segments(
            &bdir.join("bootstraps.gz"),
            res.bootstraps.len(),
            |s, buf| {
                buf.reserve(rows.len() * 8);
                for &i in &rows {
                    buf.extend_from_slice(&res.bootstraps[s][i].to_le_bytes());
                }
            },
        )?;
    }

    // libParams/flenDist.txt, logs/salmon_quant.log, and the aux dumps (shared).
    std::fs::create_dir_all(dir.join("libParams")).context("creating libParams")?;
    salmon_model::dumps::write_flen_dist(
        &dir.join("libParams").join("flenDist.txt"),
        &res.frag_len_dist,
    )
    .context("writing flenDist.txt")?;
    salmon_model::dumps::write_fld_dump(&dir.join("aux_info").join("fld.gz"), &res.frag_len_dist)
        .context("writing fld.gz")?;
    salmon_model::dumps::write_aux_bias_dumps(&dir.join("aux_info"), &res.bias_dump)
        .context("writing aux bias dumps")?;
    // `--dumpBiasModels`: the observed+expected bias models, from this pass's
    // own models — the shared writer makes it work on every
    // RAD/alignment driver (#1140: it silently vanished when deterministic
    // became the default, because only salmon-quant had a write site).
    if opts.dump_bias_models {
        salmon_model::dumps::dump_bias_models_to_file(&dir.join("bias_models.txt"), &res.bias_dump)
            .context("writing bias_models.txt")?;
    }
    std::fs::create_dir_all(dir.join("logs")).context("creating logs")?;
    // In alignment mode `num_processed` is every fragment in the input, aligned or
    // not, and `num_mapped` is those with a strand-compatible placement that were
    // quantified; label them accordingly.
    // This writer serves every RAD-based driver, so the log must describe the run
    // the *user* asked for. It used to hard-code "alignment mode" and alignment
    // wording, so a default reads run — the 2.6.0 default — was labelled as
    // something it was not, and its unmapped fragments were reported as
    // "strand-incompatible" even with num_incompatible_fragments = 0. The
    // machine-readable meta_info.json was always right; only this human-facing
    // summary lied (#1140 readiness sweep).
    //
    // `preserve_cmd_info` marks the integrated reads flow, the same signal the
    // invocation record and the fragment-length provenance use.
    let reads_run = opts.preserve_cmd_info;
    let mode = if reads_run {
        "reads mode"
    } else {
        "alignment mode"
    };
    let mapping_type = if reads_run { "mapping" } else { "alignment" };
    // In reads mode a fragment that did not map simply did not map; only
    // alignment input can distinguish "aligned but strand-incompatible".
    let (observed_label, mapped_label, rate_label) = if reads_run {
        ("observed fragments", "mapped fragments", "mapping rate")
    } else {
        (
            "observed fragments",
            "aligned and strand-compatible (quantified)",
            "compatible rate",
        )
    };
    let log = format!(
        "salmon (rust port, {mode}) v{SALMON_VERSION}\nstart: {}\nend:   {}\n\
         library type: {}\nmapping type: {}\n{observed_label}: {}\n{mapped_label}: {}\n\
         {rate_label}: {pct:.4}%\n\
         number of equivalence classes: {}\nfragment length mean (sd): {:.2} ({:.2})\n",
        res.start_time,
        asctime_now(),
        opts.lib_type,
        mapping_type,
        res.num_processed,
        res.num_mapped,
        res.num_eq_classes,
        res.frag_len_mean,
        res.frag_len_sd,
    );
    std::fs::write(dir.join("logs").join("salmon_quant.log"), log)
        .context("writing salmon_quant.log")?;

    // aux_info/ambig_info.tsv — same row set as quant.sf (decoys excluded), so
    // the two stay row-aligned as the format spec promises.
    {
        let (uniq, ambig) = &res.ambig;
        let mut w = std::io::BufWriter::new(std::fs::File::create(
            dir.join("aux_info").join("ambig_info.tsv"),
        )?);
        writeln!(w, "UniqueCount\tAmbigCount")?;
        for &i in &rows {
            writeln!(w, "{}\t{}", uniq[i], ambig[i])?;
        }
        w.flush()?;
    }

    // cmd_info.json — the invocation record.
    #[derive(Serialize)]
    struct CmdInfo {
        salmon_version: String,
        alignments: String,
        targets: String,
        #[serde(rename = "libType")]
        lib_type: String,
        output: String,
        #[serde(rename = "auxDir")]
        aux_dir: String,
    }
    let cmd = CmdInfo {
        salmon_version: SALMON_VERSION.to_string(),
        alignments: opts.bam.display().to_string(),
        targets: opts
            .transcripts
            .as_ref()
            .map(|p| p.display().to_string())
            .unwrap_or_default(),
        lib_type: opts.lib_type.clone(),
        output: opts.output_dir.display().to_string(),
        aux_dir: "aux_info".to_string(),
    };
    // Phase 2 of a two-phase driver keeps phase 1's record of the real
    // invocation rather than replacing it with a RAD-centric one (#1140).
    if !opts.preserve_cmd_info {
        std::fs::write(
            dir.join("cmd_info.json"),
            serde_json::to_string_pretty(&cmd)?,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use salmon_core::is_compatible;

    #[test]
    /// The `AS` tag can be stored in any of BAM's integer widths, so all of them
    /// must decode; a missed case would silently score placements as zero.
    fn alignment_score_value_decoding() {
        assert_eq!(value_as_i32(&Value::Int32(-12)), Some(-12));
        assert_eq!(value_as_i32(&Value::Int8(0)), Some(0));
        assert_eq!(value_as_i32(&Value::UInt16(300)), Some(300));
        // a non-integer tag value (e.g. a character) has no integer reading
        assert_eq!(value_as_i32(&Value::Character(b'A')), None);
    }

    /// Minimal `FragRecord` for classification tests (only the flag-derived
    /// fields matter to `frag_format`'s single-record branch).
    fn frag(is_reverse: bool, is_read1: bool, is_paired: bool) -> FragRecord {
        FragRecord {
            tid: 0,
            pos: 0,
            read_2bit: Vec::new(),
            ops: Vec::new(),
            score: 0,
            has_score: false,
            frag_len: 0,
            is_reverse,
            is_read1,
            is_paired,
            ref_span: 1,
            mate_tid: None,
            mate_pos: None,
            hi: None,
        }
    }

    /// Bases come from a record of the same read, reverse-complemented when the
    /// strands differ, and never from the other mate or a read of another length.
    #[test]
    fn missing_bases_come_from_the_same_read() {
        let with = |bases: &[u8], is_reverse: bool, is_read1: bool| FragRecord {
            read_2bit: bases.to_vec(),
            ops: vec![(AlnOp::Match, 3)],
            ..frag(is_reverse, is_read1, true)
        };
        let mut frags = vec![
            with(&[0, 0, 1], false, true),
            with(&[], false, true),
            with(&[], true, true),
            with(&[], false, false),
            FragRecord {
                ops: vec![(AlnOp::Match, 2)],
                ..with(&[], false, true)
            },
        ];
        fill_missing_bases(&mut frags);
        assert_eq!(frags[1].read_2bit, [0, 0, 1], "same strand");
        assert_eq!(frags[2].read_2bit, [2, 3, 3], "other strand");
        assert!(frags[3].read_2bit.is_empty(), "other mate");
        assert!(frags[4].read_2bit.is_empty(), "other length");
    }

    /// A genuine single-end read (BAM `0x1` unset) is classified `SingleEnd`, so
    /// the strand filters accept it per its own orientation — regression for
    /// issue #1057 (single-end BAM dropped under `SF`/`SR`).
    #[test]
    /// A single-end record must be classified as such and still be subject to the
    /// library-type strand filter.
    fn single_end_record_is_single_end_and_strand_filters_apply() {
        use salmon_core::LibraryFormat;
        let sf = LibraryFormat::parse("SF").unwrap();
        let sr = LibraryFormat::parse("SR").unwrap();
        let u = LibraryFormat::parse("U").unwrap();

        // forward single-end read (flags 0x0)
        let fwd = vec![frag(false, false, false)];
        let (obs, is_fw, status) = frag_format(&fwd, &[0]);
        assert_eq!(status, MateStatus::SingleEnd);
        assert!(is_fw);
        assert!(obs.is_none());
        assert!(is_compatible(sf, obs, is_fw, status)); // SF accepts forward
        assert!(!is_compatible(sr, obs, is_fw, status)); // SR rejects forward
        assert!(is_compatible(u, obs, is_fw, status)); // U accepts

        // reverse single-end read (flags 0x10) — the reporter's case
        let rev = vec![frag(true, false, false)];
        let (obs, is_fw, status) = frag_format(&rev, &[0]);
        assert_eq!(status, MateStatus::SingleEnd);
        assert!(!is_fw);
        assert!(!is_compatible(sf, obs, is_fw, status)); // SF rejects reverse
        assert!(is_compatible(sr, obs, is_fw, status)); // SR accepts reverse
        assert!(is_compatible(u, obs, is_fw, status)); // U accepts
    }

    /// A paired-end read reported alone (BAM `0x1` set) remains an orphan,
    /// classified left/right by the first-segment (`0x40`) flag — unchanged.
    #[test]
    /// An orphan must record *which* mate placed: read 1 and read 2 have mirrored
    /// strand expectations, so losing that would misapply the filter.
    fn paired_orphan_keeps_left_right_status() {
        let read1 = vec![frag(false, true, true)];
        let (_, _, status) = frag_format(&read1, &[0]);
        assert_eq!(status, MateStatus::PairedEndLeft);

        let read2 = vec![frag(true, false, true)];
        let (_, _, status) = frag_format(&read2, &[0]);
        assert_eq!(status, MateStatus::PairedEndRight);
    }

    /// The warning in #1062 must name exactly the flags the user typed, so it
    /// never claims a value was ignored that the user never supplied.
    #[test]
    /// The warning about ignored fragment-length arguments must name only the flags
    /// the user actually passed.
    fn explicit_fld_args_report_only_supplied_flags() {
        assert!(!ExplicitFldArgs::default().any());
        assert_eq!(ExplicitFldArgs::default().names(), "");

        let mean_only = ExplicitFldArgs {
            mean: true,
            ..Default::default()
        };
        assert!(mean_only.any());
        assert_eq!(mean_only.names(), "--fldMean");

        let all = ExplicitFldArgs {
            mean: true,
            sd: true,
            max: true,
        };
        assert_eq!(all.names(), "--fldMean/--fldSD/--fldMax");

        let sd_max = ExplicitFldArgs {
            mean: false,
            sd: true,
            max: true,
        };
        assert_eq!(sd_max.names(), "--fldSD/--fldMax");
    }

    /// Baked is the default, so an ordinary requant keeps reproducing the run
    /// that wrote the RAD unless the user opts out.
    #[test]
    /// The default must reproduce the run that wrote the RAD, rather than
    /// re-deriving a distribution and quietly changing the answer.
    fn fld_policy_defaults_to_baked() {
        assert_eq!(FldPolicy::default(), FldPolicy::Baked);
        let opts = AlignQuantOptions::new("x.rad".into(), "out".into());
        assert_eq!(opts.fld_policy, FldPolicy::Baked);
        assert!(!opts.explicit_fld_args.any());
    }

    /// `-t` is only read when the error or bias models need it, so a file with no
    /// records in it is a mistake. Unchecked it is silent: every reference
    /// resolves to an empty sequence and the models contribute nothing, while the
    /// run reports success.
    #[test]
    fn empty_reference_fasta_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let names = vec!["t0".to_string(), "t1".to_string()];

        for (label, contents) in [
            ("empty file", ""),
            ("no FASTA records", "just some text\nwith no header line\n"),
        ] {
            let p = dir.path().join("ref.fa");
            std::fs::write(&p, contents).unwrap();
            let err = load_ref_bytes(&p, &names, &[]).unwrap_err();
            let msg = format!("{err:#}");
            assert!(
                msg.contains("no sequences found") && msg.contains("ref.fa"),
                "{label}: the error should say what is wrong and name the file, got: {msg}"
            );
        }
    }

    /// A read failure surfaces mid-stream rather than at open, so the path has to
    /// be attached there too or a truncated file reports a bare
    /// "unexpected end of file".
    #[test]
    fn read_errors_name_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("truncated.fa.gz");
        // A gzip header and the start of a member, cut off before the end.
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        use std::io::Write;
        enc.write_all(b">t0\nACGTACGTACGTACGTACGTACGTACGT\n")
            .unwrap();
        let full = enc.finish().unwrap();
        std::fs::write(&p, &full[..full.len() - 6]).unwrap();

        let err = load_ref_bytes(&p, &["t0".to_string()], &[]).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("truncated.fa.gz"),
            "a mid-stream failure should still name the file, got: {msg}"
        );
    }

    #[test]
    fn a_wrong_length_is_a_mismatch_but_absent_or_a_clipped_tail_is_not() {
        let names: Vec<String> = ["ok", "absent", "tail", "short", "long"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let seqs: Vec<Vec<u8>> = vec![
            b"ACGT".to_vec(),
            vec![],
            b"ACGTAAaa".to_vec(),
            b"AC".to_vec(),
            b"ACGTC".to_vec(),
        ];
        let (wrong, first) = length_mismatches(&names, &seqs, &[4, 4, 4, 4, 4]);
        assert_eq!(wrong, 2);
        assert_eq!(first, Some(("short", 2, 4)));
    }

    /// Partial coverage is legitimate — a transcriptome need not cover every
    /// reference in the alignment header — so it must not be an error.
    #[test]
    fn partial_reference_coverage_is_accepted() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ref.fa");
        std::fs::write(&p, ">t0\nACGT\n").unwrap();
        let got = load_ref_bytes(&p, &["t0".to_string(), "absent".to_string()], &[]).unwrap();
        assert_eq!(got[0], b"ACGT");
        assert!(
            got[1].is_empty(),
            "an absent reference yields an empty sequence"
        );
    }
}

#[cfg(test)]
mod logsumexp_tests {
    use super::{logsumexp, logsumexp_iter, logsumexp_one};

    /// The single-value shortcut has to be the general form, bit for bit, on
    /// every input — signed zeros and non-finite values included.
    #[test]
    fn one_value_shortcut_agrees_with_the_general_form() {
        let mut cases = vec![
            0.0,
            -0.0,
            1.0,
            -1.0,
            -1e-300,
            f64::MIN_POSITIVE,
            -f64::MIN_POSITIVE / 4.0,
            f64::MAX,
            f64::MIN,
            -745.2,
            -23.998_158_637_57,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ];
        let mut x = 0x9E37_79B9_7F4A_7C15u64;
        for _ in 0..10_000 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            cases.push(f64::from_bits(x));
            cases.push(-((x % 1_000_000) as f64) / 997.0);
        }
        for x in cases {
            let general = logsumexp_iter(std::iter::once(x));
            let one = logsumexp_one(x);
            if general.is_nan() {
                assert!(one.is_nan(), "x={x:e}: general NaN, shortcut {one}");
            } else {
                assert_eq!(general.to_bits(), one.to_bits(), "x={x:e}");
            }
        }
    }

    /// `logsumexp_iter` exists so a caller holding a run of records can avoid
    /// copying their weights into a scratch slice. It is only safe to use in
    /// place of [`logsumexp`] if it agrees *bit for bit* — the eq-class weights
    /// it produces feed range factorization, where a last-bit difference can
    /// move a fragment across a bin boundary and change the class set.
    #[test]
    fn iter_form_agrees_bit_for_bit_with_the_slice_form() {
        let cases: &[&[f64]] = &[
            &[-1.0],
            &[-1.0, -1.0],
            &[0.0, -1.0, -2.0, -3.0],
            &[-1e-16, -1.0, -1e-16],
            &[-700.0, -700.5, -701.25],
            &[f64::NEG_INFINITY, -2.0],
            &[f64::NEG_INFINITY, f64::NEG_INFINITY],
            &[-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7],
        ];
        for xs in cases {
            let a = logsumexp(xs);
            let b = logsumexp_iter(xs.iter().copied());
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "slice={a} iter={b} for {xs:?} — the two forms must not diverge"
            );
        }
    }

    /// Order matters to the result, which is why the grouping sorts stably; the
    /// two forms must agree on that too, not merely on sorted input.
    #[test]
    fn iter_form_matches_on_reordered_input() {
        let xs = [-1e-16, -1.0, -1e-16, -2.0];
        let mut rev = xs;
        rev.reverse();
        assert_eq!(
            logsumexp(&rev).to_bits(),
            logsumexp_iter(rev.iter().copied()).to_bits()
        );
    }
}
