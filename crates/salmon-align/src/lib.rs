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
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_sam::Header;
use rayon::prelude::*;
use serde::Serialize;

use error_model::{AlignmentModel, AlnOp, SharedAlignmentModel};
use salmon_core::{
    is_compatible, LibraryFormat, MateStatus, PhaseTimer, ReadOrientation, ReadStrandedness,
    ReadType,
};
use salmon_eqclass::{range_factorize_bins, EquivalenceClassBuilder, TranscriptGroup};
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
    /// online-phase forgetting factor (`--forgettingFactor`, salmon default 0.65)
    pub forgetting_factor: f64,
    /// initialize the EM uniformly instead of with the online-estimate-blended
    /// warm start (`--initUniform`). Toggles behavior only in online alignment
    /// mode ([`quantify_alignments`]), the one path with a warm start to
    /// replace; the RAD quantifier always starts uniform, so there the flag is
    /// honoured by construction.
    pub init_uniform: bool,
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
    /// `--noLengthCorrection`: use the raw reference length as the effective
    /// length, instead of the fragment-length-derived one.
    ///
    /// For protocols where a fragment does not sample a transcript's interior
    /// uniformly (3' tagged-end libraries above all), the usual correction
    /// describes something the data never did, and the raw length is the honest
    /// divisor. Bias correction is disabled with it, for the same reason: its
    /// whole job is to adjust an effective length this mode is not computing.
    pub no_length_correction: bool,
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
    /// path. The deprecated online alignment path always behaves as though
    /// this were false (see the FLD note in `weigh_fragment`), so this reaches
    /// only the RAD quantifier — which is every default path (#1140, D14).
    pub model_single_frag_prob: bool,
    /// diagnostics born in an earlier phase the driver ran (e.g. the
    /// deterministic alignment pass's unusable-`AS` verdict), appended to this
    /// pass's own in `meta_info.json` — so a pipeline that only keeps the
    /// output directory sees them, not just whoever read the log (#1140)
    pub extra_diagnostics: Vec<salmon_core::Diagnostic>,
    /// `--dumpBiasModels` (hidden): write the observed+expected seq/GC/pos
    /// bias models to `bias_models.txt` for debugging / C++-parity comparison.
    /// Written by the shared output writer from the pass's [`BiasDump`], so it
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
    /// online-phase auxiliary-model training window (`--numAuxModelSamples`,
    /// default 5,000,000)
    pub num_aux_model_samples: u64,
    /// disable the lower barrier on bias-corrected effective lengths
    /// (`--noBiasLengthThreshold`)
    pub no_bias_length_threshold: bool,
    /// GC bias model bin counts (`--numGCBins` × `--conditionalGCBins`, 25×3)
    pub gc_bins: usize,
    pub cond_gc_bins: usize,
    /// skip abundance estimation + `quant.sf`, emitting only eq-classes,
    /// library type, and metadata (salmon's `--skipQuant`)
    pub skip_quant: bool,
    /// fragments processed before the FLD aux model is applied
    /// (salmon's `--numPreAuxModelSamples`; prior hardcoded value 5,000)
    pub num_pre_aux_model_samples: u64,
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
            forgetting_factor: 0.65,
            init_uniform: false,
            prior_seconds: 0.0,
            external_start_time: None,
            preserve_cmd_info: false,
            input_is_alignments: false,
            read_files: Vec::new(),
            dump_eq: false,
            dump_eq_weights: false,
            no_length_correction: false,
            no_frag_length_dist: false,
            hard_filter: false,
            model_single_frag_prob: true,
            extra_diagnostics: Vec::new(),
            dump_bias_models: false,
            sig_digits: 3,
            num_error_bins: 4,
            discard_orphans: false,
            bias_speed_samp: 5,
            num_aux_model_samples: 5_000_000,
            no_bias_length_threshold: false,
            gc_bins: salmon_model::gcbias::DEFAULT_GC_BINS,
            cond_gc_bins: salmon_model::gcbias::DEFAULT_COND_BINS,
            skip_quant: false,
            num_pre_aux_model_samples: 5_000,
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
    pub ambig: (Vec<u32>, Vec<u32>),
    /// posterior samples (bootstrap or Gibbs), one abundance vector each; empty
    /// when neither was requested. Length is `num_refs`, matching `quant.sf` rows.
    pub bootstraps: Vec<Vec<f64>>,
    /// EM/VBEM convergence: iterations run and whether the relative-difference
    /// criterion was met before `max_iter`
    pub em_iters: u32,
    pub em_converged: bool,
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

/// Current local time as an asctime-style string, matching salmon's timestamps.
pub fn asctime_now() -> String {
    jiff::Zoned::now()
        .strftime("%a %b %e %H:%M:%S %Y")
        .to_string()
}

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
fn load_ref_bytes(path: &Path, names: &[String]) -> Result<Vec<Vec<u8>>> {
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

    let seqs: Vec<Vec<u8>> = names
        .iter()
        .map(|n| by_name.remove(n).unwrap_or_default())
        .collect();
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
    Ok(seqs)
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

impl FragRecord {
    /// The read's 5' reference position: its leftmost if forward, its rightmost
    /// if reverse-complemented (salmon's `startPos`).
    #[inline]
    fn five_prime(&self) -> usize {
        if self.is_reverse {
            self.pos + self.ref_span.saturating_sub(1)
        } else {
            self.pos
        }
    }
}

/// The canonical read name: a trailing `/1` or `/2` mate suffix stripped (as
/// salmon's `getPairedNameLen` does), so the two mates of a fragment group
/// together even when the aligner kept the suffix in the QNAME.
#[inline]
fn canonical_name(name: &[u8]) -> &[u8] {
    let l = name.len();
    if l > 2 && name[l - 2] == b'/' {
        &name[..l - 2]
    } else {
        name
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

/// Does this path name a SAM file rather than a BAM, whatever it is compressed
/// with?
///
/// SAM-vs-BAM is a name question, and compression is a content question, so any
/// compression suffix is stripped before asking: `aln.sam.zst` is a SAM file
/// that happens to be zstd, and routing it to the BAM reader produced a report
/// about an invalid BGZF header for a file that was never BGZF.
fn is_sam_path(p: &Path) -> bool {
    let s = p.to_string_lossy().to_ascii_lowercase();
    salmon_core::compress::strip_compression_extension(&s).ends_with(".sam")
}

/// Open a SAM file (transparently decompressing it) as a noodles reader over a
/// boxed `BufRead`, so every input shares one concrete type.
///
/// Compression is detected from content, so the stream is decoded correctly
/// however it was compressed. Which files reach here at all is a separate,
/// name-based question, answered by [`is_sam_path`].
fn open_sam_reader(path: &Path) -> Result<sam::io::Reader<Box<dyn BufRead + Send>>> {
    let inner = salmon_core::compress::open_maybe_compressed(path)
        .with_context(|| format!("opening {}", path.display()))?;
    Ok(sam::io::Reader::new(inner))
}

/// Read just the header from a SAM/BAM input (chosen by extension).
fn read_alignment_header(path: &Path) -> Result<Header> {
    if is_sam_path(path) {
        open_sam_reader(path)?
            .read_header()
            .context("reading SAM header")
    } else {
        let mut reader = bam::io::Reader::new(
            std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?,
        );
        reader.read_header().context("reading BAM header")
    }
}

/// Convert one mapped alignment record into a `FragRecord`, returning its
/// canonical read name alongside. `None` skips the record (unmapped / unnamed /
/// no reference). Works over any noodles `alignment::Record`, so the same logic
/// serves BAM and SAM input.
fn record_to_frag<R: sam::alignment::Record>(
    record: &R,
    header: &Header,
    need_seq: bool,
) -> Option<(Vec<u8>, FragRecord)> {
    let flags = record.flags().ok()?;
    if flags.is_unmapped() {
        return None;
    }
    let name = record.name()?;
    let cname = canonical_name(name.as_ref()).to_vec();

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
    let ops: Vec<(AlnOp, usize)> = record
        .cigar()
        .iter()
        .filter_map(|r| r.ok())
        .map(|op| (kind_to_op(op.kind()), op.len()))
        .collect();
    let ref_span: usize = ops
        .iter()
        .filter(|(o, _)| {
            matches!(
                o,
                AlnOp::Match | AlnOp::SeqMatch | AlnOp::SeqMismatch | AlnOp::Del | AlnOp::RefSkip
            )
        })
        .map(|(_, l)| l)
        .sum();
    let ops = if need_seq { ops } else { Vec::new() };
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
    Some((
        cname,
        FragRecord {
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
        },
    ))
}

/// Read-only configuration + shared (thread-safe) sinks for the online pass.
/// Held across the whole pass; the worker threads borrow it immutably.
struct PassCfg<'a> {
    online: Option<&'a salmon_infer::OnlineInference>,
    fld: &'a FragmentLengthDistribution,
    eq_builder: &'a EquivalenceClassBuilder,
    ref_bytes: &'a salmon_core::RefSeqs,
    lengths: &'a [u32],
    gc_store: salmon_model::GcStore<'a>,
    length_class: Option<&'a [usize]>,
    expected_format: Option<LibraryFormat>,
    /// shared library-type detector under `-l A` (all-atomic, fed by every
    /// worker); `None` with an explicit library type
    detector: Option<&'a salmon_model::LibraryTypeDetector>,
    ignore_incompat: bool,
    incompat_prior: f64,
    paired_lib: bool,
    /// drop single-mate (orphan) placements in a paired library (`--discardOrphans`)
    discard_orphans: bool,
    /// fragments processed before the FLD aux model is applied (`--numPreAuxModelSamples`)
    pre_burnin: u64,
    range_factorization_bins: u32,
    use_error_model: bool,
    /// number of read-position bins in the alignment error model
    /// (salmon's `--numErrorBins`)
    error_bins: usize,
    seq_bias: bool,
    gc_bias: bool,
    /// GC bias model bin counts (`--conditionalGCBins` × `--numGCBins`)
    cond_gc_bins: usize,
    gc_bins: usize,
    pos_bias: bool,
    need_seq: bool,
    minibatch: usize,
    nthreads: usize,
    /// optional shared live-progress counters (updated per batch)
    progress: Option<&'a salmon_core::ProgressCounters>,
}

/// Outputs of the online pass: the observed bias models (merged from the
/// per-worker accumulators) and the processed/mapped fragment counts. The error
/// model is developed in a shared atomic structure internal to the pass and is
/// not needed afterward.
struct PassAccum {
    seq_obs: Option<(salmon_model::SBModel, salmon_model::SBModel)>,
    gc_obs: Option<salmon_model::GcFragModel>,
    pos_obs: Option<(
        Vec<salmon_model::SimplePosBias>,
        Vec<salmon_model::SimplePosBias>,
    )>,
    num_processed: u64,
    num_mapped: u64,
    /// mapped fragments where every surviving placement placed only one mate
    num_orphan: u64,
    /// fragments judged against a known expected format with ≥1 strand-compatible
    /// placement / with none (#1130); see [`AlignQuantResult`]'s fields
    num_frags_compat: u64,
    num_frags_incompat: u64,
    /// per-observed-format fragment counts; see [`AlignQuantResult`]'s field
    lib_format_counts: salmon_core::LibFormatCountsArray,
}

/// The online pass over an alignment record stream, structured as a persistent
/// producer/worker pool (mirroring salmon's parse-threads → queue →
/// quant-threads model).
///
/// One reader thread does only the cheap work — deserialize each record and
/// group consecutive mapped records by read name into raw fragment groups — and
/// pushes minibatches onto a bounded work queue (measured serial floor for this
/// is ~5% of the pass). `nthreads` persistent worker threads pull minibatches
/// continuously (no per-batch barrier; the MPMC queue load-balances), and each
/// does the expensive `record_to_frag` (2-bit encoding, CIGAR/op extraction,
/// tag parsing, allocation) *and* the per-fragment weighting. Workers read the
/// shared atomic error model for `basis` and flush their own model delta into it
/// after each minibatch; bias accumulators are per-worker and merged once at the
/// end (they are only read after the pass).
fn run_online_pass<R, I>(
    records: I,
    header: &Header,
    cfg: &PassCfg,
    acc: &mut PassAccum,
) -> Result<()>
where
    R: sam::alignment::Record + Send + Sync,
    I: Iterator<Item = std::io::Result<R>> + Send,
{
    use crossbeam_channel::bounded;

    // Fragments a worker processes between flushing its error-model delta into
    // the shared model. Small = fresher shared model (closer to salmon's live
    // per-transition training, better parity) but more atomic contention on the
    // hot match cell; 1 = per-fragment.
    const FLUSH_INTERVAL: usize = 64;

    // Shared atomic error model: read concurrently for `basis`, updated by the
    // workers flushing their per-thread deltas into it between minibatches.
    let shared_model = cfg
        .use_error_model
        .then(|| SharedAlignmentModel::new(1.0, cfg.error_bins));
    let shared_model_ref = shared_model.as_ref();
    let minibatch = cfg.minibatch;
    let need_seq = cfg.need_seq;
    // Each batch carries its forgetting-mass step, assigned by the reader in
    // batch order (so the online schedule is tied to the minibatch index, as in
    // salmon, rather than to nondeterministic worker-pull order).
    let (tx, rx) = bounded::<(f64, Vec<Vec<R>>)>(2 * cfg.nthreads + 1);
    let online_reader = cfg.online;

    std::thread::scope(|scope| -> Result<()> {
        // Reader/parser thread.
        let reader = scope.spawn(move || -> Result<()> {
            let mut cur_name: Vec<u8> = Vec::new();
            let mut have = false;
            let mut group: Vec<R> = Vec::new();
            let mut batch: Vec<Vec<R>> = Vec::with_capacity(minibatch);
            for result in records {
                let rec = result.context("reading alignment record")?;
                if rec.flags().map(|f| f.is_unmapped()).unwrap_or(true) {
                    continue;
                }
                let Some(name) = rec.name() else { continue };
                let cname = canonical_name(name.as_ref()).to_vec();
                if !have {
                    cur_name = cname;
                    have = true;
                } else if cname != cur_name {
                    batch.push(std::mem::take(&mut group));
                    cur_name = cname;
                    if batch.len() >= minibatch {
                        let fm = online_reader.map(|o| o.next_log_fm()).unwrap_or(0.0);
                        if tx.send((fm, std::mem::take(&mut batch))).is_err() {
                            return Ok(()); // all workers gone
                        }
                        batch = Vec::with_capacity(minibatch);
                    }
                }
                group.push(rec);
            }
            if have && !group.is_empty() {
                batch.push(group);
            }
            if !batch.is_empty() {
                let fm = online_reader.map(|o| o.next_log_fm()).unwrap_or(0.0);
                let _ = tx.send((fm, batch));
            }
            Ok(())
        });

        // Persistent worker pool: each pulls minibatches until the queue closes.
        let mut workers = Vec::with_capacity(cfg.nthreads);
        for _ in 0..cfg.nthreads {
            let rx = rx.clone();
            workers.push(scope.spawn(
                move || -> (
                    Local,
                    u64,
                    u64,
                    u64,
                    u64,
                    u64,
                    salmon_core::LibFormatCountsArray,
                ) {
                let mut local = Local::new(
                    cfg.use_error_model,
                    cfg.error_bins,
                    cfg.seq_bias,
                    cfg.gc_bias,
                    cfg.cond_gc_bins,
                    cfg.gc_bins,
                    cfg.pos_bias,
                );
                let mut count = 0u64;
                let mut mapped = 0u64;
                let mut orphan = 0u64;
                let mut compat = 0u64;
                let mut incompat = 0u64;
                let mut lib_fmt: salmon_core::LibFormatCountsArray =
                    [0; salmon_core::NUM_LIB_FORMATS];
                let mut frags: Vec<FragRecord> = Vec::new();
                while let Ok((log_fm, raw_batch)) = rx.recv() {
                    let ctx = FragCtx {
                        model: shared_model_ref,
                        online: cfg.online,
                        fld: cfg.fld,
                        eq_builder: cfg.eq_builder,
                        ref_bytes: cfg.ref_bytes,
                        lengths: cfg.lengths,
                        gc_store: cfg.gc_store,
                        length_class: cfg.length_class,
                        expected_format: cfg.expected_format,
                        detector: cfg.detector,
                        ignore_incompat: cfg.ignore_incompat,
                        incompat_prior: cfg.incompat_prior,
                        paired_lib: cfg.paired_lib,
                        discard_orphans: cfg.discard_orphans,
                        pre_burnin: cfg.pre_burnin,
                        range_factorization_bins: cfg.range_factorization_bins,
                        log_fm,
                    };
                    let mut since_flush = 0usize;
                    let mut batch_mapped = 0u64;
                    // Plain locals, summed once per worker at the end: the
                    // per-fragment cost is an increment, with no atomic traffic.
                    let mut batch_orphan = 0u64;
                    let mut batch_compat = 0u64;
                    let mut batch_incompat = 0u64;
                    for raw_group in &raw_batch {
                        frags.clear();
                        for r in raw_group {
                            if let Some((_, f)) = record_to_frag(r, header, need_seq) {
                                frags.push(f);
                            }
                        }
                        let (outcome, strand_judged, fmt_mask) =
                            process_fragment(&frags, &ctx, &mut local);
                        match strand_judged {
                            Some(true) => batch_compat += 1,
                            Some(false) => batch_incompat += 1,
                            None => {}
                        }
                        // Observed-format histogram (one count per distinct
                        // format among the fragment's placements).
                        let mut bits = fmt_mask;
                        while bits != 0 {
                            let i = bits.trailing_zeros() as usize;
                            lib_fmt[i] += 1;
                            bits &= bits - 1;
                        }
                        match outcome {
                            FragmentOutcome::Unmapped => {}
                            FragmentOutcome::Mapped => batch_mapped += 1,
                            FragmentOutcome::MappedOrphan => {
                                batch_mapped += 1;
                                batch_orphan += 1;
                            }
                        }
                        // Publish the error-model delta into the shared model
                        // every FLUSH_INTERVAL fragments so other workers' `basis`
                        // sees fresh-enough training. The update granularity is
                        // what governs parity with salmon's live per-transition
                        // model: ~per-fragment freshness recovers it, while the
                        // batched flush keeps hot-cell atomic contention low.
                        since_flush += 1;
                        if since_flush >= FLUSH_INTERVAL {
                            if let (Some(sm), Some(d)) = (shared_model_ref, local.model.as_mut()) {
                                sm.flush_from(d);
                                d.clear();
                            }
                            since_flush = 0;
                        }
                    }
                    if since_flush > 0 {
                        if let (Some(sm), Some(d)) = (shared_model_ref, local.model.as_mut()) {
                            sm.flush_from(d);
                            d.clear();
                        }
                    }
                    count += raw_batch.len() as u64;
                    mapped += batch_mapped;
                    orphan += batch_orphan;
                    compat += batch_compat;
                    incompat += batch_incompat;
                    // Live progress: every fragment in the BAM is "processed";
                    // only fragments with a surviving strand-compatible placement
                    // are "mapped" (so percent_mapped is correct on stranded data).
                    if let Some(p) = cfg.progress {
                        p.processed.fetch_add(
                            raw_batch.len() as u64,
                            std::sync::atomic::Ordering::Relaxed,
                        );
                        p.mapped
                            .fetch_add(batch_mapped, std::sync::atomic::Ordering::Relaxed);
                    }
                }
                (local, count, mapped, orphan, compat, incompat, lib_fmt)
            }));
        }
        drop(rx); // workers hold their own clones; lets the queue disconnect

        reader
            .join()
            .map_err(|_| anyhow::anyhow!("alignment reader thread panicked"))??;

        // Merge the per-worker bias accumulators + counts.
        let mut merged = Local::new(
            cfg.use_error_model,
            cfg.error_bins,
            cfg.seq_bias,
            cfg.gc_bias,
            cfg.cond_gc_bins,
            cfg.gc_bins,
            cfg.pos_bias,
        );
        let mut total = 0u64;
        let mut total_mapped = 0u64;
        let mut total_orphan = 0u64;
        let mut total_compat = 0u64;
        let mut total_incompat = 0u64;
        let mut total_lib_fmt: salmon_core::LibFormatCountsArray =
            [0; salmon_core::NUM_LIB_FORMATS];
        for w in workers {
            let (local, count, mapped, orphan, compat, incompat, lib_fmt) = w
                .join()
                .map_err(|_| anyhow::anyhow!("alignment worker thread panicked"))?;
            merged = merged.merge(local);
            total += count;
            total_mapped += mapped;
            total_orphan += orphan;
            total_compat += compat;
            total_incompat += incompat;
            for (t, v) in total_lib_fmt.iter_mut().zip(lib_fmt.iter()) {
                *t += v;
            }
        }
        acc.seq_obs = merged.seq_obs;
        acc.gc_obs = merged.gc_obs;
        acc.pos_obs = merged.pos_obs;
        // num_processed = every aligned fragment in the BAM; num_mapped = those
        // with a surviving strand-compatible placement (assigned to an eq-class).
        // They differ only for stranded libraries; matches reads-mode (#1025).
        acc.num_processed = total;
        acc.num_mapped = total_mapped;
        acc.num_orphan = total_orphan;
        acc.num_frags_compat = total_compat;
        acc.num_frags_incompat = total_incompat;
        acc.lib_format_counts = total_lib_fmt;
        Ok(())
    })
}

/// Dispatch the online pass over a SAM/BAM file (chosen by extension). BAM is
/// decoded on a small BGZF worker pool (the framing/parse stays serial on the
/// reader thread, but the decompression overlaps it).
/// Whether the file's first record belongs to a paired-end template (C++
/// salmon's `peekBAMIsPaired`): under `-l A` the user has not said, and the
/// library-type detector needs to know which formats are candidates. `None`
/// for a file with no records.
fn peek_alignment_is_paired(path: &Path) -> Result<Option<bool>> {
    fn first_flag<R: sam::alignment::Record, E>(
        records: impl Iterator<Item = std::result::Result<R, E>>,
    ) -> Option<bool>
    where
        E: std::fmt::Debug,
    {
        for rec in records.flatten() {
            if let Ok(flags) = rec.flags() {
                return Some(flags.is_segmented());
            }
        }
        None
    }
    if is_sam_path(path) {
        let mut reader = open_sam_reader(path)?;
        let _header = reader.read_header().context("reading SAM header")?;
        Ok(first_flag(reader.records()))
    } else {
        let file =
            std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut reader = bam::io::Reader::new(file);
        let _header = reader.read_header().context("reading BAM header")?;
        Ok(first_flag(reader.records()))
    }
}

fn stream_online_pass(bam_path: &Path, cfg: &PassCfg, acc: &mut PassAccum) -> Result<()> {
    if is_sam_path(bam_path) {
        let mut reader = open_sam_reader(bam_path)?;
        let header = reader.read_header().context("reading SAM header")?;
        run_online_pass(reader.records(), &header, cfg, acc)
    } else {
        let file = std::fs::File::open(bam_path)
            .with_context(|| format!("opening {}", bam_path.display()))?;
        let workers = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4)
            .clamp(1, 4);
        let workers = std::num::NonZeroUsize::new(workers).unwrap();
        {
            let decoder = noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, file);
            let mut reader = bam::io::Reader::from(decoder);
            let header = reader.read_header().context("reading BAM header")?;
            run_online_pass(reader.records(), &header, cfg, acc)
        }
    }
}

/// Read-only shared state for processing one fragment. Held by `&` so a whole
/// minibatch of fragments can be processed in parallel against the model and
/// abundance state as of the *previous* batch (salmon's minibatch staleness).
struct FragCtx<'a> {
    /// shared error model, read concurrently for `basis` (workers flush their
    /// own deltas into it between minibatches)
    model: Option<&'a SharedAlignmentModel>,
    online: Option<&'a salmon_infer::OnlineInference>,
    fld: &'a FragmentLengthDistribution,
    eq_builder: &'a EquivalenceClassBuilder,
    ref_bytes: &'a salmon_core::RefSeqs,
    lengths: &'a [u32],
    gc_store: salmon_model::GcStore<'a>,
    length_class: Option<&'a [usize]>,
    expected_format: Option<LibraryFormat>,
    /// see [`PassCfg::detector`]
    detector: Option<&'a salmon_model::LibraryTypeDetector>,
    ignore_incompat: bool,
    incompat_prior: f64,
    paired_lib: bool,
    /// drop single-mate (orphan) placements in a paired library (`--discardOrphans`)
    discard_orphans: bool,
    /// fragments processed before the FLD aux model is applied (`--numPreAuxModelSamples`)
    pre_burnin: u64,
    range_factorization_bins: u32,
    /// this batch's forgetting mass (online phase)
    log_fm: f64,
}

/// Per-thread accumulators for the error model and bias models. Each worker
/// folds its fragments into a `Local` (the error-model matrices seeded empty so
/// the pseudocount lives only in the global); the per-batch reduction merges
/// these and the result is folded into the globals between minibatches.
struct Local {
    model: Option<AlignmentModel>,
    seq_obs: Option<(salmon_model::SBModel, salmon_model::SBModel)>,
    gc_obs: Option<salmon_model::GcFragModel>,
    pos_obs: Option<(
        Vec<salmon_model::SimplePosBias>,
        Vec<salmon_model::SimplePosBias>,
    )>,
}

impl Local {
    #[allow(clippy::too_many_arguments)]
    fn new(
        error_model: bool,
        error_bins: usize,
        seq_bias: bool,
        gc_bias: bool,
        cond_gc_bins: usize,
        gc_bins: usize,
        pos_bias: bool,
    ) -> Self {
        let mk = || {
            (0..salmon_model::NUM_LENGTH_CLASSES)
                .map(|_| salmon_model::SimplePosBias::default())
                .collect::<Vec<_>>()
        };
        Self {
            model: error_model.then(|| AlignmentModel::empty(error_bins)),
            seq_obs: seq_bias.then(|| (salmon_model::SBModel::new(), salmon_model::SBModel::new())),
            gc_obs: gc_bias.then(|| salmon_model::GcFragModel::new(cond_gc_bins, gc_bins)),
            pos_obs: pos_bias.then(|| (mk(), mk())),
        }
    }

    fn merge(mut self, other: Self) -> Self {
        if let (Some(a), Some(b)) = (self.model.as_mut(), other.model.as_ref()) {
            a.combine(b);
        }
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
        self
    }
}

/// Process one fragment (a group of records sharing a read name): pair its
/// records into reported placements, compute each placement's conditional
/// log-weight, build the equivalence class, develop online abundances, and
/// (during burn-in) accumulate the error-model and bias deltas into `local`.
/// Pure with respect to shared state except for the concurrency-safe sinks
/// (`fld`, `online`, `eq_builder`), so it is safe to run in parallel.
/// Returns `true` if the fragment was assigned (had at least one surviving,
/// strand-compatible placement and joined an equivalence class), `false` if it
/// was dropped (every reported alignment incompatible / orphan-discarded). The
/// caller counts a fragment as *mapped* only when this returns `true`, so a
/// stranded library does not over-report `num_mapped` for fragments whose every
/// alignment is strand-incompatible (the alignment-mode analog of #1025).
/// What became of one fragment.
///
/// Distinguishes the two mapped cases so `-a` can report the concordant/orphan
/// split the reads path does, without a second pass over the placements.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum FragmentOutcome {
    /// no reported alignment survived compatibility filtering
    Unmapped,
    /// assigned, with at least one placement pairing both mates
    Mapped,
    /// assigned in a paired library, but every surviving placement placed only
    /// one mate. Matches the RAD fragment-level rule, which takes the *most
    /// complete* status across hits: a fragment that is a proper pair on one
    /// transcript and an orphan on another counts as a pair.
    MappedOrphan,
}

fn process_fragment(
    recs: &[FragRecord],
    ctx: &FragCtx,
    local: &mut Local,
) -> (FragmentOutcome, Option<bool>, u16) {
    use salmon_model::seqbias::{CONTEXT_LEFT, CONTEXT_LENGTH, CONTEXT_RIGHT};
    // salmon's LOG_EPSILON = log(0.375e-10): the orphan / implausible-length penalty.
    const LOG_EPSILON: f64 = -23.998_158_637_57;

    // Pair records into the placements the aligner reported (proper pairs +
    // orphans), NOT a cross-product of every read1/read2 on a transcript.
    let placements = pair_records(recs);
    // FLD training on this path diverges from every other, deliberately and
    // permanently (#1140, audit F17): the observation is every fragment's
    // longest reported `frag_len` at weight 1, taken *before* strand-
    // compatibility filtering and with no posterior-acceptance sampling, where
    // the one-pass reads path samples by mapping posterior and the RAD paths
    // bake an order-independent distribution. The weighting below likewise uses
    // a flat `LOG_EPSILON` orphan penalty rather than the CMF length-
    // conditioning the other paths apply. This is documented, not fixed: the
    // online alignment path is deprecated (removal planned 2.7.0), and parity
    // work on a floor being torn out would only risk the estimates it still
    // produces. `-a` without `--online` takes the deterministic path, which has
    // none of these divergences.
    let frag_len = recs.iter().map(|r| r.frag_len).max().unwrap_or(0);
    if frag_len > 0 {
        ctx.fld.add_val(frag_len as usize, 0.0);
    }

    // `-l A`: sample this fragment's observed format for the shared detector
    // while it is still collecting (preferring a proper pair's format — orphan
    // observations are single-end, which the detector ignores for a paired
    // library), then filter against the explicit type or, once the detector
    // locks in, the detected one — the same prefix-detect-then-apply the
    // reads-mode one-pass uses.
    if let Some(det) = ctx.detector {
        if det.is_active() {
            let mut sample: Option<LibraryFormat> = None;
            for pl in &placements {
                let (obs, is_fw, status) = frag_format(recs, &pl.idxs);
                let f = obs.unwrap_or_else(|| salmon_core::observed_single_format(is_fw));
                if status == MateStatus::PairedEndPaired {
                    sample = Some(f);
                    break;
                }
                if sample.is_none() {
                    sample = Some(f);
                }
            }
            if let Some(f) = sample {
                det.add_sample(f);
            }
        }
    }
    let expected = ctx
        .expected_format
        .or_else(|| ctx.detector.and_then(|d| d.resolved_format()));
    let use_aux = ctx
        .online
        .is_none_or(|o| o.num_assigned() >= ctx.pre_burnin);

    // Per surviving placement (one reported alignment): conditional log-weight
    // (eq-class) + online log-aux + fragment geometry.
    let mut sp_tid: Vec<u32> = Vec::with_capacity(placements.len());
    let mut sp_eq: Vec<f64> = Vec::with_capacity(placements.len());
    let mut sp_online: Vec<f64> = Vec::with_capacity(placements.len());
    let mut sp_geom: Vec<(usize, usize, bool)> = Vec::with_capacity(placements.len());
    let mut sp_pl: Vec<usize> = Vec::with_capacity(placements.len());
    // Orphan bookkeeping, folded into the loop below so it costs two register
    // ORs per surviving placement and no second pass.
    //
    // Deliberately *not* keyed on `proper`, which is `idxs.len() >= 2 && flen > 0`
    // and so would call a both-mates placement with no reported fragment length an
    // orphan. This mirrors `frag_format`, the canonical status derivation, which
    // looks only at how many mates were placed and at the `0x1` multi-segment flag
    // — the same signal #1057 turned on.
    let mut any_both_mates = false;
    let mut any_orphan = false;
    // Strand judgment for `lib_format_counts.json` (#1130): `None` until a
    // placement is actually compared against an expected format (a
    // `--discardOrphans`-dropped placement never is), then whether any placement
    // was compatible. Judged over all compared placements — including ones
    // `ignore_incompat` then drops — because a fragment whose every reported
    // alignment is wrong-strand is exactly what the count exists to report.
    let mut strand_judged: Option<bool> = None;
    // Observed-format set across this fragment's placements (bit per format
    // id), for the per-format histogram — recorded from every placement that
    // reaches consideration, before the strand filter, since it reports what
    // was observed rather than what was kept.
    let mut fmt_mask: u16 = 0;
    for (pi, pl) in placements.iter().enumerate() {
        let tid = pl.tid;
        let idxs = &pl.idxs;
        // --discardOrphans: in a paired library, drop single-mate placements
        // entirely instead of fragment-length-penalizing them below.
        if ctx.discard_orphans && ctx.paired_lib && idxs.len() < 2 {
            continue;
        }
        let refseq = ctx.ref_bytes.get(tid as usize).unwrap_or(&[]);
        // Conditional log-weight basis = salmon's `errLike` (Σ(fg−bg) over the
        // mate(s) under the error model; uniform 0.0 when it is disabled).
        let basis = if let Some(m) = ctx.model {
            if refseq.is_empty() {
                0.0
            } else {
                let mut ll = 0.0;
                for (rank, &i) in idxs.iter().enumerate() {
                    let r = &recs[i];
                    let (fg, bg) = m.log_likelihood(&r.read_2bit, refseq, r.pos, &r.ops, rank == 0);
                    ll += fg - bg;
                }
                ll
            }
        } else {
            0.0
        };
        let rl = ctx.lengths[tid as usize] as i32;
        let flen = idxs.iter().map(|&i| recs[i].frag_len).max().unwrap_or(0);
        let proper = idxs.len() >= 2 && flen > 0;
        let frag_start = recs[idxs[0]].pos;
        let frag_end = frag_start + (flen.max(1) as usize) - 1;
        let start_pos = if proper && flen <= rl {
            -(((rl - flen + 1) as f64).ln())
        } else {
            -((rl.max(1) as f64).ln())
        };
        let log_frag_prob = if proper {
            if use_aux {
                ctx.fld.pmf(flen as usize)
            } else {
                0.0
            }
        } else if ctx.paired_lib {
            LOG_EPSILON
        } else {
            0.0
        };
        let mut aux = basis + log_frag_prob;
        let (obs, is_fw, status) = frag_format(recs, idxs);
        fmt_mask |= 1
            << obs
                .unwrap_or_else(|| salmon_core::observed_single_format(is_fw))
                .format_id();
        if let Some(exp) = expected {
            let compat = is_compatible(exp, obs, is_fw, status);
            strand_judged = Some(strand_judged.unwrap_or(false) || compat);
            if !compat {
                if ctx.ignore_incompat {
                    continue; // this placement contributes nothing
                }
                aux += ctx.incompat_prior.ln();
            }
        }
        if idxs.len() >= 2 {
            any_both_mates = true;
        } else {
            // A lone record is an orphan only if it belongs to a multi-segment
            // template; a genuine single-end read has no mate to be missing.
            any_orphan |= recs[idxs[0]].is_paired;
        }
        sp_tid.push(tid);
        sp_eq.push(aux);
        sp_online.push(aux + start_pos);
        sp_geom.push((frag_start, frag_end, proper));
        sp_pl.push(pi);
    }
    // a fragment whose every reported alignment was incompatible is a
    // zero-probability fragment: it is not assigned and joins no eq-class —
    // but the strand judgment stands, and is the very thing
    // `num_incompatible_fragments` exists to count.
    if sp_tid.is_empty() {
        return (FragmentOutcome::Unmapped, strand_judged, fmt_mask);
    }
    // The fragment counts as an orphan only when nothing paired both mates, so a
    // fragment that is a proper pair on one transcript and an orphan on another
    // counts as a pair — the same "most complete status wins" rule the RAD
    // fragment-level type uses.
    let outcome = if any_orphan && !any_both_mates {
        FragmentOutcome::MappedOrphan
    } else {
        FragmentOutcome::Mapped
    };

    // Aggregate surviving placements by distinct transcript id (sorted).
    let mut agg: std::collections::BTreeMap<u32, Vec<usize>> = std::collections::BTreeMap::new();
    for (k, &t) in sp_tid.iter().enumerate() {
        agg.entry(t).or_default().push(k);
    }
    let tids: Vec<u32> = agg.keys().cloned().collect();
    let eq_log: Vec<f64> = agg
        .values()
        .map(|ks| logsumexp(&ks.iter().map(|&k| sp_eq[k]).collect::<Vec<_>>()))
        .collect();
    let online_log: Vec<f64> = agg
        .values()
        .map(|ks| logsumexp(&ks.iter().map(|&k| sp_online[k]).collect::<Vec<_>>()))
        .collect();

    // eq-class weights = softmax(eq_log)
    let maxe = eq_log.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut weights: Vec<f64> = eq_log.iter().map(|&l| (l - maxe).exp()).collect();
    let wsum: f64 = weights.iter().sum();
    if wsum > 0.0 {
        for w in &mut weights {
            *w /= wsum;
        }
    }

    // abundance-aware posteriors (online), per distinct transcript
    let post: Vec<f64> = match ctx.online {
        Some(o) => {
            let maps: Vec<(u32, f64)> = tids
                .iter()
                .cloned()
                .zip(online_log.iter().cloned())
                .collect();
            o.assign_fragment(&maps, ctx.log_fm)
        }
        None => weights.clone(),
    };

    // train the error model + collect bias models, weighted by posteriors
    let collecting = ctx.online.is_none_or(|o| o.collecting());
    if collecting {
        for (ti, (tid, ks)) in agg.iter().enumerate() {
            let p_tid = post[ti];
            if p_tid <= 0.0 {
                continue;
            }
            let online_log_tid = online_log[ti];
            let refseq = ctx.ref_bytes.get(*tid as usize).unwrap_or(&[]);
            for &k in ks {
                let p = p_tid * (sp_online[k] - online_log_tid).exp();
                if p <= 0.0 {
                    continue;
                }
                let idxs = &placements[sp_pl[k]].idxs;
                if let Some(m) = local.model.as_mut() {
                    if !refseq.is_empty() {
                        let lw = ctx.log_fm + p.ln();
                        for (rank, &i) in idxs.iter().enumerate() {
                            let r = &recs[i];
                            m.update(&r.read_2bit, refseq, r.pos, &r.ops, rank == 0, lw);
                        }
                    }
                }
                if refseq.is_empty() {
                    continue;
                }
                let (fs, fe, proper) = sp_geom[k];
                let rl = refseq.len();
                // Per-read 5' positions, by each read's actual strand.
                let (mut fwd_five, mut rev_five): (Option<usize>, Option<usize>) = (None, None);
                if idxs.len() == 1 {
                    let r = &recs[idxs[0]];
                    if r.is_reverse {
                        rev_five = Some(r.five_prime());
                    } else {
                        fwd_five = Some(r.five_prime());
                    }
                } else {
                    let fwd = idxs.iter().map(|&i| &recs[i]).find(|r| !r.is_reverse);
                    let rev = idxs.iter().map(|&i| &recs[i]).find(|r| r.is_reverse);
                    if let (Some(fr), Some(rr)) = (fwd, rev) {
                        let (fp, rp) = (fr.five_prime(), rr.five_prime());
                        if fp < rp {
                            fwd_five = Some(fp);
                            rev_five = Some(rp);
                        }
                    }
                }
                if let Some(obs) = local.seq_obs.as_mut() {
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
                if let (Some(gc), true) = (local.gc_obs.as_mut(), proper && fe < rl) {
                    let view = ctx.gc_store.view(*tid as usize);
                    if let Some((ff, cf)) = salmon_model::gc_desc(&view, fs as i32, fe as i32) {
                        gc.inc(ff, cf, p);
                    }
                }
                if let Some(pos) = local.pos_obs.as_mut() {
                    let lc = ctx.length_class.unwrap()[*tid as usize];
                    if let Some(five) = fwd_five {
                        pos.0[lc].add_mass(five as i32, rl as i32, p.ln());
                    }
                    if let Some(five) = rev_five {
                        pos.1[lc].add_mass(five as i32, rl as i32, p.ln());
                    }
                }
            }
        }
    }

    let group = if ctx.range_factorization_bins > 0 {
        let bins = range_factorize_bins(&weights, ctx.range_factorization_bins);
        TranscriptGroup::with_bins(tids, bins)
    } else {
        TranscriptGroup::from_sorted(tids)
    };
    ctx.eq_builder.add_group(group, weights, 1);
    (outcome, strand_judged, fmt_mask)
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
/// Shared by alignment mode and RAD mode: both collect the same per-fragment
/// observations (weighted by abundance-aware posteriors) and then correct
/// effective lengths identically — only the *source* of the observations and
/// abundances differs (online posteriors for BAM; fixed baked/derived abundances
/// for RAD). The caller runs the subsequent re-EM with the corrected lengths.
#[allow(clippy::too_many_arguments)]
fn apply_bias_correction(
    num_refs: usize,
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
            salmon_model::build_expected(num_refs, refseq_of, alphas, eff_lengths, &fld_cdf);
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
    let gc_ratio_model = if let Some(mut obs) = gc_obs {
        let mut exp = salmon_model::build_expected_gc(
            num_refs,
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
            num_refs,
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
            if alphas[tid] < 1e-8 {
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
                seq: seq_tab.as_ref().map(|(f, r)| (f, r)),
                gc: gc_ratio_model.as_ref().map(|g| (g, gc_store.view(tid))),
                pos: pos_vecs
                    .as_ref()
                    .map(|(pf, pr)| (pf.as_slice(), pr.as_slice())),
            };
            *eff_length = salmon_model::corrected_effective_length_full(
                s,
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
/// online path, the deterministic RAD writer, and genome projection — because
/// the caveat is a property of the input, not of which estimator reads it. It
/// used to live only on the online path, so the 2.6.0 default detected and said
/// nothing (COMBINE-lab/salmon#1140, audit C12).
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

pub fn quantify_alignments(opts: &AlignQuantOptions) -> Result<AlignQuantResult> {
    let start_time = opts.external_start_time.clone().unwrap_or_else(asctime_now);
    let run_timer = std::time::Instant::now();
    let mut timer = PhaseTimer::new();
    let header = read_alignment_header(&opts.bam)?;

    // Reject coordinate-sorted input up front (header-only check, no per-record cost):
    // alignment-mode requires all records of a read/pair to be adjacent (grouped by
    // read name), which a coordinate-sorted file violates.
    anyhow::ensure!(
        !coordinate_sorted_unusable(&header),
        "the input BAM/SAM appears to be coordinate-sorted (@HD SO:coordinate) and is not \
         grouped by read name (GO:query). Alignment-mode quantification requires that all \
         alignment records of a read (or read pair) are adjacent in the file. Please collate \
         it by read name first, e.g. `samtools collate` or `samtools sort -n`."
    );

    // References (transcripts) in @SQ order define the transcript ids.
    let names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();
    let lengths: Vec<u32> = header
        .reference_sequences()
        .values()
        .map(|rs| rs.length().get() as u32)
        .collect();
    let num_refs = names.len();
    anyhow::ensure!(num_refs > 0, "BAM header has no reference sequences");

    let eq_builder = EquivalenceClassBuilder::new();
    let mut fld =
        FragmentLengthDistribution::new(1.0, opts.fld_max, opts.fld_mean, opts.fld_sd, 4, 0.5, 1);

    // The error model and bias models need the transcriptome (salmon requires `-t`).
    let use_error_model = opts.transcripts.is_some() && !opts.no_error_model;
    let bias_on = opts.seq_bias || opts.gc_bias || opts.pos_bias;
    anyhow::ensure!(
        !bias_on || opts.transcripts.is_some(),
        "--seqBias/--gcBias/--posBias in alignment mode require -t/--targets (the transcriptome FASTA)"
    );
    let ref_bytes: salmon_core::RefSeqs = if use_error_model || bias_on {
        salmon_core::RefSeqs::from_sequences(load_ref_bytes(
            opts.transcripts.as_ref().unwrap(),
            &names,
        )?)
    } else {
        salmon_core::RefSeqs::default()
    };

    // Online (dual-phase) abundances: develop running estimates so the error
    // model and bias models are trained/collected with abundance-aware posteriors
    // in a single streaming pass (salmon's online phase), rather than two passes.
    let ref_lens_u64: Vec<u64> = lengths.iter().map(|&l| l as u64).collect();
    let online = (use_error_model || bias_on).then(|| {
        salmon_infer::OnlineInference::new(
            &ref_lens_u64,
            0.05,
            opts.forgetting_factor,
            opts.num_aux_model_samples,
        )
    });

    // Per-transcript inputs the bias collection needs (the observed bias models
    // themselves are accumulated per-worker inside the pass).
    // GC cumulative-count backing: one rank bitvector over the concatenated
    // references (salmon's `--reduceGCMemory`, now the default — faster and ~2x
    // leaner than dense per-transcript prefixes, identical results). `gc_store`
    // presents per-transcript `GcView`s.
    let (gc_rank, gc_offsets): (Option<salmon_model::GcRank>, Vec<u64>) = if opts.gc_bias {
        // `RefSeqs` already stores the references concatenated with their
        // endpoints, which is what the rank bitvector needs; rebuilding both here
        // copied every base a second time.
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

    // ---- online pass (reader/worker pipeline) ------------------------------
    // A dedicated reader thread groups raw records by name; the worker pool
    // converts (record_to_frag -- the bulk of per-record cost) and weights them
    // in parallel. Per-thread error-model/bias deltas merge into the globals
    // between minibatches (salmon's minibatch staleness).
    const MINIBATCH: usize = 1000;
    // `-l A`: infer the library type from the alignments themselves, with the
    // same prefix-sampling detector the reads-mode one-pass uses — sample
    // observed formats until the budget is spent, lock in, and filter the rest
    // of the run against the detected type. (C++ salmon's alignment mode never
    // detected strandedness — it peeked pairedness and fell back to IU/U — but
    // the rewrite's deterministic `-a` path detects via its RAD pass, and the
    // online path should agree with it.) Read type comes from the first
    // record's flags, C++'s `peekBAMIsPaired`.
    let auto_detect = LibraryFormat::is_auto(&opts.lib_type);
    let peeked_paired = if auto_detect {
        peek_alignment_is_paired(&opts.bam)?
    } else {
        None
    };
    let detector = auto_detect.then(|| {
        warn_auto_detect_from_alignments();
        salmon_model::LibraryTypeDetector::new(if peeked_paired.unwrap_or(true) {
            salmon_core::ReadType::PairedEnd
        } else {
            salmon_core::ReadType::SingleEnd
        })
    });
    // A paired library expects two mates; a single mate to a transcript is then
    // an "unexpected orphan" and is fragment-length-penalized. (Single-end libs
    // aren't.) Under `-l A`, pairedness comes from the peek above.
    let paired_lib = match peeked_paired {
        Some(p) => p,
        None => !matches!(opts.lib_type.as_str(), "U" | "SF" | "SR" | "S"),
    };
    // Orientation-compatibility filtering (salmon): drop alignments whose
    // orientation is incompatible with the expected library type. `None` under
    // auto (`A`), where the detector above supplies the format mid-run.
    let expected_format = if auto_detect {
        None
    } else {
        LibraryFormat::parse(&opts.lib_type).ok()
    };
    let ignore_incompat = opts.incompat_prior <= 0.0;
    let nthreads = rayon::current_num_threads().max(1);

    // The bias accumulators + counters are developed inside the pass; the model
    // is trained there too but not needed afterward. Scope `cfg`/`acc` so their
    // borrows (fld, eq_builder, online, ...) are released before the post-pass
    // effective-length / EM work below.
    let (
        seq_obs,
        gc_obs,
        pos_obs,
        num_processed,
        num_mapped,
        num_orphan,
        num_compat,
        num_incompat,
        lib_format_counts,
    ) = {
        let cfg = PassCfg {
            online: online.as_ref(),
            fld: &fld,
            eq_builder: &eq_builder,
            ref_bytes: &ref_bytes,
            lengths: &lengths,
            gc_store,
            length_class: length_class.as_deref(),
            expected_format,
            detector: detector.as_ref(),
            ignore_incompat,
            incompat_prior: opts.incompat_prior,
            paired_lib,
            discard_orphans: opts.discard_orphans,
            pre_burnin: opts.num_pre_aux_model_samples,
            range_factorization_bins: opts.range_factorization_bins,
            use_error_model,
            error_bins: opts.num_error_bins,
            seq_bias: opts.seq_bias,
            gc_bias: opts.gc_bias,
            cond_gc_bins: opts.cond_gc_bins,
            gc_bins: opts.gc_bins,
            pos_bias: opts.pos_bias,
            need_seq: use_error_model || bias_on,
            minibatch: MINIBATCH,
            nthreads,
            progress: opts.progress.as_deref(),
        };
        let mut acc = PassAccum {
            seq_obs: None,
            gc_obs: None,
            pos_obs: None,
            num_processed: 0,
            num_mapped: 0,
            num_orphan: 0,
            num_frags_compat: 0,
            num_frags_incompat: 0,
            lib_format_counts: [0; salmon_core::NUM_LIB_FORMATS],
        };
        stream_online_pass(&opts.bam, &cfg, &mut acc)?;
        (
            acc.seq_obs,
            acc.gc_obs,
            acc.pos_obs,
            acc.num_processed,
            acc.num_mapped,
            acc.num_orphan,
            acc.num_frags_compat,
            acc.num_frags_incompat,
            acc.lib_format_counts,
        )
    };

    // `-l A`: the end-of-run resolved type — the locked-in format when the
    // sample budget was reached mid-run (in which case it also filtered the
    // rest of the run), else the best guess from whatever samples a smaller
    // run collected. Reported, and used as `lib_format_counts.json`'s
    // expected format.
    let detected_library_type = detector.as_ref().map(|d| {
        let f = d.final_format();
        tracing::info!(
            "library format under `-l A`: {} (detected from alignments)",
            f.canonical()
        );
        f.canonical().to_string()
    });

    // ---- base effective lengths --------------------------------------------
    fld.cache();
    let cond_means = fld.conditional_means();
    let mut eff_lengths = vec![0f64; num_refs];
    for (tid, len) in lengths.iter().enumerate() {
        eff_lengths[tid] = salmon_model::smoothed_effective_length(&cond_means, *len as usize);
    }

    let mut collapsed = eq_builder.finish();
    collapsed.update_eff_lengths(&eff_lengths);
    let num_eq_classes = collapsed.len();
    // BAM read + online equivalence-class build.
    timer.mark("mapping");

    // Count-blended EM initialization (salmon's `CollapsedEMOptimizer::optimize`):
    // seed abundances with a linear combination of the online-phase abundance
    // estimates (`projectedCounts`) and the uniform distribution, weighted by the
    // fraction of the burn-in target observed. A warm start near the solution
    // cuts the number of EM iterations to convergence.
    let init_alphas: Option<Vec<f64>> = online.as_ref().map(|o| {
        let masses: Vec<f64> = (0..num_refs).map(|t| o.mass_log(t).exp()).collect();
        let mass_sum: f64 = masses.iter().sum();
        let total_reads = num_mapped as f64;
        // online estimates scaled to a proper count distribution (sum = reads)
        let projected: Vec<f64> = if mass_sum > 0.0 {
            masses.iter().map(|&m| m / mass_sum * total_reads).collect()
        } else {
            vec![0.0; num_refs]
        };
        let uniform_prior = if num_refs > 0 {
            total_reads / num_refs as f64
        } else {
            0.0
        };
        // salmon's numRequiredFragments default (the online-phase target)
        const NUM_REQUIRED_FRAGMENTS: f64 = 50_000_000.0;
        let frac_observed = (total_reads / NUM_REQUIRED_FRAGMENTS).min(0.999);
        projected
            .iter()
            .map(|&p| p * frac_observed + uniform_prior * (1.0 - frac_observed))
            .collect()
    });
    // The packed CSR layout is built once and reused by the EM, the post-bias
    // re-run and bootstrap / Gibbs (see the reads-mode driver for the rationale).
    // `weights` is Gibbs-only, so it is populated only when Gibbs will run.
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
    // `--initUniform` forces the plain uniform EM start; otherwise warm-start
    // from the online-estimate-blended init.
    let mut em = if opts.skip_quant {
        // `--skipQuant`: emit eq-classes + library type + metadata, skip the
        // optimizer and quant.sf. Abundances left at zero.
        salmon_infer::EmResult {
            alphas: vec![0.0; num_refs],
            iters: 0,
            converged: true,
            dropped_mass: 0.0,
        }
    } else if opts.init_uniform {
        salmon_infer::optimize_packed_with_init(&packed, &opts.em, true, None, Some(&eff_lengths))
    } else {
        salmon_infer::optimize_packed_with_init(
            &packed,
            &opts.em,
            true,
            init_alphas.as_deref(),
            Some(&eff_lengths),
        )
    };

    // ---- bias-corrected effective lengths (shared with RAD mode) ------------
    let mut bias_dump = salmon_model::dumps::BiasDump::default();
    if bias_on && !opts.skip_quant {
        bias_dump = apply_bias_correction(
            num_refs,
            &ref_bytes,
            &gc_store,
            &lengths,
            length_class.as_deref(),
            length_quantiles.as_deref(),
            &fld,
            &em.alphas,
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
        // Only the combined weights changed; patch them in place.
        packed.refresh_combined(&collapsed);
        em = salmon_infer::optimize_packed_with_init(
            &packed,
            &opts.em,
            true,
            None,
            Some(&eff_lengths),
        );
    }
    let inference_truncated_mass = em.dropped_mass;
    let (em_iters, em_converged) = (em.iters, em.converged);
    let counts = em.alphas;
    // Bias correction + EM/VBEM point estimate.
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

    // Run diagnostics from end-of-run aggregates. The aggregate observed-format
    // estimate is the detector's answer under `-l A`, else inferred from the
    // end-of-run observed-format histogram with the same decision rule — which
    // is what lets the library-type mismatch warning fire here, as it does in
    // the `--rad` path, when an explicit `-l` disagrees with the data.
    let observed_estimate = detected_library_type.clone().or_else(|| {
        (num_mapped > 0).then(|| {
            salmon_model::infer_format_from_counts(
                &lib_format_counts,
                if paired_lib {
                    salmon_core::ReadType::PairedEnd
                } else {
                    salmon_core::ReadType::SingleEnd
                },
            )
            .canonical()
            .to_string()
        })
    });
    let requested_lib = LibraryFormat::parse(&opts.lib_type)
        .map(|f| f.canonical().to_string())
        .unwrap_or_else(|_| opts.lib_type.clone());
    let diagnostics = salmon_core::input_diagnostics(
        num_processed,
        num_mapped,
        LibraryFormat::is_auto(&opts.lib_type),
        &requested_lib,
        observed_estimate.as_deref(),
    );
    for d in &diagnostics {
        if d.severity == "error" {
            tracing::error!("{}", d.message);
        } else {
            tracing::warn!("{}", d.message);
        }
    }

    // Posterior uncertainty (bootstrap / Gibbs) + ambiguity. The packed CSR
    // layout is identical to reads mode — alignment mode simply feeds it from the
    // BAM-derived eq-classes — so the same `salmon_infer` samplers apply. Bootstrap
    // takes precedence over Gibbs, matching reads mode.
    let ambig = salmon_infer::ambiguity_counts(&packed);
    let bootstraps: Vec<Vec<f64>> = if opts.skip_quant {
        Vec::new()
    } else if opts.num_bootstraps > 0 {
        salmon_infer::bootstrap(
            &packed,
            &opts.em,
            Some(&eff_lengths),
            opts.num_bootstraps,
            0x5A13_0000,
        )
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
            // Honour --perNucleotidePrior here as the EM does; hardcoding a
            // per-transcript prior made the flag a no-op for Gibbs on every
            // path (#1140, audit D13).
            per_transcript_prior: !opts.em.per_nucleotide_prior,
        };
        salmon_infer::gibbs_sample(&packed, &eff_lengths, &counts, &gopts, 0x6217_0000)
    } else {
        Vec::new()
    };
    // Posterior sampling (bootstrap / Gibbs); empty when neither is requested.
    timer.mark("posterior");

    let result = AlignQuantResult {
        // Tallied while streaming the BAM, so `-a` reports the same
        // concordant/orphan split the reads path does.
        num_orphan: Some(num_orphan),
        num_compatible_fragments: num_compat,
        num_incompatible_fragments: num_incompat,
        lib_format_counts,
        provenance: crate::rad::RadProvenance::default(),
        // Straight from the BAM being quantified, no RAD in between.
        source_programs: source_program_lines(&header),
        source: FragmentSource::Bam,
        names,
        lengths,
        eff_lengths,
        tpm,
        counts,
        num_processed,
        num_mapped,
        inference_truncated_mass,
        num_eq_classes,
        frag_len_mean: fld.mean(),
        frag_len_sd: fld.sd(),
        frag_len_source: salmon_model::FragLengthSource::Alignments,
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
    write_outputs(opts, &result)?;
    timer.mark("output");
    Ok(result)
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
    let (first_decoy, num_decoys) = res
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
    let rows: Vec<usize> =
        salmon_core::quant_row_indices(res.names.len(), first_decoy, num_decoys).collect();

    // quant.sf (EffectiveLength + NumReads at --sigDigits decimals; TPM at 6).
    // Skipped under --skipQuant (no abundances), like salmon.
    if !opts.skip_quant {
        let sd = opts.sig_digits as usize;
        let mut w = std::io::BufWriter::new(std::fs::File::create(dir.join("quant.sf"))?);
        writeln!(w, "Name\tLength\tEffectiveLength\tTPM\tNumReads")?;
        for &i in &rows {
            writeln!(
                w,
                "{}\t{}\t{:.*}\t{:.6}\t{:.*}",
                res.names[i], res.lengths[i], sd, res.eff_lengths[i], res.tpm[i], sd, res.counts[i]
            )?;
        }
        w.flush()?;
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
        detected_library_type: Option<String>,
        /// which inference path produced these results: `deterministic` (RAD
        /// input — including `--deterministic`'s phase 2) or `online` (the
        /// streaming BAM pass). Makes the online-to-deterministic transition
        /// auditable from existing output (#1140).
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
    // every quantify_rad driver and the online alignment path, so the 2.6.0
    // default surfaces the condition it computes instead of only recording it
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
        gc_bias_correct: opts.gc_bias,
        pos_bias_correct: opts.pos_bias,
        // Not the running bias models' bin counts: C++ writes
        // `readBias(FORWARD).counts.size()` (GZipWriter.cpp), the legacy fixed
        // 4096-bin k-mer counter behind `observed_bias.gz`, which this port
        // replaced with SBModel/GC/pos models and stubs out in the aux dumps.
        // There is no such bin count to report here, and the models that did
        // run publish their own sizes in their dumps, so 0 stands (matching the
        // reads writer) rather than a number of a different model's bins.
        num_bias_bins: 0,
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
        // concern (#1140). Appended here, in the writer both `quantify_rad`
        // and `quantify_alignments` share, so the field's promise holds on
        // every path rather than only the RAD one.
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

        let f = std::fs::File::create(bdir.join("bootstraps.gz"))?;
        let mut enc = GzEncoder::new(f, Compression::new(6));
        for sample in &res.bootstraps {
            // Restricted to the same rows as quant.sf so the two stay aligned
            // positionally, which is what tximport assumes.
            for &i in &rows {
                enc.write_all(&sample[i].to_le_bytes())?;
            }
        }
        enc.finish()?;
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
    // `--dumpBiasModels`: the same dump the one-pass reads path writes, from
    // this pass's own models — the shared writer makes it work on every
    // RAD/alignment driver (#1140: it silently vanished when deterministic
    // became the default, because only salmon-quant had a write site).
    if opts.dump_bias_models {
        salmon_model::dumps::dump_bias_models_to_file(&dir.join("bias_models.txt"), &res.bias_dump)
            .context("writing bias_models.txt")?;
    }
    std::fs::create_dir_all(dir.join("logs")).context("creating logs")?;
    // In alignment mode the input records are already aligned, so `num_processed`
    // is the count of aligned fragments and `num_mapped` is those with a strand-
    // compatible placement that were quantified; label them accordingly so a
    // <100% rate on a stranded library is not read as lost alignments.
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
            "aligned fragments",
            "strand-compatible (quantified)",
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
            let err = load_ref_bytes(&p, &names).unwrap_err();
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

        let err = load_ref_bytes(&p, &["t0".to_string()]).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("truncated.fa.gz"),
            "a mid-stream failure should still name the file, got: {msg}"
        );
    }

    /// Partial coverage is legitimate — a transcriptome need not cover every
    /// reference in the alignment header — so it must not be an error.
    #[test]
    fn partial_reference_coverage_is_accepted() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ref.fa");
        std::fs::write(&p, ">t0\nACGT\n").unwrap();
        let got = load_ref_bytes(&p, &["t0".to_string(), "absent".to_string()]).unwrap();
        assert_eq!(got[0], b"ACGT");
        assert!(
            got[1].is_empty(),
            "an absent reference yields an empty sequence"
        );
    }
}

#[cfg(test)]
mod logsumexp_tests {
    use super::{logsumexp, logsumexp_iter};

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
