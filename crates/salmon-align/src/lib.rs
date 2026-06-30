//! `salmon-align`: alignment-based (BAM) quantification.
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
mod rad;

pub use bam_rad::{write_alignment_rad, AlignRadSummary};
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
use serde::Serialize;

use error_model::{AlignmentModel, AlnOp, SharedAlignmentModel};
use salmon_core::{
    is_compatible, LibraryFormat, MateStatus, ReadOrientation, ReadStrandedness, ReadType,
};
use salmon_eqclass::{range_factorize_bins, EquivalenceClassBuilder, TranscriptGroup};
use salmon_infer::{optimize, optimize_with_init, EmOptions};
use salmon_model::FragmentLengthDistribution;

const SALMON_VERSION: &str = env!("CARGO_PKG_VERSION");

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
    pub ref_seqs: Option<Vec<Vec<u8>>>,
    /// disable the alignment error model (salmon's `--noErrorModel`)
    pub no_error_model: bool,
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
    /// online-phase forgetting factor (`--forgettingFactor`, salmon default 0.65)
    pub forgetting_factor: f64,
    /// initialize the EM uniformly instead of with the online-estimate-blended
    /// warm start (`--initUniform`)
    pub init_uniform: bool,
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
            seq_bias: false,
            gc_bias: false,
            pos_bias: false,
            bias_seed_em_iters: 11,
            incompat_prior: 0.0,
            fld_mean: 250.0,
            fld_sd: 25.0,
            fld_max: 1000,
            forgetting_factor: 0.65,
            init_uniform: false,
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
#[derive(Debug, Clone)]
pub struct AlignQuantResult {
    pub names: Vec<String>,
    pub lengths: Vec<u32>,
    pub eff_lengths: Vec<f64>,
    pub tpm: Vec<f64>,
    pub counts: Vec<f64>,
    pub num_processed: u64,
    pub num_mapped: u64,
    /// fragment mass dropped in the final min-alpha redistribution (every member
    /// of an equivalence class truncated); normally 0.
    pub inference_truncated_mass: f64,
    pub num_eq_classes: usize,
    pub frag_len_mean: f64,
    pub frag_len_sd: f64,
    pub length_classes: Vec<u32>,
    pub frag_len_dist: Vec<f64>,
    pub start_time: String,
    pub bias_dump: salmon_model::dumps::BiasDump,
    /// per-transcript (unique, ambiguous) fragment counts for `ambig_info.tsv`
    pub ambig: (Vec<u32>, Vec<u32>),
    /// posterior samples (bootstrap or Gibbs), one abundance vector each; empty
    /// when neither was requested. Length is `num_refs`, matching `quant.sf` rows.
    pub bootstraps: Vec<Vec<f64>>,
}

/// Current local time as an asctime-style string, matching salmon's timestamps.
fn asctime_now() -> String {
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
        let status = if r.is_read1 {
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

/// Load a (optionally gzip'd) transcriptome FASTA and return the (ASCII) base
/// sequences aligned to the BAM's reference order (`names`); a name absent from
/// the FASTA yields an empty sequence (its model contributions are skipped). The
/// same bytes feed both the error model (2-bit on the fly) and the bias models.
fn load_ref_bytes(path: &Path, names: &[String]) -> Result<Vec<Vec<u8>>> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut magic = [0u8; 2];
    let is_gz = {
        use std::io::Read;
        let mut f = std::fs::File::open(path)?;
        f.read_exact(&mut magic).is_ok() && magic == [0x1f, 0x8b]
    };
    let reader: Box<dyn BufRead> = if is_gz {
        Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(
            file,
        )))
    } else {
        Box::new(std::io::BufReader::new(file))
    };

    let mut by_name: HashMap<String, Vec<u8>> = HashMap::new();
    let mut cur_name: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line?;
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
    Ok(names
        .iter()
        .map(|n| by_name.remove(n).unwrap_or_default())
        .collect())
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
    frag_len: i32,
    /// reverse-strand alignment (BAM `0x10` flag)
    is_reverse: bool,
    /// first mate of the pair (BAM `0x40` flag)
    is_read1: bool,
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

/// Does this path name a SAM file (plain or gzipped) rather than a BAM?
fn is_sam_path(p: &Path) -> bool {
    let s = p.to_string_lossy().to_ascii_lowercase();
    s.ends_with(".sam") || s.ends_with(".sam.gz")
}

/// Open a SAM file (transparently gunzipping `.sam.gz`) as a noodles reader over
/// a boxed `BufRead`, so plain and gzipped inputs share one concrete type.
fn open_sam_reader(path: &Path) -> Result<sam::io::Reader<Box<dyn BufRead + Send>>> {
    let file = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let lower = path.to_string_lossy().to_ascii_lowercase();
    let inner: Box<dyn BufRead + Send> = if lower.ends_with(".gz") {
        Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(
            file,
        )))
    } else {
        Box::new(std::io::BufReader::new(file))
    };
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
    let score = record
        .data()
        .get(&Tag::ALIGNMENT_SCORE)
        .and_then(|r| r.ok())
        .and_then(|v| value_as_i32(&v))
        .unwrap_or(0);
    let frag_len = record.template_length().map(|t| t.abs()).unwrap_or(0);
    let is_reverse = flags.is_reverse_complemented();
    let is_read1 = flags.is_first_segment();
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
            frag_len,
            is_reverse,
            is_read1,
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
    ref_bytes: &'a [Vec<u8>],
    lengths: &'a [u32],
    gc_store: salmon_model::GcStore<'a>,
    length_class: Option<&'a [usize]>,
    expected_format: Option<LibraryFormat>,
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
            workers.push(scope.spawn(move || -> (Local, u64, u64) {
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
                    for raw_group in &raw_batch {
                        frags.clear();
                        for r in raw_group {
                            if let Some((_, f)) = record_to_frag(r, header, need_seq) {
                                frags.push(f);
                            }
                        }
                        if process_fragment(&frags, &ctx, &mut local) {
                            batch_mapped += 1;
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
                (local, count, mapped)
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
        for w in workers {
            let (local, count, mapped) = w
                .join()
                .map_err(|_| anyhow::anyhow!("alignment worker thread panicked"))?;
            merged = merged.merge(local);
            total += count;
            total_mapped += mapped;
        }
        acc.seq_obs = merged.seq_obs;
        acc.gc_obs = merged.gc_obs;
        acc.pos_obs = merged.pos_obs;
        // num_processed = every aligned fragment in the BAM; num_mapped = those
        // with a surviving strand-compatible placement (assigned to an eq-class).
        // They differ only for stranded libraries; matches reads-mode (#1025).
        acc.num_processed = total;
        acc.num_mapped = total_mapped;
        Ok(())
    })
}

/// Dispatch the online pass over a SAM/BAM file (chosen by extension). BAM is
/// decoded on a small BGZF worker pool (the framing/parse stays serial on the
/// reader thread, but the decompression overlaps it).
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
        let decoder = noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, file);
        let mut reader = bam::io::Reader::from(decoder);
        let header = reader.read_header().context("reading BAM header")?;
        run_online_pass(reader.records(), &header, cfg, acc)
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
    ref_bytes: &'a [Vec<u8>],
    lengths: &'a [u32],
    gc_store: salmon_model::GcStore<'a>,
    length_class: Option<&'a [usize]>,
    expected_format: Option<LibraryFormat>,
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
fn process_fragment(recs: &[FragRecord], ctx: &FragCtx, local: &mut Local) -> bool {
    use salmon_model::seqbias::{CONTEXT_LEFT, CONTEXT_LENGTH, CONTEXT_RIGHT};
    // salmon's LOG_EPSILON = log(0.375e-10): the orphan / implausible-length penalty.
    const LOG_EPSILON: f64 = -23.998_158_637_57;

    // Pair records into the placements the aligner reported (proper pairs +
    // orphans), NOT a cross-product of every read1/read2 on a transcript.
    let placements = pair_records(recs);
    let frag_len = recs.iter().map(|r| r.frag_len).max().unwrap_or(0);
    if frag_len > 0 {
        ctx.fld.add_val(frag_len as usize, 0.0);
    }
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
    for (pi, pl) in placements.iter().enumerate() {
        let tid = pl.tid;
        let idxs = &pl.idxs;
        // --discardOrphans: in a paired library, drop single-mate placements
        // entirely instead of fragment-length-penalizing them below.
        if ctx.discard_orphans && ctx.paired_lib && idxs.len() < 2 {
            continue;
        }
        let refseq = ctx
            .ref_bytes
            .get(tid as usize)
            .map(|v| v.as_slice())
            .unwrap_or(&[]);
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
        if let Some(exp) = ctx.expected_format {
            let (obs, is_fw, status) = frag_format(recs, idxs);
            if !is_compatible(exp, obs, is_fw, status) {
                if ctx.ignore_incompat {
                    continue; // this placement contributes nothing
                }
                aux += ctx.incompat_prior.ln();
            }
        }
        sp_tid.push(tid);
        sp_eq.push(aux);
        sp_online.push(aux + start_pos);
        sp_geom.push((frag_start, frag_end, proper));
        sp_pl.push(pi);
    }
    // a fragment whose every reported alignment was incompatible is a
    // zero-probability fragment: it is not assigned and joins no eq-class.
    if sp_tid.is_empty() {
        return false;
    }

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
            let refseq = ctx
                .ref_bytes
                .get(*tid as usize)
                .map(|v| v.as_slice())
                .unwrap_or(&[]);
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
    true
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
    ref_bytes: &[Vec<u8>],
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
    let refseq_of = |t: usize| ref_bytes[t].as_slice();

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

    for tid in 0..num_refs {
        if alphas[tid] < 1e-8 {
            continue;
        }
        let s = ref_bytes[tid].as_slice();
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
        eff_lengths[tid] = salmon_model::corrected_effective_length_full(
            s,
            &fld_cdf,
            fld_low,
            fld_high,
            &bias,
            eff_lengths[tid],
            bias_speed_samp,
            no_bias_length_threshold,
        );
    }
    bias_dump
}

/// Run alignment-based quantification end-to-end.
pub fn quantify_alignments(opts: &AlignQuantOptions) -> Result<AlignQuantResult> {
    let start_time = asctime_now();
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
    let ref_bytes: Vec<Vec<u8>> = if use_error_model || bias_on {
        load_ref_bytes(opts.transcripts.as_ref().unwrap(), &names)?
    } else {
        Vec::new()
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
        let mut concat: Vec<u8> = Vec::new();
        let mut offs: Vec<u64> = Vec::with_capacity(ref_bytes.len() + 1);
        offs.push(0);
        for s in &ref_bytes {
            concat.extend_from_slice(s);
            offs.push(concat.len() as u64);
        }
        (Some(salmon_model::GcRank::new(&concat)), offs)
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
    // A paired library expects two mates; a single mate to a transcript is then
    // an "unexpected orphan" and is fragment-length-penalized. (Single-end libs
    // aren't.)
    let paired_lib = !matches!(opts.lib_type.as_str(), "U" | "SF" | "SR" | "S");
    // Orientation-compatibility filtering (salmon): drop alignments whose
    // orientation is incompatible with the expected library type. Skipped under
    // auto (`A`) library type.
    let expected_format = match opts.lib_type.as_str() {
        "A" => None,
        s => LibraryFormat::parse(s).ok(),
    };
    let ignore_incompat = opts.incompat_prior <= 0.0;
    let nthreads = rayon::current_num_threads().max(1);

    // The bias accumulators + counters are developed inside the pass; the model
    // is trained there too but not needed afterward. Scope `cfg`/`acc` so their
    // borrows (fld, eq_builder, online, ...) are released before the post-pass
    // effective-length / EM work below.
    let (seq_obs, gc_obs, pos_obs, num_processed, num_mapped) = {
        let cfg = PassCfg {
            online: online.as_ref(),
            fld: &fld,
            eq_builder: &eq_builder,
            ref_bytes: &ref_bytes,
            lengths: &lengths,
            gc_store,
            length_class: length_class.as_deref(),
            expected_format,
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
        };
        stream_online_pass(&opts.bam, &cfg, &mut acc)?;
        (
            acc.seq_obs,
            acc.gc_obs,
            acc.pos_obs,
            acc.num_processed,
            acc.num_mapped,
        )
    };

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
        optimize(&collapsed, num_refs, &opts.em, Some(&eff_lengths))
    } else {
        optimize_with_init(
            &collapsed,
            num_refs,
            &opts.em,
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
        em = optimize(&collapsed, num_refs, &opts.em, Some(&eff_lengths));
    }
    let inference_truncated_mass = em.dropped_mass;
    let counts = em.alphas;

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

    // Posterior uncertainty (bootstrap / Gibbs) + ambiguity. The packed CSR
    // layout is identical to reads mode — alignment mode simply feeds it from the
    // BAM-derived eq-classes — so the same `salmon_infer` samplers apply. Bootstrap
    // takes precedence over Gibbs, matching reads mode.
    let packed = salmon_infer::PackedEqClasses::from_collapsed(&collapsed, num_refs);
    let ambig = salmon_infer::ambiguity_counts(&packed);
    let bootstraps: Vec<Vec<f64>> = if opts.skip_quant {
        Vec::new()
    } else if opts.num_bootstraps > 0 {
        salmon_infer::bootstrap(&packed, &opts.em, opts.num_bootstraps, 0x5A13_0000)
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
            per_transcript_prior: true,
        };
        salmon_infer::gibbs_sample(&packed, &eff_lengths, &counts, &gopts, 0x6217_0000)
    } else {
        Vec::new()
    };

    let result = AlignQuantResult {
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
        length_classes,
        frag_len_dist: fld.log_pmf().iter().map(|lp| lp.exp()).collect(),
        start_time,
        bias_dump,
        ambig,
        bootstraps,
    };
    write_outputs(opts, &result)?;
    Ok(result)
}

fn write_outputs(opts: &AlignQuantOptions, res: &AlignQuantResult) -> Result<()> {
    let dir = &opts.output_dir;
    std::fs::create_dir_all(dir.join("aux_info")).context("creating output dirs")?;

    // quant.sf (EffectiveLength + NumReads at --sigDigits decimals; TPM at 6).
    // Skipped under --skipQuant (no abundances), like salmon.
    if !opts.skip_quant {
        let sd = opts.sig_digits as usize;
        let mut w = std::io::BufWriter::new(std::fs::File::create(dir.join("quant.sf"))?);
        writeln!(w, "Name\tLength\tEffectiveLength\tTPM\tNumReads")?;
        for i in 0..res.names.len() {
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
        seq_bias_correct: bool,
        gc_bias_correct: bool,
        num_bias_bins: usize,
        mapping_type: String,
        // index hashes are empty in alignment mode (no salmon index)
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
    let meta = MetaInfo {
        salmon_version: SALMON_VERSION.to_string(),
        samp_type: "none".to_string(),
        opt_type: if opts.em.use_vbem { "vb" } else { "em" }.to_string(),
        quant_errors: vec![],
        num_libraries: 1,
        library_types: vec![opts.lib_type.clone()],
        frag_dist_length: res.frag_len_dist.len(),
        frag_length_mean: res.frag_len_mean,
        frag_length_sd: res.frag_len_sd,
        seq_bias_correct: opts.seq_bias,
        gc_bias_correct: opts.gc_bias,
        num_bias_bins: 0,
        mapping_type: "alignment".to_string(),
        index_seq_hash: String::new(),
        index_name_hash: String::new(),
        index_seq_hash512: String::new(),
        index_name_hash512: String::new(),
        index_decoy_seq_hash: String::new(),
        index_decoy_name_hash: String::new(),
        num_valid_targets: res.names.len(),
        num_decoy_targets: 0,
        num_eq_classes: res.num_eq_classes,
        serialized_eq_classes: false,
        eq_class_properties,
        length_classes: res.length_classes.clone(),
        num_processed: res.num_processed,
        num_mapped: res.num_mapped,
        num_dovetail_fragments: 0,
        num_fragments_filtered_vm: 0,
        num_alignments_below_threshold_for_mapped_fragments_vm: 0,
        percent_mapped: pct,
        num_decoy_fragments: 0,
        inference_truncated_mass: res.inference_truncated_mass,
        num_bootstraps: res.bootstraps.len() as u32,
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
        for (i, name) in res.names.iter().enumerate() {
            if i > 0 {
                enc.write_all(b"\t")?;
            }
            enc.write_all(name.as_bytes())?;
        }
        enc.write_all(b"\n")?;
        enc.finish()?;

        let f = std::fs::File::create(bdir.join("bootstraps.gz"))?;
        let mut enc = GzEncoder::new(f, Compression::new(6));
        for sample in &res.bootstraps {
            for v in sample {
                enc.write_all(&v.to_le_bytes())?;
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
    std::fs::create_dir_all(dir.join("logs")).context("creating logs")?;
    // In alignment mode the input records are already aligned, so `num_processed`
    // is the count of aligned fragments and `num_mapped` is those with a strand-
    // compatible placement that were quantified; label them accordingly so a
    // <100% rate on a stranded library is not read as lost alignments.
    let log = format!(
        "salmon (rust port, alignment mode) v{SALMON_VERSION}\nstart: {}\nend:   {}\n\
         library type: {}\naligned fragments: {}\nstrand-compatible (quantified): {}\n\
         compatible rate: {pct:.4}%\n\
         number of equivalence classes: {}\nfragment length mean (sd): {:.2} ({:.2})\n",
        res.start_time,
        asctime_now(),
        opts.lib_type,
        res.num_processed,
        res.num_mapped,
        res.num_eq_classes,
        res.frag_len_mean,
        res.frag_len_sd,
    );
    std::fs::write(dir.join("logs").join("salmon_quant.log"), log)
        .context("writing salmon_quant.log")?;

    // aux_info/ambig_info.tsv
    {
        let (uniq, ambig) = &res.ambig;
        let mut w = std::io::BufWriter::new(std::fs::File::create(
            dir.join("aux_info").join("ambig_info.tsv"),
        )?);
        writeln!(w, "UniqueCount\tAmbigCount")?;
        for i in 0..uniq.len() {
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
    std::fs::write(
        dir.join("cmd_info.json"),
        serde_json::to_string_pretty(&cmd)?,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn alignment_score_value_decoding() {
        assert_eq!(value_as_i32(&Value::Int32(-12)), Some(-12));
        assert_eq!(value_as_i32(&Value::Int8(0)), Some(0));
        assert_eq!(value_as_i32(&Value::UInt16(300)), Some(300));
        // a non-integer tag value (e.g. a character) has no integer reading
        assert_eq!(value_as_i32(&Value::Character(b'A')), None);
    }
}
