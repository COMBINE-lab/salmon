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

mod error_model;

use std::collections::HashMap;
use std::io::{BufRead, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use serde::Serialize;

use error_model::{AlignmentModel, AlnOp};
use salmon_core::{
    is_compatible, LibraryFormat, MateStatus, ReadOrientation, ReadStrandedness, ReadType,
};
use salmon_eqclass::{range_factorize_bins, EquivalenceClassBuilder, TranscriptGroup};
use salmon_infer::{optimize, EmOptions};
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
    /// disable the alignment error model (salmon's `--noErrorModel`)
    pub no_error_model: bool,
    /// enable sequence-specific bias correction (`--seqBias`)
    pub seq_bias: bool,
    /// enable fragment-GC bias correction (`--gcBias`)
    pub gc_bias: bool,
    /// enable positional bias correction (`--posBias`)
    pub pos_bias: bool,
    /// weight multiplier for orientation-incompatible alignments; `0` drops them
    /// (salmon's default `ignoreIncompat` behavior)
    pub incompat_prior: f64,
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
            no_error_model: false,
            seq_bias: false,
            gc_bias: false,
            pos_bias: false,
            incompat_prior: 0.0,
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
    pub num_eq_classes: usize,
    pub frag_len_mean: f64,
    pub frag_len_sd: f64,
    pub length_classes: Vec<u32>,
    pub frag_len_dist: Vec<f64>,
    pub start_time: String,
    pub bias_dump: salmon_model::dumps::BiasDump,
    /// per-transcript (unique, ambiguous) fragment counts for `ambig_info.tsv`
    pub ambig: (Vec<u32>, Vec<u32>),
}

/// Current local time as an asctime-style string, matching salmon's timestamps.
fn asctime_now() -> String {
    jiff::Zoned::now().strftime("%a %b %e %H:%M:%S %Y").to_string()
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
        let (Some(mtid), Some(mpos)) = (r1.mate_tid, r1.mate_pos) else { continue };
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
            placements.push(Placement { tid: recs[i].tid, idxs: vec![i] });
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
        let r1 = idxs.iter().map(|&i| &recs[i]).find(|r| r.is_read1).unwrap_or(&recs[idxs[0]]);
        let r2 = idxs.iter().map(|&i| &recs[i]).find(|r| !r.is_read1).unwrap_or(&recs[idxs[1]]);
        let (r1_fw, r2_fw) = (!r1.is_reverse, !r2.is_reverse);
        let (orientation, strandedness) = if r1_fw != r2_fw {
            let (fw, rc) = if r1_fw { (r1, r2) } else { (r2, r1) };
            let orientation = if fw.pos <= rc.pos { ReadOrientation::Toward } else { ReadOrientation::Away };
            let strandedness = if r1_fw { ReadStrandedness::SA } else { ReadStrandedness::AS };
            (orientation, strandedness)
        } else {
            (ReadOrientation::Same, if r1_fw { ReadStrandedness::S } else { ReadStrandedness::A })
        };
        (Some(LibraryFormat::new(ReadType::PairedEnd, orientation, strandedness)), r1_fw, MateStatus::PairedEndPaired)
    } else {
        let r = &recs[idxs[0]];
        let status = if r.is_read1 { MateStatus::PairedEndLeft } else { MateStatus::PairedEndRight };
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
        Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(file)))
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

/// Stream a BAM, grouping consecutive mapped records that share a read name into
/// fragments, and invoke `f` once per fragment. (Used twice: train the error
/// model, then build equivalence classes — avoids holding the whole BAM in memory.)
fn for_each_fragment<F>(bam_path: &Path, need_seq: bool, mut f: F) -> Result<()>
where
    F: FnMut(&[FragRecord]),
{
    let mut reader = bam::io::Reader::new(
        std::fs::File::open(bam_path).with_context(|| format!("opening {}", bam_path.display()))?,
    );
    let _header = reader.read_header().context("reading BAM header")?;

    let mut cur_name: Vec<u8> = Vec::new();
    let mut have_group = false;
    let mut group: Vec<FragRecord> = Vec::new();

    for result in reader.records() {
        let record = result.context("reading BAM record")?;
        if record.flags().is_unmapped() {
            continue;
        }
        let Some(name) = record.name() else { continue };
        let name_bytes: &[u8] = name.as_ref();
        let cname = canonical_name(name_bytes);
        if !have_group {
            cur_name = cname.to_vec();
            have_group = true;
        } else if cname != cur_name.as_slice() {
            f(&group);
            group.clear();
            cur_name.clear();
            cur_name.extend_from_slice(cname);
        }

        let Some(Ok(tid)) = record.reference_sequence_id() else { continue };
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
            .filter(|(o, _)| matches!(o, AlnOp::Match | AlnOp::SeqMatch | AlnOp::SeqMismatch | AlnOp::Del | AlnOp::RefSkip))
            .map(|(_, l)| l)
            .sum();
        let ops = if need_seq { ops } else { Vec::new() };
        let score = record
            .data()
            .get(&Tag::ALIGNMENT_SCORE)
            .and_then(|r| r.ok())
            .and_then(|v| value_as_i32(&v))
            .unwrap_or(0);
        let frag_len = record.template_length().abs();
        let flags = record.flags();
        let is_reverse = flags.is_reverse_complemented();
        let is_read1 = flags.is_first_segment();
        // Mate linkage as the aligner recorded it (RNEXT/PNEXT); a mate that is
        // unmapped or absent leaves these `None`, making the record an orphan.
        let mate_tid = (!flags.is_mate_unmapped())
            .then(|| record.mate_reference_sequence_id().and_then(|r| r.ok()))
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
        group.push(FragRecord {
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
        });
    }
    if have_group {
        f(&group);
    }
    Ok(())
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
    let Some(hd) = header.header() else { return false };
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

/// Run alignment-based quantification end-to-end.
pub fn quantify_alignments(opts: &AlignQuantOptions) -> Result<AlignQuantResult> {
    let start_time = asctime_now();
    let mut reader = bam::io::Reader::new(
        std::fs::File::open(&opts.bam).with_context(|| format!("opening {}", opts.bam.display()))?,
    );
    let header = reader.read_header().context("reading BAM header")?;

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
    let mut fld = FragmentLengthDistribution::new(1.0, 1000, 250.0, 25.0, 4, 0.5, 1);

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

    let mut model = use_error_model.then(|| AlignmentModel::new(1.0, 4));
    // Online (dual-phase) abundances: develop running estimates so the error
    // model and bias models are trained/collected with abundance-aware posteriors
    // in a single streaming pass (salmon's online phase), rather than two passes.
    let ref_lens_u64: Vec<u64> = lengths.iter().map(|&l| l as u64).collect();
    let online = (use_error_model || bias_on)
        .then(|| salmon_infer::OnlineInference::new(&ref_lens_u64, 0.05, 0.65, 5_000_000));

    // Observed bias accumulators (single-threaded → owned).
    use salmon_model::seqbias::{CONTEXT_LEFT, CONTEXT_LENGTH, CONTEXT_RIGHT};
    let mut seq_obs = opts
        .seq_bias
        .then(|| (salmon_model::SBModel::new(), salmon_model::SBModel::new()));
    let mut gc_obs = opts.gc_bias.then(salmon_model::GcFragModel::default_model);
    let gc_prefix: Vec<Vec<u32>> = if opts.gc_bias {
        ref_bytes.iter().map(|s| salmon_model::gc_prefix(s)).collect()
    } else {
        Vec::new()
    };
    let length_quantiles: Option<Vec<u32>> = opts.pos_bias.then(|| {
        salmon_model::compute_length_quantiles(&lengths, salmon_model::NUM_LENGTH_CLASSES)
    });
    let length_class: Option<Vec<usize>> = length_quantiles.as_ref().map(|q| {
        lengths.iter().map(|&l| salmon_model::length_class_index(q, l)).collect()
    });
    let mut pos_obs = opts.pos_bias.then(|| {
        let mk = || {
            (0..salmon_model::NUM_LENGTH_CLASSES)
                .map(|_| salmon_model::SimplePosBias::default())
                .collect::<Vec<_>>()
        };
        (mk(), mk())
    });

    // ---- single online pass: train error model + collect bias + eq-classes ---
    const MINIBATCH: u64 = 1000;
    // salmon's LOG_EPSILON = log(0.375e-10): the orphan / implausible-length penalty.
    const LOG_EPSILON: f64 = -23.998_158_637_57;
    // A paired library expects two mates; a single mate to a transcript is then an
    // "unexpected orphan" and is fragment-length-penalized. (Single-end libs aren't.)
    let paired_lib = !matches!(opts.lib_type.as_str(), "U" | "SF" | "SR" | "S");
    // Orientation-compatibility filtering (salmon): drop alignments whose
    // orientation is incompatible with the expected library type (and drop a
    // fragment entirely if all its alignments are incompatible — these become
    // zero-probability fragments). Skipped under auto (`A`) library type.
    let expected_format = match opts.lib_type.as_str() {
        "A" => None,
        s => LibraryFormat::parse(s).ok(),
    };
    let ignore_incompat = opts.incompat_prior <= 0.0;
    let mut num_processed = 0u64;
    let mut num_mapped = 0u64;
    let mut frag_count = 0u64;
    let mut log_fm = 0.0;
    for_each_fragment(&opts.bam, use_error_model || bias_on, |recs| {
        num_processed += 1;
        if recs.is_empty() {
            return;
        }
        num_mapped += 1;
        // one forgetting-mass timestep per minibatch of fragments
        if let Some(o) = online.as_ref() {
            if frag_count % MINIBATCH == 0 {
                log_fm = o.next_log_fm();
            }
        }
        frag_count += 1;

        // Pair records into the placements the aligner reported (proper pairs +
        // orphans), NOT a cross-product of every read1/read2 on a transcript.
        let placements = pair_records(recs);
        let frag_len = recs.iter().map(|r| r.frag_len).max().unwrap_or(0);
        if frag_len > 0 {
            fld.add_val(frag_len as usize, 0.0);
        }
        let use_aux = online.as_ref().is_none_or(|o| o.num_assigned() >= 5000);

        // Per surviving placement (one reported alignment): conditional log-weight
        // (eq-class) + online log-aux + fragment geometry. A placement that fails
        // the orientation-compatibility filter is dropped here.
        let mut sp_tid: Vec<u32> = Vec::with_capacity(placements.len());
        let mut sp_eq: Vec<f64> = Vec::with_capacity(placements.len());
        let mut sp_online: Vec<f64> = Vec::with_capacity(placements.len());
        let mut sp_geom: Vec<(usize, usize, bool)> = Vec::with_capacity(placements.len());
        let mut sp_pl: Vec<usize> = Vec::with_capacity(placements.len());
        for (pi, pl) in placements.iter().enumerate() {
            let tid = pl.tid;
            let idxs = &pl.idxs;
            let refseq = ref_bytes.get(tid as usize).map(|v| v.as_slice()).unwrap_or(&[]);
            // Conditional log-weight basis = salmon's `errLike`. With the error model
            // it is Σ(fg−bg) over the mate(s). Without the error model salmon leaves
            // errLike = LOG_1 for CIGAR-bearing aligners (it only uses the alignment
            // score for its own CIGAR-less aligners, pufferfish/rapmap, via
            // `useASWithoutCIGAR`). To match, a disabled error model yields a uniform
            // (0.0) basis rather than an AS-derived one.
            let basis = if let Some(m) = model.as_ref() {
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
            // fragment geometry + length-normalized online aux (start-pos + FLD)
            let rl = lengths[tid as usize] as i32;
            let flen = idxs.iter().map(|&i| recs[i].frag_len).max().unwrap_or(0);
            let proper = idxs.len() >= 2 && flen > 0;
            let frag_start = recs[idxs[0]].pos;
            let frag_end = frag_start + (flen.max(1) as usize) - 1;
            let start_pos = if proper && flen <= rl {
                -(((rl - flen + 1) as f64).ln())
            } else {
                -((rl.max(1) as f64).ln())
            };
            // Fragment-length probability (salmon's logFragProb): the FLD pmf for
            // a proper pair, or LOG_EPSILON for an unexpected orphan in a paired
            // library. This is part of the eq-class conditional weight (salmon's
            // auxProb = logFragProb + errLike + compat), NOT just the abundance
            // posterior — it down-weights placements whose implied insert size is
            // implausible (critical for permissive aligners that report many
            // discordant/wrong-length multimappers).
            let log_frag_prob = if proper {
                if use_aux { fld.pmf(flen as usize) } else { 0.0 }
            } else if paired_lib {
                LOG_EPSILON
            } else {
                0.0
            };
            let mut aux = basis + log_frag_prob;
            // orientation compatibility (salmon's logAlignCompatProb): drop the
            // alignment if incompatible and we're ignoring incompatible fragments,
            // else down-weight it by the incompatibility prior.
            if let Some(exp) = expected_format {
                let (obs, is_fw, status) = frag_format(recs, idxs);
                if !is_compatible(exp, obs, is_fw, status) {
                    if ignore_incompat {
                        continue; // this placement contributes nothing
                    }
                    aux += opts.incompat_prior.ln();
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
            return;
        }

        // Aggregate surviving placements by distinct transcript id (sorted): a
        // transcript that the fragment hits in several reported placements appears
        // once in the eq-class, its weight the logsumexp over those placements.
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

        // abundance-aware posteriors (online), per distinct transcript, drive
        // model/bias training
        let post: Vec<f64> = match online.as_ref() {
            Some(o) => {
                let maps: Vec<(u32, f64)> = tids.iter().cloned().zip(online_log.iter().cloned()).collect();
                o.assign_fragment(&maps, log_fm)
            }
            None => weights.clone(),
        };

        // train the error model + collect bias models, weighted by posteriors;
        // a transcript's posterior is split across its placements by their
        // within-transcript online softmax weight.
        let collecting = online.as_ref().is_none_or(|o| o.collecting());
        if collecting {
            for (ti, (tid, ks)) in agg.iter().enumerate() {
                let p_tid = post[ti];
                if p_tid <= 0.0 {
                    continue;
                }
                let online_log_tid = online_log[ti];
                let refseq = ref_bytes.get(*tid as usize).map(|v| v.as_slice()).unwrap_or(&[]);
                for &k in ks {
                    let p = p_tid * (sp_online[k] - online_log_tid).exp();
                    if p <= 0.0 {
                        continue;
                    }
                    let idxs = &placements[sp_pl[k]].idxs;
                    if let Some(m) = model.as_mut() {
                        if !refseq.is_empty() {
                            let lw = log_fm + p.ln();
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

                // Per-read 5' positions, by each read's *actual* strand (salmon's
                // logic): the forward read's 5' feeds the FW model, the reverse
                // read's 5' (its rightmost ref coordinate, RC'd) feeds the RC
                // model. For a concordant pair only the inward case (fwPos < rcPos)
                // is collected. This is correct for any library orientation, not
                // just FR.
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

                // sequence bias: forward 5' -> FW model (window [5'-3, 5'+6));
                // reverse 5' -> RC model (window [5'-5, 5'+4), reverse-complemented).
                if let Some(obs) = seq_obs.as_mut() {
                    if let Some(five) = fwd_five {
                        let s = five as i32 - CONTEXT_LEFT as i32;
                        if s >= 0 && (s as usize + CONTEXT_LENGTH) <= rl {
                            obs.0.add_context(&refseq[s as usize..s as usize + CONTEXT_LENGTH], false, p);
                        }
                    }
                    if let Some(five) = rev_five {
                        let s = five as i32 - CONTEXT_RIGHT as i32;
                        if s >= 0 && (s as usize + CONTEXT_LENGTH) <= rl {
                            obs.1.add_context(&refseq[s as usize..s as usize + CONTEXT_LENGTH], true, p);
                        }
                    }
                }
                // fragment-GC bias is span-based (orientation-independent).
                if let (Some(gc), true) = (gc_obs.as_mut(), proper && fe < rl) {
                    if let Some((ff, cf)) =
                        salmon_model::gc_desc(&gc_prefix[*tid as usize], rl as i32, fs as i32, fe as i32)
                    {
                        gc.inc(ff, cf, p);
                    }
                }
                // positional bias: forward 5' -> fwd model, reverse 5' -> RC model.
                if let Some(pos) = pos_obs.as_mut() {
                    let lc = length_class.as_ref().unwrap()[*tid as usize];
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

        let group = if opts.range_factorization_bins > 0 {
            let bins = range_factorize_bins(&weights, opts.range_factorization_bins);
            TranscriptGroup::with_bins(tids, bins)
        } else {
            TranscriptGroup::from_sorted(tids)
        };
        eq_builder.add_group(group, weights, 1);
    })?;

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
    let mut em = optimize(&collapsed, num_refs, &opts.em);

    // ---- bias-corrected effective lengths (shared with reads mode) ----------
    let mut bias_dump = salmon_model::dumps::BiasDump::default();
    if bias_on {
        let log_pmf = fld.log_pmf();
        let pmf_lin: Vec<f64> = log_pmf.iter().map(|lp| lp.exp()).collect();
        let (fld_cdf, fld_low, fld_high) = salmon_model::seqbias::fld_cdf_and_bounds(&pmf_lin);
        let k = if opts.seq_bias { CONTEXT_LENGTH } else { 1 };
        let refseq_of = |t: usize| ref_bytes[t].as_slice();

        let seq = seq_obs.map(|(mut of, mut or)| {
            of.normalize();
            or.normalize();
            let (ef, er) =
                salmon_model::build_expected(num_refs, refseq_of, &em.alphas, &eff_lengths, &fld_cdf);
            (of, or, ef, er)
        });
        if let Some((of, or, ef, er)) = seq.as_ref() {
            bias_dump.obs5_seq = of.dump().to_vec();
            bias_dump.obs3_seq = or.dump().to_vec();
            bias_dump.exp5_seq = ef.dump().to_vec();
            bias_dump.exp3_seq = er.dump().to_vec();
        }
        let gc_ratio_model = if let Some(mut obs) = gc_obs {
            let mut exp = salmon_model::build_expected_gc(
                num_refs,
                refseq_of,
                |t| gc_prefix[t].as_slice(),
                &em.alphas,
                &eff_lengths,
                &fld_cdf,
                fld_low,
                fld_high,
                salmon_model::gcbias::DEFAULT_COND_BINS,
                salmon_model::gcbias::DEFAULT_GC_BINS,
                k,
                salmon_model::GC_SAMP_STRIDE,
            );
            obs.normalize();
            exp.normalize();
            bias_dump.obs_gc = obs.dump().to_vec();
            bias_dump.exp_gc = exp.dump().to_vec();
            Some(salmon_model::gc_ratio(&mut obs, &mut exp, salmon_model::gcbias::GC_MAX_RATIO))
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
                &em.alphas,
                &eff_lengths,
                &fld_cdf,
                length_quantiles.as_ref().unwrap(),
                k,
            );
            (of, or, ef, er)
        });
        if let Some((ofw, orc, efw, erc)) = pos_models.as_ref() {
            let masses = |v: &[salmon_model::SimplePosBias]| v.iter().map(|m| m.masses().to_vec()).collect();
            bias_dump.obs5_pos = masses(ofw);
            bias_dump.obs3_pos = masses(orc);
            bias_dump.exp5_pos = masses(efw);
            bias_dump.exp3_pos = masses(erc);
        }

        for tid in 0..num_refs {
            if em.alphas[tid] < 1e-8 {
                continue;
            }
            let s = ref_bytes[tid].as_slice();
            let pos_vecs: Option<(Vec<f64>, Vec<f64>)> = pos_models.as_ref().map(|(ofw, orc, efw, erc)| {
                let lc = length_class.as_ref().unwrap()[tid];
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
                seq: seq.as_ref().map(|(of, or, ef, er)| (of, ef, or, er)),
                gc: gc_ratio_model.as_ref().map(|g| (g, gc_prefix[tid].as_slice())),
                pos: pos_vecs.as_ref().map(|(pf, pr)| (pf.as_slice(), pr.as_slice())),
            };
            eff_lengths[tid] = salmon_model::corrected_effective_length_full(
                s,
                &fld_cdf,
                fld_low,
                fld_high,
                &bias,
                eff_lengths[tid],
                salmon_model::GC_SAMP_STRIDE,
            );
        }
        collapsed.update_eff_lengths(&eff_lengths);
        em = optimize(&collapsed, num_refs, &opts.em);
    }
    let counts = em.alphas;

    let rates: Vec<f64> = (0..num_refs)
        .map(|i| if eff_lengths[i] > 0.0 { counts[i] / eff_lengths[i] } else { 0.0 })
        .collect();
    let rate_sum: f64 = rates.iter().sum();
    let tpm: Vec<f64> = rates
        .iter()
        .map(|r| if rate_sum > 0.0 { r / rate_sum * 1e6 } else { 0.0 })
        .collect();

    let length_classes =
        salmon_model::compute_length_quantiles(&lengths, salmon_model::NUM_LENGTH_CLASSES);
    let ambig = salmon_infer::ambiguity_counts(&salmon_infer::PackedEqClasses::from_collapsed(
        &collapsed, num_refs,
    ));
    let result = AlignQuantResult {
        names,
        lengths,
        eff_lengths,
        tpm,
        counts,
        num_processed,
        num_mapped,
        num_eq_classes,
        frag_len_mean: fld.mean(),
        frag_len_sd: fld.sd(),
        length_classes,
        frag_len_dist: fld.log_pmf().iter().map(|lp| lp.exp()).collect(),
        start_time,
        bias_dump,
        ambig,
    };
    write_outputs(opts, &result)?;
    Ok(result)
}

fn write_outputs(opts: &AlignQuantOptions, res: &AlignQuantResult) -> Result<()> {
    let dir = &opts.output_dir;
    std::fs::create_dir_all(dir.join("aux_info")).context("creating output dirs")?;

    // quant.sf
    let mut w = std::io::BufWriter::new(std::fs::File::create(dir.join("quant.sf"))?);
    writeln!(w, "Name\tLength\tEffectiveLength\tTPM\tNumReads")?;
    for i in 0..res.names.len() {
        writeln!(
            w,
            "{}\t{}\t{:.3}\t{:.6}\t{:.3}",
            res.names[i], res.lengths[i], res.eff_lengths[i], res.tpm[i], res.counts[i]
        )?;
    }
    w.flush()?;

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
        num_bootstraps: 0,
        call: "quant".to_string(),
        start_time: res.start_time.clone(),
        end_time: asctime_now(),
    };
    std::fs::write(
        dir.join("aux_info").join("meta_info.json"),
        serde_json::to_string_pretty(&meta)?,
    )?;

    // libParams/flenDist.txt, logs/salmon_quant.log, and the aux dumps (shared).
    std::fs::create_dir_all(dir.join("libParams")).context("creating libParams")?;
    salmon_model::dumps::write_flen_dist(&dir.join("libParams").join("flenDist.txt"), &res.frag_len_dist)
        .context("writing flenDist.txt")?;
    salmon_model::dumps::write_fld_dump(&dir.join("aux_info").join("fld.gz"), &res.frag_len_dist)
        .context("writing fld.gz")?;
    salmon_model::dumps::write_aux_bias_dumps(&dir.join("aux_info"), &res.bias_dump)
        .context("writing aux bias dumps")?;
    std::fs::create_dir_all(dir.join("logs")).context("creating logs")?;
    let log = format!(
        "salmon (rust port, alignment mode) v{SALMON_VERSION}\nstart: {}\nend:   {}\n\
         library type: {}\nobserved fragments: {}\nmapped fragments:   {}\nmapping rate: {pct:.4}%\n\
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
    std::fs::write(dir.join("logs").join("salmon_quant.log"), log).context("writing salmon_quant.log")?;

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
    std::fs::write(dir.join("cmd_info.json"), serde_json::to_string_pretty(&cmd)?)?;
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
