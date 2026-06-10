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

/// Group a fragment's records by transcript id (ids sorted ascending; each id's
/// record indices sorted by reference position, so index 0 is the left mate).
fn group_by_tid(recs: &[FragRecord]) -> Vec<(u32, Vec<usize>)> {
    let mut m: HashMap<u32, Vec<usize>> = HashMap::new();
    for (i, r) in recs.iter().enumerate() {
        m.entry(r.tid).or_default().push(i);
    }
    let mut v: Vec<(u32, Vec<usize>)> = m.into_iter().collect();
    v.sort_by_key(|(t, _)| *t);
    for (_, idxs) in &mut v {
        idxs.sort_by_key(|&i| recs[i].pos);
    }
    v
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

/// Load a (optionally gzip'd) transcriptome FASTA and return 2-bit-encoded
/// sequences aligned to the BAM's reference order (`names`); a name absent from
/// the FASTA yields an empty sequence (its error-model contributions are skipped).
fn load_ref_2bit(path: &Path, names: &[String]) -> Result<Vec<Vec<u8>>> {
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
            cur_seq.extend(line.trim_end().bytes().map(base_2bit));
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
        if !have_group {
            cur_name = name_bytes.to_vec();
            have_group = true;
        } else if name_bytes != cur_name.as_slice() {
            f(&group);
            group.clear();
            cur_name.clear();
            cur_name.extend_from_slice(name_bytes);
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
        let ops: Vec<(AlnOp, usize)> = if need_seq {
            record
                .cigar()
                .iter()
                .filter_map(|r| r.ok())
                .map(|op| (kind_to_op(op.kind()), op.len()))
                .collect()
        } else {
            Vec::new()
        };
        let score = record
            .data()
            .get(&Tag::ALIGNMENT_SCORE)
            .and_then(|r| r.ok())
            .and_then(|v| value_as_i32(&v))
            .unwrap_or(0);
        let frag_len = record.template_length().abs();
        group.push(FragRecord { tid: tid as u32, pos, read_2bit, ops, score, frag_len });
    }
    if have_group {
        f(&group);
    }
    Ok(())
}

/// Run alignment-based quantification end-to-end.
pub fn quantify_alignments(opts: &AlignQuantOptions) -> Result<AlignQuantResult> {
    let mut reader = bam::io::Reader::new(
        std::fs::File::open(&opts.bam).with_context(|| format!("opening {}", opts.bam.display()))?,
    );
    let header = reader.read_header().context("reading BAM header")?;

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

    // The alignment error model needs the transcriptome (salmon requires `-t`).
    let use_error_model = opts.transcripts.is_some() && !opts.no_error_model;
    let ref_2bit: Vec<Vec<u8>> = if use_error_model {
        load_ref_2bit(opts.transcripts.as_ref().unwrap(), &names)?
    } else {
        Vec::new()
    };
    let mut model = use_error_model.then(|| AlignmentModel::new(1.0, 4));

    // ---- pass 1: train the error model (uniform per-fragment posterior) -----
    // salmon trains the model online weighted by the evolving posterior; here we
    // train once with a per-fragment-uniform posterior (each fragment contributes
    // total weight 1), which yields essentially the same converged error profile
    // since it is dominated by confidently/uniquely assigned reads.
    if let Some(model) = model.as_mut() {
        for_each_fragment(&opts.bam, true, |recs| {
            if recs.is_empty() {
                return;
            }
            let by_tid = group_by_tid(recs);
            let lw = (1.0 / by_tid.len() as f64).ln();
            for (tid, idxs) in &by_tid {
                let refseq = &ref_2bit[*tid as usize];
                if refseq.is_empty() {
                    continue;
                }
                for (rank, &i) in idxs.iter().enumerate() {
                    let r = &recs[i];
                    // left mate = the alignment with the smaller reference position
                    model.update(&r.read_2bit, refseq, r.pos, &r.ops, rank == 0, lw);
                }
            }
        })?;
    }

    // ---- pass 2: build the (weighted) equivalence classes -------------------
    let mut num_processed = 0u64;
    let mut num_mapped = 0u64;
    for_each_fragment(&opts.bam, use_error_model, |recs| {
        num_processed += 1;
        if recs.is_empty() {
            return;
        }
        num_mapped += 1;
        let by_tid = group_by_tid(recs);
        let frag_len = recs.iter().map(|r| r.frag_len).max().unwrap_or(0);
        if frag_len > 0 {
            fld.add_val(frag_len as usize, 0.0);
        }

        let mut tids: Vec<u32> = Vec::with_capacity(by_tid.len());
        let mut logw: Vec<f64> = Vec::with_capacity(by_tid.len());
        match model.as_ref() {
            // error model: conditional log-weight = Σ_mates (fg − bg)
            Some(model) => {
                for (tid, idxs) in &by_tid {
                    let refseq = &ref_2bit[*tid as usize];
                    let mut ll = 0.0;
                    for (rank, &i) in idxs.iter().enumerate() {
                        let r = &recs[i];
                        if !refseq.is_empty() {
                            let (fg, bg) =
                                model.log_likelihood(&r.read_2bit, refseq, r.pos, &r.ops, rank == 0);
                            ll += fg - bg;
                        }
                    }
                    tids.push(*tid);
                    logw.push(ll);
                }
            }
            // AS fallback (no transcriptome / error model disabled)
            None => {
                let best = by_tid
                    .iter()
                    .map(|(_, idxs)| idxs.iter().map(|&i| recs[i].score).sum::<i32>())
                    .max()
                    .unwrap_or(0);
                for (tid, idxs) in &by_tid {
                    let as_sum: i32 = idxs.iter().map(|&i| recs[i].score).sum();
                    tids.push(*tid);
                    logw.push(-opts.score_exp * (best - as_sum) as f64);
                }
            }
        }

        // softmax → per-fragment conditional probabilities (sum to 1)
        let maxlw = logw.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mut weights: Vec<f64> = logw.iter().map(|&l| (l - maxlw).exp()).collect();
        let wsum: f64 = weights.iter().sum();
        if wsum > 0.0 {
            for w in &mut weights {
                *w /= wsum;
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

    // ---- effective lengths, EM, TPM (shared logic with reads mode) ----------
    // Base effective length = `refLen − E[L | L ≤ refLen]` (salmon's
    // `computeSmoothedEffectiveLengths`); see the reads-mode driver for the
    // rationale over the truncated-PMF estimate.
    fld.cache();
    let cond_means = fld.conditional_means();
    let mut eff_lengths = vec![0f64; num_refs];
    for (tid, len) in lengths.iter().enumerate() {
        eff_lengths[tid] = salmon_model::smoothed_effective_length(&cond_means, *len as usize);
    }

    let mut collapsed = eq_builder.finish();
    collapsed.update_eff_lengths(&eff_lengths);
    let num_eq_classes = collapsed.len();
    let counts = optimize(&collapsed, num_refs, &opts.em).alphas;

    let rates: Vec<f64> = (0..num_refs)
        .map(|i| if eff_lengths[i] > 0.0 { counts[i] / eff_lengths[i] } else { 0.0 })
        .collect();
    let rate_sum: f64 = rates.iter().sum();
    let tpm: Vec<f64> = rates
        .iter()
        .map(|r| if rate_sum > 0.0 { r / rate_sum * 1e6 } else { 0.0 })
        .collect();

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

    // aux_info/meta_info.json
    #[derive(Serialize)]
    struct MetaInfo {
        salmon_version: String,
        mapping_type: String,
        library_types: Vec<String>,
        num_valid_targets: usize,
        num_eq_classes: usize,
        num_processed: u64,
        num_mapped: u64,
        percent_mapped: f64,
        frag_length_mean: f64,
        call: String,
    }
    let pct = if res.num_processed > 0 {
        100.0 * res.num_mapped as f64 / res.num_processed as f64
    } else {
        0.0
    };
    let meta = MetaInfo {
        salmon_version: SALMON_VERSION.to_string(),
        mapping_type: "alignment".to_string(),
        library_types: vec![opts.lib_type.clone()],
        num_valid_targets: res.names.len(),
        num_eq_classes: res.num_eq_classes,
        num_processed: res.num_processed,
        num_mapped: res.num_mapped,
        percent_mapped: pct,
        frag_length_mean: res.frag_len_mean,
        call: "quant-alignment".to_string(),
    };
    std::fs::write(
        dir.join("aux_info").join("meta_info.json"),
        serde_json::to_string_pretty(&meta)?,
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
