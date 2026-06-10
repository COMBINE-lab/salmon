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

use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result};
use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use serde::Serialize;

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

/// Per-transcript accumulator for the current fragment.
#[derive(Default, Clone, Copy)]
struct TidAccum {
    score: i32,
    frag_len: i32,
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
    let mut num_processed = 0u64;
    let mut num_mapped = 0u64;

    // Group consecutive records by read name; accumulate per-transcript scores.
    let mut cur_name: Vec<u8> = Vec::new();
    let mut have_group = false;
    let mut accum: HashMap<u32, TidAccum> = HashMap::new();

    let mut flush = |accum: &mut HashMap<u32, TidAccum>,
                     fld: &FragmentLengthDistribution,
                     num_processed: &mut u64,
                     num_mapped: &mut u64| {
        *num_processed += 1;
        if accum.is_empty() {
            return;
        }
        *num_mapped += 1;

        // best per-transcript score and a representative fragment length
        let best = accum.values().map(|a| a.score).max().unwrap_or(0);
        let frag_len = accum.values().map(|a| a.frag_len).max().unwrap_or(0);
        if frag_len > 0 {
            fld.add_val(frag_len as usize, 0.0);
        }

        let mut entries: Vec<(u32, f64)> = accum
            .iter()
            .map(|(&tid, a)| (tid, (-opts.score_exp * (best - a.score) as f64).exp()))
            .collect();
        entries.sort_by_key(|e| e.0);
        let wsum: f64 = entries.iter().map(|e| e.1).sum();
        let tids: Vec<u32> = entries.iter().map(|e| e.0).collect();
        let mut weights: Vec<f64> = entries.iter().map(|e| e.1).collect();
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
        accum.clear();
    };

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
            flush(&mut accum, &fld, &mut num_processed, &mut num_mapped);
            cur_name.clear();
            cur_name.extend_from_slice(name_bytes);
        }

        let Some(Ok(tid)) = record.reference_sequence_id() else {
            continue;
        };
        let score = record
            .data()
            .get(&Tag::ALIGNMENT_SCORE)
            .and_then(|r| r.ok())
            .and_then(|v| value_as_i32(&v))
            .unwrap_or(0);
        let frag_len = record.template_length().abs();

        let e = accum.entry(tid as u32).or_default();
        e.score += score;
        if frag_len > e.frag_len {
            e.frag_len = frag_len;
        }
    }
    if have_group {
        flush(&mut accum, &fld, &mut num_processed, &mut num_mapped);
    }

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
