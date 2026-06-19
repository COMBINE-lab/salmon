//! Output writers: `quant.sf`, `cmd_info.json`, `lib_format_counts.json`,
//! `aux_info/meta_info.json`, `libParams/flenDist.txt`, and `logs/salmon_quant.log`.
//! Functional parity with salmon's column/field layout (the binary aux bias
//! dumps and the index-hash meta fields are deferred to a later step).

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;

use crate::{QuantOptions, QuantResult};

pub(crate) const SALMON_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Write the standard salmon output set into `opts.output_dir`.
pub fn write_outputs(opts: &QuantOptions, res: &QuantResult) -> Result<()> {
    let dir = &opts.output_dir;
    std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    std::fs::create_dir_all(dir.join("aux_info")).context("creating aux_info")?;

    // Under `--skipQuant` there are no abundances, so (like salmon) no quant.sf;
    // the eq-classes, metadata, and library counts are still written.
    if !opts.skip_quant {
        write_quant_sf(&dir.join("quant.sf"), res, opts.sig_digits as usize)?;
    }
    write_cmd_info(&dir.join("cmd_info.json"), opts)?;
    write_lib_counts(&dir.join("lib_format_counts.json"), opts, res)?;
    write_meta_info(&dir.join("aux_info").join("meta_info.json"), opts, res)?;
    write_ambig_info(&dir.join("aux_info").join("ambig_info.tsv"), res)?;
    std::fs::create_dir_all(dir.join("libParams")).context("creating libParams")?;
    salmon_model::dumps::write_flen_dist(
        &dir.join("libParams").join("flenDist.txt"),
        &res.frag_len_dist,
    )
    .context("writing flenDist.txt")?;
    write_quant_log(&dir.join("logs").join("salmon_quant.log"), opts, res)?;
    salmon_model::dumps::write_fld_dump(&dir.join("aux_info").join("fld.gz"), &res.frag_len_dist)
        .context("writing fld.gz")?;
    salmon_model::dumps::write_aux_bias_dumps(&dir.join("aux_info"), &res.bias_dump)
        .context("writing aux bias dumps")?;
    if !res.bootstraps.is_empty() {
        write_bootstraps(&dir.join("aux_info").join("bootstrap"), res)?;
    }
    Ok(())
}
/// `logs/salmon_quant.log`: a concise run summary. salmon's log is more verbose,
/// but downstream tools key off the JSON; this records the salient counters.
fn write_quant_log(path: &Path, opts: &QuantOptions, res: &QuantResult) -> Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).context("creating logs")?;
    }
    let pct = if res.num_processed > 0 {
        100.0 * res.num_mapped as f64 / res.num_processed as f64
    } else {
        0.0
    };
    let end_time = crate::asctime_now();
    let mut s = String::new();
    s.push_str(&format!("salmon (rust port) v{SALMON_VERSION}\n"));
    s.push_str(&format!("start: {}\n", res.start_time));
    s.push_str(&format!("end:   {end_time}\n"));
    s.push_str(&format!("library type: {}\n", res.library_type));
    s.push_str(&format!(
        "mapping type: {}\n",
        if opts.sketch { "pseudo" } else { "mapping" }
    ));
    s.push_str(&format!("observed fragments: {}\n", res.num_processed));
    s.push_str(&format!("mapped fragments:   {}\n", res.num_mapped));
    s.push_str(&format!("mapping rate: {pct:.4}%\n"));
    s.push_str(&format!(
        "number of equivalence classes: {}\n",
        res.num_eq_classes
    ));
    s.push_str(&format!(
        "fragment length mean (sd): {:.2} ({:.2})\n",
        res.frag_len_mean, res.frag_len_sd
    ));
    std::fs::write(path, s).with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

/// `aux_info/ambig_info.tsv`: `UniqueCount\tAmbigCount` per transcript (id order).
fn write_ambig_info(path: &Path, res: &QuantResult) -> Result<()> {
    let (uniq, ambig) = &res.ambig;
    let mut w = std::io::BufWriter::new(
        std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?,
    );
    writeln!(w, "UniqueCount\tAmbigCount")?;
    for i in 0..uniq.len() {
        writeln!(w, "{}\t{}", uniq[i], ambig[i])?;
    }
    Ok(())
}

/// `aux_info/bootstrap/{names.tsv.gz, bootstraps.gz}`: salmon's posterior-sample
/// layout — `names.tsv.gz` lists the transcript names (tab-separated, id order),
/// `bootstraps.gz` is the gzip of raw little-endian `f64`s with each sample's
/// `num_txps` values written contiguously (tximport reads it directly).
fn write_bootstraps(dir: &Path, res: &QuantResult) -> Result<()> {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;

    let f = std::fs::File::create(dir.join("names.tsv.gz"))?;
    let mut enc = GzEncoder::new(f, Compression::new(6));
    for (i, name) in res.names.iter().enumerate() {
        if i > 0 {
            enc.write_all(b"\t")?;
        }
        enc.write_all(name.as_bytes())?;
    }
    enc.write_all(b"\n")?;
    enc.finish()?;

    let f = std::fs::File::create(dir.join("bootstraps.gz"))?;
    let mut enc = GzEncoder::new(f, Compression::new(6));
    for sample in &res.bootstraps {
        for v in sample {
            enc.write_all(&v.to_le_bytes())?;
        }
    }
    enc.finish()?;
    Ok(())
}

/// Reference indices to emit in `quant.sf`, in order, given the decoy layout.
///
/// The index numbering is `[transcripts][decoys][short transcripts]`: transcripts
/// occupy `[0, first_decoy_index)`, the decoy block `[first_decoy_index,
/// first_decoy_index + num_decoys)` is skipped, and sub-`k` "short" transcripts
/// (no k-mers, never seeded, always 0 reads) occupy the tail
/// `[first_decoy_index + num_decoys, total)` and are still reported.
///
/// `first_decoy_index == None` means no decoys (emit everything). A legacy index
/// that recorded decoys but not `num_decoys` (== 0) keeps the old behavior:
/// everything from `first_decoy_index` on is treated as decoy and dropped.
fn quant_row_indices(
    total: usize,
    first_decoy_index: Option<usize>,
    num_decoys: usize,
) -> impl Iterator<Item = usize> {
    let fdi = first_decoy_index.unwrap_or(total).min(total);
    let decoy_end = match first_decoy_index {
        Some(_) if num_decoys > 0 => (fdi + num_decoys).min(total),
        _ => total,
    };
    (0..fdi).chain(decoy_end..total)
}

/// `quant.sf`: `Name  Length  EffectiveLength  TPM  NumReads`.
fn write_quant_sf(path: &Path, res: &QuantResult, sig_digits: usize) -> Result<()> {
    let mut w = std::io::BufWriter::new(
        std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?,
    );
    writeln!(w, "Name\tLength\tEffectiveLength\tTPM\tNumReads")?;
    // The reference numbering is `[transcripts][decoys][short transcripts]`:
    // decoys occupy `[first_decoy_index, first_decoy_index + num_decoys)` and are
    // never quantified (excluded here, matching salmon), while sub-`k` "short"
    // transcripts (no k-mers, never seeded) are appended after the decoy block by
    // the index builder. We still report the shorts — always 0 reads / 0 TPM — so
    // `quant.sf` lists every input transcript.
    // salmon writes EffectiveLength and NumReads with `--sigDigits` decimals
    // (default 3) and TPM with the fixed `{:f}` (6-decimal) format.
    for i in quant_row_indices(res.names.len(), res.first_decoy_index, res.num_decoys) {
        writeln!(
            w,
            "{}\t{}\t{:.*}\t{:.6}\t{:.*}",
            res.names[i],
            res.lengths[i],
            sig_digits,
            res.eff_lengths[i],
            res.tpm[i],
            sig_digits,
            res.counts[i]
        )?;
    }
    w.flush()?;
    Ok(())
}

#[derive(Serialize)]
struct CmdInfo {
    salmon_version: String,
    index: String,
    #[serde(rename = "libType")]
    lib_type: String,
    output: String,
    mates1: Vec<String>,
    mates2: Vec<String>,
    #[serde(rename = "unmatedReads")]
    unmated_reads: Vec<String>,
    threads: usize,
    sketch: bool,
}

fn paths_to_strings(paths: &[std::path::PathBuf]) -> Vec<String> {
    paths.iter().map(|p| p.display().to_string()).collect()
}

fn write_cmd_info(path: &Path, opts: &QuantOptions) -> Result<()> {
    let info = CmdInfo {
        salmon_version: SALMON_VERSION.to_string(),
        index: opts.index_dir.display().to_string(),
        lib_type: opts.lib_type.clone(),
        output: opts.output_dir.display().to_string(),
        mates1: paths_to_strings(&opts.mates1),
        mates2: paths_to_strings(&opts.mates2),
        unmated_reads: paths_to_strings(&opts.unmated),
        threads: opts.num_threads,
        sketch: opts.sketch,
    };
    write_json(path, &info)
}

#[derive(Serialize)]
struct LibCounts {
    read_files: Vec<String>,
    expected_format: String,
    compatible_fragment_ratio: f64,
    num_compatible_fragments: u64,
    num_assigned_fragments: u64,
    num_frags_with_concordant_consistent_mappings: u64,
    num_frags_with_inconsistent_or_orphan_mappings: u64,
    strand_mapping_bias: f64,
}

fn write_lib_counts(path: &Path, opts: &QuantOptions, res: &QuantResult) -> Result<()> {
    // No strand-compatibility filtering yet: every mapped fragment is treated
    // as compatible with the declared library type.
    let mut read_files = paths_to_strings(&opts.mates1);
    read_files.extend(paths_to_strings(&opts.mates2));
    read_files.extend(paths_to_strings(&opts.unmated));
    let counts = LibCounts {
        read_files,
        expected_format: res.library_type.clone(),
        compatible_fragment_ratio: 1.0,
        num_compatible_fragments: res.num_mapped,
        num_assigned_fragments: res.num_mapped,
        num_frags_with_concordant_consistent_mappings: res.num_mapped - res.num_orphan,
        num_frags_with_inconsistent_or_orphan_mappings: res.num_orphan,
        strand_mapping_bias: 0.0,
    };
    write_json(path, &counts)
}

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
    keep_duplicates: bool,
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

fn write_meta_info(path: &Path, opts: &QuantOptions, res: &QuantResult) -> Result<()> {
    let percent_mapped = if res.num_processed > 0 {
        100.0 * res.num_mapped as f64 / res.num_processed as f64
    } else {
        0.0
    };
    let samp_type = if opts.num_bootstraps > 0 {
        "bootstrap"
    } else if opts.num_gibbs_samples > 0 {
        "gibbs"
    } else {
        "none"
    };
    // eq-class properties match salmon's default in-memory representation.
    let eq_class_properties = if opts.range_factorization_bins > 0 {
        vec!["range_factorized".to_string(), "gzipped".to_string()]
    } else {
        vec!["gzipped".to_string()]
    };
    let meta = MetaInfo {
        salmon_version: SALMON_VERSION.to_string(),
        samp_type: samp_type.to_string(),
        opt_type: if opts.em.use_vbem { "vb" } else { "em" }.to_string(),
        quant_errors: vec![],
        num_libraries: 1,
        library_types: vec![res.library_type.clone()],
        frag_dist_length: res.frag_len_dist.len(),
        frag_length_mean: res.frag_len_mean,
        frag_length_sd: res.frag_len_sd,
        seq_bias_correct: opts.seq_bias,
        gc_bias_correct: opts.gc_bias,
        num_bias_bins: 0,
        mapping_type: if opts.sketch { "pseudo" } else { "mapping" }.to_string(),
        keep_duplicates: res.keep_duplicates,
        index_seq_hash: res.index_seq_hash.clone(),
        index_name_hash: res.index_name_hash.clone(),
        index_seq_hash512: res.index_seq_hash512.clone(),
        index_name_hash512: res.index_name_hash512.clone(),
        index_decoy_seq_hash: res.index_decoy_seq_hash.clone(),
        index_decoy_name_hash: res.index_decoy_name_hash.clone(),
        // Valid targets = all non-decoy references (transcripts, including the
        // sub-k shorts we report with 0 counts); decoys are the explicit block.
        // Legacy indices (decoys present, `num_decoys` unrecorded) fall back to
        // the old suffix interpretation.
        num_valid_targets: if res.num_decoys > 0 {
            res.names.len().saturating_sub(res.num_decoys)
        } else {
            res.first_decoy_index.unwrap_or(res.names.len())
        },
        num_decoy_targets: if res.num_decoys > 0 {
            res.num_decoys
        } else {
            res.first_decoy_index
                .map_or(0, |fdi| res.names.len().saturating_sub(fdi))
        },
        num_eq_classes: res.num_eq_classes,
        serialized_eq_classes: opts.dump_eq || opts.dump_eq_weights,
        eq_class_properties,
        length_classes: res.length_classes.clone(),
        num_processed: res.num_processed,
        num_mapped: res.num_mapped,
        num_dovetail_fragments: res.num_dovetail_fragments,
        num_fragments_filtered_vm: res.num_fragments_filtered_vm,
        num_alignments_below_threshold_for_mapped_fragments_vm: res
            .num_alignments_below_threshold_for_mapped_fragments_vm,
        percent_mapped,
        num_decoy_fragments: res.num_decoy_fragments,
        num_bootstraps: res.bootstraps.len() as u32,
        call: "quant".to_string(),
        start_time: res.start_time.clone(),
        end_time: crate::asctime_now(),
    };
    write_json(path, &meta)
}

fn write_json<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let json = serde_json::to_string_pretty(value)?;
    std::fs::write(path, json).with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::quant_row_indices;

    fn rows(total: usize, fdi: Option<usize>, nd: usize) -> Vec<usize> {
        quant_row_indices(total, fdi, nd).collect()
    }

    #[test]
    fn no_decoys_emits_all() {
        assert_eq!(rows(5, None, 0), vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn decoys_then_shorts_skips_decoys_keeps_shorts() {
        // [0,3) transcripts, [3,5) decoys, [5,7) short transcripts.
        assert_eq!(rows(7, Some(3), 2), vec![0, 1, 2, 5, 6]);
    }

    #[test]
    fn decoys_as_clean_suffix_no_shorts() {
        assert_eq!(rows(5, Some(3), 2), vec![0, 1, 2]);
    }

    #[test]
    fn legacy_index_without_num_decoys_drops_suffix() {
        // num_decoys unrecorded (0) but decoys present -> old behavior: drop
        // everything from first_decoy_index on (no spurious decoy/short rows).
        assert_eq!(rows(7, Some(3), 0), vec![0, 1, 2]);
    }

    #[test]
    fn clamps_out_of_range_decoy_block() {
        // num_decoys overshooting total must not panic or emit phantom rows.
        assert_eq!(rows(5, Some(3), 99), vec![0, 1, 2]);
    }
}
