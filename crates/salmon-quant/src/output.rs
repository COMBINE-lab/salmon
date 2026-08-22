//! Output writers: `quant.sf`, `cmd_info.json`, `lib_format_counts.json`,
//! `aux_info/meta_info.json`, `libParams/flenDist.txt`, and `logs/salmon_quant.log`.
//! Functional parity with salmon's column/field layout (the binary aux bias
//! dumps and the index-hash meta fields are deferred to a later step).
//!
//! # What a salmon output directory contains
//!
//! * `quant.sf` — the result: one row per transcript with its length, effective
//!   length, TPM and estimated fragment count. This is what downstream analysis
//!   reads.
//! * `cmd_info.json` — exactly how the run was invoked, so a result can be
//!   reproduced.
//! * `aux_info/meta_info.json` — how the run went: mapping rate, fragment-length
//!   statistics, index hashes, timing, diagnostics. Tools like tximport and
//!   MultiQC key off this.
//! * `aux_info/bootstrap/` — posterior samples, when requested.
//! * `libParams/`, `logs/`, and the binary bias dumps — diagnostics.
//!
//! The exact filenames, column order and JSON field names are part of salmon's
//! contract with a large ecosystem of downstream tools, so they are matched
//! deliberately rather than chosen.

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;

use crate::{QuantOptions, QuantResult};

/// This build's version, read from Cargo.toml at compile time and recorded in
/// every output file for provenance.
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
///
/// A cheap confidence indicator: a transcript whose evidence is mostly unique is
/// well determined, one whose evidence is mostly shared is not.
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
///
/// The emitted set is restricted to the same rows as `quant.sf` via
/// `quant_row_indices` — decoys are excluded and the sub-`k` "short" transcripts
/// (always 0) are included — so both files align positionally and by name with
/// `quant.sf`, which is what tximport and friends assume when they read
/// `bootstraps.gz` against the `quant.sf` rows. (The per-sample vectors are
/// indexed by reference id over the full `[transcripts][decoys][shorts]` range;
/// without this filter they would carry the decoy columns and misalign.)
fn write_bootstraps(dir: &Path, res: &QuantResult) -> Result<()> {
    use flate2::write::GzEncoder;
    use flate2::Compression;
    std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;

    let row_indices: Vec<usize> =
        quant_row_indices(res.names.len(), res.first_decoy_index, res.num_decoys).collect();

    let f = std::fs::File::create(dir.join("names.tsv.gz"))?;
    let mut enc = GzEncoder::new(f, Compression::new(6));
    enc.write_all(&select_names_tsv(&res.names, &row_indices))?;
    enc.finish()?;

    let f = std::fs::File::create(dir.join("bootstraps.gz"))?;
    let mut enc = GzEncoder::new(f, Compression::new(6));
    for sample in &res.bootstraps {
        enc.write_all(&select_sample_bytes(sample, &row_indices))?;
    }
    enc.finish()?;
    Ok(())
}

/// Tab-separated, newline-terminated `names.tsv.gz` body for the selected rows.
fn select_names_tsv(names: &[String], row_indices: &[usize]) -> Vec<u8> {
    let mut out = Vec::new();
    for (j, &i) in row_indices.iter().enumerate() {
        if j > 0 {
            out.push(b'\t');
        }
        out.extend_from_slice(names[i].as_bytes());
    }
    out.push(b'\n');
    out
}

/// Little-endian `f64` bytes for one posterior sample, restricted to the selected
/// rows (in `row_indices` order) so the stream aligns positionally with `quant.sf`.
fn select_sample_bytes(sample: &[f64], row_indices: &[usize]) -> Vec<u8> {
    let mut out = Vec::with_capacity(row_indices.len() * 8);
    for &i in row_indices {
        out.extend_from_slice(&sample[i].to_le_bytes());
    }
    out
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
/// everything from `first_decoy_index` on is treated as decoy and dropped. The
/// build guarantees decoys are one contiguous block, so the emitted set is simply
/// the two non-decoy ranges below.
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
///
/// The four numbers per transcript: its annotated length, the effective length
/// used as the divisor, TPM (a within-sample relative abundance that already
/// accounts for length, summing to one million), and the estimated fragment
/// count — fractional, because ambiguous fragments are split between
/// transcripts.
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

/// `cmd_info.json`: the run's inputs and key settings.
///
/// `#[derive(Serialize)]` generates the JSON encoding from the struct, so the
/// field names below *are* the JSON keys — which is why several carry an explicit
/// `rename` to match salmon's camelCase spelling.
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

/// Render paths for the JSON outputs (they are recorded as text, not as
/// platform-specific path objects).
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

/// `lib_format_counts.json`: how the mapped fragments' observed library formats
/// relate to the declared (or detected) library type — the first thing to check
/// when a mapping rate looks wrong. The struct, its C++-matching field
/// semantics, and its `num_incompatible_fragments` extension (#1130) live in
/// [`salmon_core::LibFormatCountsFile`], shared with the alignment/RAD writer.
fn write_lib_counts(path: &Path, opts: &QuantOptions, res: &QuantResult) -> Result<()> {
    let mut read_files = paths_to_strings(&opts.mates1);
    read_files.extend(paths_to_strings(&opts.mates2));
    read_files.extend(paths_to_strings(&opts.unmated));
    let counts = salmon_core::LibFormatCountsFile::new(
        read_files,
        res.library_type.clone(),
        res.num_mapped,
        res.num_compatible_fragments,
        res.num_incompatible_fragments,
        &res.lib_format_counts,
    );
    counts.log_warnings();
    write_json(path, &counts)
}

/// `aux_info/meta_info.json`: everything a downstream tool or a human needs to
/// judge the run — provenance, settings, counters, timing and diagnostics.
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
    /// provenance of the two fields above: observed data or the
    /// `--fldMean`/`--fldSD` prior (see `FragLengthSource`)
    frag_length_source: &'static str,
    seq_bias_correct: bool,
    gc_bias_correct: bool,
    pos_bias_correct: bool,
    num_bias_bins: usize,
    mapping_type: String,
    /// how the index was built; omitted when the index predates recording it,
    /// rather than reported as a `false` the run cannot actually vouch for
    #[serde(skip_serializing_if = "Option::is_none")]
    keep_duplicates: Option<bool>,
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
    /// fragment mass dropped in the final min-alpha redistribution (every member
    /// of an equivalence class truncated); normally 0.
    inference_truncated_mass: f64,
    num_bootstraps: u32,
    /// mapped fragments placed as orphans (only one mate mapped)
    num_orphan: u64,
    /// whether every field here is a true observation of the run.
    ///
    /// Mapping the reads directly supplies almost everything, but not quite
    /// everything: an index built before it recorded how it was built cannot say,
    /// and guessing would be worse than saying so. Each gap is listed below with
    /// the upstream responsible, so it is actionable.
    meta_info_complete: bool,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    incomplete_meta_info_fields: Vec<crate::MissingMetaField>,
    /// range-factorization bins used for equivalence classes (0 = disabled)
    range_factorization_bins: u32,
    /// EM/VBEM convergence
    num_em_iterations: u32,
    em_converged: bool,
    /// library format observed by the auto-detector (null if not observed)
    detected_library_type: Option<String>,
    /// which inference path produced these results: `online` (the one-pass
    /// streaming quant), or `none` for a mapping-only `--skipQuant` pass (a
    /// `--deterministic` run's phase 2 rewrites this file with
    /// `deterministic`). Makes the online-to-deterministic transition
    /// auditable from existing output (#1140).
    inference_path: &'static str,
    /// total wall-clock seconds and peak resident set size (KiB, Linux)
    total_time_seconds: f64,
    peak_rss_kb: u64,
    /// machine-readable run diagnostics / bad-input warnings
    diagnostics: Vec<crate::Diagnostic>,
    call: String,
    start_time: String,
    end_time: String,
}

fn write_meta_info(path: &Path, opts: &QuantOptions, res: &QuantResult) -> Result<()> {
    // Mapping the reads answers nearly everything directly; what it cannot
    // answer comes from the index, which may predate recording it.
    let mut missing_meta_info_fields = Vec::new();
    if res.keep_duplicates.is_none() {
        missing_meta_info_fields.push(crate::MissingMetaField::index(
            "keep_duplicates",
            "this index predates recording how it was built; rebuild the index \
             to record it",
        ));
    }
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
        frag_length_source: res.frag_len_source.as_str(),
        seq_bias_correct: opts.seq_bias,
        gc_bias_correct: opts.gc_bias,
        pos_bias_correct: opts.pos_bias,
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
        inference_truncated_mass: res.inference_truncated_mass,
        num_bootstraps: res.bootstraps.len() as u32,
        num_orphan: res.num_orphan,
        meta_info_complete: missing_meta_info_fields.is_empty(),
        incomplete_meta_info_fields: missing_meta_info_fields,
        range_factorization_bins: opts.range_factorization_bins,
        num_em_iterations: res.em_iters,
        em_converged: res.em_converged,
        detected_library_type: res.detected_library_type.clone(),
        inference_path: if opts.skip_quant { "none" } else { "online" },
        total_time_seconds: res.total_seconds,
        peak_rss_kb: res.peak_rss_kb,
        diagnostics: res.diagnostics.clone(),
        call: "quant".to_string(),
        start_time: res.start_time.clone(),
        end_time: crate::asctime_now(),
    };
    write_json(path, &meta)
}

/// Write any serializable value as pretty-printed JSON. Pretty rather than
/// compact because these files are read by people as often as by tools.
fn write_json<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let json = serde_json::to_string_pretty(value)?;
    std::fs::write(path, json).with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{quant_row_indices, select_names_tsv, select_sample_bytes};

    /// Collect the emitted row indices for a given reference layout.
    fn rows(total: usize, fdi: Option<usize>, nd: usize) -> Vec<usize> {
        quant_row_indices(total, fdi, nd).collect()
    }

    /// The ordinary case: no decoys, so every reference is reported.
    #[test]
    fn no_decoys_emits_all() {
        assert_eq!(rows(5, None, 0), vec![0, 1, 2, 3, 4]);
    }

    /// The case that motivates the whole function: the trailing short
    /// transcripts sit *after* the decoy block, so a naive "drop everything from
    /// the first decoy" rule would silently lose real transcripts.
    #[test]
    fn decoys_then_shorts_skips_decoys_keeps_shorts() {
        // [0,3) transcripts, [3,5) decoys, [5,7) short transcripts.
        assert_eq!(rows(7, Some(3), 2), vec![0, 1, 2, 5, 6]);
    }

    /// With no short transcripts the decoy block really is the suffix, and both
    /// interpretations agree.
    #[test]
    fn decoys_as_clean_suffix_no_shorts() {
        assert_eq!(rows(5, Some(3), 2), vec![0, 1, 2]);
    }

    /// An index built before `num_decoys` was recorded must keep its original
    /// interpretation rather than being reinterpreted and emitting decoy rows.
    #[test]
    fn legacy_index_without_num_decoys_drops_suffix() {
        // num_decoys unrecorded (0) but decoys present -> old behavior: drop
        // everything from first_decoy_index on (no spurious decoy/short rows).
        assert_eq!(rows(7, Some(3), 0), vec![0, 1, 2]);
    }

    /// A corrupt or inconsistent header must not panic or invent rows.
    #[test]
    fn clamps_out_of_range_decoy_block() {
        // num_decoys overshooting total must not panic or emit phantom rows.
        assert_eq!(rows(5, Some(3), 99), vec![0, 1, 2]);
    }

    /// tximport reads `bootstraps.gz` positionally against the `quant.sf` rows,
    /// so a misalignment here would silently attribute every posterior sample to
    /// the wrong transcript.
    #[test]
    fn bootstrap_output_aligns_with_quant_rows() {
        // Reference layout: [t0 t1 t2][d0 d1][s0] — 3 transcripts, 2 decoys, 1
        // short. The per-sample vector is indexed by reference id over the full
        // range (decoys carry residual mass, the short is always 0). The emitted
        // bootstrap rows must drop the decoys, keep the short at 0, and stay in
        // quant.sf order so names.tsv.gz/bootstraps.gz align positionally.
        let names: Vec<String> = ["t0", "t1", "t2", "d0", "d1", "s0"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let sample = vec![10.0_f64, 20.0, 30.0, 7.0, 9.0, 0.0];
        let row_indices: Vec<usize> = quant_row_indices(names.len(), Some(3), 2).collect();
        assert_eq!(row_indices, vec![0, 1, 2, 5]);

        let names_tsv = select_names_tsv(&names, &row_indices);
        assert_eq!(names_tsv, b"t0\tt1\tt2\ts0\n");

        let bytes = select_sample_bytes(&sample, &row_indices);
        let got: Vec<f64> = bytes
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        // decoy values 7.0/9.0 excluded; short s0 present at 0.0.
        assert_eq!(got, vec![10.0, 20.0, 30.0, 0.0]);
    }
}
