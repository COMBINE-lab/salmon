//! Output writers: `quant.sf`, `cmd_info.json`, `lib_format_counts.json`, and
//! `aux_info/meta_info.json`. Functional parity with salmon's column/field
//! layout (the binary aux bias/FLD files are deferred to a later phase).

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;

use crate::{QuantOptions, QuantResult};

const SALMON_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Write the standard salmon output set into `opts.output_dir`.
pub fn write_outputs(opts: &QuantOptions, res: &QuantResult) -> Result<()> {
    let dir = &opts.output_dir;
    std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    std::fs::create_dir_all(dir.join("aux_info")).context("creating aux_info")?;

    write_quant_sf(&dir.join("quant.sf"), res)?;
    write_cmd_info(&dir.join("cmd_info.json"), opts)?;
    write_lib_counts(&dir.join("lib_format_counts.json"), opts, res)?;
    write_meta_info(&dir.join("aux_info").join("meta_info.json"), opts, res)?;
    Ok(())
}

/// `quant.sf`: `Name  Length  EffectiveLength  TPM  NumReads`.
fn write_quant_sf(path: &Path, res: &QuantResult) -> Result<()> {
    let mut w = std::io::BufWriter::new(
        std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?,
    );
    writeln!(w, "Name\tLength\tEffectiveLength\tTPM\tNumReads")?;
    for i in 0..res.names.len() {
        writeln!(
            w,
            "{}\t{}\t{:.3}\t{:.6}\t{:.3}",
            res.names[i], res.lengths[i], res.eff_lengths[i], res.tpm[i], res.counts[i]
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
        num_frags_with_concordant_consistent_mappings: res.num_mapped,
        num_frags_with_inconsistent_or_orphan_mappings: 0,
        strand_mapping_bias: 0.0,
    };
    write_json(path, &counts)
}

#[derive(Serialize)]
struct MetaInfo {
    salmon_version: String,
    samp_type: String,
    opt_type: String,
    num_libraries: usize,
    library_types: Vec<String>,
    frag_dist_length: usize,
    frag_length_mean: f64,
    seq_bias_correct: bool,
    gc_bias_correct: bool,
    num_bias_bins: usize,
    mapping_type: String,
    num_valid_targets: usize,
    num_decoy_targets: usize,
    num_eq_classes: usize,
    num_processed: u64,
    num_mapped: u64,
    percent_mapped: f64,
    num_decoy_fragments: u64,
    call: String,
}

fn write_meta_info(path: &Path, opts: &QuantOptions, res: &QuantResult) -> Result<()> {
    let percent_mapped = if res.num_processed > 0 {
        100.0 * res.num_mapped as f64 / res.num_processed as f64
    } else {
        0.0
    };
    let meta = MetaInfo {
        salmon_version: SALMON_VERSION.to_string(),
        samp_type: "none".to_string(),
        opt_type: if opts.em.use_vbem { "vb" } else { "em" }.to_string(),
        num_libraries: 1,
        library_types: vec![res.library_type.clone()],
        frag_dist_length: 1001,
        frag_length_mean: res.frag_len_mean,
        seq_bias_correct: false,
        gc_bias_correct: false,
        num_bias_bins: 0,
        mapping_type: if opts.sketch { "pseudo" } else { "mapping" }.to_string(),
        num_valid_targets: res.names.len(),
        num_decoy_targets: 0,
        num_eq_classes: res.num_eq_classes,
        num_processed: res.num_processed,
        num_mapped: res.num_mapped,
        percent_mapped,
        num_decoy_fragments: 0,
        call: "quant".to_string(),
    };
    write_json(path, &meta)
}

fn write_json<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    let json = serde_json::to_string_pretty(value)?;
    std::fs::write(path, json).with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}
