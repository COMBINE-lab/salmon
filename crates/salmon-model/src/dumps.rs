//! Aux-output dump helpers shared by reads-mode and alignment-mode quant.
//!
//! salmon writes several `aux_info` files in raw little-endian binary, gzipped:
//! `fld.gz` (fragment-length sample histogram), the legacy simple-count seq-bias
//! model (`observed_bias`/`observed_bias_3p`/`expected_bias`), and — under bias
//! correction — the per-model observed/expected tables. The Rust port computes
//! the `SBModel`/GC/positional models (dumped here as a documented Rust format:
//! gzip of raw LE arrays; positional files carry a small header) but not salmon's
//! legacy simple-count model, which is written as a documented stub for
//! file-presence parity. `libParams/flenDist.txt` is the text PMF.

use std::io::Write;
use std::path::Path;

/// Flattened observed/expected bias-model tables captured for the dump files.
/// Each group is empty when its correction was not enabled. Seq tables are the
/// [`SBModel`](crate::seqbias::SBModel) transition tables; GC the `cond×gc`
/// matrices; pos the per-length-class bin masses.
#[derive(Debug, Clone, Default)]
pub struct BiasDump {
    pub obs5_seq: Vec<f64>,
    pub obs3_seq: Vec<f64>,
    pub exp5_seq: Vec<f64>,
    pub exp3_seq: Vec<f64>,
    pub obs_gc: Vec<f64>,
    pub exp_gc: Vec<f64>,
    pub obs5_pos: Vec<Vec<f64>>,
    pub obs3_pos: Vec<Vec<f64>>,
    pub exp5_pos: Vec<Vec<f64>>,
    pub exp3_pos: Vec<Vec<f64>>,
}

/// gzip a raw byte buffer to `path` (level 6, matching salmon's aux dumps).
pub fn gz_write(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    let f = std::fs::File::create(path)?;
    let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::new(6));
    enc.write_all(bytes)?;
    enc.finish()?;
    Ok(())
}

/// gzip of a raw little-endian `f64` array.
pub fn write_f64_gz(path: &Path, vals: &[f64]) -> std::io::Result<()> {
    let mut b = Vec::with_capacity(vals.len() * 8);
    for v in vals {
        b.extend_from_slice(&v.to_le_bytes());
    }
    gz_write(path, &b)
}

/// gzip of a raw little-endian `i32` array.
pub fn write_i32_gz(path: &Path, vals: &[i32]) -> std::io::Result<()> {
    let mut b = Vec::with_capacity(vals.len() * 4);
    for v in vals {
        b.extend_from_slice(&v.to_le_bytes());
    }
    gz_write(path, &b)
}

/// gzip of a per-length-class positional model: header `[u32 num_models][u32
/// bins_per_model]` then the models' bin masses as `f64` LE, row-major.
pub fn write_pos_gz(path: &Path, models: &[Vec<f64>]) -> std::io::Result<()> {
    let bins = models.first().map(|m| m.len()).unwrap_or(0) as u32;
    let mut b = Vec::new();
    b.extend_from_slice(&(models.len() as u32).to_le_bytes());
    b.extend_from_slice(&bins.to_le_bytes());
    for m in models {
        for v in m {
            b.extend_from_slice(&v.to_le_bytes());
        }
    }
    gz_write(path, &b)
}

/// `aux_info/fld.gz`: per-length sample histogram. salmon draws 10,000 samples
/// from the log-PMF and writes the per-length `i32` counts; we write the
/// deterministic expected histogram `round(10000 * pmf[len])` (same type/layout).
pub fn write_fld_dump(path: &Path, pmf: &[f64]) -> std::io::Result<()> {
    const N_SAMPLES: f64 = 10000.0;
    let hist: Vec<i32> = pmf
        .iter()
        .map(|&p| (p * N_SAMPLES).round() as i32)
        .collect();
    write_i32_gz(path, &hist)
}

/// `libParams/flenDist.txt`: the normalized fragment-length PMF as a single line
/// of tab-separated scientific-notation values (salmon's format).
pub fn write_flen_dist(path: &Path, pmf: &[f64]) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut s = String::with_capacity(pmf.len() * 14);
    for (i, p) in pmf.iter().enumerate() {
        if i > 0 {
            s.push('\t');
        }
        s.push_str(&format!("{p:e}"));
    }
    s.push('\n');
    std::fs::write(path, s)
}

/// Write the `aux_info` bias dumps into `aux_dir`: documented stubs for salmon's
/// legacy simple-count model (not implemented in the port), plus the computed
/// seq/GC/pos observed+expected tables present in `dump`.
pub fn write_aux_bias_dumps(aux_dir: &Path, dump: &BiasDump) -> std::io::Result<()> {
    // Legacy simple-count seq-bias model (the port uses SBModel instead): stubs.
    write_i32_gz(&aux_dir.join("observed_bias.gz"), &[0])?;
    write_i32_gz(&aux_dir.join("observed_bias_3p.gz"), &[0])?;
    write_f64_gz(&aux_dir.join("expected_bias.gz"), &[1.0])?;

    if !dump.obs5_seq.is_empty() {
        write_f64_gz(&aux_dir.join("obs5_seq.gz"), &dump.obs5_seq)?;
        write_f64_gz(&aux_dir.join("obs3_seq.gz"), &dump.obs3_seq)?;
        write_f64_gz(&aux_dir.join("exp5_seq.gz"), &dump.exp5_seq)?;
        write_f64_gz(&aux_dir.join("exp3_seq.gz"), &dump.exp3_seq)?;
    }
    if !dump.obs_gc.is_empty() {
        write_f64_gz(&aux_dir.join("obs_gc.gz"), &dump.obs_gc)?;
        write_f64_gz(&aux_dir.join("exp_gc.gz"), &dump.exp_gc)?;
    }
    if !dump.obs5_pos.is_empty() {
        write_pos_gz(&aux_dir.join("obs5_pos.gz"), &dump.obs5_pos)?;
        write_pos_gz(&aux_dir.join("obs3_pos.gz"), &dump.obs3_pos)?;
        write_pos_gz(&aux_dir.join("exp5_pos.gz"), &dump.exp5_pos)?;
        write_pos_gz(&aux_dir.join("exp3_pos.gz"), &dump.exp3_pos)?;
    }
    Ok(())
}
