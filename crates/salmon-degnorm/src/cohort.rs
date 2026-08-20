//! Fitting the degradation model across a cohort of samples, and the tables
//! that come out of it.
//!
//! [`crate::nmf`] handles one transcript's coverage matrix. This module is the
//! bookkeeping around it: line up several samples' coverage dumps, decide which
//! transcripts are worth fitting, run the fits in parallel, and write the
//! degradation-index and adjusted-count matrices.
//!
//! # Why samples must come from the same index
//!
//! Position bin `b` of transcript `t` only means the same thing in two samples
//! if `t` has the same length in both, and rows only line up if the transcript
//! order matches. Both hold automatically when the samples were quantified
//! against one index, and nothing salvages the fit when they were not — so a
//! mismatch is a hard error naming the first target that disagrees, not a
//! best-effort intersection.
//!
//! # Sequencing depth
//!
//! A deeper sample has higher coverage everywhere, which the rank-one model
//! would happily absorb into its abundance term `a`. It is still worth dividing
//! it out first: the fit minimises a squared error over the whole matrix, so
//! without scaling the deepest sample dominates the shared envelope. Coverage
//! is therefore scaled by the run's mapped-fragment count relative to the
//! cohort mean before anything else happens.

use std::collections::HashMap;
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;

use crate::coverage::CoverageProfiles;
use crate::nmf::{self, FitOptions};

/// One sample's inputs: the coverage dump, and optionally the `quant.sf` whose
/// counts get adjusted.
#[derive(Debug, Clone)]
pub struct Sample {
    /// Display name used for the output matrix columns.
    pub name: String,
    /// `aux_info/coverage.gz` from that sample's `salmon quant`.
    pub coverage: PathBuf,
    /// That sample's `quant.sf`. Without it the run still reports degradation
    /// indices, just no adjusted counts.
    pub quant: Option<PathBuf>,
}

/// Cohort-level knobs.
#[derive(Debug, Clone)]
pub struct CohortOptions {
    /// Per-transcript model settings.
    pub fit: FitOptions,
    /// Skip transcripts shorter than this. Short transcripts have too few
    /// distinct positions for a coverage *shape* to mean anything, and they are
    /// where spurious indices come from.
    pub min_length: u32,
    /// Skip transcripts whose mean depth, averaged over samples and positions,
    /// is below this. A transcript nobody sequenced has no recoverable shape.
    pub min_mean_coverage: f64,
    /// Ceiling on the index used to inflate counts. The reported DI is
    /// untouched; this only bounds `1 / (1 − DI)` so a nearly-empty transcript
    /// cannot turn a handful of reads into a huge adjusted count.
    pub max_adjust_di: f64,
}

impl Default for CohortOptions {
    fn default() -> Self {
        Self {
            fit: FitOptions::default(),
            min_length: 200,
            min_mean_coverage: 0.1,
            max_adjust_di: 0.9,
        }
    }
}

/// What a cohort run produced, in memory. [`write_tables`] puts it on disk.
#[derive(Debug, Clone)]
pub struct CohortResult {
    pub sample_names: Vec<String>,
    pub target_names: Vec<String>,
    /// Row-major `targets × samples` degradation indices. Targets that were not
    /// fitted are all-zero rows.
    pub di: Vec<f64>,
    /// Whether each target was actually fitted (rather than filtered out).
    pub fitted: Vec<bool>,
    /// Raw `NumReads` per target per sample, row-major; empty when no `quant.sf`
    /// was supplied.
    pub raw_counts: Vec<f64>,
    /// `raw / (1 − min(DI, max_adjust_di))`; empty when `raw_counts` is.
    pub adjusted_counts: Vec<f64>,
    /// Per-sample depth scale factor applied to the coverage before fitting.
    pub depth_factors: Vec<f64>,
    /// Mapped fragments per sample, as recorded in each coverage dump.
    pub num_mapped: Vec<u64>,
}

impl CohortResult {
    pub fn num_samples(&self) -> usize {
        self.sample_names.len()
    }

    pub fn num_targets(&self) -> usize {
        self.target_names.len()
    }

    /// Mean degradation index per sample over the fitted targets. The headline
    /// number: "how degraded is this library".
    pub fn mean_di(&self) -> Vec<f64> {
        let m = self.num_samples();
        let mut sums = vec![0.0; m];
        let mut n = 0usize;
        for (t, &ok) in self.fitted.iter().enumerate() {
            if !ok {
                continue;
            }
            n += 1;
            for s in 0..m {
                sums[s] += self.di[t * m + s];
            }
        }
        if n == 0 {
            return vec![0.0; m];
        }
        sums.iter().map(|x| x / n as f64).collect()
    }
}

/// Load every sample, fit every eligible transcript, and return the cohort
/// tables.
pub fn run(samples: &[Sample], opts: &CohortOptions) -> Result<CohortResult> {
    if samples.len() < 2 {
        bail!(
            "degradation normalization needs at least 2 samples (got {}); \
             a single coverage curve is rank one by construction and its \
             degradation index is zero by definition",
            samples.len()
        );
    }

    let profiles: Vec<CoverageProfiles> = samples
        .iter()
        .map(|s| {
            CoverageProfiles::read(&s.coverage)
                .with_context(|| format!("reading coverage for sample '{}'", s.name))
        })
        .collect::<Result<_>>()?;
    check_comparable(samples, &profiles)?;

    let m = samples.len();
    let n_bins = profiles[0].num_bins;
    let num_targets = profiles[0].num_targets();
    let target_names = profiles[0].names.clone();
    let ref_lens = profiles[0].ref_lens.clone();

    // Depth scaling; see the module docs. A dump that recorded no mapped
    // fragments (an empty or truncated run) gets a factor of 1 rather than an
    // infinity.
    let num_mapped: Vec<u64> = profiles.iter().map(|p| p.num_mapped).collect();
    let mean_mapped = num_mapped.iter().sum::<u64>() as f64 / m as f64;
    let depth_factors: Vec<f64> = num_mapped
        .iter()
        .map(|&d| {
            if d == 0 || mean_mapped <= 0.0 {
                1.0
            } else {
                mean_mapped / d as f64
            }
        })
        .collect();

    // One independent fit per transcript. `map_init` keeps the `m × n_bins`
    // scratch matrix alive per rayon worker instead of reallocating it for
    // every one of a few hundred thousand transcripts.
    let per_target: Vec<(bool, Vec<f64>)> = (0..num_targets)
        .into_par_iter()
        .map_init(
            || vec![0.0f64; m * n_bins],
            |k, t| {
                for (i, p) in profiles.iter().enumerate() {
                    let row = p.row(t);
                    for j in 0..n_bins {
                        k[i * n_bins + j] = f64::from(row[j]) * depth_factors[i];
                    }
                }
                let mean_cov = k.iter().sum::<f64>() / (m * n_bins) as f64;
                if ref_lens[t] < opts.min_length || mean_cov < opts.min_mean_coverage {
                    return (false, vec![0.0; m]);
                }
                (true, nmf::fit(k, m, n_bins, &opts.fit).di)
            },
        )
        .collect();

    let mut di = Vec::with_capacity(num_targets * m);
    let mut fitted = Vec::with_capacity(num_targets);
    for (ok, row) in per_target {
        fitted.push(ok);
        di.extend_from_slice(&row);
    }

    // Counts, when the caller pointed at the quant files.
    let (raw_counts, adjusted_counts) = if samples.iter().all(|s| s.quant.is_some()) {
        let mut raw = vec![0.0f64; num_targets * m];
        for (i, s) in samples.iter().enumerate() {
            let path = s.quant.as_ref().expect("checked above");
            let counts = read_quant_counts(path)
                .with_context(|| format!("reading counts for sample '{}'", s.name))?;
            for (t, name) in target_names.iter().enumerate() {
                if let Some(&c) = counts.get(name) {
                    raw[t * m + i] = c;
                }
            }
        }
        let adjusted = raw
            .iter()
            .zip(&di)
            .map(|(&c, &d)| c / (1.0 - d.min(opts.max_adjust_di)))
            .collect();
        (raw, adjusted)
    } else {
        (Vec::new(), Vec::new())
    };

    Ok(CohortResult {
        sample_names: samples.iter().map(|s| s.name.clone()).collect(),
        target_names,
        di,
        fitted,
        raw_counts,
        adjusted_counts,
        depth_factors,
        num_mapped,
    })
}

/// Reject cohorts whose coverage dumps cannot be compared position by position.
fn check_comparable(samples: &[Sample], profiles: &[CoverageProfiles]) -> Result<()> {
    let first = &profiles[0];
    for (s, p) in samples.iter().zip(profiles).skip(1) {
        if p.num_bins != first.num_bins {
            bail!(
                "sample '{}' has {} coverage bins but '{}' has {}; re-run quant with the \
                 same --degCovBins for every sample",
                s.name,
                p.num_bins,
                samples[0].name,
                first.num_bins
            );
        }
        if p.num_targets() != first.num_targets() {
            bail!(
                "sample '{}' covers {} transcripts but '{}' covers {}; all samples must be \
                 quantified against the same index",
                s.name,
                p.num_targets(),
                samples[0].name,
                first.num_targets()
            );
        }
        // Name *and* length: same names with different lengths means different
        // annotation versions, which silently misaligns every bin.
        for (t, ((n1, l1), (n2, l2))) in first
            .names
            .iter()
            .zip(&first.ref_lens)
            .zip(p.names.iter().zip(&p.ref_lens))
            .enumerate()
        {
            if n1 != n2 || l1 != l2 {
                bail!(
                    "sample '{}' disagrees with '{}' at transcript {t}: {n2} (len {l2}) vs \
                     {n1} (len {l1}); all samples must be quantified against the same index",
                    s.name,
                    samples[0].name
                );
            }
        }
    }
    Ok(())
}

/// `Name -> NumReads` from a `quant.sf`.
fn read_quant_counts(path: &Path) -> Result<HashMap<String, f64>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut out = HashMap::new();
    for (lineno, line) in std::io::BufReader::new(f).lines().enumerate() {
        let line = line?;
        if lineno == 0 || line.trim().is_empty() {
            continue; // header / trailing blank
        }
        // salmon's fixed layout: Name, Length, EffectiveLength, TPM, NumReads.
        let toks: Vec<&str> = line.split('\t').collect();
        if toks.len() < 5 {
            continue;
        }
        let v: f64 = toks[4]
            .parse()
            .with_context(|| format!("{}:{}: unparseable NumReads", path.display(), lineno + 1))?;
        out.insert(toks[0].to_string(), v);
    }
    Ok(out)
}

/// Write the cohort tables into `dir`: `degradation_index.tsv`, and — when
/// counts were available — `counts_raw.tsv` and `counts_adjusted.tsv`, plus a
/// `degnorm_summary.json` with the run's parameters and headline numbers.
pub fn write_tables(dir: &Path, res: &CohortResult, opts: &CohortOptions) -> Result<()> {
    std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    write_matrix(
        &dir.join("degradation_index.tsv"),
        res,
        &res.di,
        // A fitted DI of 0 and a skipped transcript are different statements,
        // and a downstream `1/(1-DI)` would silently treat the second as the
        // first. `NA` forces the distinction to be handled.
        true,
    )?;
    if !res.raw_counts.is_empty() {
        write_matrix(&dir.join("counts_raw.tsv"), res, &res.raw_counts, false)?;
        write_matrix(
            &dir.join("counts_adjusted.tsv"),
            res,
            &res.adjusted_counts,
            false,
        )?;
    }

    let summary = serde_json::json!({
        "samples": res.sample_names,
        "num_targets": res.num_targets(),
        "num_targets_fitted": res.fitted.iter().filter(|&&f| f).count(),
        "num_mapped": res.num_mapped,
        "depth_factors": res.depth_factors,
        "mean_degradation_index": res.mean_di(),
        "parameters": {
            "min_length": opts.min_length,
            "min_mean_coverage": opts.min_mean_coverage,
            "max_adjust_di": opts.max_adjust_di,
            "max_iter": opts.fit.max_iter,
            "tol": opts.fit.tol,
            "als_iter": opts.fit.als_iter,
            "mask_frac": opts.fit.mask_frac,
        },
    });
    let path = dir.join("degnorm_summary.json");
    std::fs::write(&path, serde_json::to_vec_pretty(&summary)?)
        .with_context(|| format!("writing {}", path.display()))?;
    Ok(())
}

/// One `Name`-plus-a-column-per-sample TSV. `na_for_unfitted` writes `NA` for
/// targets the model skipped instead of their placeholder zeros.
fn write_matrix(
    path: &Path,
    res: &CohortResult,
    values: &[f64],
    na_for_unfitted: bool,
) -> Result<()> {
    let f = std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let mut w = BufWriter::new(f);
    write!(w, "Name")?;
    for s in &res.sample_names {
        write!(w, "\t{s}")?;
    }
    writeln!(w)?;
    let m = res.num_samples();
    for (t, name) in res.target_names.iter().enumerate() {
        write!(w, "{name}")?;
        for s in 0..m {
            if na_for_unfitted && !res.fitted[t] {
                write!(w, "\tNA")?;
            } else {
                write!(w, "\t{:.6}", values[t * m + s])?;
            }
        }
        writeln!(w)?;
    }
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coverage::CoverageAccumulator;

    /// Write a coverage dump whose single transcript has the given per-bin
    /// depths, and return the sample pointing at it.
    fn sample(dir: &Path, name: &str, rows: &[&[f32]], lens: &[u32], num_mapped: u64) -> Sample {
        let n_bins = rows[0].len();
        let acc = CoverageAccumulator::new(lens.to_vec(), n_bins);
        // Build the profile directly rather than through `add_fragment`: this
        // is testing the cohort layer, not the binning.
        let mut profiles = acc.finish(
            (0..lens.len()).map(|i| format!("tx{i}")).collect(),
            num_mapped,
        );
        profiles.values = rows.iter().flat_map(|r| r.iter().copied()).collect();
        let path = dir.join(format!("{name}.cov.gz"));
        profiles.write(&path).unwrap();
        Sample {
            name: name.to_string(),
            coverage: path,
            quant: None,
        }
    }

    fn flat(v: f32, n: usize) -> Vec<f32> {
        vec![v; n]
    }

    #[test]
    fn a_truncated_sample_is_flagged_and_the_others_are_not() {
        let dir = tempfile::tempdir().unwrap();
        let n = 20;
        let intact = flat(10.0, n);
        let mut truncated = flat(10.0, n);
        truncated[..n / 2].fill(0.0);
        let samples = vec![
            sample(dir.path(), "a", &[&intact], &[1000], 1_000_000),
            sample(dir.path(), "b", &[&truncated], &[1000], 1_000_000),
            sample(dir.path(), "c", &[&intact], &[1000], 1_000_000),
        ];
        let res = run(&samples, &CohortOptions::default()).unwrap();
        assert!(res.fitted[0]);
        assert!(res.di[0] < 0.05, "intact: {}", res.di[0]);
        assert!((res.di[1] - 0.5).abs() < 0.05, "truncated: {}", res.di[1]);
        assert!(res.di[2] < 0.05, "intact: {}", res.di[2]);
    }

    #[test]
    fn short_and_shallow_transcripts_are_skipped_rather_than_fitted() {
        let dir = tempfile::tempdir().unwrap();
        let n = 20;
        let deep = flat(10.0, n);
        let shallow = flat(0.001, n);
        // tx0: long enough but nearly empty. tx1: well covered but too short.
        let rows: Vec<&[f32]> = vec![&shallow, &deep];
        let samples = vec![
            sample(dir.path(), "a", &rows, &[1000, 50], 1_000_000),
            sample(dir.path(), "b", &rows, &[1000, 50], 1_000_000),
        ];
        let res = run(&samples, &CohortOptions::default()).unwrap();
        assert_eq!(res.fitted, vec![false, false]);
    }

    #[test]
    fn one_sample_is_refused_with_an_explanation() {
        let dir = tempfile::tempdir().unwrap();
        let s = sample(dir.path(), "a", &[&flat(5.0, 8)], &[1000], 100);
        let err = run(std::slice::from_ref(&s), &CohortOptions::default()).unwrap_err();
        assert!(
            err.to_string().contains("at least 2 samples"),
            "unhelpful error: {err}"
        );
    }

    #[test]
    fn mismatched_indices_are_refused() {
        let dir = tempfile::tempdir().unwrap();
        let a = sample(dir.path(), "a", &[&flat(5.0, 8)], &[1000], 100);
        let b = sample(dir.path(), "b", &[&flat(5.0, 8)], &[900], 100);
        let err = run(&[a, b], &CohortOptions::default()).unwrap_err();
        assert!(
            err.to_string().contains("same index"),
            "unhelpful error: {err}"
        );
    }

    #[test]
    fn adjusted_counts_undo_the_measured_loss() {
        let dir = tempfile::tempdir().unwrap();
        let n = 20;
        let intact = flat(10.0, n);
        let mut truncated = flat(10.0, n);
        truncated[..n / 2].fill(0.0);
        let mut samples = vec![
            sample(dir.path(), "a", &[&intact], &[1000], 1_000_000),
            sample(dir.path(), "b", &[&truncated], &[1000], 1_000_000),
        ];
        // A degraded sample reports half the reads of the intact one; the
        // adjustment should bring the two back into line.
        for (s, reads) in samples.iter_mut().zip([1000.0, 500.0]) {
            let q = dir.path().join(format!("{}.quant.sf", s.name));
            std::fs::write(
                &q,
                format!(
                    "Name\tLength\tEffectiveLength\tTPM\tNumReads\ntx0\t1000\t800\t1.0\t{reads}\n"
                ),
            )
            .unwrap();
            s.quant = Some(q);
        }
        let res = run(&samples, &CohortOptions::default()).unwrap();
        assert_eq!(res.raw_counts, vec![1000.0, 500.0]);
        let ratio = res.adjusted_counts[1] / res.adjusted_counts[0];
        assert!(
            (ratio - 1.0).abs() < 0.15,
            "adjustment left a {ratio:.2}x gap: {:?}",
            res.adjusted_counts
        );
    }

    #[test]
    fn tables_land_on_disk_with_na_for_skipped_targets() {
        let dir = tempfile::tempdir().unwrap();
        let rows: Vec<&[f32]> = vec![&[0.0; 12], &[8.0; 12]];
        let samples = vec![
            sample(dir.path(), "a", &rows, &[1000, 1000], 500),
            sample(dir.path(), "b", &rows, &[1000, 1000], 500),
        ];
        let opts = CohortOptions::default();
        let res = run(&samples, &opts).unwrap();
        let out = dir.path().join("out");
        write_tables(&out, &res, &opts).unwrap();
        let di = std::fs::read_to_string(out.join("degradation_index.tsv")).unwrap();
        let lines: Vec<&str> = di.lines().collect();
        assert_eq!(lines[0], "Name\ta\tb");
        assert_eq!(lines[1], "tx0\tNA\tNA", "empty transcript should be NA");
        assert!(lines[2].starts_with("tx1\t0."));
        let summary = std::fs::read_to_string(out.join("degnorm_summary.json")).unwrap();
        assert!(summary.contains("\"num_targets_fitted\": 1"));
        // No quant files were supplied, so no count tables.
        assert!(!out.join("counts_raw.tsv").exists());
    }
}
