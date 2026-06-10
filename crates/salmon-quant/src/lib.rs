//! `salmon-quant`: the reads-mode quantification driver.
//!
//! Ties the pipeline together for FASTQ input: parse reads in parallel
//! (paraseq) → map each fragment ([`salmon_map`]) → accumulate weighted
//! equivalence classes ([`salmon_eqclass`]) and fragment lengths
//! ([`salmon_model`]) → compute effective lengths → estimate abundances by EM
//! ([`salmon_infer`]) → write `quant.sf` and the JSON metadata
//! ([`output`]).
//!
//! Selective alignment is the default; set [`QuantOptions::sketch`] for the
//! alignment-free pseudoalignment path.

mod output;
mod processor;

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;

use piscem_rs::io::fastx::{reader_with_batch_size, Collection, CollectionType};
use piscem_rs::mapping::hit_searcher::SkippingStrategy;

use salmon_core::{LibraryFormat, ReadType, Transcript};
use salmon_eqclass::EquivalenceClassBuilder;
use salmon_index::SalmonIndex;
use salmon_infer::{optimize, EmOptions};
use salmon_map::MapConfig;
use salmon_model::{FragmentLengthDistribution, LibraryTypeDetector, SBModel};

use processor::{QuantProcessor, Shared};

pub use output::write_outputs;

/// Options for a reads-mode quantification run.
#[derive(Debug, Clone)]
pub struct QuantOptions {
    /// directory of a salmon index (produced by `salmon-index`)
    pub index_dir: PathBuf,
    /// paired-end mate-1 files (parallel to `mates2`)
    pub mates1: Vec<PathBuf>,
    /// paired-end mate-2 files
    pub mates2: Vec<PathBuf>,
    /// single-end read files
    pub unmated: Vec<PathBuf>,
    /// output directory
    pub output_dir: PathBuf,
    /// library type string (e.g. `"IU"`, `"ISR"`, `"A"` for auto)
    pub lib_type: String,
    /// worker threads (`0` = all cores)
    pub num_threads: usize,
    /// use the alignment-free pseudoalignment path
    pub sketch: bool,
    /// selective-alignment mapping configuration
    pub map_config: MapConfig,
    /// EM/VBEM inference options
    pub em: EmOptions,
    /// range-factorization bins for equivalence classes (0 = disabled, salmon default 4)
    pub range_factorization_bins: u32,
    /// prior weight for strand-incompatible mappings; `0` drops them (salmon default)
    pub incompat_prior: f64,
    /// write `aux_info/eq_classes.txt` (naive transcript-set classes)
    pub dump_eq: bool,
    /// enable sequence-specific bias correction (`--seqBias`)
    pub seq_bias: bool,
}

impl QuantOptions {
    /// Options with salmon-like defaults for the given index and output dir.
    pub fn new(index_dir: PathBuf, output_dir: PathBuf) -> Self {
        Self {
            index_dir,
            mates1: Vec::new(),
            mates2: Vec::new(),
            unmated: Vec::new(),
            output_dir,
            lib_type: "A".to_string(),
            num_threads: 0,
            sketch: false,
            map_config: MapConfig::default(),
            em: EmOptions::default(),
            range_factorization_bins: 4,
            incompat_prior: 0.0,
            dump_eq: false,
            seq_bias: false,
        }
    }

    fn is_paired(&self) -> bool {
        !self.mates1.is_empty()
    }
}

/// Quantification results (also written to disk by [`write_outputs`]).
#[derive(Debug, Clone)]
pub struct QuantResult {
    pub names: Vec<String>,
    pub lengths: Vec<u32>,
    pub eff_lengths: Vec<f64>,
    pub tpm: Vec<f64>,
    /// estimated read counts per transcript
    pub counts: Vec<f64>,
    pub num_processed: u64,
    pub num_mapped: u64,
    pub num_eq_classes: usize,
    pub frag_len_mean: f64,
    /// the library type used: the detected format when `-l A`, else the
    /// user-specified one
    pub library_type: String,
}

/// Run quantification end-to-end, writing outputs and returning the results.
pub fn quantify(opts: &QuantOptions) -> Result<QuantResult> {
    let salmon = SalmonIndex::load(&opts.index_dir)
        .with_context(|| format!("loading index {}", opts.index_dir.display()))?;
    let num_refs = salmon.num_refs();

    let eq_builder = EquivalenceClassBuilder::new();
    // Gaussian-prior fragment length distribution (mean 250, sd 25); empirical
    // observations from concordant pairs refine it.
    let mut fld = FragmentLengthDistribution::new(1.0, 1000, 250.0, 25.0, 4, 0.5, 1);
    let num_processed = AtomicU64::new(0);
    let num_mapped = AtomicU64::new(0);
    let nthreads = if opts.num_threads == 0 {
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
    } else {
        opts.num_threads
    };

    // Auto-detect the library type when requested (`-l A`).
    let auto_detect = LibraryFormat::is_auto(&opts.lib_type);
    let read_type = if opts.is_paired() {
        ReadType::PairedEnd
    } else {
        ReadType::SingleEnd
    };
    let detector = auto_detect.then(|| LibraryTypeDetector::new(read_type));

    // Strand-compatibility filtering applies only with an explicit library
    // type; in auto mode the type is unknown during the pass.
    let expected_format = if auto_detect {
        None
    } else {
        Some(
            LibraryFormat::parse(&opts.lib_type)
                .with_context(|| format!("invalid library type '{}'", opts.lib_type))?,
        )
    };
    let ignore_incompat = opts.incompat_prior <= 0.0;

    // Shared observed sequence-bias models (merged per-thread) for `--seqBias`.
    let seqbias_obs = opts
        .seq_bias
        .then(|| std::sync::Mutex::new((SBModel::new(), SBModel::new())));

    // ---- parallel mapping pass (borrows the accumulators) -------------------
    {
        let shared = Shared {
            salmon: &salmon,
            eq: &eq_builder,
            fld: &fld,
            detector: detector.as_ref(),
            map_cfg: &opts.map_config,
            sketch: opts.sketch,
            skip: SkippingStrategy::Permissive,
            range_factorization_bins: opts.range_factorization_bins,
            expected_format,
            incompat_prior: opts.incompat_prior,
            ignore_incompat,
            collect_seqbias: opts.seq_bias,
            seqbias_obs: seqbias_obs.as_ref(),
            num_processed: &num_processed,
            num_mapped: &num_mapped,
        };
        let mut proc = QuantProcessor::new(shared);
        if opts.is_paired() {
            run_paired(&opts.mates1, &opts.mates2, &mut proc, nthreads)?;
        } else {
            run_single(&opts.unmated, &mut proc, nthreads)?;
        }
    }

    // Resolve the library type: the detected format (when auto), else the
    // user-specified string. Fall back to a sensible default if detection saw
    // no usable samples.
    let library_type = if let Some(det) = &detector {
        det.most_likely_type()
            .map(|f| f.canonical().to_string())
            .unwrap_or_else(|| {
                if opts.is_paired() { "IU" } else { "U" }.to_string()
            })
    } else {
        opts.lib_type.clone()
    };

    // ---- effective lengths from the (now frozen) FLD ------------------------
    fld.cache();
    let log_pmf = fld.log_pmf();
    let mut eff_lengths = vec![0f64; num_refs];
    for tid in 0..num_refs {
        let t = Transcript::new(tid as u32, salmon.ref_name(tid), salmon.ref_len(tid) as u32);
        eff_lengths[tid] = t.compute_log_effective_length(log_pmf, 0).exp();
    }

    // ---- inference ----------------------------------------------------------
    let mut collapsed = eq_builder.finish();
    collapsed.update_eff_lengths(&eff_lengths);
    let num_eq_classes = collapsed.len();

    if opts.dump_eq {
        dump_eq_classes(&opts.output_dir, &salmon, &collapsed)
            .context("writing eq_classes.txt")?;
    }
    let mut em = optimize(&collapsed, num_refs, &opts.em);

    // ---- sequence-bias correction: re-estimate effective lengths -----------
    // Build expected models from the initial abundance estimate, recompute
    // bias-corrected effective lengths, then re-run inference.
    if let Some(obs_mtx) = seqbias_obs {
        use rayon::prelude::*;
        let (mut obs_fw, mut obs_rc) = obs_mtx.into_inner().unwrap();
        obs_fw.normalize();
        obs_rc.normalize();
        let pmf_lin: Vec<f64> = log_pmf.iter().map(|lp| lp.exp()).collect();
        let (fld_cdf, fld_low, fld_high) = salmon_model::seqbias::fld_cdf_and_bounds(&pmf_lin);
        let (exp_fw, exp_rc) = salmon_model::build_expected(
            num_refs,
            |t| salmon.ref_seq(t as u32),
            &em.alphas,
            &eff_lengths,
            &fld_cdf,
        );

        // Both forward and RC observed models are populated by read orientation
        // (forward reads' 5' contexts -> fw, reverse reads' 5' contexts -> rc),
        // so both sides of the correction are informative for single- and
        // paired-end data alike.
        let corrected: Vec<f64> = (0..num_refs)
            .into_par_iter()
            .map(|tid| {
                if em.alphas[tid] < 1e-8 {
                    eff_lengths[tid]
                } else {
                    salmon_model::corrected_effective_length(
                        salmon.ref_seq(tid as u32),
                        &fld_cdf,
                        fld_low,
                        fld_high,
                        &obs_fw,
                        &exp_fw,
                        &obs_rc,
                        &exp_rc,
                        eff_lengths[tid],
                        salmon_model::seqbias::FLD_SAMP_STRIDE,
                    )
                }
            })
            .collect();
        eff_lengths = corrected;
        collapsed.update_eff_lengths(&eff_lengths);
        em = optimize(&collapsed, num_refs, &opts.em);
    }
    let counts = em.alphas;

    // ---- TPM ----------------------------------------------------------------
    let rates: Vec<f64> = (0..num_refs)
        .map(|i| if eff_lengths[i] > 0.0 { counts[i] / eff_lengths[i] } else { 0.0 })
        .collect();
    let rate_sum: f64 = rates.iter().sum();
    let tpm: Vec<f64> = rates
        .iter()
        .map(|r| if rate_sum > 0.0 { r / rate_sum * 1e6 } else { 0.0 })
        .collect();

    let result = QuantResult {
        names: (0..num_refs).map(|i| salmon.ref_name(i).to_string()).collect(),
        lengths: (0..num_refs).map(|i| salmon.ref_len(i) as u32).collect(),
        eff_lengths,
        tpm,
        counts,
        num_processed: num_processed.load(Ordering::Relaxed),
        num_mapped: num_mapped.load(Ordering::Relaxed),
        num_eq_classes,
        frag_len_mean: fld.mean(),
        library_type,
    };

    write_outputs(opts, &result)?;
    Ok(result)
}

/// Write the naive equivalence classes (collapsing any range-factorized
/// sub-classes back to their transcript set) for comparison/diagnostics.
fn dump_eq_classes(
    dir: &Path,
    salmon: &SalmonIndex,
    collapsed: &salmon_eqclass::CollapsedEqClasses,
) -> Result<()> {
    use std::collections::BTreeMap;
    // Accumulate count and summed combined-weights per transcript set, so the
    // dump shows the average per-fragment conditional probability per transcript.
    let mut naive: BTreeMap<Vec<u32>, (u64, Vec<f64>)> = BTreeMap::new();
    for (group, value) in &collapsed.classes {
        let entry = naive
            .entry(group.txps.clone())
            .or_insert_with(|| (0, vec![0.0; group.txps.len()]));
        entry.0 += value.count;
        for (i, w) in value.combined_weights.iter().enumerate() {
            entry.1[i] += w * value.count as f64;
        }
    }
    std::fs::create_dir_all(dir.join("aux_info"))?;
    let mut w = std::io::BufWriter::new(std::fs::File::create(dir.join("aux_info").join("eq_classes.txt"))?);
    use std::io::Write as _;
    for (txps, (count, wsum)) in &naive {
        let names: Vec<&str> = txps.iter().map(|&t| salmon.ref_name(t as usize)).collect();
        let avg: Vec<String> = wsum
            .iter()
            .map(|s| format!("{:.4}", s / *count as f64))
            .collect();
        writeln!(w, "{}\t{}\t{}\t{}", txps.len(), names.join(","), avg.join(","), count)?;
    }
    w.flush()?;
    Ok(())
}

/// Open a (optionally gzipped) read file as a boxed reader.
fn open_reader(path: &Path) -> Result<Box<dyn std::io::Read + Send>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Ok(Box::new(MultiGzDecoder::new(f)))
    } else {
        Ok(Box::new(f))
    }
}

fn run_single(paths: &[PathBuf], proc: &mut QuantProcessor, nthreads: usize) -> Result<()> {
    let mut readers = Vec::with_capacity(paths.len());
    for p in paths {
        readers.push(
            reader_with_batch_size(open_reader(p)?)
                .map_err(|e| anyhow::anyhow!("opening {}: {e}", p.display()))?,
        );
    }
    let collection = Collection::new(readers, CollectionType::Single)
        .map_err(|e| anyhow::anyhow!("building read collection: {e}"))?;
    collection
        .process_parallel(proc, nthreads, None)
        .map_err(|e| anyhow::anyhow!("mapping failed: {e}"))?;
    Ok(())
}

fn run_paired(
    mates1: &[PathBuf],
    mates2: &[PathBuf],
    proc: &mut QuantProcessor,
    nthreads: usize,
) -> Result<()> {
    anyhow::ensure!(
        mates1.len() == mates2.len(),
        "mates1 and mates2 must have the same number of files"
    );
    let mut readers = Vec::with_capacity(mates1.len() * 2);
    for (a, b) in mates1.iter().zip(mates2.iter()) {
        readers.push(
            reader_with_batch_size(open_reader(a)?)
                .map_err(|e| anyhow::anyhow!("opening {}: {e}", a.display()))?,
        );
        readers.push(
            reader_with_batch_size(open_reader(b)?)
                .map_err(|e| anyhow::anyhow!("opening {}: {e}", b.display()))?,
        );
    }
    let collection = Collection::new(readers, CollectionType::Paired)
        .map_err(|e| anyhow::anyhow!("building paired read collection: {e}"))?;
    collection
        .process_parallel_paired(proc, nthreads, None)
        .map_err(|e| anyhow::anyhow!("mapping failed: {e}"))?;
    Ok(())
}
