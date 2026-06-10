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

use salmon_core::{LibraryFormat, ReadType};
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
    /// enable fragment-GC bias correction (`--gcBias`)
    pub gc_bias: bool,
    /// enable positional bias correction (`--posBias`)
    pub pos_bias: bool,
    /// number of bootstrap replicates (`--numBootstraps`); 0 = off
    pub num_bootstraps: u32,
    /// number of Gibbs posterior samples (`--numGibbsSamples`); 0 = off
    pub num_gibbs_samples: u32,
    /// Gibbs thinning factor (`--thinningFactor`, salmon default 16)
    pub thinning_factor: u32,
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
            gc_bias: false,
            pos_bias: false,
            num_bootstraps: 0,
            num_gibbs_samples: 0,
            thinning_factor: 16,
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
    /// posterior samples (bootstrap or Gibbs), one abundance vector each; empty
    /// when neither was requested
    pub bootstraps: Vec<Vec<f64>>,
    /// per-transcript (unique, ambiguous) fragment counts for `ambig_info.tsv`
    pub ambig: (Vec<u32>, Vec<u32>),
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

    // For `--gcBias`: per-transcript cumulative G+C prefix sums (salmon's
    // `Transcript::GCCount_`) and the shared observed fragment-GC model.
    let gc_prefix: Option<Vec<Vec<u32>>> = opts.gc_bias.then(|| {
        (0..num_refs)
            .map(|tid| salmon_model::gc_prefix(salmon.ref_seq(tid as u32)))
            .collect()
    });
    let gcbias_obs = opts
        .gc_bias
        .then(|| std::sync::Mutex::new(salmon_model::GcFragModel::default_model()));

    // For `--posBias`: transcript length-class quantiles + per-transcript class,
    // and the shared observed (5', 3') positional-bias models per length class.
    let length_quantiles: Option<Vec<u32>> = opts.pos_bias.then(|| {
        let lens: Vec<u32> = (0..num_refs).map(|t| salmon.ref_len(t) as u32).collect();
        salmon_model::compute_length_quantiles(&lens, salmon_model::NUM_LENGTH_CLASSES)
    });
    let length_class: Option<Vec<usize>> = length_quantiles.as_ref().map(|q| {
        (0..num_refs)
            .map(|t| salmon_model::length_class_index(q, salmon.ref_len(t) as u32))
            .collect()
    });
    let posbias_obs = opts.pos_bias.then(|| {
        let mk = || {
            (0..salmon_model::NUM_LENGTH_CLASSES)
                .map(|_| salmon_model::SimplePosBias::default())
                .collect::<Vec<_>>()
        };
        std::sync::Mutex::new((mk(), mk()))
    });

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
            collect_gcbias: opts.gc_bias,
            gc_prefix: gc_prefix.as_deref(),
            gcbias_obs: gcbias_obs.as_ref(),
            collect_posbias: opts.pos_bias,
            length_class: length_class.as_deref(),
            posbias_obs: posbias_obs.as_ref(),
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
    // salmon's base effective length is `refLen − E[L | L ≤ refLen]`
    // (`computeSmoothedEffectiveLengths` via `correctionFactorsFromMass`), NOT the
    // truncated-PMF `Σ pmf(l)·(refLen−l+1)` estimate, which falls back to the raw
    // reference length for any transcript shorter than the FLD mean.
    fld.cache();
    let log_pmf = fld.log_pmf();
    let cond_means = fld.conditional_means();
    let mut eff_lengths = vec![0f64; num_refs];
    for tid in 0..num_refs {
        eff_lengths[tid] =
            salmon_model::smoothed_effective_length(&cond_means, salmon.ref_len(tid) as usize);
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

    // ---- bias correction: re-estimate effective lengths --------------------
    // Build the expected (background) bias models from the initial abundance
    // estimate, recompute bias-corrected effective lengths, then re-run
    // inference. Composes `--seqBias`, `--gcBias`, and `--posBias` in any
    // combination via the unified convolution (salmon's updateEffectiveLengths).
    if opts.seq_bias || opts.gc_bias || opts.pos_bias {
        use rayon::prelude::*;
        let pmf_lin: Vec<f64> = log_pmf.iter().map(|lp| lp.exp()).collect();
        let (fld_cdf, fld_low, fld_high) = salmon_model::seqbias::fld_cdf_and_bounds(&pmf_lin);
        // K excludes the leading sequence context only when seq-correcting.
        let k = if opts.seq_bias {
            salmon_model::seqbias::CONTEXT_LENGTH
        } else {
            1
        };

        // Sequence-bias observed + expected models.
        let seq = seqbias_obs.map(|m| {
            let (mut obs_fw, mut obs_rc) = m.into_inner().unwrap();
            obs_fw.normalize();
            obs_rc.normalize();
            let (exp_fw, exp_rc) = salmon_model::build_expected(
                num_refs,
                |t| salmon.ref_seq(t as u32),
                &em.alphas,
                &eff_lengths,
                &fld_cdf,
            );
            (obs_fw, obs_rc, exp_fw, exp_rc)
        });

        // Fragment-GC observed + expected models -> clamped ratio model.
        let prefixes = gc_prefix.as_ref();
        let gc_ratio_model = gcbias_obs.map(|m| {
            let mut obs = m.into_inner().unwrap();
            let mut exp = salmon_model::build_expected_gc(
                num_refs,
                |t| salmon.ref_seq(t as u32),
                |t| prefixes.unwrap()[t].as_slice(),
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
            salmon_model::gc_ratio(&mut obs, &mut exp, salmon_model::gcbias::GC_MAX_RATIO)
        });

        // Positional-bias observed (finalized) + expected (5'/3' per length class).
        let pos_models = posbias_obs.map(|m| {
            let (mut obs_fw, mut obs_rc) = m.into_inner().unwrap();
            for x in obs_fw.iter_mut().chain(obs_rc.iter_mut()) {
                x.finalize();
            }
            let (exp_fw, exp_rc) = salmon_model::build_expected_pos(
                num_refs,
                |t| salmon.ref_len(t) as usize,
                &em.alphas,
                &eff_lengths,
                &fld_cdf,
                length_quantiles.as_ref().unwrap(),
                k,
            );
            // Debug: dump finalized pos models for parity comparison with salmon.
            if std::env::var("SALMON_DEBUG_POS").is_ok() {
                use std::io::Write;
                let mut f = std::fs::File::create(
                    opts.output_dir.join("rust_pos_models.txt"),
                )
                .unwrap();
                for (name, models) in [
                    ("obs5", &obs_fw),
                    ("obs3", &obs_rc),
                    ("exp5", &exp_fw),
                    ("exp3", &exp_rc),
                ] {
                    for (lc, m) in models.iter().enumerate() {
                        write!(f, "{name} {lc}").unwrap();
                        for v in m.masses() {
                            write!(f, " {v:.6}").unwrap();
                        }
                        writeln!(f).unwrap();
                    }
                }
            }
            (obs_fw, obs_rc, exp_fw, exp_rc)
        });

        let corrected: Vec<f64> = (0..num_refs)
            .into_par_iter()
            .map(|tid| {
                if em.alphas[tid] < 1e-8 {
                    return eff_lengths[tid];
                }
                let s = salmon.ref_seq(tid as u32);
                let ref_len = s.len();
                // Per-position 5'/3' positional-bias factors (obs/exp projected).
                let pos_vecs: Option<(Vec<f64>, Vec<f64>)> =
                    pos_models.as_ref().map(|(ofw, orc, efw, erc)| {
                        let lc = length_class.as_ref().unwrap()[tid];
                        let (mut o5, mut e5) = (vec![0.0; ref_len], vec![0.0; ref_len]);
                        let (mut o3, mut e3) = (vec![0.0; ref_len], vec![0.0; ref_len]);
                        ofw[lc].project_weights(&mut o5);
                        efw[lc].project_weights(&mut e5);
                        orc[lc].project_weights(&mut o3);
                        erc[lc].project_weights(&mut e3);
                        let pf: Vec<f64> = o5.iter().zip(&e5).map(|(o, e)| o / e).collect();
                        let pr: Vec<f64> = o3.iter().zip(&e3).map(|(o, e)| o / e).collect();
                        (pf, pr)
                    });
                let bias = salmon_model::BiasInputs {
                    seq: seq.as_ref().map(|(of, or, ef, er)| (of, ef, or, er)),
                    gc: gc_ratio_model
                        .as_ref()
                        .map(|g| (g, prefixes.unwrap()[tid].as_slice())),
                    pos: pos_vecs.as_ref().map(|(pf, pr)| (pf.as_slice(), pr.as_slice())),
                };
                salmon_model::corrected_effective_length_full(
                    s,
                    &fld_cdf,
                    fld_low,
                    fld_high,
                    &bias,
                    eff_lengths[tid],
                    salmon_model::GC_SAMP_STRIDE,
                )
            })
            .collect();
        eff_lengths = corrected;
        collapsed.update_eff_lengths(&eff_lengths);
        em = optimize(&collapsed, num_refs, &opts.em);
    }
    let counts = em.alphas;

    // ---- posterior uncertainty (bootstrap / Gibbs) + ambiguity --------------
    // The packed CSR layout (piscem-infer style) makes these parallel-friendly.
    let packed = salmon_infer::PackedEqClasses::from_collapsed(&collapsed, num_refs);
    let ambig = salmon_infer::ambiguity_counts(&packed);
    let num_mapped_frags = num_mapped.load(Ordering::Relaxed);
    let bootstraps: Vec<Vec<f64>> = if opts.num_bootstraps > 0 {
        salmon_infer::bootstrap(
            &packed,
            &opts.em,
            opts.num_bootstraps,
            num_mapped_frags,
            true, // useScaledCounts (selective-alignment, no orphan-only quasi)
            0x5A13_0000,
        )
    } else if opts.num_gibbs_samples > 0 {
        // Gibbs prior follows the main optimizer (salmon): with VBEM and a
        // per-transcript prior it is `max(1.0, vbPrior)`; with plain EM it is
        // 1e-3 per transcript. The rust VBEM uses a constant per-transcript prior.
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
        salmon_infer::gibbs_sample(&packed, &eff_lengths, &counts, &gopts, num_mapped_frags, 0x6217_0000)
    } else {
        Vec::new()
    };

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
        bootstraps,
        ambig,
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
