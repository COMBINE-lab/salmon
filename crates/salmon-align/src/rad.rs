//! Parallel RAD-format input quantification (`salmon quant --rad`).
//!
//! Quantifies from a RAD file of mappings instead of FASTQ/BAM, taking a path
//! analogous to alignment mode: each fragment's placements are turned into a
//! weighted equivalence class, then the shared inference (EM/VBEM + bootstrap/
//! Gibbs) and output (`write_outputs`) run exactly as for BAM input.
//!
//! Accepts:
//!   * **piscem `map-bulk`** RAD (`bulk_with_pos`) — sketch-style, no scores, all
//!     placements equally-best (uniform eq-class weights);
//!   * **salmon-generated** RAD — sketch (uniform) or selective-alignment (per-hit
//!     alignment scores → soft-weighted like internal selective alignments).
//!
//! The RAD chunks are already collated by read, so we read them in parallel with
//! libradicl's [`ParallelRadReader`] (a producer feeding a bounded work queue)
//! and a pool of worker threads, mirroring piscem-infer's `process_parallel`.
//!
//! v1 scope: no sequence/GC/positional bias correction (RAD carries no read
//! sequences; bias would need `-t` and the expected-model machinery) and no decoy
//! filtering for piscem RAD (salmon RAD already excludes decoys at write time).
//! These are tracked follow-ups.

use std::fs::File;
use std::io::BufReader;
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::{Context, Result};
use libradicl::header::RadPrelude;
use libradicl::rad_types::{MappedFragmentOrientation, TagMap};
use libradicl::readers::{ParallelRadReader, EMPTY_METACHUNK_CALLBACK};
use libradicl::record::{MappedRecord, PiscemBulkReadRecord, RecordContext};

use salmon_core::{
    is_compatible, LibraryFormat, MateStatus, ReadOrientation, ReadStrandedness, ReadType,
};
use salmon_eqclass::{range_factorize_bins, EquivalenceClassBuilder, TranscriptGroup};
use salmon_model::FragmentLengthDistribution;
use salmon_rad::{RadInputProfile, SalmonBulkRecord};

use crate::{asctime_now, logsumexp, AlignQuantOptions, AlignQuantResult};

/// salmon's `LOG_EPSILON = log(0.375e-10)` — the orphan / implausible-length penalty.
const LOG_EPSILON: f64 = -23.998_158_637_57;
/// fragments to sample (single-threaded prefix scan) when estimating the FLD.
const FLD_SAMPLE_TARGET: u64 = 100_000;

// ---------------------------------------------------------------------------
// Placement adapter
// ---------------------------------------------------------------------------

/// One placement parsed from a RAD record.
struct RadPlacement {
    tid: u32,
    /// alignment score (meaningful only for the salmon selective-alignment profile)
    score: i32,
    is_fw: bool,
    mate_fw: bool,
    status: MateStatus,
    /// fragment length, or [`salmon_rad::FRAG_LEN_UNPAIRED`] for orphan/single-end
    frag_len: u16,
}

/// A RAD record that can be turned into placements.
trait RadFragment {
    fn placements(&self) -> Vec<RadPlacement>;
}

fn mapping_type_to_status(frag_type: u8) -> MateStatus {
    use salmon_rad::frag_map_type::*;
    match frag_type {
        MAPPED_PAIR => MateStatus::PairedEndPaired,
        LEFT_ORPHAN => MateStatus::PairedEndLeft,
        RIGHT_ORPHAN => MateStatus::PairedEndRight,
        _ => MateStatus::SingleEnd,
    }
}

/// piscem `MappedFragmentOrientation` → `(read1_fw, mate_fw)`.
fn orientation_to_strands(o: MappedFragmentOrientation) -> (bool, bool) {
    match o {
        MappedFragmentOrientation::Forward => (true, false),
        MappedFragmentOrientation::Reverse => (false, false),
        MappedFragmentOrientation::ForwardReverse => (true, false),
        MappedFragmentOrientation::ReverseForward => (false, true),
        MappedFragmentOrientation::ForwardForward => (true, true),
        MappedFragmentOrientation::ReverseReverse => (false, false),
        MappedFragmentOrientation::Unknown => (true, false),
    }
}

impl RadFragment for PiscemBulkReadRecord {
    fn placements(&self) -> Vec<RadPlacement> {
        let status = mapping_type_to_status(self.frag_type);
        (0..self.refs.len())
            .map(|i| {
                let (is_fw, mate_fw) = orientation_to_strands(self.dirs[i]);
                RadPlacement {
                    tid: self.refs[i],
                    score: 0,
                    is_fw,
                    mate_fw,
                    status,
                    frag_len: self.frag_lengths[i],
                }
            })
            .collect()
    }
}

impl RadFragment for SalmonBulkRecord {
    fn placements(&self) -> Vec<RadPlacement> {
        let status = mapping_type_to_status(self.frag_type);
        self.hits
            .iter()
            .map(|h| RadPlacement {
                tid: h.tid,
                score: h.score,
                is_fw: h.is_fw,
                mate_fw: h.mate_fw,
                status,
                frag_len: h.frag_len,
            })
            .collect()
    }
}

/// Best-effort observed library format for a RAD placement (proper pairs only).
/// Mirrors `frag_format` for the inward case; orientation defaults to `Toward`
/// (RAD bulk records carry only the leftmost fragment position, so outward pairs
/// cannot be distinguished — acceptable for the common inward libraries).
fn rad_frag_format(p: &RadPlacement) -> (Option<LibraryFormat>, bool, MateStatus) {
    if p.status == MateStatus::PairedEndPaired {
        let (orientation, strandedness) = if p.is_fw != p.mate_fw {
            (
                ReadOrientation::Toward,
                if p.is_fw {
                    ReadStrandedness::SA
                } else {
                    ReadStrandedness::AS
                },
            )
        } else {
            (
                ReadOrientation::Same,
                if p.is_fw {
                    ReadStrandedness::S
                } else {
                    ReadStrandedness::A
                },
            )
        };
        (
            Some(LibraryFormat::new(
                ReadType::PairedEnd,
                orientation,
                strandedness,
            )),
            p.is_fw,
            MateStatus::PairedEndPaired,
        )
    } else {
        (None, p.is_fw, p.status)
    }
}

// ---------------------------------------------------------------------------
// Per-fragment equivalence-class assembly
// ---------------------------------------------------------------------------

struct FragCfg<'a> {
    /// salmon selective-alignment profile → soft-weight by score; else uniform
    scored: bool,
    score_exp: f64,
    fld: &'a FragmentLengthDistribution,
    lengths: &'a [u32],
    expected_format: Option<LibraryFormat>,
    ignore_incompat: bool,
    incompat_prior: f64,
    paired_lib: bool,
    range_factorization_bins: u32,
    eq_builder: &'a EquivalenceClassBuilder,
}

/// Weight one fragment's placements and add the equivalence class. Returns
/// `true` if the fragment was assigned (had a surviving placement). A trimmed
/// analog of [`crate::process_fragment`]: no error model, bias, or online phase.
fn process_rad_fragment(placements: &[RadPlacement], cfg: &FragCfg) -> bool {
    if placements.is_empty() {
        return false;
    }
    let best = placements.iter().map(|p| p.score).max().unwrap_or(0);

    let mut sp_tid: Vec<u32> = Vec::with_capacity(placements.len());
    let mut sp_eq: Vec<f64> = Vec::with_capacity(placements.len());
    for p in placements {
        // eq-class log-weight basis: soft-weight by score (SA) or uniform (sketch).
        let basis = if cfg.scored {
            -cfg.score_exp * (best - p.score) as f64
        } else {
            0.0
        };
        let rl = cfg.lengths.get(p.tid as usize).copied().unwrap_or(0) as i32;
        let proper = p.status == MateStatus::PairedEndPaired
            && p.frag_len != salmon_rad::FRAG_LEN_UNPAIRED
            && p.frag_len > 0;
        let log_frag_prob = if proper {
            // length-conditioned proper-pair probability (FLD PMF); cmf(rl) ~ 0
            // for transcripts longer than the FLD support, matching align mode.
            let flen = (p.frag_len as i32).min(rl.max(1)) as usize;
            cfg.fld.pmf(flen)
        } else if cfg.paired_lib {
            LOG_EPSILON
        } else {
            0.0
        };
        let mut aux = basis + log_frag_prob;
        if let Some(exp) = cfg.expected_format {
            let (obs, is_fw, status) = rad_frag_format(p);
            if !is_compatible(exp, obs, is_fw, status) {
                if cfg.ignore_incompat {
                    continue;
                }
                aux += cfg.incompat_prior.ln();
            }
        }
        sp_tid.push(p.tid);
        sp_eq.push(aux);
    }
    if sp_tid.is_empty() {
        return false;
    }

    // Aggregate by distinct transcript id (sorted), softmax over the logsumexp'd
    // per-transcript weights — identical to alignment mode's eq-class formation.
    let mut agg: std::collections::BTreeMap<u32, Vec<usize>> = std::collections::BTreeMap::new();
    for (k, &t) in sp_tid.iter().enumerate() {
        agg.entry(t).or_default().push(k);
    }
    let tids: Vec<u32> = agg.keys().cloned().collect();
    let eq_log: Vec<f64> = agg
        .values()
        .map(|ks| logsumexp(&ks.iter().map(|&k| sp_eq[k]).collect::<Vec<_>>()))
        .collect();
    let maxe = eq_log.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut weights: Vec<f64> = eq_log.iter().map(|&l| (l - maxe).exp()).collect();
    let wsum: f64 = weights.iter().sum();
    if wsum > 0.0 {
        for w in &mut weights {
            *w /= wsum;
        }
    }

    let group = if cfg.range_factorization_bins > 0 {
        let bins = range_factorize_bins(&weights, cfg.range_factorization_bins);
        TranscriptGroup::with_bins(tids, bins)
    } else {
        TranscriptGroup::from_sorted(tids)
    };
    cfg.eq_builder.add_group(group, weights, 1);
    true
}

// ---------------------------------------------------------------------------
// RAD file helpers
// ---------------------------------------------------------------------------

/// Open the file and read its prelude + file-tag values, returning the reader
/// positioned at the first chunk.
fn open_prelude(path: &Path) -> Result<(BufReader<File>, RadPrelude, TagMap)> {
    let mut reader = BufReader::new(
        File::open(path).with_context(|| format!("opening RAD file {}", path.display()))?,
    );
    let prelude = RadPrelude::from_bytes(&mut reader).context("parsing RAD prelude")?;
    let file_tag_map = prelude
        .file_tags
        .parse_tags_from_bytes(&mut reader)
        .context("parsing RAD file tags")?;
    Ok((reader, prelude, file_tag_map))
}

/// Reference lengths from the RAD `ref_lengths` file tag (one per reference).
fn ref_lengths_from_tags(prelude: &RadPrelude, file_tag_map: &TagMap) -> Result<Vec<u32>> {
    use libradicl::rad_types::TagValue;
    let v = match file_tag_map.get("ref_lengths") {
        Some(TagValue::ArrayU32(v)) => v.clone(),
        _ => anyhow::bail!(
            "RAD file is missing a `ref_lengths` (Array<u32>) file tag; \
             cannot compute effective lengths"
        ),
    };
    debug_assert_eq!(v.len(), prelude.hdr.ref_names.len());
    Ok(v)
}

// ---------------------------------------------------------------------------
// FLD estimation (single-threaded prefix scan)
// ---------------------------------------------------------------------------

/// Train the FLD from a prefix of the file's proper-pair fragment lengths.
fn estimate_fld<R>(path: &Path, fld: &FragmentLengthDistribution) -> Result<()>
where
    R: MappedRecord + RadFragment,
    R::ParsingContext: RecordContext,
{
    use libradicl::chunk::Chunk;
    let (mut reader, prelude, _ftm) = open_prelude(path)?;
    let ctx: R::ParsingContext = prelude.get_record_context()?;
    let mut sampled: u64 = 0;
    for _ in 0..prelude.hdr.num_chunks {
        if sampled >= FLD_SAMPLE_TARGET {
            break;
        }
        let chunk = Chunk::<R>::from_bytes(&mut reader, &ctx);
        for rec in &chunk.reads {
            for p in rec.placements() {
                if p.status == MateStatus::PairedEndPaired
                    && p.frag_len != salmon_rad::FRAG_LEN_UNPAIRED
                    && p.frag_len > 0
                {
                    fld.add_val(p.frag_len as usize, 0.0);
                    sampled += 1;
                }
            }
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Parallel weighting pass
// ---------------------------------------------------------------------------

/// Run the parallel RAD pass: `ParallelRadReader` feeds MetaChunks to a worker
/// pool that weights each fragment into the shared equivalence-class builder.
#[allow(clippy::too_many_arguments)]
fn run_rad_pass<R>(
    path: &Path,
    nthreads: usize,
    cfg: &FragCfg,
    num_processed: &AtomicU64,
    num_mapped: &AtomicU64,
) -> Result<()>
where
    R: MappedRecord + RadFragment + Send,
    R::ParsingContext: RecordContext + Clone + Send + Sync,
{
    let (reader, prelude, file_tag_map) = open_prelude(path)?;
    let nz = NonZeroUsize::new(nthreads.max(1)).unwrap();
    let mut prr =
        ParallelRadReader::<R, _>::from_prelude_and_file_tag_map(reader, prelude, file_tag_map, nz);
    let queue = prr.get_queue();
    let done = prr.is_done();

    std::thread::scope(|s| -> Result<()> {
        for _ in 0..nthreads.max(1) {
            let queue = queue.clone();
            let done = done.clone();
            s.spawn(move || loop {
                if let Some(mc) = queue.pop() {
                    for chunk in mc.iter() {
                        for rec in &chunk.reads {
                            let assigned = process_rad_fragment(&rec.placements(), cfg);
                            num_processed.fetch_add(1, Ordering::Relaxed);
                            if assigned {
                                num_mapped.fetch_add(1, Ordering::Relaxed);
                            }
                        }
                    }
                } else if done.load(Ordering::Acquire) {
                    break;
                } else {
                    std::hint::spin_loop();
                }
            });
        }
        // Producer: fill the queue until the file is exhausted (blocks).
        prr.start_chunk_parsing(EMPTY_METACHUNK_CALLBACK)
            .context("parsing RAD chunks")
    })
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Quantify from a RAD file (`salmon quant --rad`).
pub fn quantify_rad(opts: &AlignQuantOptions, rad_path: &Path) -> Result<AlignQuantResult> {
    let start_time = asctime_now();

    anyhow::ensure!(
        !(opts.seq_bias || opts.gc_bias || opts.pos_bias),
        "--seqBias/--gcBias/--posBias are not yet supported with --rad input"
    );

    // Header: profile, reference names + lengths (from the RAD file tags).
    let (_reader, prelude, file_tag_map) = open_prelude(rad_path)?;
    let profile = salmon_rad::detect_input_profile(&prelude)?;
    let names: Vec<String> = prelude.hdr.ref_names.clone();
    let lengths = ref_lengths_from_tags(&prelude, &file_tag_map)?;
    let num_refs = names.len();
    anyhow::ensure!(num_refs > 0, "RAD file has no reference sequences");
    anyhow::ensure!(
        lengths.len() == num_refs,
        "RAD ref_lengths ({}) does not match ref_names ({num_refs})",
        lengths.len()
    );
    let scored = matches!(
        profile,
        RadInputProfile::Salmon(salmon_rad::RadProfile::SelectiveAlignment)
    );
    tracing::info!(
        "quantifying from RAD ({}), {num_refs} references",
        match profile {
            RadInputProfile::PiscemBulk => "piscem bulk",
            RadInputProfile::Salmon(salmon_rad::RadProfile::Sketch) => "salmon sketch",
            RadInputProfile::Salmon(salmon_rad::RadProfile::SelectiveAlignment) => "salmon SA",
        }
    );

    // FLD: estimate from a sample prefix, then freeze.
    let mut fld =
        FragmentLengthDistribution::new(1.0, opts.fld_max, opts.fld_mean, opts.fld_sd, 4, 0.5, 1);
    match profile {
        RadInputProfile::PiscemBulk => estimate_fld::<PiscemBulkReadRecord>(rad_path, &fld)?,
        RadInputProfile::Salmon(_) => estimate_fld::<SalmonBulkRecord>(rad_path, &fld)?,
    }
    fld.cache();

    // Library-type handling (same as alignment mode).
    let paired_lib = !matches!(opts.lib_type.as_str(), "U" | "SF" | "SR" | "S");
    let expected_format = match opts.lib_type.as_str() {
        "A" => None,
        s => LibraryFormat::parse(s).ok(),
    };
    let ignore_incompat = opts.incompat_prior <= 0.0;
    let nthreads = rayon::current_num_threads().max(1);

    let eq_builder = EquivalenceClassBuilder::new();
    let num_processed = AtomicU64::new(0);
    let num_mapped = AtomicU64::new(0);
    {
        let cfg = FragCfg {
            scored,
            score_exp: opts.score_exp,
            fld: &fld,
            lengths: &lengths,
            expected_format,
            ignore_incompat,
            incompat_prior: opts.incompat_prior,
            paired_lib,
            range_factorization_bins: opts.range_factorization_bins,
            eq_builder: &eq_builder,
        };
        match profile {
            RadInputProfile::PiscemBulk => run_rad_pass::<PiscemBulkReadRecord>(
                rad_path,
                nthreads,
                &cfg,
                &num_processed,
                &num_mapped,
            )?,
            RadInputProfile::Salmon(_) => run_rad_pass::<SalmonBulkRecord>(
                rad_path,
                nthreads,
                &cfg,
                &num_processed,
                &num_mapped,
            )?,
        }
    }
    let num_processed = num_processed.into_inner();
    let num_mapped = num_mapped.into_inner();
    tracing::info!("mapped {num_mapped} / {num_processed} fragments from RAD");

    // ---- effective lengths + inference (no bias correction in v1) ----------
    let cond_means = fld.conditional_means();
    let mut eff_lengths = vec![0f64; num_refs];
    for (tid, &len) in lengths.iter().enumerate() {
        eff_lengths[tid] = salmon_model::smoothed_effective_length(&cond_means, len as usize);
    }

    let mut collapsed = eq_builder.finish();
    collapsed.update_eff_lengths(&eff_lengths);
    let num_eq_classes = collapsed.len();

    let em = if opts.skip_quant {
        salmon_infer::EmResult {
            alphas: vec![0.0; num_refs],
            iters: 0,
            converged: true,
            dropped_mass: 0.0,
        }
    } else {
        salmon_infer::optimize(&collapsed, num_refs, &opts.em, Some(&eff_lengths))
    };
    let inference_truncated_mass = em.dropped_mass;
    let counts = em.alphas;

    let rates: Vec<f64> = (0..num_refs)
        .map(|i| {
            if eff_lengths[i] > 0.0 {
                counts[i] / eff_lengths[i]
            } else {
                0.0
            }
        })
        .collect();
    let rate_sum: f64 = rates.iter().sum();
    let tpm: Vec<f64> = rates
        .iter()
        .map(|r| {
            if rate_sum > 0.0 {
                r / rate_sum * 1e6
            } else {
                0.0
            }
        })
        .collect();

    let length_classes =
        salmon_model::compute_length_quantiles(&lengths, salmon_model::NUM_LENGTH_CLASSES);

    let packed = salmon_infer::PackedEqClasses::from_collapsed(&collapsed, num_refs);
    let ambig = salmon_infer::ambiguity_counts(&packed);
    let bootstraps: Vec<Vec<f64>> = if opts.skip_quant {
        Vec::new()
    } else if opts.num_bootstraps > 0 {
        salmon_infer::bootstrap(&packed, &opts.em, opts.num_bootstraps, 0x5A13_0000)
    } else if opts.num_gibbs_samples > 0 {
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
        salmon_infer::gibbs_sample(&packed, &eff_lengths, &counts, &gopts, 0x6217_0000)
    } else {
        Vec::new()
    };

    let result = AlignQuantResult {
        names,
        lengths,
        eff_lengths,
        tpm,
        counts,
        num_processed,
        num_mapped,
        inference_truncated_mass,
        num_eq_classes,
        frag_len_mean: fld.mean(),
        frag_len_sd: fld.sd(),
        length_classes,
        frag_len_dist: fld.log_pmf().iter().map(|lp| lp.exp()).collect(),
        start_time,
        bias_dump: salmon_model::dumps::BiasDump::default(),
        ambig,
        bootstraps,
    };
    crate::write_outputs(opts, &result)?;
    Ok(result)
}
