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
    /// enable sequence-specific bias correction (`--seqBias`)
    pub seq_bias: bool,
    /// enable fragment-GC bias correction (`--gcBias`)
    pub gc_bias: bool,
    /// enable positional bias correction (`--posBias`)
    pub pos_bias: bool,
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
            seq_bias: false,
            gc_bias: false,
            pos_bias: false,
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

/// Load a (optionally gzip'd) transcriptome FASTA and return the (ASCII) base
/// sequences aligned to the BAM's reference order (`names`); a name absent from
/// the FASTA yields an empty sequence (its model contributions are skipped). The
/// same bytes feed both the error model (2-bit on the fly) and the bias models.
fn load_ref_bytes(path: &Path, names: &[String]) -> Result<Vec<Vec<u8>>> {
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
            cur_seq.extend(line.trim_end().bytes());
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

    // The error model and bias models need the transcriptome (salmon requires `-t`).
    let use_error_model = opts.transcripts.is_some() && !opts.no_error_model;
    let bias_on = opts.seq_bias || opts.gc_bias || opts.pos_bias;
    anyhow::ensure!(
        !bias_on || opts.transcripts.is_some(),
        "--seqBias/--gcBias/--posBias in alignment mode require -t/--targets (the transcriptome FASTA)"
    );
    let ref_bytes: Vec<Vec<u8>> = if use_error_model || bias_on {
        load_ref_bytes(opts.transcripts.as_ref().unwrap(), &names)?
    } else {
        Vec::new()
    };

    let mut model = use_error_model.then(|| AlignmentModel::new(1.0, 4));
    // Online (dual-phase) abundances: develop running estimates so the error
    // model and bias models are trained/collected with abundance-aware posteriors
    // in a single streaming pass (salmon's online phase), rather than two passes.
    let ref_lens_u64: Vec<u64> = lengths.iter().map(|&l| l as u64).collect();
    let online = (use_error_model || bias_on)
        .then(|| salmon_infer::OnlineInference::new(&ref_lens_u64, 0.05, 0.65, 5_000_000));

    // Observed bias accumulators (single-threaded → owned).
    use salmon_model::seqbias::{CONTEXT_LEFT, CONTEXT_LENGTH, CONTEXT_RIGHT};
    let mut seq_obs = opts
        .seq_bias
        .then(|| (salmon_model::SBModel::new(), salmon_model::SBModel::new()));
    let mut gc_obs = opts.gc_bias.then(salmon_model::GcFragModel::default_model);
    let gc_prefix: Vec<Vec<u32>> = if opts.gc_bias {
        ref_bytes.iter().map(|s| salmon_model::gc_prefix(s)).collect()
    } else {
        Vec::new()
    };
    let length_quantiles: Option<Vec<u32>> = opts.pos_bias.then(|| {
        salmon_model::compute_length_quantiles(&lengths, salmon_model::NUM_LENGTH_CLASSES)
    });
    let length_class: Option<Vec<usize>> = length_quantiles.as_ref().map(|q| {
        lengths.iter().map(|&l| salmon_model::length_class_index(q, l)).collect()
    });
    let mut pos_obs = opts.pos_bias.then(|| {
        let mk = || {
            (0..salmon_model::NUM_LENGTH_CLASSES)
                .map(|_| salmon_model::SimplePosBias::default())
                .collect::<Vec<_>>()
        };
        (mk(), mk())
    });

    // ---- single online pass: train error model + collect bias + eq-classes ---
    const MINIBATCH: u64 = 1000;
    let mut num_processed = 0u64;
    let mut num_mapped = 0u64;
    let mut frag_count = 0u64;
    let mut log_fm = 0.0;
    for_each_fragment(&opts.bam, use_error_model || bias_on, |recs| {
        num_processed += 1;
        if recs.is_empty() {
            return;
        }
        num_mapped += 1;
        // one forgetting-mass timestep per minibatch of fragments
        if let Some(o) = online.as_ref() {
            if frag_count % MINIBATCH == 0 {
                log_fm = o.next_log_fm();
            }
        }
        frag_count += 1;

        let by_tid = group_by_tid(recs);
        let n = by_tid.len();
        let frag_len = recs.iter().map(|r| r.frag_len).max().unwrap_or(0);
        if frag_len > 0 {
            fld.add_val(frag_len as usize, 0.0);
        }
        let use_aux = online.as_ref().is_none_or(|o| o.num_assigned() >= 5000);

        // Per-tid: alignment conditional log-weight (eq-class) + online log-aux.
        let mut tids: Vec<u32> = Vec::with_capacity(n);
        let mut eq_log: Vec<f64> = Vec::with_capacity(n);
        let mut online_log: Vec<f64> = Vec::with_capacity(n);
        // geometry for bias (per tid): (frag_start, frag_end, proper)
        let mut geom: Vec<(usize, usize, bool)> = Vec::with_capacity(n);
        let best_as = by_tid
            .iter()
            .map(|(_, idxs)| idxs.iter().map(|&i| recs[i].score).sum::<i32>())
            .max()
            .unwrap_or(0);
        for (tid, idxs) in &by_tid {
            let refseq = ref_bytes.get(*tid as usize).map(|v| v.as_slice()).unwrap_or(&[]);
            // conditional log-weight: error model Σ(fg-bg), else AS-based
            let basis = if let Some(m) = model.as_ref() {
                if refseq.is_empty() {
                    0.0
                } else {
                    let mut ll = 0.0;
                    for (rank, &i) in idxs.iter().enumerate() {
                        let r = &recs[i];
                        let (fg, bg) = m.log_likelihood(&r.read_2bit, refseq, r.pos, &r.ops, rank == 0);
                        ll += fg - bg;
                    }
                    ll
                }
            } else {
                let as_sum: i32 = idxs.iter().map(|&i| recs[i].score).sum();
                -opts.score_exp * (best_as - as_sum) as f64
            };
            // fragment geometry + length-normalized online aux (start-pos + FLD)
            let rl = lengths[*tid as usize] as i32;
            let flen = idxs.iter().map(|&i| recs[i].frag_len).max().unwrap_or(0);
            let proper = idxs.len() >= 2 && flen > 0;
            let frag_start = recs[idxs[0]].pos;
            let frag_end = frag_start + (flen.max(1) as usize) - 1;
            let start_pos = if proper && flen <= rl {
                -(((rl - flen + 1) as f64).ln())
            } else {
                -((rl.max(1) as f64).ln())
            };
            let fld_term = if proper && use_aux { fld.pmf(flen as usize) } else { 0.0 };
            tids.push(*tid);
            eq_log.push(basis);
            online_log.push(basis + start_pos + fld_term);
            geom.push((frag_start, frag_end, proper));
        }

        // eq-class weights = softmax(eq_log)
        let maxe = eq_log.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mut weights: Vec<f64> = eq_log.iter().map(|&l| (l - maxe).exp()).collect();
        let wsum: f64 = weights.iter().sum();
        if wsum > 0.0 {
            for w in &mut weights {
                *w /= wsum;
            }
        }

        // abundance-aware posteriors (online) drive model/bias training
        let post: Vec<f64> = match online.as_ref() {
            Some(o) => {
                let maps: Vec<(u32, f64)> = tids.iter().cloned().zip(online_log.iter().cloned()).collect();
                o.assign_fragment(&maps, log_fm)
            }
            None => weights.clone(),
        };

        // train the error model + collect bias models, weighted by posteriors
        let collecting = online.as_ref().is_none_or(|o| o.collecting());
        if collecting {
            for (ti, (tid, idxs)) in by_tid.iter().enumerate() {
                let p = post[ti];
                if p <= 0.0 {
                    continue;
                }
                let refseq = ref_bytes.get(*tid as usize).map(|v| v.as_slice()).unwrap_or(&[]);
                if let Some(m) = model.as_mut() {
                    if !refseq.is_empty() {
                        let lw = log_fm + p.ln();
                        for (rank, &i) in idxs.iter().enumerate() {
                            let r = &recs[i];
                            m.update(&r.read_2bit, refseq, r.pos, &r.ops, rank == 0, lw);
                        }
                    }
                }
                if refseq.is_empty() {
                    continue;
                }
                let (fs, fe, proper) = geom[ti];
                let rl = refseq.len();
                // sequence bias: 5' context (fwd) at frag start, 3' (rc) at frag end
                if let Some(obs) = seq_obs.as_mut() {
                    let s5 = fs as i32 - CONTEXT_LEFT as i32;
                    if s5 >= 0 && (s5 as usize + CONTEXT_LENGTH) <= rl {
                        obs.0.add_context(&refseq[s5 as usize..s5 as usize + CONTEXT_LENGTH], false, p);
                    }
                    if proper {
                        let s3 = fe as i32 - CONTEXT_RIGHT as i32;
                        if s3 >= 0 && (s3 as usize + CONTEXT_LENGTH) <= rl {
                            obs.1.add_context(&refseq[s3 as usize..s3 as usize + CONTEXT_LENGTH], true, p);
                        }
                    }
                }
                // fragment-GC bias (paired fragments only)
                if let (Some(gc), true) = (gc_obs.as_mut(), proper && fe < rl) {
                    if let Some((ff, cf)) =
                        salmon_model::gc_desc(&gc_prefix[*tid as usize], rl as i32, fs as i32, fe as i32)
                    {
                        gc.inc(ff, cf, p);
                    }
                }
                // positional bias
                if let Some(pos) = pos_obs.as_mut() {
                    let lc = length_class.as_ref().unwrap()[*tid as usize];
                    pos.0[lc].add_mass(fs as i32, rl as i32, p.ln());
                    if proper {
                        pos.1[lc].add_mass(fe as i32, rl as i32, p.ln());
                    }
                }
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

    // ---- base effective lengths --------------------------------------------
    fld.cache();
    let cond_means = fld.conditional_means();
    let mut eff_lengths = vec![0f64; num_refs];
    for (tid, len) in lengths.iter().enumerate() {
        eff_lengths[tid] = salmon_model::smoothed_effective_length(&cond_means, *len as usize);
    }

    let mut collapsed = eq_builder.finish();
    collapsed.update_eff_lengths(&eff_lengths);
    let num_eq_classes = collapsed.len();
    let mut em = optimize(&collapsed, num_refs, &opts.em);

    // ---- bias-corrected effective lengths (shared with reads mode) ----------
    if bias_on {
        let log_pmf = fld.log_pmf();
        let pmf_lin: Vec<f64> = log_pmf.iter().map(|lp| lp.exp()).collect();
        let (fld_cdf, fld_low, fld_high) = salmon_model::seqbias::fld_cdf_and_bounds(&pmf_lin);
        let k = if opts.seq_bias { CONTEXT_LENGTH } else { 1 };
        let refseq_of = |t: usize| ref_bytes[t].as_slice();

        let seq = seq_obs.map(|(mut of, mut or)| {
            of.normalize();
            or.normalize();
            let (ef, er) =
                salmon_model::build_expected(num_refs, refseq_of, &em.alphas, &eff_lengths, &fld_cdf);
            (of, or, ef, er)
        });
        let gc_ratio_model = gc_obs.map(|mut obs| {
            let mut exp = salmon_model::build_expected_gc(
                num_refs,
                refseq_of,
                |t| gc_prefix[t].as_slice(),
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
        let pos_models = pos_obs.map(|(mut of, mut or)| {
            for x in of.iter_mut().chain(or.iter_mut()) {
                x.finalize();
            }
            let (ef, er) = salmon_model::build_expected_pos(
                num_refs,
                |t| lengths[t] as usize,
                &em.alphas,
                &eff_lengths,
                &fld_cdf,
                length_quantiles.as_ref().unwrap(),
                k,
            );
            (of, or, ef, er)
        });

        for tid in 0..num_refs {
            if em.alphas[tid] < 1e-8 {
                continue;
            }
            let s = ref_bytes[tid].as_slice();
            let pos_vecs: Option<(Vec<f64>, Vec<f64>)> = pos_models.as_ref().map(|(ofw, orc, efw, erc)| {
                let lc = length_class.as_ref().unwrap()[tid];
                let rl = s.len();
                let (mut o5, mut e5) = (vec![0.0; rl], vec![0.0; rl]);
                let (mut o3, mut e3) = (vec![0.0; rl], vec![0.0; rl]);
                ofw[lc].project_weights(&mut o5);
                efw[lc].project_weights(&mut e5);
                orc[lc].project_weights(&mut o3);
                erc[lc].project_weights(&mut e3);
                (
                    salmon_model::positional_factor(&o5, &e5),
                    salmon_model::positional_factor(&o3, &e3),
                )
            });
            let bias = salmon_model::BiasInputs {
                seq: seq.as_ref().map(|(of, or, ef, er)| (of, ef, or, er)),
                gc: gc_ratio_model.as_ref().map(|g| (g, gc_prefix[tid].as_slice())),
                pos: pos_vecs.as_ref().map(|(pf, pr)| (pf.as_slice(), pr.as_slice())),
            };
            eff_lengths[tid] = salmon_model::corrected_effective_length_full(
                s,
                &fld_cdf,
                fld_low,
                fld_high,
                &bias,
                eff_lengths[tid],
                salmon_model::GC_SAMP_STRIDE,
            );
        }
        collapsed.update_eff_lengths(&eff_lengths);
        em = optimize(&collapsed, num_refs, &opts.em);
    }
    let counts = em.alphas;

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
