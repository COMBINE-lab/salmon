//! Deterministic alignment-mode RAD producer.
//!
//! `salmon quant -a <bam> --deterministic` runs this single pass over a
//! name-grouped transcriptome BAM: it pairs each fragment's records (reusing the
//! alignment-mode [`crate::pair_records`] / [`crate::frag_format`]), writes the
//! placements as a salmon RAD (selective-alignment profile, per-hit score = the
//! BAM `AS` tag), and accumulates everything the RAD reader needs baked into the
//! header so requant is a single pass:
//!   * an order-independent [`DiscreteFld`] (fragment-length distribution +
//!     library-format detection), and
//!   * when bias correction is requested, a [`NaiveEqBuilder`] feeding a rough
//!     seed EM whose abundances are baked as `initial_abundances`.
//!
//! The caller (CLI `run_deterministic_align`) then hands the RAD to the existing
//! [`crate::quantify_rad`], which takes its fully-baked single-pass fused branch.
//!
//! This is score-based: it carries the BAM `AS` tag and does not run salmon's
//! online alignment error model (that path stays in [`crate::quantify_alignments`]).
//! A fast-follow adds an order-independent error model for fidelity.

use std::path::Path;

use anyhow::{bail, Context, Result};
use crossbeam_channel::bounded;
use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::Header;

use salmon_core::{LibraryFormat, MateStatus};
use salmon_eqclass::{NaiveEqBuilder, NaivePlacement, NAIVE_NO_FMT};
use salmon_model::{smoothed_effective_length, DiscreteFld};
use salmon_rad::{
    frag_map_type, ChunkCodec, FragmentChunkBuf, RadHit, RadOutputWriter, RadProfile,
    SalmonBulkRecord,
};

use crate::{
    coordinate_sorted_unusable, frag_format, is_sam_path, open_sam_reader, pair_records,
    read_alignment_header, record_to_frag, AlignQuantOptions, FragRecord,
};

/// Flush a worker's chunk buffer once it reaches this many record bytes.
const CHUNK_FLUSH_BYTES: usize = 64 * 1024;
/// Fragment groups per reader→worker minibatch (matches the online pass).
const MINIBATCH: usize = 1000;

/// Outcome of the deterministic BAM→RAD pass.
pub struct AlignRadSummary {
    pub num_processed: u64,
    pub num_mapped: u64,
}

/// Project a name-grouped transcriptome BAM into a fully-baked salmon RAD at
/// `rad_path`. Reference ids/names/lengths come from the BAM `@SQ` header.
pub fn write_alignment_rad(
    opts: &AlignQuantOptions,
    rad_path: &Path,
    codec: ChunkCodec,
) -> Result<AlignRadSummary> {
    let header = read_alignment_header(&opts.bam)?;
    if coordinate_sorted_unusable(&header) {
        bail!(
            "alignment input looks coordinate-sorted; salmon needs read-name-grouped \
             alignments (e.g. `samtools collate` or `samtools sort -n`)"
        );
    }
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

    let bias_on = opts.seq_bias || opts.gc_bias || opts.pos_bias;
    // Paired unless the library type is explicitly single-end (matches the
    // alignment-mode `paired_lib` derivation in `quantify_alignments`).
    let is_paired = !matches!(opts.lib_type.as_str(), "U" | "SF" | "SR" | "S");
    // The user-fixed expected format (`None` for auto / unstranded), resolved
    // against the FLD-detected format at end of pass.
    let expected_format = match opts.lib_type.as_str() {
        "A" | "IU" | "U" => None,
        s => LibraryFormat::parse(s).ok(),
    };

    let name_refs: Vec<&str> = names.iter().map(String::as_str).collect();
    let mut writer = RadOutputWriter::create(
        rad_path,
        &name_refs,
        &lengths,
        is_paired,
        RadProfile::SelectiveAlignment,
        opts.fld_max + 1,
        codec,
    )?;

    let fld = DiscreteFld::new(opts.fld_max);
    let naive = bias_on.then(NaiveEqBuilder::new);

    let nthreads = rayon::current_num_threads().max(1);
    let summary = stream_grouped_pass(&opts.bam, nthreads, &writer, &fld, naive.as_ref())?;

    // ---- end-of-pass bake (one RAD pass) ---------------------------------
    let (frag_dist, det_fmt) = fld.finish(opts.fld_mean, opts.fld_sd);
    writer.set_frag_length_dist(frag_dist.log_pmf());
    let resolved_fmt = expected_format.or(det_fmt);
    if let Some(f) = resolved_fmt {
        writer.set_library_format(f.format_id());
    }

    if let Some(nb) = naive.as_ref() {
        // Rough seed EM on the naive (uniform-weight) eq-classes, with
        // library-incompatible placements dropped against the resolved format,
        // over the fixed FLD's effective lengths — deliberately under-converged,
        // deterministic, and baked so `quantify_rad` needs no extra pass.
        let cond_means = frag_dist.conditional_means();
        let eff_lengths: Vec<f64> = lengths
            .iter()
            .map(|&len| smoothed_effective_length(&cond_means, len as usize))
            .collect();
        let mut collapsed = nb.finish(resolved_fmt);
        collapsed.update_eff_lengths(&eff_lengths);
        let mut ro = opts.em.clone();
        ro.min_iter = opts.bias_seed_em_iters;
        ro.max_iter = opts.bias_seed_em_iters;
        ro.min_alpha = 0.0;
        let em = salmon_infer::optimize(&collapsed, num_refs, &ro, Some(&eff_lengths));
        writer.set_initial_abundances(&em.alphas);
    }

    writer.finalize()?;
    Ok(summary)
}

/// Dispatch the grouped pass over SAM or BAM (chosen by extension), mirroring
/// [`crate::stream_online_pass`]'s decoder setup.
fn stream_grouped_pass(
    bam_path: &Path,
    nthreads: usize,
    writer: &RadOutputWriter,
    fld: &DiscreteFld,
    naive: Option<&NaiveEqBuilder>,
) -> Result<AlignRadSummary> {
    if is_sam_path(bam_path) {
        let mut reader = open_sam_reader(bam_path)?;
        let header = reader.read_header().context("reading SAM header")?;
        run_grouped_pass(reader.records(), &header, nthreads, writer, fld, naive)
    } else {
        let file = std::fs::File::open(bam_path)
            .with_context(|| format!("opening {}", bam_path.display()))?;
        let bgzf_workers = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4)
            .clamp(1, 4);
        let bgzf_workers = std::num::NonZeroUsize::new(bgzf_workers).unwrap();
        let decoder = noodles_bgzf::io::MultithreadedReader::with_worker_count(bgzf_workers, file);
        let mut reader = bam::io::Reader::from(decoder);
        let header = reader.read_header().context("reading BAM header")?;
        run_grouped_pass(reader.records(), &header, nthreads, writer, fld, naive)
    }
}

/// One reader thread groups records by canonical read name into minibatches; a
/// worker pool turns each group into a RAD record + FLD/naive-eq tallies.
fn run_grouped_pass<R, I>(
    records: I,
    header: &Header,
    nthreads: usize,
    writer: &RadOutputWriter,
    fld: &DiscreteFld,
    naive: Option<&NaiveEqBuilder>,
) -> Result<AlignRadSummary>
where
    R: sam::alignment::Record + Send + Sync,
    I: Iterator<Item = std::io::Result<R>> + Send,
{
    let (tx, rx) = bounded::<Vec<Vec<R>>>(2 * nthreads + 1);
    let codec = writer.codec();

    std::thread::scope(|scope| -> Result<AlignRadSummary> {
        // Reader/parser thread: group consecutive same-canonical-name records.
        let reader = scope.spawn(move || -> Result<()> {
            let mut cur_name: Vec<u8> = Vec::new();
            let mut have = false;
            let mut group: Vec<R> = Vec::new();
            let mut batch: Vec<Vec<R>> = Vec::with_capacity(MINIBATCH);
            for result in records {
                let rec = result.context("reading alignment record")?;
                if rec.flags().map(|f| f.is_unmapped()).unwrap_or(true) {
                    continue;
                }
                let Some(name) = rec.name() else { continue };
                let cname = crate::canonical_name(name.as_ref()).to_vec();
                if !have {
                    cur_name = cname;
                    have = true;
                } else if cname != cur_name {
                    batch.push(std::mem::take(&mut group));
                    cur_name = cname;
                    if batch.len() >= MINIBATCH {
                        if tx.send(std::mem::take(&mut batch)).is_err() {
                            return Ok(()); // all workers gone
                        }
                        batch = Vec::with_capacity(MINIBATCH);
                    }
                }
                group.push(rec);
            }
            if have && !group.is_empty() {
                batch.push(group);
            }
            if !batch.is_empty() {
                let _ = tx.send(batch);
            }
            Ok(())
        });

        // Worker pool.
        let mut workers = Vec::with_capacity(nthreads);
        for _ in 0..nthreads {
            let rx = rx.clone();
            workers.push(scope.spawn(move || -> Result<(u64, u64)> {
                let mut buf = FragmentChunkBuf::with_capacity_codec(CHUNK_FLUSH_BYTES, codec);
                let mut frags: Vec<FragRecord> = Vec::new();
                let mut processed = 0u64;
                let mut mapped = 0u64;
                while let Ok(raw_batch) = rx.recv() {
                    for raw_group in &raw_batch {
                        processed += 1;
                        frags.clear();
                        for r in raw_group {
                            if let Some((_, f)) = record_to_frag(r, header, false) {
                                frags.push(f);
                            }
                        }
                        if frags.is_empty() {
                            continue;
                        }
                        let placements = pair_records(&frags);
                        if placements.is_empty() {
                            continue;
                        }
                        let single = placements.len() == 1;
                        let mut hits = Vec::with_capacity(placements.len());
                        let mut statuses = Vec::with_capacity(placements.len());
                        let mut naive_sig: Vec<NaivePlacement> = Vec::new();
                        // (frag_len, is_fw, mate_fw) of a unique proper pair → FLD.
                        let mut fld_sample: Option<(usize, bool, bool)> = None;
                        for pl in &placements {
                            let (fmt, is_fw, status) = frag_format(&frags, &pl.idxs);
                            let paired = status == MateStatus::PairedEndPaired;
                            // idxs are sorted by ref pos, so [0] is the leftmost.
                            let pos = frags[pl.idxs[0]].pos as i32;
                            let frag_len = pl
                                .idxs
                                .iter()
                                .map(|&i| frags[i].frag_len)
                                .max()
                                .unwrap_or(0);
                            // Orphan/SE read-length proxy for the bounded-CMF weight:
                            // the alignment's reference span (no read_len in the BAM
                            // FragRecord).
                            let read_len = frags[pl.idxs[0]].ref_span as i32;
                            // Fragment alignment score = sum of the mates' AS tags.
                            let score: i32 = pl.idxs.iter().map(|&i| frags[i].score).sum();
                            let mate_fw = paired
                                && pl
                                    .idxs
                                    .iter()
                                    .map(|&i| &frags[i])
                                    .find(|r| !r.is_read1)
                                    .map(|r2| !r2.is_reverse)
                                    .unwrap_or(false);
                            hits.push(RadHit::for_placement(
                                pl.tid, is_fw, mate_fw, pos, paired, frag_len, read_len, score,
                            ));
                            statuses.push(status);
                            if single && paired && frag_len > 0 {
                                fld_sample = Some((frag_len as usize, is_fw, mate_fw));
                            }
                            if naive.is_some() {
                                let fmt_id = fmt.map(|f| f.format_id()).unwrap_or(NAIVE_NO_FMT);
                                naive_sig.push(NaivePlacement {
                                    tid: pl.tid,
                                    fmt_id,
                                    is_fw,
                                    status,
                                });
                            }
                        }
                        let frag_type = frag_map_type::fragment_level(statuses.iter().copied());
                        let rec = SalmonBulkRecord::new(frag_type, hits);
                        buf.write(&rec, writer.context())?;
                        if buf.byte_len() >= CHUNK_FLUSH_BYTES {
                            let bytes = buf.take_bytes()?;
                            writer.append_chunk_bytes(&bytes)?;
                        }
                        if let Some((len, is_fw, mate_fw)) = fld_sample {
                            fld.add(len, is_fw, mate_fw);
                        }
                        if let Some(nb) = naive {
                            nb.add(naive_sig);
                        }
                        mapped += 1;
                    }
                }
                if buf.nrec() > 0 {
                    let bytes = buf.take_bytes()?;
                    writer.append_chunk_bytes(&bytes)?;
                }
                Ok((processed, mapped))
            }));
        }
        drop(rx);

        reader
            .join()
            .map_err(|_| anyhow::anyhow!("alignment reader thread panicked"))??;

        let mut num_processed = 0u64;
        let mut num_mapped = 0u64;
        for w in workers {
            let (p, m) = w
                .join()
                .map_err(|_| anyhow::anyhow!("alignment worker thread panicked"))??;
            num_processed += p;
            num_mapped += m;
        }
        Ok(AlignRadSummary {
            num_processed,
            num_mapped,
        })
    })
}
