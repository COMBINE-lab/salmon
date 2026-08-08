//! Deterministic alignment-mode RAD producer.
//!
//! # What "deterministic" means here
//!
//! An ordinary run learns its models *while* mapping, from whichever fragments
//! each thread happened to see first, so tiny differences accumulate and the
//! output depends on the thread count. `--deterministic` splits the run in two:
//! a first pass that only reads alignments and accumulates order-independent
//! integer tallies, and a second that quantifies from those fixed models. This
//! module is the first pass — it converts a BAM into a RAD with everything the
//! second pass needs already baked into the header.
//!
//! `salmon quant -a <bam> --deterministic` runs over a name-grouped
//! transcriptome BAM, pairs each fragment's records (reusing the alignment-mode
//! [`crate::pair_records`] / [`crate::frag_format`]), writes the placements as a
//! salmon RAD (selective-alignment profile), and bakes everything the RAD reader
//! needs into the header so requant is a single pass:
//!   * an order-independent [`DiscreteFld`] (fragment-length distribution +
//!     library-format detection), and
//!   * when bias correction is requested, a [`NaiveEqBuilder`] feeding a rough
//!     seed EM whose abundances are baked as `initial_abundances`.
//!
//! The caller (CLI `run_deterministic_align`) then hands the RAD to the existing
//! [`crate::quantify_rad`], which takes its fully-baked single-pass fused branch.
//!
//! **Per-hit score.** Two interpretations, chosen by whether a transcriptome is
//! available (`-t`, or index sequences via `ref_seqs`):
//!   * **No `-t` (or `--noErrorModel`)** — one BAM pass; the per-hit score is the
//!     BAM `AS` tag (`score_kind = AS`), soft-weighted at quant time. This is the
//!     original MVP.
//!   * **With `-t`** — an **order-independent alignment error model**. Pass A
//!     trains a [`CountingAlignmentModel`] (integer transition counts, merged by
//!     integer add across threads ⇒ partition-independent) alongside the FLD /
//!     naive-eq tally; the counts are normalized once into a fixed log-space
//!     model; pass B scores each placement `Σ(fg−bg)` against that fixed model
//!     and stores the quantized log-weight (`score_kind = LOGWEIGHT`). This
//!     restores the fidelity of salmon's online error model while staying
//!     byte-reproducible. Training is uniform (one count per placement-mate, no
//!     posterior weighting) — a deliberate, documented divergence from the online
//!     model's posterior-weighted training (which is order-dependent by design).

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

use crate::error_model::{AlignmentModel, CountingAlignmentModel};
use crate::{
    coordinate_sorted_unusable, frag_format, is_sam_path, load_ref_bytes, open_sam_reader,
    pair_records, read_alignment_header, record_to_frag, AlignQuantOptions, FragRecord,
};

/// Flush a worker's chunk buffer once it reaches this many record bytes.
const CHUNK_FLUSH_BYTES: usize = 64 * 1024;
/// Fragment groups per reader→worker minibatch (matches the online pass).
const MINIBATCH: usize = 1000;
/// Laplace pseudocount per error-model transition cell (salmon's default alpha).
const ERROR_MODEL_ALPHA: f64 = 1.0;

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

    // Reference sequences for the error model: supplied index sequences take
    // precedence, else load the `-t` FASTA (mirrors `quantify_rad`). Absent ⇒ no
    // error model (the MVP AS-score path).
    let ref_bytes: Option<salmon_core::RefSeqs> = if let Some(rs) = &opts.ref_seqs {
        Some(rs.clone())
    } else if let Some(t) = &opts.transcripts {
        Some(salmon_core::RefSeqs::from_sequences(load_ref_bytes(
            t, &names,
        )?))
    } else {
        None
    };
    // The error model is OPT-IN in deterministic mode (`--errorModel`). It is off
    // by default: across every truth-bearing benchmark (uniform + realistic
    // Illumina errors, 50/76 bp) AS scoring is at least as accurate, and it runs a
    // single BAM pass instead of two. See the deterministic-error-model notes.
    if opts.deterministic_error_model && ref_bytes.is_none() {
        tracing::warn!(
            "--errorModel needs a transcriptome (-t / index sequences) to train against; \
             none provided, scoring by the BAM AS tag instead"
        );
    }
    let use_error_model =
        opts.deterministic_error_model && ref_bytes.is_some() && !opts.no_error_model;

    let name_refs: Vec<&str> = names.iter().map(String::as_str).collect();
    let mut writer = RadOutputWriter::create(
        rad_path,
        &name_refs,
        &lengths,
        is_paired,
        RadProfile::SelectiveAlignment,
        opts.fld_max + 1,
        codec,
        // Fragments came from an external aligner, and there is no salmon index
        // behind them; a reader reads the absent index as "unknown".
        &salmon_rad::WriterProvenance {
            mapping_type: salmon_rad::MappingType::Alignment,
            index: None,
        },
    )?;

    let fld = DiscreteFld::new(opts.fld_max);
    let naive = bias_on.then(NaiveEqBuilder::new);
    let nthreads = rayon::current_num_threads().max(1);

    let summary = if use_error_model {
        let ref_bytes = ref_bytes.as_ref().unwrap();
        let error_bins = opts.num_error_bins.max(1);
        // Pass A — train the counting error model + FLD + naive-eq (records not
        // written yet; scores aren't final). Per-worker counting models are
        // merged by integer add, so the total is partition-independent.
        let models = stream_grouped(&opts.bam, nthreads, || {
            TrainWorker::new(error_bins, &fld, naive.as_ref(), ref_bytes)
        })?;
        let mut counting = CountingAlignmentModel::new(error_bins);
        for m in &models {
            counting.combine(m);
        }
        let fixed = counting.finalize(ERROR_MODEL_ALPHA);
        // Pass B — score each placement against the fixed model and write records.
        let sums = stream_grouped(&opts.bam, nthreads, || {
            ScoreWriteWorker::new(&writer, ref_bytes, &fixed)
        })?;
        writer.set_score_kind(salmon_rad::SCORE_KIND_LOGWEIGHT);
        reduce_summary(sums)
    } else {
        // One pass — write records carrying the BAM `AS` score + FLD + naive-eq.
        let sums = stream_grouped(&opts.bam, nthreads, || {
            WriteWorker::new(&writer, &fld, naive.as_ref())
        })?;
        reduce_summary(sums)
    };

    // ---- end-of-pass bake (single quant pass) ----------------------------
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

fn reduce_summary(parts: Vec<(u64, u64)>) -> AlignRadSummary {
    let (mut num_processed, mut num_mapped) = (0u64, 0u64);
    for (p, m) in parts {
        num_processed += p;
        num_mapped += m;
    }
    AlignRadSummary {
        num_processed,
        num_mapped,
    }
}

// ---------------------------------------------------------------------------
// Shared per-group analysis
// ---------------------------------------------------------------------------

/// Per-placement geometry every worker needs (RAD record + FLD + naive-eq),
/// independent of how the per-hit score is computed.
struct PlacementInfo {
    tid: u32,
    /// fragment record indices for this placement (sorted by reference pos), used
    /// to compute the per-hit score (AS sum, or error-model `Σ(fg−bg)`).
    idxs: Vec<usize>,
    is_fw: bool,
    mate_fw: bool,
    pos: i32,
    paired: bool,
    frag_len: i32,
    read_len: i32,
    status: MateStatus,
    /// observed library format id for the naive-eq signature (`NAIVE_NO_FMT` if none).
    fmt_id: u8,
}

/// One read group's paired records + per-placement geometry + the unique-proper
/// FLD sample. `None` when the group produced no usable placement.
struct GroupData {
    frags: Vec<FragRecord>,
    infos: Vec<PlacementInfo>,
    /// `(frag_len, is_fw, mate_fw)` of a unique proper pair → the FLD.
    fld_sample: Option<(usize, bool, bool)>,
}

/// Pair a read group's records and compute the shared per-placement geometry.
/// `need_seq` requests the read's 2-bit bases (needed by the error model).
fn analyze_group<R: sam::alignment::Record>(
    group: &[R],
    header: &Header,
    need_seq: bool,
    scratch: &mut Vec<FragRecord>,
) -> Option<GroupData> {
    scratch.clear();
    for r in group {
        if let Some((_, f)) = record_to_frag(r, header, need_seq) {
            scratch.push(f);
        }
    }
    if scratch.is_empty() {
        return None;
    }
    let frags = std::mem::take(scratch);
    let placements = pair_records(&frags);
    if placements.is_empty() {
        *scratch = frags;
        return None;
    }
    let single = placements.len() == 1;
    let mut infos = Vec::with_capacity(placements.len());
    let mut fld_sample = None;
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
        // Orphan/SE read-length proxy for the bounded-CMF weight: the alignment's
        // reference span (no read_len in the BAM FragRecord).
        let read_len = frags[pl.idxs[0]].ref_span as i32;
        let mate_fw = paired
            && pl
                .idxs
                .iter()
                .map(|&i| &frags[i])
                .find(|r| !r.is_read1)
                .map(|r2| !r2.is_reverse)
                .unwrap_or(false);
        if single && paired && frag_len > 0 {
            fld_sample = Some((frag_len as usize, is_fw, mate_fw));
        }
        infos.push(PlacementInfo {
            tid: pl.tid,
            idxs: pl.idxs.clone(),
            is_fw,
            mate_fw,
            pos,
            paired,
            frag_len,
            read_len,
            status,
            fmt_id: fmt.map(|f| f.format_id()).unwrap_or(NAIVE_NO_FMT),
        });
    }
    Some(GroupData {
        frags,
        infos,
        fld_sample,
    })
}

/// Feed a group's placements to the shared FLD + naive-eq accumulators.
fn tally_fld_naive(data: &GroupData, fld: &DiscreteFld, naive: Option<&NaiveEqBuilder>) {
    if let Some((len, is_fw, mate_fw)) = data.fld_sample {
        fld.add(len, is_fw, mate_fw);
    }
    if let Some(nb) = naive {
        let sig: Vec<NaivePlacement> = data
            .infos
            .iter()
            .map(|i| NaivePlacement {
                tid: i.tid,
                fmt_id: i.fmt_id,
                is_fw: i.is_fw,
                status: i.status,
            })
            .collect();
        nb.add(sig);
    }
}

/// Serialize a group's placements to a RAD record with per-hit `scores` (one per
/// placement, same order as `data.infos`) and append to the worker buffer.
fn write_record(
    data: &GroupData,
    scores: &[i32],
    buf: &mut FragmentChunkBuf,
    writer: &RadOutputWriter,
) -> Result<()> {
    let hits: Vec<RadHit> = data
        .infos
        .iter()
        .zip(scores)
        .map(|(i, &score)| {
            RadHit::for_placement(
                i.tid, i.is_fw, i.mate_fw, i.pos, i.paired, i.frag_len, i.read_len, score,
            )
        })
        .collect();
    let frag_type = frag_map_type::fragment_level(data.infos.iter().map(|i| i.status));
    let rec = SalmonBulkRecord::new(frag_type, hits);
    buf.write(&rec, writer.context())?;
    if buf.byte_len() >= CHUNK_FLUSH_BYTES {
        let bytes = buf.take_bytes()?;
        writer.append_chunk_bytes(&bytes)?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Per-group workers
// ---------------------------------------------------------------------------

/// One reader thread groups records by canonical read name into minibatches; a
/// pool of these workers turns each group into whatever the pass needs. The
/// per-`R` `handle` keeps the worker state (buffers, models) `R`-independent so
/// the same worker serves both the SAM and BAM decoders.
trait GroupWorker: Send {
    type Output: Send;
    fn handle<R: sam::alignment::Record>(&mut self, header: &Header, group: &[R]) -> Result<()>;
    fn finish(self) -> Result<Self::Output>;
}

/// MVP pass: write records carrying the BAM `AS` score + tally FLD / naive-eq.
struct WriteWorker<'a> {
    buf: FragmentChunkBuf,
    writer: &'a RadOutputWriter,
    fld: &'a DiscreteFld,
    naive: Option<&'a NaiveEqBuilder>,
    scratch: Vec<FragRecord>,
    processed: u64,
    mapped: u64,
}

impl<'a> WriteWorker<'a> {
    fn new(
        writer: &'a RadOutputWriter,
        fld: &'a DiscreteFld,
        naive: Option<&'a NaiveEqBuilder>,
    ) -> Self {
        Self {
            buf: FragmentChunkBuf::with_capacity_codec(CHUNK_FLUSH_BYTES, writer.codec()),
            writer,
            fld,
            naive,
            scratch: Vec::new(),
            processed: 0,
            mapped: 0,
        }
    }
}

impl GroupWorker for WriteWorker<'_> {
    type Output = (u64, u64);
    fn handle<R: sam::alignment::Record>(&mut self, header: &Header, group: &[R]) -> Result<()> {
        self.processed += 1;
        let Some(data) = analyze_group(group, header, false, &mut self.scratch) else {
            return Ok(());
        };
        // Fragment alignment score = sum of the mates' AS tags.
        let scores: Vec<i32> = data
            .infos
            .iter()
            .map(|i| i.idxs.iter().map(|&j| data.frags[j].score).sum())
            .collect();
        write_record(&data, &scores, &mut self.buf, self.writer)?;
        tally_fld_naive(&data, self.fld, self.naive);
        self.mapped += 1;
        self.scratch = data.frags;
        Ok(())
    }
    fn finish(mut self) -> Result<(u64, u64)> {
        if self.buf.nrec() > 0 {
            let bytes = self.buf.take_bytes()?;
            self.writer.append_chunk_bytes(&bytes)?;
        }
        Ok((self.processed, self.mapped))
    }
}

/// Error-model pass A: tally the counting model + FLD / naive-eq (no writing).
struct TrainWorker<'a> {
    model: CountingAlignmentModel,
    fld: &'a DiscreteFld,
    naive: Option<&'a NaiveEqBuilder>,
    ref_bytes: &'a salmon_core::RefSeqs,
    scratch: Vec<FragRecord>,
}

impl<'a> TrainWorker<'a> {
    fn new(
        error_bins: usize,
        fld: &'a DiscreteFld,
        naive: Option<&'a NaiveEqBuilder>,
        ref_bytes: &'a salmon_core::RefSeqs,
    ) -> Self {
        Self {
            model: CountingAlignmentModel::new(error_bins),
            fld,
            naive,
            ref_bytes,
            scratch: Vec::new(),
        }
    }
}

impl GroupWorker for TrainWorker<'_> {
    type Output = CountingAlignmentModel;
    fn handle<R: sam::alignment::Record>(&mut self, header: &Header, group: &[R]) -> Result<()> {
        let Some(data) = analyze_group(group, header, true, &mut self.scratch) else {
            return Ok(());
        };
        // Uniform per-placement training: one count per placement-mate. Integer
        // increments make the merged model independent of thread partitioning.
        for info in &data.infos {
            let refseq = self.ref_bytes.get(info.tid as usize);
            let Some(refseq) = refseq.filter(|s| !s.is_empty()) else {
                continue;
            };
            for (rank, &i) in info.idxs.iter().enumerate() {
                let r = &data.frags[i];
                self.model
                    .count(&r.read_2bit, refseq, r.pos, &r.ops, rank == 0);
            }
        }
        tally_fld_naive(&data, self.fld, self.naive);
        self.scratch = data.frags;
        Ok(())
    }
    fn finish(self) -> Result<CountingAlignmentModel> {
        Ok(self.model)
    }
}

/// Error-model pass B: score each placement against the fixed model, write records.
struct ScoreWriteWorker<'a> {
    buf: FragmentChunkBuf,
    writer: &'a RadOutputWriter,
    ref_bytes: &'a salmon_core::RefSeqs,
    model: &'a AlignmentModel,
    scratch: Vec<FragRecord>,
    processed: u64,
    mapped: u64,
}

impl<'a> ScoreWriteWorker<'a> {
    fn new(
        writer: &'a RadOutputWriter,
        ref_bytes: &'a salmon_core::RefSeqs,
        model: &'a AlignmentModel,
    ) -> Self {
        Self {
            buf: FragmentChunkBuf::with_capacity_codec(CHUNK_FLUSH_BYTES, writer.codec()),
            writer,
            ref_bytes,
            model,
            scratch: Vec::new(),
            processed: 0,
            mapped: 0,
        }
    }
}

impl GroupWorker for ScoreWriteWorker<'_> {
    type Output = (u64, u64);
    fn handle<R: sam::alignment::Record>(&mut self, header: &Header, group: &[R]) -> Result<()> {
        self.processed += 1;
        let Some(data) = analyze_group(group, header, true, &mut self.scratch) else {
            return Ok(());
        };
        // Per placement: Σ(fg−bg) over its mate(s) under the fixed model (0 when
        // the reference is missing, matching the online path), quantized to i32.
        let scores: Vec<i32> = data
            .infos
            .iter()
            .map(|info| {
                let refseq = self
                    .ref_bytes
                    .get(info.tid as usize)
                    .filter(|s| !s.is_empty());
                let ll = match refseq {
                    None => 0.0,
                    Some(refseq) => info
                        .idxs
                        .iter()
                        .enumerate()
                        .map(|(rank, &i)| {
                            let r = &data.frags[i];
                            let (fg, bg) = self.model.log_likelihood(
                                &r.read_2bit,
                                refseq,
                                r.pos,
                                &r.ops,
                                rank == 0,
                            );
                            fg - bg
                        })
                        .sum::<f64>(),
                };
                (ll * salmon_rad::SCORE_LOG_SCALE)
                    .round()
                    .clamp(i32::MIN as f64, i32::MAX as f64) as i32
            })
            .collect();
        write_record(&data, &scores, &mut self.buf, self.writer)?;
        self.mapped += 1;
        self.scratch = data.frags;
        Ok(())
    }
    fn finish(mut self) -> Result<(u64, u64)> {
        if self.buf.nrec() > 0 {
            let bytes = self.buf.take_bytes()?;
            self.writer.append_chunk_bytes(&bytes)?;
        }
        Ok((self.processed, self.mapped))
    }
}

// ---------------------------------------------------------------------------
// Reader / worker orchestration
// ---------------------------------------------------------------------------

/// Dispatch the grouped pass over SAM or BAM (chosen by extension), mirroring
/// [`crate::stream_online_pass`]'s decoder setup, driving one `make_worker()`
/// per thread.
fn stream_grouped<MK, W>(
    bam_path: &Path,
    nthreads: usize,
    make_worker: MK,
) -> Result<Vec<W::Output>>
where
    MK: Fn() -> W + Sync + Send,
    W: GroupWorker,
{
    if is_sam_path(bam_path) {
        let mut reader = open_sam_reader(bam_path)?;
        let header = reader.read_header().context("reading SAM header")?;
        run_grouped(reader.records(), &header, nthreads, make_worker)
    } else {
        let file = std::fs::File::open(bam_path)
            .with_context(|| format!("opening {}", bam_path.display()))?;
        let bgzf_workers = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4)
            .clamp(1, 4);
        let bgzf_workers = std::num::NonZeroUsize::new(bgzf_workers).unwrap();
        {
            let decoder =
                noodles_bgzf::io::MultithreadedReader::with_worker_count(bgzf_workers, file);
            let mut reader = bam::io::Reader::from(decoder);
            let header = reader.read_header().context("reading BAM header")?;
            run_grouped(reader.records(), &header, nthreads, make_worker)
        }
    }
}

/// One reader thread groups consecutive same-canonical-name records into
/// minibatches; a worker pool consumes them. Returns each worker's `finish()`.
fn run_grouped<R, I, MK, W>(
    records: I,
    header: &Header,
    nthreads: usize,
    make_worker: MK,
) -> Result<Vec<W::Output>>
where
    R: sam::alignment::Record + Send + Sync,
    I: Iterator<Item = std::io::Result<R>> + Send,
    MK: Fn() -> W + Sync,
    W: GroupWorker,
{
    let (tx, rx) = bounded::<Vec<Vec<R>>>(2 * nthreads + 1);

    std::thread::scope(|scope| -> Result<Vec<W::Output>> {
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
        let make_worker = &make_worker;
        let mut workers = Vec::with_capacity(nthreads);
        for _ in 0..nthreads {
            let rx = rx.clone();
            workers.push(scope.spawn(move || -> Result<W::Output> {
                let mut w = make_worker();
                while let Ok(batch) = rx.recv() {
                    for group in &batch {
                        w.handle(header, group)?;
                    }
                }
                w.finish()
            }));
        }
        drop(rx);

        reader
            .join()
            .map_err(|_| anyhow::anyhow!("alignment reader thread panicked"))??;

        let mut outs = Vec::with_capacity(nthreads);
        for w in workers {
            outs.push(
                w.join()
                    .map_err(|_| anyhow::anyhow!("alignment worker thread panicked"))??,
            );
        }
        Ok(outs)
    })
}
