//! Genome-alignment quantification via bramble projection.
//!
//! `salmon quant -a <genome.bam> --annotation <gtf|gff>` reads a name-collated
//! genome BAM, projects each read group into transcriptome coordinates with
//! `bramble-rs`, and writes a salmon RAD (selective-alignment profile) that the
//! existing [`crate::quantify_rad`] then consumes — so the whole inference / bias
//! / EM / output stack is reused unchanged and the result is deterministic.
//!
//! The transcriptome reference set (names + lengths) comes entirely from the
//! annotation via the bramble `G2TTree`, so no transcriptome FASTA is needed for
//! the ref set; a genome FASTA (`--genome`) is needed only to reconstruct
//! transcript sequences when bias correction is requested.
//!
//! bramble exposes no projected CIGAR, so salmon's alignment error model cannot
//! run here; each placement's quality is bramble's normalized `similarity_score`.
//! As with `--deterministic`, this is a single RAD pass: the FLD, library format,
//! and (when bias is on) rough abundances are baked during projection.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use crossbeam_channel::bounded;
use noodles_bam as bam;
use noodles_sam as sam;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_sam::Header;

use bramble_rs::annotation::{load_transcripts, Transcript};
use bramble_rs::fasta::FastaDb;
use bramble_rs::g2t::build_g2t_from_refnames;
use bramble_rs::{
    project_group_with, G2TTree, GenomicAlignment, ProjectedAlignment, ProjectionConfig,
    ProjectionContext,
};

use salmon_core::{observed_paired_format, LibraryFormat, MateStatus};
use salmon_eqclass::{NaiveEqBuilder, NaivePlacement, NAIVE_NO_FMT};
use salmon_model::{smoothed_effective_length, DiscreteFld};
use salmon_rad::{
    frag_map_type, ChunkCodec, FragmentChunkBuf, RadHit, RadOutputWriter, RadProfile,
    SalmonBulkRecord,
};

use crate::{coordinate_sorted_unusable, is_sam_path, open_sam_reader, read_alignment_header};

const CHUNK_FLUSH_BYTES: usize = 64 * 1024;
const MINIBATCH: usize = 1000;

/// Options for the genome→transcriptome projection pass.
pub struct GenomeProjectOptions {
    /// name-collated genome BAM (spliced alignments)
    pub bam: PathBuf,
    /// GTF/GFF transcript annotation
    pub annotation: PathBuf,
    /// genome FASTA, required only to reconstruct transcript sequences for bias
    pub genome_fasta: Option<PathBuf>,
    /// salmon library type string (resolved against the FLD-detected format)
    pub lib_type: String,
    pub junc_miss_discount: Option<f64>,
    /// any of seq/GC/positional bias requested (drives the rough-EM + ref_seqs)
    pub bias: bool,
    pub fld_mean: f64,
    pub fld_sd: f64,
    pub fld_max: usize,
    pub bias_seed_em_iters: u32,
    pub em: salmon_infer::EmOptions,
    pub rad_codec: ChunkCodec,
}

/// What the projection pass produced: the transcriptome ref set + (for bias) the
/// reconstructed transcript sequences, handed to `quantify_rad`.
pub struct ProjectionArtifacts {
    pub names: Vec<String>,
    pub lengths: Vec<u32>,
    pub ref_seqs: Option<Vec<Vec<u8>>>,
    pub num_processed: u64,
    pub num_mapped: u64,
}

/// Project a genome BAM into a fully-baked salmon RAD at `rad_path`.
pub fn project_genome_bam_to_rad(
    opts: &GenomeProjectOptions,
    rad_path: &Path,
) -> Result<ProjectionArtifacts> {
    let header = read_alignment_header(&opts.bam)?;
    if coordinate_sorted_unusable(&header) {
        bail!(
            "genome alignment input looks coordinate-sorted; projection needs read-name-grouped \
             alignments (e.g. `samtools collate` or `samtools sort -n`)"
        );
    }
    // @SQ order defines the reference-id space bramble's `ref_id` indexes into.
    let refnames: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();
    anyhow::ensure!(
        !refnames.is_empty(),
        "genome BAM header has no reference sequences"
    );

    let txs = load_transcripts(&opts.annotation)
        .with_context(|| format!("loading annotation {}", opts.annotation.display()))?;
    anyhow::ensure!(
        !txs.is_empty(),
        "annotation {} has no transcripts",
        opts.annotation.display()
    );
    let g2t = build_g2t_from_refnames(&txs, &refnames, None)
        .context("building the genome->transcriptome index from the annotation")?;
    let names: Vec<String> = g2t.transcript_names().to_vec();
    let lengths: Vec<u32> = g2t.transcript_lengths().to_vec();
    let num_refs = names.len();

    // Transcript sequences for bias (seq/GC/pos): reconstruct from the genome.
    let ref_seqs = if opts.bias {
        let gfa = opts.genome_fasta.as_ref().ok_or_else(|| {
            anyhow::anyhow!(
                "bias correction in genome-projection mode needs --genome <FASTA> to reconstruct \
                 transcript sequences"
            )
        })?;
        let fdb = FastaDb::load(gfa)
            .with_context(|| format!("loading genome FASTA {}", gfa.display()))?;
        Some(reconstruct_ref_seqs(&txs, &fdb))
    } else {
        None
    };

    let is_paired = !matches!(opts.lib_type.as_str(), "U" | "SF" | "SR" | "S");
    let name_refs: Vec<&str> = names.iter().map(String::as_str).collect();
    let mut writer = RadOutputWriter::create(
        rad_path,
        &name_refs,
        &lengths,
        is_paired,
        RadProfile::SelectiveAlignment,
        opts.fld_max + 1,
        opts.rad_codec,
    )?;

    // Short-read projection only. Long-read genome quantification is out of scope
    // for salmon — use oarfish (https://github.com/COMBINE-lab/oarfish) instead.
    let mut cfg = ProjectionConfig::short_read();
    if let Some(d) = opts.junc_miss_discount {
        cfg.junc_miss_discount = d;
    }

    let fld = DiscreteFld::new(opts.fld_max);
    let naive = opts.bias.then(NaiveEqBuilder::new);

    // Transcript strand per bramble tid (tid order == `txs` order; `build_g2t`
    // assigns ids in that order), for the genomic->transcript orientation flip.
    let tx_is_minus: Vec<bool> = txs.iter().map(|t| t.strand == '-').collect();

    let nthreads = rayon::current_num_threads().max(1);
    let is_short = true;
    let (num_processed, num_mapped) = stream_project_pass(
        &opts.bam,
        nthreads,
        &g2t,
        &cfg,
        &writer,
        &fld,
        naive.as_ref(),
        &tx_is_minus,
        is_short,
    )?;

    // ---- single-pass bake -------------------------------------------------
    let (frag_dist, det_fmt) = fld.finish(opts.fld_mean, opts.fld_sd);
    writer.set_frag_length_dist(frag_dist.log_pmf());
    let expected_format = match opts.lib_type.as_str() {
        "A" | "IU" | "U" => None,
        s => LibraryFormat::parse(s).ok(),
    };
    let resolved_fmt = expected_format.or(det_fmt);
    if let Some(f) = resolved_fmt {
        writer.set_library_format(f.format_id());
    }
    if let Some(nb) = naive.as_ref() {
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

    Ok(ProjectionArtifacts {
        names,
        lengths,
        ref_seqs,
        num_processed,
        num_mapped,
    })
}

/// Reconstruct each transcript's 5'→3' sequence (tid order = `txs` slice order,
/// which is how `build_g2t` assigns ids) by concatenating its exon slices in
/// ascending genomic order and reverse-complementing `-`-strand transcripts.
fn reconstruct_ref_seqs(txs: &[Transcript], fdb: &FastaDb) -> Vec<Vec<u8>> {
    txs.iter()
        .map(|tx| {
            let mut exons = tx.exons.clone();
            exons.sort_by_key(|e| e.start);
            let mut seq = Vec::new();
            for e in &exons {
                if let Some(s) = fdb.get_slice(&tx.seqname, e.start, e.end) {
                    seq.extend_from_slice(&s);
                }
            }
            if tx.strand == '-' {
                revcomp(&mut seq);
            }
            seq
        })
        .collect()
}

/// In-place reverse complement (ACGT; other bases → `N`).
fn revcomp(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = match *b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        };
    }
}

/// Map a noodles CIGAR op kind to its SAM op code (the form bramble expects).
fn kind_to_sam_code(k: sam::alignment::record::cigar::op::Kind) -> u8 {
    use sam::alignment::record::cigar::op::Kind;
    match k {
        Kind::Match => 0,
        Kind::Insertion => 1,
        Kind::Deletion => 2,
        Kind::Skip => 3,
        Kind::SoftClip => 4,
        Kind::HardClip => 5,
        Kind::Pad => 6,
        Kind::SequenceMatch => 7,
        Kind::SequenceMismatch => 8,
    }
}

/// Read an `A`-type (character) tag as a `+`/`-` strand char.
fn strand_tag<R: sam::alignment::Record>(rec: &R, tag: [u8; 2]) -> Option<char> {
    match rec.data().get(&Tag::from(tag)).and_then(|r| r.ok())? {
        Value::Character(c) => Some(c as char),
        _ => None,
    }
}

/// The original aligner (e.g. STAR) `AS` alignment score of a genome record.
/// bramble's C++ passes this through unchanged as the projected short-read score
/// (it only recomputes AS for long reads), so we carry it alongside each
/// `GenomicAlignment` and use it as the placement score. Missing/absent AS -> 0.
fn genome_alignment_score<R: sam::alignment::Record>(rec: &R) -> i32 {
    rec.data()
        .get(&Tag::ALIGNMENT_SCORE)
        .and_then(|r| r.ok())
        .and_then(|v| match v {
            Value::Int8(x) => Some(x as i32),
            Value::UInt8(x) => Some(x as i32),
            Value::Int16(x) => Some(x as i32),
            Value::UInt16(x) => Some(x as i32),
            Value::Int32(x) => Some(x),
            Value::UInt32(x) => Some(x as i32),
            _ => None,
        })
        .unwrap_or(0)
}

/// Build a bramble `GenomicAlignment` from one mapped noodles record. `query_name`
/// is the (shared) canonical group name. Returns `None` for unmapped/unnamed/
/// no-reference records.
fn record_to_genomic_alignment<R: sam::alignment::Record>(
    rec: &R,
    header: &Header,
    query_name: &str,
) -> Option<GenomicAlignment> {
    let flags = rec.flags().ok()?;
    if flags.is_unmapped() {
        return None;
    }
    let ref_id = rec.reference_sequence_id(header)?.ok()? as i32;
    let ref_start = rec.alignment_start().and_then(|r| r.ok())?.get() as i64; // 1-based SAM POS
    let cigar: Vec<(u32, u8)> = rec
        .cigar()
        .iter()
        .filter_map(|r| r.ok())
        .map(|op| (op.len() as u32, kind_to_sam_code(op.kind())))
        .collect();
    // Query length from the query-consuming CIGAR ops (M/I/S/=/X), robust to a
    // missing SEQ field.
    let read_len: usize = cigar
        .iter()
        .filter(|(_, c)| matches!(c, 0 | 1 | 4 | 7 | 8))
        .map(|(l, _)| *l as usize)
        .sum();
    let hit_index = rec
        .data()
        .get(&Tag::HIT_INDEX)
        .and_then(|r| r.ok())
        .and_then(|v| match v {
            Value::Int8(x) => Some(x as i32),
            Value::UInt8(x) => Some(x as i32),
            Value::Int16(x) => Some(x as i32),
            Value::UInt16(x) => Some(x as i32),
            Value::Int32(x) => Some(x),
            Value::UInt32(x) => Some(x as i32),
            _ => None,
        })
        .unwrap_or(0);
    let mate_unmapped = flags.is_mate_unmapped();
    let mate_ref_id = (!mate_unmapped)
        .then(|| rec.mate_reference_sequence_id(header).and_then(|r| r.ok()))
        .flatten()
        .map(|t| t as i32);
    let mate_ref_start = rec
        .mate_alignment_start()
        .and_then(|r| r.ok())
        .map(|p| p.get() as i64);

    Some(GenomicAlignment {
        query_name: query_name.to_string(),
        ref_id,
        ref_start,
        is_reverse: flags.is_reverse_complemented(),
        cigar,
        sequence: None, // soft-clip rescue (use_fasta) is off in this MVP
        is_paired: flags.is_segmented(),
        is_first_in_pair: flags.is_first_segment(),
        xs_strand: strand_tag(rec, [b'X', b'S']),
        ts_strand: strand_tag(rec, [b't', b's']),
        hit_index,
        mate_ref_id,
        mate_ref_start,
        mate_is_unmapped: mate_unmapped,
        read_len,
    })
}

/// Collapse one read group's projected alignments into a salmon RAD record.
/// Returns the record, the (frag_len,is_fw,mate_fw) FLD sample for a fragment
/// that maps as a proper pair to exactly one transcript, and the naive-eq
/// signature. `None` if nothing projected.
fn build_projected_record(
    gas: &[GenomicAlignment],
    gn_as: &[i32],
    proj: &[ProjectedAlignment],
    tx_is_minus: &[bool],
    is_short: bool,
) -> Option<(
    SalmonBulkRecord,
    Option<(usize, bool, bool)>,
    Vec<NaivePlacement>,
)> {
    if proj.is_empty() {
        return None;
    }
    let mut by_tid: BTreeMap<u32, Vec<&ProjectedAlignment>> = BTreeMap::new();
    for p in proj {
        by_tid.entry(p.transcript_id).or_default().push(p);
    }
    let single_tid = by_tid.len() == 1;
    let mut hits = Vec::with_capacity(by_tid.len());
    let mut statuses = Vec::with_capacity(by_tid.len());
    let mut naive_sig = Vec::with_capacity(by_tid.len());
    let mut fld_sample = None;

    let is_read1 =
        |e: &ProjectedAlignment| gas.get(e.input_index).is_none_or(|g| g.is_first_in_pair);
    // Fragment orientation on the transcript. bramble's *public* projected
    // `is_reverse` is `matched_strand != read_splice_strand`, which is useless
    // for unstranded reads (no XS/ts tag → read strand `.` → always reverse).
    // bramble-cli instead sets each mate's output strand as
    //   output_reverse = genomic_reverse XOR (transcript_strand == '-')
    // (the genomic BAM strand, flipped when the transcript is on the minus
    // strand). We replicate that exactly: the genomic strand comes from the
    // source record and the transcript strand from the annotation (tid order ==
    // `txs` order, as `build_g2t` assigns ids). This yields correct inward pairs
    // and a correct `-l A` format (IU on unstranded input, not a spurious ISF).
    let genomic_fw = |e: &ProjectedAlignment| {
        let genomic_rev = gas
            .get(e.input_index)
            .map_or(e.is_reverse, |g| g.is_reverse);
        let tx_minus = tx_is_minus
            .get(e.transcript_id as usize)
            .copied()
            .unwrap_or(false);
        !(genomic_rev ^ tx_minus)
    };

    // Per-placement score, matching C++ bramble. Short reads: pass the original
    // aligner (STAR) AS through unchanged — C++ never rewrites short-read AS and
    // leaves `similarity_score` unweighted, so junc_hits does NOT skew the split
    // across a read's isoforms. Long reads: C++ rewrites AS = (genome_AS + clip) *
    // similarity_score; clip rescue isn't exposed on ProjectedAlignment, so use
    // genome_AS * similarity_score (clip ~ 0 without soft-clip rescue).
    let placement_score = |e: &ProjectedAlignment| -> i32 {
        let g = gn_as.get(e.input_index).copied().unwrap_or(0);
        if is_short {
            g
        } else {
            (g as f64 * e.similarity_score).round() as i32
        }
    };

    for (&tid, entries) in &by_tid {
        let proper = entries.len() >= 2
            && entries
                .iter()
                .any(|e| e.same_transcript_as_mate && e.insert_size != 0);
        if proper {
            let r1 = entries
                .iter()
                .copied()
                .find(|e| is_read1(e))
                .unwrap_or(entries[0]);
            let r2 = entries
                .iter()
                .copied()
                .find(|e| !is_read1(e))
                .unwrap_or(entries[1]);
            let pos = r1
                .transcript_start
                .min(r2.transcript_start)
                .saturating_sub(1) as i32;
            let is_fw = genomic_fw(r1);
            let mate_fw = genomic_fw(r2);
            let frag_len = entries
                .iter()
                .map(|e| e.insert_size.unsigned_abs())
                .max()
                .unwrap_or(0) as i32;
            let score = placement_score(r1) + placement_score(r2);
            hits.push(RadHit::for_placement(
                tid, is_fw, mate_fw, pos, true, frag_len, 0, score,
            ));
            statuses.push(MateStatus::PairedEndPaired);
            naive_sig.push(NaivePlacement {
                tid,
                fmt_id: observed_paired_format(is_fw, mate_fw).format_id(),
                is_fw,
                status: MateStatus::PairedEndPaired,
            });
            if single_tid && frag_len > 0 {
                fld_sample = Some((frag_len as usize, is_fw, mate_fw));
            }
        } else {
            // Orphan(s) / single-end on this transcript: one hit per projected
            // alignment, with the read length in the frag_len slot.
            for e in entries {
                let paired = gas.get(e.input_index).is_some_and(|g| g.is_paired);
                let status = if !paired {
                    MateStatus::SingleEnd
                } else if is_read1(e) {
                    MateStatus::PairedEndLeft
                } else {
                    MateStatus::PairedEndRight
                };
                let is_fw = genomic_fw(e);
                let pos = e.transcript_start.saturating_sub(1) as i32;
                let read_len = e.query_aligned_len as i32;
                let score = placement_score(e);
                hits.push(RadHit::for_placement(
                    tid, is_fw, false, pos, false, 0, read_len, score,
                ));
                statuses.push(status);
                naive_sig.push(NaivePlacement {
                    tid,
                    fmt_id: NAIVE_NO_FMT,
                    is_fw,
                    status,
                });
            }
        }
    }
    let frag_type = frag_map_type::fragment_level(statuses.iter().copied());
    Some((
        SalmonBulkRecord::new(frag_type, hits),
        fld_sample,
        naive_sig,
    ))
}

/// Dispatch the projection pass over SAM or BAM (chosen by extension).
#[allow(clippy::too_many_arguments)]
fn stream_project_pass(
    bam_path: &Path,
    nthreads: usize,
    g2t: &G2TTree,
    cfg: &ProjectionConfig,
    writer: &RadOutputWriter,
    fld: &DiscreteFld,
    naive: Option<&NaiveEqBuilder>,
    tx_is_minus: &[bool],
    is_short: bool,
) -> Result<(u64, u64)> {
    if is_sam_path(bam_path) {
        let mut reader = open_sam_reader(bam_path)?;
        let header = reader.read_header().context("reading SAM header")?;
        run_project_pass(
            reader.records(),
            &header,
            nthreads,
            g2t,
            cfg,
            writer,
            fld,
            naive,
            tx_is_minus,
            is_short,
        )
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
            run_project_pass(
                reader.records(),
                &header,
                nthreads,
                g2t,
                cfg,
                writer,
                fld,
                naive,
                tx_is_minus,
                is_short,
            )
        }
    }
}

/// Reader thread groups records by canonical name; workers project each group via
/// a reused `ProjectionContext` and write the resulting RAD records.
#[allow(clippy::too_many_arguments)]
fn run_project_pass<R, I>(
    records: I,
    header: &Header,
    nthreads: usize,
    g2t: &G2TTree,
    cfg: &ProjectionConfig,
    writer: &RadOutputWriter,
    fld: &DiscreteFld,
    naive: Option<&NaiveEqBuilder>,
    tx_is_minus: &[bool],
    is_short: bool,
) -> Result<(u64, u64)>
where
    R: sam::alignment::Record + Send + Sync,
    I: Iterator<Item = std::io::Result<R>> + Send,
{
    let (tx, rx) = bounded::<Vec<Vec<R>>>(2 * nthreads + 1);
    let codec = writer.codec();

    std::thread::scope(|scope| -> Result<(u64, u64)> {
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
                            return Ok(());
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

        let mut workers = Vec::with_capacity(nthreads);
        for _ in 0..nthreads {
            let rx = rx.clone();
            workers.push(scope.spawn(move || -> Result<(u64, u64)> {
                let mut pctx = ProjectionContext::new();
                let mut buf = FragmentChunkBuf::with_capacity_codec(CHUNK_FLUSH_BYTES, codec);
                let mut gas: Vec<GenomicAlignment> = Vec::new();
                // Original aligner AS per genome record, parallel to `gas` (same
                // index space as `ProjectedAlignment::input_index`).
                let mut gn_as: Vec<i32> = Vec::new();
                let mut processed = 0u64;
                let mut mapped = 0u64;
                while let Ok(raw_batch) = rx.recv() {
                    for raw_group in &raw_batch {
                        processed += 1;
                        gas.clear();
                        gn_as.clear();
                        // The whole group shares one canonical query name.
                        let qname = raw_group
                            .iter()
                            .find_map(|r| {
                                r.name().map(|n| {
                                    String::from_utf8_lossy(crate::canonical_name(n.as_ref()))
                                        .into_owned()
                                })
                            })
                            .unwrap_or_default();
                        for r in raw_group {
                            if let Some(ga) = record_to_genomic_alignment(r, header, &qname) {
                                gas.push(ga);
                                gn_as.push(genome_alignment_score(r));
                            }
                        }
                        if gas.is_empty() {
                            continue;
                        }
                        let proj = project_group_with(&gas, g2t, cfg, &mut pctx);
                        let Some((rec, fld_sample, naive_sig)) =
                            build_projected_record(&gas, &gn_as, &proj, tx_is_minus, is_short)
                        else {
                            continue; // nothing passed the similarity filter
                        };
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
            .map_err(|_| anyhow::anyhow!("genome-BAM reader thread panicked"))??;

        let mut num_processed = 0u64;
        let mut num_mapped = 0u64;
        for w in workers {
            let (p, m) = w
                .join()
                .map_err(|_| anyhow::anyhow!("projection worker thread panicked"))??;
            num_processed += p;
            num_mapped += m;
        }
        Ok((num_processed, num_mapped))
    })
}
