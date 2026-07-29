//! Chunked `--writeBam` output.
//!
//! # Ordering contract
//!
//! Exactly one guarantee: **all records for a fragment are contiguous**. A
//! worker encodes a whole fragment before testing the chunk-size threshold, and
//! chunks are written whole, so a fragment can never be split across chunks and
//! interleaved with another worker's output.
//!
//! Nothing else is ordered. Fragments appear in whatever order workers finish
//! chunks, which varies with thread count and scheduling. That is deliberate:
//! any stronger guarantee would mean making workers wait for each other.
//!
//! # Where the pipeline is (and is not) serial
//!
//! - Record encoding: fully parallel, into per-worker buffers, no shared state.
//! - Handoff: a bounded channel; sending moves only the `Vec` descriptor, never
//!   the bytes.
//! - This module's writer thread: copies chunk bytes into BGZF-sized blocks.
//!   Serial, but a `memcpy` — ~4% utilized at the ~410 MiB/s peak record rate
//!   measured here.
//! - Deflate: fully parallel across [`compression_worker_count`] workers. This
//!   is the expensive part, and it is the part that scales.
//! - noodles' frame writer: emits compressed blocks in submission order. Serial,
//!   but only a `memcpy` plus the file write.
//!
//! So the only serialized work is proportional to bytes moved, not to
//! compression effort, which is the most concurrency an ordered BGZF stream
//! admits. Nothing waits on a *particular* worker at any point.

use std::io::{self, Write};
use std::num::NonZeroUsize;
use std::path::Path;
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

use crossbeam_channel::{Receiver, Sender};
use noodles_bgzf as bgzf;
use salmon_index::SalmonIndex;
use salmon_map::ScoredMapping;

use crate::mapping_record::{self, AlignmentRecord, CigarKind};

const CHUNK_TARGET: usize = 2 * 1024 * 1024;
const INITIAL_CAPACITY: usize = 64 * 1024;
/// A chunk only overshoots [`CHUNK_TARGET`] by one fragment's records, so this
/// ceiling exists to drop pathological growth, not as a routine limit.
const RETAIN_CAPACITY: usize = 2 * CHUNK_TARGET;

/// Uncompressed BAM bytes one BGZF deflate worker consumes per second.
///
/// Measured with `noodles_bgzf::io::MultithreadedWriter` (zlib-rs backend,
/// default level) over 1 GiB of real salmon BAM records: 168.5, 167.3, 167.6,
/// 167.0, 165.1 and 153.9 MiB/s per worker at 1/2/4/8/16/32 workers. Flat, so
/// deflate scales linearly and one constant describes it.
const DEFLATE_MIB_PER_SEC: f64 = 165.0;

/// Peak uncompressed BAM bytes one mapping thread produces per second.
///
/// Measured by sweeping the deflate worker count until compression stopped
/// being the limit: record output plateaus at ~413 MiB/s across 8 mapping
/// threads, i.e. ~52 MiB/s each. This is an intentionally pessimistic number —
/// it comes from a tiny-transcriptome fixture where mapping is nearly free and
/// almost every read yields several records, so real workloads (expensive
/// mapping, same record volume) produce *less* per thread and need fewer
/// deflate workers than this predicts.
const RECORD_MIB_PER_SEC_PER_MAPPER: f64 = 52.0;

/// Ceiling on deflate workers. `DEFLATE_MIB_PER_SEC * MAX_COMPRESSION_WORKERS`
/// is ~1.3 GiB/s of compression, three times the fastest record production
/// measured at any thread count, so past this point extra workers cannot be the
/// thing standing between salmon and full throughput — they would only take
/// cores from mapping.
const MAX_COMPRESSION_WORKERS: usize = 8;

/// How many BGZF deflate workers to run alongside `mapping_threads`.
///
/// Deflate is the only stage whose cost scales with output volume, and each
/// worker is a core not mapping reads, so we provision the smallest number that
/// stops compression being the bottleneck and no more. Production scales with
/// the mapping threads and consumption with the deflate workers, so the balance
/// point is
///
/// ```text
/// workers = ceil(mapping_threads * 52 MiB/s / 165 MiB/s)  ~=  mapping_threads / 3
/// ```
///
/// The sweep this comes from (2.4M fragments, mapping-phase time relative to a
/// no-BAM-output run, best of 5):
///
/// ```text
///  -p 8    C=1 2.32x   C=2 1.18x   C=3 1.09x   C=4 1.06x   C=6 1.00x   C=8 0.98x
/// ```
///
/// At `-p 16` and `-p 32` only the `C=1` point reproduced cleanly (2.08x and
/// 2.20x); every `C >= 2` point there fell within run-to-run noise on a shared
/// machine, so the shape above is what the constants are fitted to.
///
/// The curve is sharply asymmetric, and that asymmetry is what this heuristic is
/// shaped around: one worker too few costs >2x — record output pins to exactly
/// one worker's 165 MiB/s — while extra workers past the balance point cost
/// nothing measurable, because they simply block on an empty queue. So round up,
/// and cap at [`MAX_COMPRESSION_WORKERS`].
///
/// The residual ~9% at `-p 8, C=3` is deliberately left on the table. Buying it
/// back means spending cores that, on a workload where mapping is not nearly
/// free, contribute more to mapping than to compression; `--bamCompressThreads`
/// is there for anyone who wants to make that trade explicitly.
///
/// A low estimate also degrades gracefully rather than catastrophically: the
/// chunk queue is bounded, so mapping workers block instead of growing memory.
fn compression_worker_count(mapping_threads: usize) -> NonZeroUsize {
    let mapping_threads = mapping_threads.max(1);
    let balanced =
        (mapping_threads as f64 * RECORD_MIB_PER_SEC_PER_MAPPER / DEFLATE_MIB_PER_SEC).ceil();
    // Never more workers than mapping threads: below that ratio the pipeline is
    // producer-limited no matter how much compression capacity is available.
    let workers = (balanced as usize).clamp(1, MAX_COMPRESSION_WORKERS.min(mapping_threads));
    NonZeroUsize::new(workers).expect("clamped to at least 1")
}

struct WriterFailure(Mutex<Option<String>>);

impl WriterFailure {
    fn record(&self, error: &io::Error) {
        let mut slot = self.0.lock().unwrap();
        if slot.is_none() {
            *slot = Some(error.to_string());
        }
    }

    fn error(&self) -> Option<io::Error> {
        self.0
            .lock()
            .unwrap()
            .as_ref()
            .map(|message| io::Error::other(message.clone()))
    }
}

pub struct BamOutput {
    sender: Option<Sender<Vec<u8>>>,
    recycled: Receiver<Vec<u8>>,
    failure: Arc<WriterFailure>,
    handle: Option<JoinHandle<io::Result<()>>>,
}

impl BamOutput {
    /// Open `path` and start the writer thread. `threads` is the mapping thread
    /// count (`-p`); `compress_threads` overrides the derived BGZF worker count
    /// (`--bamCompressThreads`) for benchmarking or unusual storage.
    pub fn create(
        path: &Path,
        salmon: &SalmonIndex,
        command: &str,
        threads: usize,
        compress_threads: Option<NonZeroUsize>,
    ) -> anyhow::Result<Self> {
        use anyhow::Context as _;
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("creating BAM output dir {}", parent.display()))?;
            }
        }
        let file = std::fs::File::create(path)
            .with_context(|| format!("creating BAM output {}", path.display()))?;
        let queue_capacity = (threads.max(1) * 2).clamp(4, 32);
        let (sender, receiver): (Sender<Vec<u8>>, Receiver<Vec<u8>>) =
            crossbeam_channel::bounded(queue_capacity);
        let (recycle_sender, recycled): (Sender<Vec<u8>>, Receiver<Vec<u8>>) =
            crossbeam_channel::bounded(queue_capacity + threads);
        let failure = Arc::new(WriterFailure(Mutex::new(None)));
        let writer_failure = Arc::clone(&failure);
        let header = encode_header(salmon, command)?;
        let compression_workers =
            compress_threads.unwrap_or_else(|| compression_worker_count(threads));
        tracing::debug!(
            mapping_threads = threads,
            bgzf_workers = compression_workers.get(),
            "BAM output pipeline"
        );
        let handle = std::thread::Builder::new()
            .name("salmon-bam-writer".into())
            .spawn(move || {
                let result = (|| {
                    let file = std::io::BufWriter::with_capacity(1 << 20, file);
                    let mut writer =
                        bgzf::io::MultithreadedWriter::with_worker_count(compression_workers, file);
                    writer.write_all(&header)?;
                    while let Ok(mut chunk) = receiver.recv() {
                        writer.write_all(&chunk)?;
                        chunk.clear();
                        if chunk.capacity() <= RETAIN_CAPACITY {
                            let _ = recycle_sender.try_send(chunk);
                        }
                    }
                    // `finish` hands back the inner writer after appending the
                    // BGZF EOF block; that block is still sitting in the
                    // `BufWriter`, and `BufWriter::drop` flushes only on a
                    // best-effort basis and *discards* any error. Flush it here
                    // so a failure at the very end of the run (ENOSPC, I/O
                    // error) surfaces instead of silently truncating the BAM.
                    let mut file = writer.finish()?;
                    file.flush()?;
                    Ok(())
                })();
                if let Err(error) = &result {
                    writer_failure.record(error);
                }
                result
            })?;
        Ok(Self {
            sender: Some(sender),
            recycled,
            failure,
            handle: Some(handle),
        })
    }

    /// A per-worker encode buffer bound to this output. Tying the borrow into
    /// the scratch makes "a scratch always has its writer" a type-level fact,
    /// so callers never have to re-derive it.
    pub fn scratch(&self) -> BamScratch<'_> {
        BamScratch {
            output: self,
            buffer: self.take_buffer(),
        }
    }

    /// A cleared chunk buffer: a recycled one if the writer has returned any,
    /// otherwise a fresh modest allocation that grows toward `CHUNK_TARGET`.
    fn take_buffer(&self) -> Vec<u8> {
        self.recycled
            .try_recv()
            .unwrap_or_else(|_| Vec::with_capacity(INITIAL_CAPACITY))
    }

    fn send(&self, chunk: Vec<u8>) -> io::Result<()> {
        if let Some(error) = self.failure.error() {
            return Err(error);
        }
        self.sender.as_ref().unwrap().send(chunk).map_err(|_| {
            self.failure
                .error()
                .unwrap_or_else(|| io::Error::other("BAM writer stopped"))
        })
    }

    pub fn finish(mut self) -> anyhow::Result<()> {
        drop(self.sender.take());
        self.handle
            .take()
            .unwrap()
            .join()
            .map_err(|_| anyhow::anyhow!("BAM writer thread panicked"))?
            .map_err(anyhow::Error::from)
    }
}

impl Drop for BamOutput {
    fn drop(&mut self) {
        drop(self.sender.take());
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

/// One mapping worker's private BAM encode buffer, bound to the output it feeds.
pub struct BamScratch<'a> {
    output: &'a BamOutput,
    buffer: Vec<u8>,
}

impl BamScratch<'_> {
    /// Encode every record for one fragment, then hand the chunk off if it has
    /// reached the target size.
    ///
    /// The size check deliberately happens *after* the whole fragment is
    /// encoded, never between its records: all records for a fragment are
    /// therefore contiguous in one chunk, and chunks are written whole. That is
    /// what keeps the output collated by fragment even though nothing else
    /// about the record order is constrained.
    pub fn write_fragment(
        &mut self,
        salmon: &SalmonIndex,
        r1_id: &[u8],
        r1_seq: &[u8],
        r2: Option<(&[u8], &[u8])>,
        maps: &[ScoredMapping],
    ) -> io::Result<()> {
        mapping_record::emit_fragment_records(salmon, r1_id, r1_seq, r2, maps, |record| {
            encode_record(&mut self.buffer, record)
        })?;
        if self.buffer.len() >= CHUNK_TARGET {
            self.flush()?;
        }
        Ok(())
    }

    /// Hand the accumulated chunk to the writer thread. Only the `Vec`
    /// descriptor moves; the backing allocation is not copied.
    pub fn flush(&mut self) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let chunk = std::mem::take(&mut self.buffer);
        self.output.send(chunk)?;
        self.buffer = self.output.take_buffer();
        Ok(())
    }
}

fn push_i32(buf: &mut Vec<u8>, value: i32) {
    buf.extend_from_slice(&value.to_le_bytes());
}

fn push_u32(buf: &mut Vec<u8>, value: u32) {
    buf.extend_from_slice(&value.to_le_bytes());
}

/// Encode the BAM header block: the shared SAM header text (byte-identical to
/// what `--writeMappings` writes) followed by the binary reference table.
fn encode_header(salmon: &SalmonIndex, command: &str) -> anyhow::Result<Vec<u8>> {
    let text = mapping_record::header_text(salmon, command);
    let mut buf = Vec::with_capacity(text.len() + salmon.num_refs() * 32);
    buf.extend_from_slice(b"BAM\x01");
    push_i32(&mut buf, i32::try_from(text.len())?);
    buf.extend_from_slice(text.as_bytes());
    push_i32(&mut buf, i32::try_from(salmon.num_refs())?);
    for tid in 0..salmon.num_refs() {
        let name = salmon.ref_name(tid).as_bytes();
        push_i32(&mut buf, i32::try_from(name.len() + 1)?);
        buf.extend_from_slice(name);
        buf.push(0);
        push_i32(
            &mut buf,
            i32::try_from(salmon.ref_len(tid))
                .map_err(|_| anyhow::anyhow!("reference {} exceeds BAM length limit", tid))?,
        );
    }
    Ok(buf)
}

fn reg2bin(start: i32, end: i32) -> u16 {
    let start = start.max(0) as u32;
    let end = end.max(start as i32 + 1) as u32 - 1;
    if start >> 14 == end >> 14 {
        (4681 + (start >> 14)) as u16
    } else if start >> 17 == end >> 17 {
        (585 + (start >> 17)) as u16
    } else if start >> 20 == end >> 20 {
        (73 + (start >> 20)) as u16
    } else if start >> 23 == end >> 23 {
        (9 + (start >> 23)) as u16
    } else if start >> 26 == end >> 26 {
        (1 + (start >> 26)) as u16
    } else {
        0
    }
}

fn base_code(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' => 8,
        _ => 15,
    }
}

fn record_base(record: &AlignmentRecord<'_>, index: usize) -> u8 {
    if record.reverse_complement {
        mapping_record::complement(record.sequence[record.sequence.len() - 1 - index])
    } else {
        record.sequence[index]
    }
}

fn push_aux_i64(buf: &mut Vec<u8>, tag: [u8; 2], value: i64) {
    buf.extend_from_slice(&tag);
    if let Ok(value) = u8::try_from(value) {
        buf.push(b'C');
        buf.push(value);
    } else if let Ok(value) = i8::try_from(value) {
        buf.push(b'c');
        buf.push(value as u8);
    } else if let Ok(value) = u16::try_from(value) {
        buf.push(b'S');
        buf.extend_from_slice(&value.to_le_bytes());
    } else if let Ok(value) = i16::try_from(value) {
        buf.push(b's');
        buf.extend_from_slice(&value.to_le_bytes());
    } else if let Ok(value) = u32::try_from(value) {
        buf.push(b'I');
        buf.extend_from_slice(&value.to_le_bytes());
    } else {
        buf.push(b'i');
        buf.extend_from_slice(&(value as i32).to_le_bytes());
    }
}

fn encode_record(buf: &mut Vec<u8>, record: &AlignmentRecord<'_>) -> io::Result<()> {
    let block_start = buf.len();
    push_u32(buf, 0);
    push_i32(buf, record.reference_id as i32);
    push_i32(buf, record.position.max(0));
    let name_len = record.name.len() + 1;
    let cigar_len = record.cigar.as_slice().len();
    let reference_span: i32 = record
        .cigar
        .as_slice()
        .iter()
        .filter(|op| matches!(op.kind, CigarKind::Match))
        .map(|op| op.len as i32)
        .sum();
    buf.push(u8::try_from(name_len).map_err(|_| io::Error::other("BAM read name too long"))?);
    buf.push(record.mapping_quality);
    buf.extend_from_slice(
        &reg2bin(record.position, record.position + reference_span).to_le_bytes(),
    );
    buf.extend_from_slice(&(cigar_len as u16).to_le_bytes());
    buf.extend_from_slice(&record.flags.to_le_bytes());
    push_i32(buf, record.sequence.len() as i32);
    push_i32(buf, record.mate_reference_id.map_or(-1, |id| id as i32));
    push_i32(buf, record.mate_position.unwrap_or(-1));
    push_i32(buf, record.template_length);
    buf.extend_from_slice(record.name);
    buf.push(0);
    for op in record.cigar.as_slice() {
        let code = match op.kind {
            CigarKind::Match => 0,
            CigarKind::SoftClip => 4,
        };
        push_u32(buf, ((op.len as u32) << 4) | code);
    }
    for i in (0..record.sequence.len()).step_by(2) {
        let high = base_code(record_base(record, i));
        let low = if i + 1 < record.sequence.len() {
            base_code(record_base(record, i + 1))
        } else {
            0
        };
        buf.push((high << 4) | low);
    }
    buf.resize(buf.len() + record.sequence.len(), 0xff);
    push_aux_i64(buf, *b"NH", record.nh as i64);
    push_aux_i64(buf, *b"HI", record.hi as i64);
    buf.extend_from_slice(b"XTA");
    buf.push(record.xt);
    push_aux_i64(buf, *b"AS", record.alignment_score as i64);
    let block_size = u32::try_from(buf.len() - block_start - 4)
        .map_err(|_| io::Error::other("BAM record too large"))?;
    buf[block_start..block_start + 4].copy_from_slice(&block_size.to_le_bytes());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compression_workers_track_the_measured_balance_point() {
        // One deflate worker per ~3 mapping threads (52 / 165 MiB/s).
        for (threads, expected) in [(1, 1), (2, 1), (3, 1), (4, 2), (8, 3), (16, 6), (24, 8)] {
            assert_eq!(
                compression_worker_count(threads).get(),
                expected,
                "-p {threads}"
            );
        }
    }

    #[test]
    fn compression_workers_are_bounded() {
        // Never zero, never past the cap, and never more than the mapping
        // threads they are meant to keep up with.
        for threads in [0usize, 1, 7, 32, 64, 256, usize::MAX] {
            let workers = compression_worker_count(threads).get();
            assert!(workers >= 1);
            assert!(workers <= MAX_COMPRESSION_WORKERS, "-p {threads}");
            assert!(workers <= threads.max(1), "-p {threads}");
        }
    }
}
