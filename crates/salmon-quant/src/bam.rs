//! Chunked `--writeBam` output.

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
const RETAIN_CAPACITY: usize = 8 * 1024 * 1024;

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
    pub fn create(
        path: &Path,
        salmon: &SalmonIndex,
        command: &str,
        threads: usize,
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
        let compression_workers = NonZeroUsize::new((threads / 4).clamp(1, 4)).unwrap();
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
                    writer.finish()?;
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

    pub fn scratch(&self) -> BamScratch {
        BamScratch {
            buffer: self
                .recycled
                .try_recv()
                .unwrap_or_else(|_| Vec::with_capacity(INITIAL_CAPACITY)),
        }
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

pub struct BamScratch {
    buffer: Vec<u8>,
}

impl BamScratch {
    pub fn write_fragment(
        &mut self,
        output: &BamOutput,
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
            self.flush(output)?;
        }
        Ok(())
    }

    pub fn flush(&mut self, output: &BamOutput) -> io::Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let chunk = std::mem::take(&mut self.buffer);
        output.send(chunk)?;
        self.buffer = output
            .recycled
            .try_recv()
            .unwrap_or_else(|_| Vec::with_capacity(INITIAL_CAPACITY));
        Ok(())
    }
}

fn push_i32(buf: &mut Vec<u8>, value: i32) {
    buf.extend_from_slice(&value.to_le_bytes());
}

fn push_u32(buf: &mut Vec<u8>, value: u32) {
    buf.extend_from_slice(&value.to_le_bytes());
}

fn encode_header(salmon: &SalmonIndex, command: &str) -> anyhow::Result<Vec<u8>> {
    use std::fmt::Write as _;
    let mut text = String::from("@HD\tVN:1.6\tSO:unknown\n");
    for tid in 0..salmon.num_refs() {
        writeln!(
            text,
            "@SQ\tSN:{}\tLN:{}",
            salmon.ref_name(tid),
            salmon.ref_len(tid)
        )?;
    }
    writeln!(
        text,
        "@PG\tID:salmon\tPN:salmon\tVN:{}\tCL:{}",
        crate::output::SALMON_VERSION,
        command
    )?;
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
