//! The salmon RAD output writer.
//!
//! A single [`RadOutputWriter`] is created per output file; it writes the
//! prelude + file-tag values up front and exposes a thread-safe
//! [`append_chunk_bytes`](RadOutputWriter::append_chunk_bytes) for worker
//! threads. Each worker accumulates records into a per-thread
//! [`FragmentChunkBuf`] and flushes it (one RAD chunk per flush) to the shared
//! writer. [`finalize`](RadOutputWriter::finalize) backpatches the chunk count.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use libradicl::chunk::ChunkBuf;
use libradicl::header::RadPrelude;
use libradicl::writers::{ConcurrentChunkWriter, RadFileWriter};

use crate::record::{SalmonBulkContext, SalmonBulkRecord};
use crate::{schema, RadProfile};

/// Shared, thread-safe RAD output sink for one file.
pub struct RadOutputWriter {
    ccw: ConcurrentChunkWriter<BufWriter<File>>,
    ctx: SalmonBulkContext,
}

impl RadOutputWriter {
    /// Create a RAD file at `path`, writing the prelude + file tags for the
    /// given profile. `ref_names`/`ref_lengths` cover the full reference table.
    pub fn create(
        path: &Path,
        ref_names: &[&str],
        ref_lengths: &[u32],
        is_paired: bool,
        profile: RadProfile,
    ) -> anyhow::Result<Self> {
        let (prelude, file_tag_map): (RadPrelude, _) =
            schema::build_prelude(profile, is_paired, ref_names, ref_lengths);
        let file = BufWriter::new(File::create(path)?);
        let fw = RadFileWriter::new(file, &prelude, &file_tag_map)?;
        Ok(Self {
            ccw: ConcurrentChunkWriter::new(fw),
            ctx: SalmonBulkContext { profile },
        })
    }

    /// The parsing/serialization context (profile) for records in this file.
    pub fn context(&self) -> &SalmonBulkContext {
        &self.ctx
    }

    /// Append a fully-framed chunk (from [`FragmentChunkBuf::take_bytes`]) to the
    /// file. Thread-safe.
    pub fn append_chunk_bytes(&self, bytes: &[u8]) -> anyhow::Result<()> {
        self.ccw.append_chunk_bytes(bytes)
    }

    /// Flush and finalize the file, backpatching the chunk count. Call after all
    /// workers have flushed and dropped their references.
    pub fn finalize(self) -> anyhow::Result<()> {
        self.ccw.finalize()?;
        Ok(())
    }
}

/// Per-thread accumulator of salmon RAD records, flushed as one chunk.
pub struct FragmentChunkBuf {
    buf: ChunkBuf,
    cap: usize,
}

impl FragmentChunkBuf {
    /// Create a per-thread buffer with the given initial byte capacity.
    pub fn with_capacity(cap: usize) -> Self {
        Self {
            buf: ChunkBuf::with_capacity(cap),
            cap,
        }
    }

    /// Serialize one fragment's record into the buffer.
    pub fn write(&mut self, rec: &SalmonBulkRecord, ctx: &SalmonBulkContext) -> anyhow::Result<()> {
        self.buf.write_record(rec, ctx)
    }

    /// Number of records accumulated so far.
    pub fn nrec(&self) -> u32 {
        self.buf.nrec()
    }

    /// Number of record bytes accumulated so far (excluding the chunk header).
    pub fn byte_len(&self) -> usize {
        self.buf.byte_len()
    }

    /// Take the accumulated records as a single framed chunk (`nbytes`+`nrec`
    /// header prepended) and reset the buffer for reuse.
    pub fn take_bytes(&mut self) -> Vec<u8> {
        std::mem::replace(&mut self.buf, ChunkBuf::with_capacity(self.cap)).into_bytes()
    }
}

impl Default for FragmentChunkBuf {
    fn default() -> Self {
        Self::with_capacity(64 * 1024)
    }
}
