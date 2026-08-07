//! The salmon RAD output writer.
//!
//! # The concurrency shape
//!
//! Mapping runs on many threads at once, and all of them produce records for the
//! same output file. Having every thread lock the file per record would serialize
//! the whole run, so the work is split in two:
//!
//! * each worker thread owns a private [`FragmentChunkBuf`] and serializes its
//!   records into it with no locking at all;
//! * when that buffer is big enough the worker turns it into one self-contained,
//!   already-compressed RAD *chunk* and hands the finished bytes to the shared
//!   [`RadOutputWriter`], which appends it under a brief lock.
//!
//! So the expensive parts (serialization and compression) are fully parallel, and
//! the shared writer only ever does a bulk append.
//!
//! A single [`RadOutputWriter`] is created per output file; it writes the
//! prelude and file-tag values up front and exposes a thread-safe
//! [`append_chunk_bytes`](RadOutputWriter::append_chunk_bytes) for worker
//! threads. [`finalize`](RadOutputWriter::finalize) backpatches the chunk count.
//!
//! Because chunks land in whatever order the threads finish, chunk order in the
//! file is not deterministic — which is fine, since a reader treats chunks as an
//! unordered bag and salmon's determinism comes from order-independent
//! aggregation downstream, not from file order.

use std::fs::File;
use std::io::{BufWriter, Cursor, Seek, SeekFrom, Write};
use std::path::Path;

use libradicl::chunk::ChunkBuf;
use libradicl::header::RadPrelude;
use libradicl::rad_types::TagValue;
use libradicl::writers::{ConcurrentChunkWriter, RadFileWriter};
use libradicl::ChunkCodec;

use crate::record::{SalmonBulkContext, SalmonBulkRecord};
use crate::{schema, RadProfile, BAKED_ABUND, BAKED_FLD, BAKED_LIBFMT, BAKED_SCORE_KIND};

/// Where a [`RadOutputWriter`]'s bytes go.
///
/// `libradicl`'s writer is already generic over `Write + Seek` (it needs `Seek`
/// for the header backpatching at finalize), so supporting an in-memory image
/// costs only this dispatch. It is an enum rather than a type parameter
/// deliberately: a generic would leak through every worker that holds a
/// `&RadOutputWriter`, for a choice made once at construction.
enum Sink {
    File(BufWriter<File>),
    Memory(Cursor<Vec<u8>>),
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            Sink::File(w) => w.write(buf),
            Sink::Memory(w) => w.write(buf),
        }
    }
    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            Sink::File(w) => w.flush(),
            Sink::Memory(w) => w.flush(),
        }
    }
}

impl Seek for Sink {
    fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        match self {
            Sink::File(w) => w.seek(pos),
            Sink::Memory(w) => w.seek(pos),
        }
    }
}

/// Shared, thread-safe RAD output sink.
pub struct RadOutputWriter {
    /// The underlying sink, wrapped so that appends from many threads are safe.
    ccw: ConcurrentChunkWriter<Sink>,
    /// Which profile records are written in (decides the per-hit layout).
    ctx: SalmonBulkContext,
    /// reserved-slot lengths, so baked values are padded/truncated to fit exactly.
    ///
    /// Backpatching overwrites bytes in place, so a baked array must be exactly
    /// as long as the placeholder reserved in the header.
    fld_len: usize,
    num_refs: usize,
    /// values to backpatch at finalize (set after the pass; `None` ⇒ placeholder).
    pending_fld: Option<Vec<f64>>,
    pending_abund: Option<Vec<f64>>,
    /// resolved library format as a `format_id` (see [`crate::LIBRARY_FORMAT_TAG`]).
    pending_libfmt: Option<u8>,
    /// non-default score interpretation (see [`crate::SCORE_KIND_TAG`]); only the
    /// scored (SA) profile reserves the slot, so set only in that profile.
    pending_score_kind: Option<u8>,
    /// chunk compression codec applied by worker buffers (advertised in the header).
    codec: ChunkCodec,
}

impl RadOutputWriter {
    /// Create a RAD file at `path`, writing the prelude + file tags for the
    /// given profile. `ref_names`/`ref_lengths` cover the full reference table.
    /// `fld_len` reserves a fixed-size fragment-length-distribution slot (log-PMF
    /// over raw lengths `[0, fld_len)`) to be backpatched at finalize.
    pub fn create(
        path: &Path,
        ref_names: &[&str],
        ref_lengths: &[u32],
        is_paired: bool,
        profile: RadProfile,
        fld_len: usize,
        codec: ChunkCodec,
    ) -> anyhow::Result<Self> {
        let (prelude, file_tag_map): (RadPrelude, _) =
            schema::build_prelude(profile, is_paired, ref_names, ref_lengths, fld_len, codec);
        // `BufWriter` batches small writes into large ones; without it every
        // chunk append would become a separate system call.
        let file = Sink::File(BufWriter::new(File::create(path)?));
        // Writing the prelude here means the file is valid (if empty) from the
        // moment it is created.
        let fw = RadFileWriter::new(file, &prelude, &file_tag_map)?;
        Ok(Self {
            ccw: ConcurrentChunkWriter::new(fw),
            ctx: SalmonBulkContext { profile },
            fld_len,
            num_refs: ref_names.len(),
            // Nothing is baked until the run finishes.
            pending_fld: None,
            pending_abund: None,
            pending_libfmt: None,
            pending_score_kind: None,
            codec,
        })
    }

    /// Create a RAD image held entirely in memory.
    ///
    /// Same prelude, same records, same backpatching as [`create`](Self::create)
    /// — only the destination differs, so the bytes are indistinguishable from a
    /// file's and any reader accepts them.
    ///
    /// This exists for the `--deterministic` pipeline, which writes an
    /// intermediate RAD, reads it back (up to five separate passes, depending on
    /// what has to be derived), then deletes it. Keeping it in memory removes a
    /// full write and every one of those re-reads. The trade is explicit and the
    /// caller makes it: the whole image is resident instead of on disk.
    #[allow(clippy::too_many_arguments)]
    pub fn create_in_memory(
        ref_names: &[&str],
        ref_lengths: &[u32],
        is_paired: bool,
        profile: RadProfile,
        fld_len: usize,
        codec: ChunkCodec,
    ) -> anyhow::Result<Self> {
        let (prelude, file_tag_map): (RadPrelude, _) =
            schema::build_prelude(profile, is_paired, ref_names, ref_lengths, fld_len, codec);
        let fw = RadFileWriter::new(
            Sink::Memory(Cursor::new(Vec::new())),
            &prelude,
            &file_tag_map,
        )?;
        Ok(Self {
            ccw: ConcurrentChunkWriter::new(fw),
            ctx: SalmonBulkContext { profile },
            fld_len,
            num_refs: ref_names.len(),
            pending_fld: None,
            pending_abund: None,
            pending_libfmt: None,
            pending_score_kind: None,
            codec,
        })
    }

    /// The chunk compression codec for this file. Worker buffers must compress
    /// with this codec (see [`FragmentChunkBuf::with_capacity_codec`]), because
    /// the header advertises exactly one codec for the whole file.
    pub fn codec(&self) -> ChunkCodec {
        self.codec
    }

    /// The parsing/serialization context (profile) for records in this file.
    /// Workers need it to serialize records the same way the header promises.
    pub fn context(&self) -> &SalmonBulkContext {
        &self.ctx
    }

    /// Bake the fragment-length distribution (a log-PMF over raw lengths) into the
    /// header at finalize, so a reader can quantify in a single pass with exact
    /// FLD parity. Padded/truncated to the reserved `fld_len`.
    pub fn set_frag_length_dist(&mut self, log_pmf: &[f64]) {
        self.pending_fld = Some(log_pmf.to_vec());
    }

    /// Bake initial per-reference abundance estimates into the header at finalize
    /// (a prior for future bias-aware requant). Padded/truncated to `num_refs`.
    pub fn set_initial_abundances(&mut self, abundances: &[f64]) {
        self.pending_abund = Some(abundances.to_vec());
    }

    /// Bake the resolved library format (a `format_id`, see
    /// [`crate::LIBRARY_FORMAT_TAG`]) into the header at finalize, so a reader can
    /// apply concordance filtering under `-l A` without re-inferring the type.
    pub fn set_library_format(&mut self, format_id: u8) {
        self.pending_libfmt = Some(format_id);
    }

    /// Bake the per-hit score interpretation (see [`crate::SCORE_KIND_TAG`]) into
    /// the header at finalize. Only valid on the scored (selective-alignment)
    /// profile, whose prelude reserves the slot; a no-op-worthy default
    /// ([`crate::SCORE_KIND_AS`]) need not be set (absent ⇒ that default).
    pub fn set_score_kind(&mut self, score_kind: u8) {
        self.pending_score_kind = Some(score_kind);
    }

    /// Append a fully-framed chunk (from [`FragmentChunkBuf::take_bytes`]) to the
    /// file. Thread-safe: takes `&self`, so every worker can call it directly.
    pub fn append_chunk_bytes(&self, bytes: &[u8]) -> anyhow::Result<()> {
        self.ccw.append_chunk_bytes(bytes)
    }

    /// Backpatch any baked FLD / abundances + the `baked_flags` marker, then flush
    /// and finalize (backpatching the chunk count).
    ///
    /// Takes `self` by value, so the type system enforces that nothing can be
    /// appended afterwards. Call after all workers have flushed and dropped their
    /// references.
    pub fn finalize(self) -> anyhow::Result<()> {
        self.finalize_sink()?;
        Ok(())
    }

    /// Finalize and take the completed image, for a writer built by
    /// [`create_in_memory`](Self::create_in_memory).
    ///
    /// Returns an error if this writer is backed by a file — the caller asked for
    /// bytes from something that put them on disk, which is a programming error
    /// rather than a runtime condition.
    pub fn finalize_into_bytes(self) -> anyhow::Result<Vec<u8>> {
        match self.finalize_sink()? {
            Sink::Memory(c) => Ok(c.into_inner()),
            Sink::File(_) => {
                anyhow::bail!("finalize_into_bytes called on a file-backed RAD writer")
            }
        }
    }

    /// Backpatch the baked slots, finalize, and hand back the sink.
    fn finalize_sink(self) -> anyhow::Result<Sink> {
        // Move the fields out before consuming `self.ccw` below.
        let fld_len = self.fld_len;
        let num_refs = self.num_refs;
        let pending_fld = self.pending_fld;
        let pending_abund = self.pending_abund;
        let pending_libfmt = self.pending_libfmt;
        let pending_score_kind = self.pending_score_kind;
        // Unwrap the concurrent wrapper back into a plain writer; this succeeds
        // only once every worker handle is gone, which is what makes the
        // in-place patching below safe.
        let mut w = self.ccw.into_writer()?;
        // Accumulates which slots we actually filled, recorded in the header so a
        // reader can tell real data from an untouched placeholder.
        let mut flags: u8 = 0;
        if let Some(mut pmf) = pending_fld {
            // pad with -inf (log 0) / truncate to the reserved slot length.
            // -inf is the log of probability zero, so padding adds impossible
            // lengths rather than uniformly likely ones.
            pmf.resize(fld_len, f64::NEG_INFINITY);
            w.backpatch_file_tag_value(crate::FRAG_LENGTH_DIST_TAG, &TagValue::ArrayF64(pmf))?;
            flags |= BAKED_FLD;
        }
        if let Some(mut abund) = pending_abund {
            abund.resize(num_refs, 0.0);
            w.backpatch_file_tag_value(crate::INITIAL_ABUNDANCES_TAG, &TagValue::ArrayF64(abund))?;
            flags |= BAKED_ABUND;
        }
        if let Some(fmt) = pending_libfmt {
            w.backpatch_file_tag_value(crate::LIBRARY_FORMAT_TAG, &TagValue::U8(fmt))?;
            flags |= BAKED_LIBFMT;
        }
        // Only bake a non-default score kind; the reserved slot already holds the
        // AS default, so leaving it untouched keeps older readers' behavior.
        if let Some(k) = pending_score_kind {
            if k != crate::SCORE_KIND_AS {
                w.backpatch_file_tag_value(crate::SCORE_KIND_TAG, &TagValue::U8(k))?;
                flags |= BAKED_SCORE_KIND;
            }
        }
        if flags != 0 {
            w.backpatch_file_tag_value(crate::BAKED_FLAGS_TAG, &TagValue::U8(flags))?;
        }
        // Writes the final chunk count into the header and flushes the sink.
        w.finalize()
    }
}

/// Per-thread accumulator of salmon RAD records, flushed as one chunk.
///
/// Private to one thread, so none of this needs locking.
pub struct FragmentChunkBuf {
    buf: ChunkBuf,
    /// Capacity to give the replacement buffer after a flush, so a long run
    /// stops reallocating after the first chunk.
    cap: usize,
    codec: ChunkCodec,
}

impl FragmentChunkBuf {
    /// Create a per-thread buffer with the given initial byte capacity and no
    /// chunk compression.
    pub fn with_capacity(cap: usize) -> Self {
        Self::with_capacity_codec(cap, ChunkCodec::None)
    }

    /// Create a per-thread buffer that compresses each flushed chunk with
    /// `codec` (must match the codec advertised by the owning [`RadOutputWriter`]).
    pub fn with_capacity_codec(cap: usize, codec: ChunkCodec) -> Self {
        Self {
            buf: ChunkBuf::with_capacity(cap),
            cap,
            codec,
        }
    }

    /// Serialize one fragment's record into the buffer.
    pub fn write(&mut self, rec: &SalmonBulkRecord, ctx: &SalmonBulkContext) -> anyhow::Result<()> {
        self.buf.write_record(rec, ctx)
    }

    /// Number of records accumulated so far. Callers use it, with
    /// [`Self::byte_len`], to decide when a chunk is worth flushing.
    pub fn nrec(&self) -> u32 {
        self.buf.nrec()
    }

    /// Number of record bytes accumulated so far (excluding the chunk header).
    pub fn byte_len(&self) -> usize {
        self.buf.byte_len()
    }

    /// Take the accumulated records as a single framed chunk (`nbytes`+`nrec`
    /// header prepended, payload compressed per the configured codec) and reset
    /// the buffer for reuse. Fails only if compression fails (e.g. a zstd codec
    /// in a libradicl built without the `zstd` feature).
    pub fn take_bytes(&mut self) -> anyhow::Result<Vec<u8>> {
        // `mem::replace` swaps a fresh buffer into place and hands back the old
        // one, so the caller gets ownership of the finished bytes without
        // copying and this buffer is immediately usable again.
        std::mem::replace(&mut self.buf, ChunkBuf::with_capacity(self.cap))
            .into_bytes_with_codec(self.codec)
    }
}

impl Default for FragmentChunkBuf {
    /// 64 KiB: large enough that chunk framing and lock overhead are negligible,
    /// small enough that many worker threads' buffers stay cheap.
    fn default() -> Self {
        Self::with_capacity(64 * 1024)
    }
}
