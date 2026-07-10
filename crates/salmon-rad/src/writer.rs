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
use libradicl::rad_types::TagValue;
use libradicl::writers::{ConcurrentChunkWriter, RadFileWriter};
use libradicl::ChunkCodec;

use crate::record::{SalmonBulkContext, SalmonBulkRecord};
use crate::{schema, RadProfile, BAKED_ABUND, BAKED_FLD, BAKED_LIBFMT, BAKED_SCORE_KIND};

/// Shared, thread-safe RAD output sink for one file.
pub struct RadOutputWriter {
    ccw: ConcurrentChunkWriter<BufWriter<File>>,
    ctx: SalmonBulkContext,
    /// reserved-slot lengths, so baked values are padded/truncated to fit exactly.
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
        let file = BufWriter::new(File::create(path)?);
        let fw = RadFileWriter::new(file, &prelude, &file_tag_map)?;
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
    /// with this codec (see [`FragmentChunkBuf::with_capacity_codec`]).
    pub fn codec(&self) -> ChunkCodec {
        self.codec
    }

    /// The parsing/serialization context (profile) for records in this file.
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
    /// file. Thread-safe.
    pub fn append_chunk_bytes(&self, bytes: &[u8]) -> anyhow::Result<()> {
        self.ccw.append_chunk_bytes(bytes)
    }

    /// Backpatch any baked FLD / abundances + the `baked_flags` marker, then flush
    /// and finalize (backpatching the chunk count). Call after all workers have
    /// flushed and dropped their references.
    pub fn finalize(self) -> anyhow::Result<()> {
        let fld_len = self.fld_len;
        let num_refs = self.num_refs;
        let pending_fld = self.pending_fld;
        let pending_abund = self.pending_abund;
        let pending_libfmt = self.pending_libfmt;
        let pending_score_kind = self.pending_score_kind;
        let mut w = self.ccw.into_writer()?;
        let mut flags: u8 = 0;
        if let Some(mut pmf) = pending_fld {
            // pad with -inf (log 0) / truncate to the reserved slot length.
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
        w.finalize()?;
        Ok(())
    }
}

/// Per-thread accumulator of salmon RAD records, flushed as one chunk.
pub struct FragmentChunkBuf {
    buf: ChunkBuf,
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

    /// Number of records accumulated so far.
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
        std::mem::replace(&mut self.buf, ChunkBuf::with_capacity(self.cap))
            .into_bytes_with_codec(self.codec)
    }
}

impl Default for FragmentChunkBuf {
    fn default() -> Self {
        Self::with_capacity(64 * 1024)
    }
}
