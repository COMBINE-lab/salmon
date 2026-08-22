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

use anyhow::Context;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use libradicl::chunk::ChunkBuf;
use libradicl::header::RadPrelude;
use libradicl::rad_types::TagValue;
use libradicl::writers::{ConcurrentChunkWriter, RadFileWriter};
use libradicl::ChunkCodec;

use crate::record::{SalmonBulkContext, SalmonBulkRecord};
use crate::{schema, RadProfile, BAKED_ABUND, BAKED_FLD, BAKED_LIBFMT, BAKED_SCORE_KIND};

/// Shared, thread-safe RAD output sink for one file.
pub struct RadOutputWriter {
    /// where the finished file should end up
    final_path: PathBuf,
    /// where bytes are actually written until `finalize` renames it
    partial_path: PathBuf,
    /// The underlying file, wrapped so that appends from many threads are safe.
    ccw: ConcurrentChunkWriter<BufWriter<File>>,
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
    /// Mapping-pass counters, known only once the pass ends.
    pending_counters: Option<crate::MapCounters>,
    /// Whether index identity was written at create time. The tags are written
    /// up front but `baked_flags` is only patched at finalize, so this carries
    /// the fact across; without it a reader would ignore the tags it can see.
    has_index_prov: bool,
    /// chunk compression codec applied by worker buffers (advertised in the header).
    codec: ChunkCodec,
    /// Bytes appended since the last free-disk-space probe (all workers).
    ///
    /// A filling disk should fail in seconds with an actionable message, not
    /// at the end of a mapping pass measured in minutes (#1140): every
    /// [`SPACE_CHECK_INTERVAL`] bytes, the writer probes the target
    /// filesystem and errors out early when free space falls under
    /// [`LOW_SPACE_BYTES`]. The counter reset races benignly between workers —
    /// it paces the probe, it does not gate correctness.
    bytes_since_space_check: std::sync::atomic::AtomicU64,
}

/// Probe the target filesystem's free space every this many appended bytes.
const SPACE_CHECK_INTERVAL: u64 = 256 * 1024 * 1024;
/// Fail the run when the target filesystem's free space falls below this while
/// chunks are still being appended.
const LOW_SPACE_BYTES: u64 = 512 * 1024 * 1024;

/// The sibling path a RAD file is written to before it is finalized.
///
/// Keeping it next to the destination (rather than in a temp dir) means the
/// final rename is within one filesystem, so it is atomic; a cross-device rename
/// would degrade to a copy and reintroduce the partial-file window this exists
/// to close. The pid keeps concurrent runs targeting the same output from
/// colliding.
fn partial_path_for(path: &Path) -> PathBuf {
    let mut name = path.file_name().unwrap_or_default().to_os_string();
    name.push(format!(".partial-{}", std::process::id()));
    path.with_file_name(name)
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
        provenance: &crate::WriterProvenance,
    ) -> anyhow::Result<Self> {
        let (prelude, file_tag_map): (RadPrelude, _) = schema::build_prelude(
            profile,
            is_paired,
            ref_names,
            ref_lengths,
            fld_len,
            codec,
            provenance,
        );
        // Write to a sibling temporary path and rename onto `path` only once the
        // run has finalized successfully. A partial RAD therefore never appears
        // under the name the user asked for, so a later `--rad` run cannot pick
        // one up believing it is complete. Rename rather than delete-on-error,
        // because rename is atomic and also covers `SIGKILL`, where no cleanup
        // code of ours would run at all.
        let partial_path = partial_path_for(path);
        // `BufWriter` batches small writes into large ones; without it every
        // chunk append would become a separate system call.
        let file = BufWriter::new(
            File::create(&partial_path)
                .with_context(|| format!("creating RAD output file {}", partial_path.display()))?,
        );
        // Writing the prelude here means the file is valid (if empty) from the
        // moment it is created.
        let fw = RadFileWriter::new(file, &prelude, &file_tag_map)?;
        Ok(Self {
            final_path: path.to_path_buf(),
            partial_path,
            ccw: ConcurrentChunkWriter::new(fw),
            ctx: SalmonBulkContext { profile },
            fld_len,
            num_refs: ref_names.len(),
            // Nothing is baked until the run finishes.
            pending_fld: None,
            pending_abund: None,
            pending_libfmt: None,
            pending_score_kind: None,
            pending_counters: None,
            has_index_prov: provenance.index.is_some(),
            codec,
            bytes_since_space_check: std::sync::atomic::AtomicU64::new(0),
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
    /// Record what the mapping pass observed, to be baked at finalize.
    ///
    /// Without this a requant cannot report a mapping rate: the file holds only
    /// the fragments that mapped, so it would report 100% by construction.
    pub fn set_map_counters(&mut self, counters: crate::MapCounters) {
        self.pending_counters = Some(counters);
    }

    pub fn set_score_kind(&mut self, score_kind: u8) {
        self.pending_score_kind = Some(score_kind);
    }

    /// Append a fully-framed chunk (from [`FragmentChunkBuf::take_bytes`]) to the
    /// file. Thread-safe: takes `&self`, so every worker can call it directly.
    pub fn append_chunk_bytes(&self, bytes: &[u8]) -> anyhow::Result<()> {
        self.ccw.append_chunk_bytes(bytes).with_context(|| {
            format!(
                "writing RAD chunk to {}\n  \
                 the incomplete file has been left there; it is NOT usable as `--rad` input",
                self.partial_path.display()
            )
        })?;
        use std::sync::atomic::Ordering;
        let since = self
            .bytes_since_space_check
            .fetch_add(bytes.len() as u64, Ordering::Relaxed)
            + bytes.len() as u64;
        if since >= SPACE_CHECK_INTERVAL {
            self.bytes_since_space_check.store(0, Ordering::Relaxed);
            let dir = self
                .partial_path
                .parent()
                .unwrap_or(std::path::Path::new("."));
            if let Some(free) = salmon_core::free_disk_bytes(dir) {
                if free < LOW_SPACE_BYTES {
                    anyhow::bail!(
                        "running out of disk space while writing {} ({} MiB free in {}); \
                         stopping now rather than at the end of the pass. \
                         Consider `--radScratchDir <dir>` to place the intermediate on a \
                         larger/faster volume, or `--radCompress zstd` (~30% smaller).",
                        self.partial_path.display(),
                        free / (1024 * 1024),
                        dir.display()
                    );
                }
            }
        }
        Ok(())
    }

    /// Backpatch any baked FLD / abundances + the `baked_flags` marker, then flush
    /// and finalize (backpatching the chunk count).
    ///
    /// Takes `self` by value, so the type system enforces that nothing can be
    /// appended afterwards. Call after all workers have flushed and dropped their
    /// references.
    pub fn finalize(self) -> anyhow::Result<()> {
        // Move the fields out before consuming `self.ccw` below.
        let fld_len = self.fld_len;
        let num_refs = self.num_refs;
        let pending_fld = self.pending_fld;
        let pending_abund = self.pending_abund;
        let pending_libfmt = self.pending_libfmt;
        let pending_score_kind = self.pending_score_kind;
        let pending_counters = self.pending_counters;
        let has_index_prov = self.has_index_prov;
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
        if let Some(c) = pending_counters {
            for (tag, value) in [
                (crate::NUM_PROCESSED_TAG, c.num_processed),
                (crate::NUM_DOVETAIL_TAG, c.num_dovetail),
                (crate::NUM_FILTERED_VM_TAG, c.num_filtered_vm),
                (crate::NUM_BELOW_THRESH_VM_TAG, c.num_below_threshold_vm),
                (crate::NUM_DECOY_FRAGMENTS_TAG, c.num_decoy_fragments),
            ] {
                w.backpatch_file_tag_value(tag, &TagValue::U64(value))?;
            }
            flags |= crate::BAKED_MAP_COUNTERS;
        }
        if has_index_prov {
            flags |= crate::BAKED_INDEX_PROV;
        }
        // Reaching here means every chunk was written and the header is about to
        // be patched, so record completeness alongside whatever else was baked.
        // See `WRITE_COMPLETE`: readers treat this as confirmatory, never required.
        flags |= crate::WRITE_COMPLETE;
        w.backpatch_file_tag_value(crate::BAKED_FLAGS_TAG, &TagValue::U8(flags))?;
        // Writes the final chunk count into the header and flushes the file.
        let buf = w
            .finalize()
            .with_context(|| format!("finalizing {}", self.partial_path.display()))?;
        // `into_inner` flushes, and unlike dropping the `BufWriter` it reports a
        // failure to do so rather than swallowing it.
        let file = buf.into_inner().map_err(|e| {
            anyhow::anyhow!(
                "flushing {}: {}",
                self.partial_path.display(),
                e.into_error()
            )
        })?;
        // Flush to the device before the rename makes this the user-visible
        // output. On NFS in particular a deferred `ENOSPC` surfaces only here, so
        // dropping the handle instead would discard exactly the error this fix
        // exists to report (salmon#1105).
        file.sync_all()
            .with_context(|| format!("flushing {} to disk", self.partial_path.display()))?;
        // Only now does the output take the name the user asked for. Any failure
        // before this point leaves the bytes under the `.partial-*` name, where a
        // later `--rad` run will not mistake them for a complete file.
        std::fs::rename(&self.partial_path, &self.final_path).with_context(|| {
            format!(
                "renaming {} to {}",
                self.partial_path.display(),
                self.final_path.display()
            )
        })?;
        Ok(())
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
