//! RAD-format mapping store for the deterministic two-pass path (`--deterministic`).
//!
//! Pass 1 maps in parallel and writes each mapped fragment's `(key, mappings)`
//! here as a RAD record; pass 2 reads them back and replays inference in a fixed
//! key-sorted order. We reuse the RAD container and tag framework (the format
//! used by piscem, piscem-infer, and alevin-fry, via `libradicl`'s type
//! primitives) rather than a bespoke binary format, so this store shares an
//! ingestion route with those tools and lays the groundwork for RAD-format
//! *input* to salmon-quant (`piscem map` -> `salmon quant`).
//!
//! ## Profile
//!
//! Standard bulk RAD keeps `(ref | orientation, position, frag_length)` per
//! alignment. salmon's selective-alignment inference also needs the per-mapping
//! coverage weight and a few more per-mapping fields, so we declare a
//! salmon-specific alignment-level tag schema that is a superset of the bulk
//! profile. A RAD that lacks the extra tags (for example one produced by
//! `piscem map`) maps onto the uniform-weight path instead; consuming that is a
//! separate, future feature.
//!
//! Only the fields the per-fragment inference ([`crate::processor`]'s `record`)
//! reads are stored; the SAM-output-only fields of [`ScoredMapping`] (`score`,
//! `r1_pos`, `r2_pos`, `r2_fw`, `r1_score`) are written during pass 1 and are not
//! needed in pass 2, so they are reconstructed at their inert defaults on read.
//!
//! ## Byte layout (RAD-conformant, little-endian)
//!
//! - prelude: `is_paired (u8)`, `ref_count (u64)`, `ref_count x (name_len u16 +
//!   name bytes)`, `num_chunks (u64)`, then the file / read / alignment
//!   [`TagSection`]s (written via `libradicl`).
//! - one or more chunks: `nbytes (u32, incl. this 8-byte header)`, `nrec (u32)`,
//!   then `nrec` records.
//! - record: `num_aln (u32)`, the read-level `frag_key (u128)`, then per
//!   alignment the alignment-level tag values in the declared order.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::Mutex;

use anyhow::{bail, Context, Result};
use libradicl::rad_types::{RadFloatId, RadIntId, RadType, TagDesc, TagSection, TagSectionLabel};

use salmon_core::{LibraryFormat, MateStatus};
use salmon_index::SalmonIndex;
use salmon_map::ScoredMapping;

/// Sentinel `format` tag value standing for `ScoredMapping::format == None`
/// (no library format determinable). Library `format_id`s are small, so the top
/// `u8` value is free to use as the "absent" marker.
const FORMAT_NONE: u8 = 0xFF;

/// The read-level tag section: the stable per-fragment sort key.
fn read_tag_section() -> TagSection {
    let mut s = TagSection::new_with_label(TagSectionLabel::ReadTags);
    s.add_tag_desc(TagDesc {
        name: "frag_key".to_string(),
        typeid: RadType::Int(RadIntId::U128),
    });
    s
}

/// The alignment-level tag section: salmon's full-fidelity per-mapping fields.
/// Signed values (`fragment_len`, `read_len`, `*_pos`, which use `-1` sentinels)
/// are stored as their `u32` two's-complement bit pattern (RAD ints are
/// unsigned) and cast back on read.
fn aln_tag_section() -> TagSection {
    let mut s = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
    let int = |name: &str, id| TagDesc {
        name: name.to_string(),
        typeid: RadType::Int(id),
    };
    s.add_tag_desc(int("tid", RadIntId::U32));
    s.add_tag_desc(TagDesc {
        name: "weight".to_string(),
        typeid: RadType::Float(RadFloatId::F64),
    });
    s.add_tag_desc(int("status", RadIntId::U8));
    s.add_tag_desc(int("is_fw", RadIntId::U8));
    s.add_tag_desc(int("fragment_len", RadIntId::U32));
    s.add_tag_desc(int("read_len", RadIntId::U32));
    s.add_tag_desc(int("ref_pos", RadIntId::U32));
    s.add_tag_desc(int("fw_pos", RadIntId::U32));
    s.add_tag_desc(int("rc_pos", RadIntId::U32));
    s.add_tag_desc(int("format", RadIntId::U8));
    s
}

fn status_to_u8(s: MateStatus) -> u8 {
    match s {
        MateStatus::PairedEndLeft => 0,
        MateStatus::PairedEndRight => 1,
        MateStatus::PairedEndPaired => 2,
        MateStatus::SingleEnd => 3,
    }
}

fn status_from_u8(v: u8) -> Result<MateStatus> {
    Ok(match v {
        0 => MateStatus::PairedEndLeft,
        1 => MateStatus::PairedEndRight,
        2 => MateStatus::PairedEndPaired,
        3 => MateStatus::SingleEnd,
        other => bail!("invalid MateStatus tag value {other}"),
    })
}

/// Append one fragment's record (`num_aln`, `frag_key`, per-alignment tags) to
/// `buf`. Hand-written little-endian bytes in the declared tag order: this is the
/// hot path, so it avoids per-record allocation (as piscem's own RAD writer does).
fn write_record(buf: &mut Vec<u8>, key: u128, maps: &[ScoredMapping]) {
    buf.extend_from_slice(&(maps.len() as u32).to_le_bytes());
    buf.extend_from_slice(&key.to_le_bytes());
    for m in maps {
        buf.extend_from_slice(&m.tid.to_le_bytes());
        buf.extend_from_slice(&m.weight.to_le_bytes());
        buf.push(status_to_u8(m.status));
        buf.push(m.is_fw as u8);
        buf.extend_from_slice(&(m.fragment_len as u32).to_le_bytes());
        buf.extend_from_slice(&(m.read_len as u32).to_le_bytes());
        buf.extend_from_slice(&(m.ref_pos as u32).to_le_bytes());
        buf.extend_from_slice(&(m.fw_pos as u32).to_le_bytes());
        buf.extend_from_slice(&(m.rc_pos as u32).to_le_bytes());
        buf.push(m.format.map_or(FORMAT_NONE, |f| f.format_id()));
    }
}

/// Write the RAD prelude: header + the three tag sections.
fn write_prelude<W: Write>(
    w: &mut W,
    salmon: &SalmonIndex,
    is_paired: bool,
    num_chunks: u64,
) -> Result<()> {
    w.write_all(&[is_paired as u8])?;
    let num_refs = salmon.num_refs();
    w.write_all(&(num_refs as u64).to_le_bytes())?;
    for t in 0..num_refs {
        let name = salmon.ref_name(t);
        let len: u16 = name
            .len()
            .try_into()
            .context("reference name longer than 65535 bytes")?;
        w.write_all(&len.to_le_bytes())?;
        w.write_all(name.as_bytes())?;
    }
    w.write_all(&num_chunks.to_le_bytes())?;
    // file tags: none (the deterministic store carries no file-level values).
    TagSection::new_with_label(TagSectionLabel::FileTags)
        .write(w)
        .context("writing RAD file tag section")?;
    read_tag_section()
        .write(w)
        .context("writing RAD read tag section")?;
    aln_tag_section()
        .write(w)
        .context("writing RAD alignment tag section")?;
    Ok(())
}

/// Serialize one chunk (`nbytes`, `nrec`, records) into a standalone byte
/// buffer. Doing this off-lock lets the parallel store hold its file mutex only
/// for the final `write_all`.
fn serialize_chunk(frags: &[(u128, Vec<ScoredMapping>)]) -> Vec<u8> {
    let mut body = Vec::new();
    for (key, maps) in frags {
        write_record(&mut body, *key, maps);
    }
    // A single chunk over 4 GiB would need splitting; per-batch chunks are far
    // below that, so saturate rather than fail (the reader streams to EOF).
    let nbytes: u32 = (body.len() + 8).try_into().unwrap_or(u32::MAX);
    let mut out = Vec::with_capacity(body.len() + 8);
    out.extend_from_slice(&nbytes.to_le_bytes());
    out.extend_from_slice(&(frags.len() as u32).to_le_bytes());
    out.extend_from_slice(&body);
    out
}

/// Write a complete single-chunk RAD file holding `frags`. Convenience for tests
/// and small in-memory stores; the parallel path uses [`RadChunkWriter`].
pub fn write_rad(
    path: &Path,
    salmon: &SalmonIndex,
    is_paired: bool,
    frags: &[(u128, Vec<ScoredMapping>)],
) -> Result<()> {
    let f = File::create(path).with_context(|| format!("creating RAD file {}", path.display()))?;
    let mut w = BufWriter::new(f);
    write_prelude(&mut w, salmon, is_paired, 1)?;
    w.write_all(&serialize_chunk(frags))?;
    w.flush()?;
    Ok(())
}

/// Concurrent, append-only RAD writer for the deterministic two-pass store. The
/// prelude (with `num_chunks = 0`, "stream to EOF") is written at construction;
/// each mapping worker then appends its batch as a chunk via [`Self::write_chunk`]
/// (serialized off-lock, so the file mutex is held only for the append). Records
/// land in thread-arrival order; pass 2 sorts them by key, so the order on disk
/// does not affect the result.
pub struct RadChunkWriter {
    inner: Mutex<BufWriter<File>>,
    /// first append error, surfaced by [`Self::flush`] (the per-chunk write path
    /// is called from `paraseq` callbacks that cannot return our error type).
    err: Mutex<Option<io::Error>>,
}

impl RadChunkWriter {
    pub fn create(path: &Path, salmon: &SalmonIndex, is_paired: bool) -> Result<Self> {
        let f =
            File::create(path).with_context(|| format!("creating RAD file {}", path.display()))?;
        let mut w = BufWriter::new(f);
        write_prelude(&mut w, salmon, is_paired, 0)?;
        Ok(Self {
            inner: Mutex::new(w),
            err: Mutex::new(None),
        })
    }

    /// Append `frags` as one chunk. No-op when empty. Errors are captured and
    /// reported later by [`Self::flush`].
    pub fn write_chunk(&self, frags: &[(u128, Vec<ScoredMapping>)]) {
        if frags.is_empty() {
            return;
        }
        let bytes = serialize_chunk(frags);
        let mut w = self.inner.lock().unwrap();
        if let Err(e) = w.write_all(&bytes) {
            let mut slot = self.err.lock().unwrap();
            if slot.is_none() {
                *slot = Some(e);
            }
        }
    }

    /// Flush the buffered writer and surface the first append error, if any.
    pub fn flush(&self) -> Result<()> {
        if let Some(e) = self.err.lock().unwrap().take() {
            return Err(anyhow::Error::from(e).context("appending to RAD mapping store"));
        }
        self.inner
            .lock()
            .unwrap()
            .flush()
            .context("flushing RAD mapping store")?;
        Ok(())
    }
}

fn read_u8<R: Read>(r: &mut R) -> Result<u8> {
    let mut b = [0u8; 1];
    r.read_exact(&mut b)?;
    Ok(b[0])
}
fn read_u16<R: Read>(r: &mut R) -> Result<u16> {
    let mut b = [0u8; 2];
    r.read_exact(&mut b)?;
    Ok(u16::from_le_bytes(b))
}
fn read_u32<R: Read>(r: &mut R) -> Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}
/// Read a `u32`, returning `None` on a clean end of file (used to detect the last
/// chunk when `num_chunks` is the "stream to EOF" sentinel `0`).
fn read_u32_opt<R: Read>(r: &mut R) -> Result<Option<u32>> {
    let mut b = [0u8; 4];
    match r.read_exact(&mut b) {
        Ok(()) => Ok(Some(u32::from_le_bytes(b))),
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e.into()),
    }
}
fn read_u64<R: Read>(r: &mut R) -> Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}
fn read_u128<R: Read>(r: &mut R) -> Result<u128> {
    let mut b = [0u8; 16];
    r.read_exact(&mut b)?;
    Ok(u128::from_le_bytes(b))
}
fn read_f64<R: Read>(r: &mut R) -> Result<f64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}

/// Advance the reader past one tag section (`num_tags u16`, then each
/// `name_len u16 + name + type u8`, plus two type ids for an `Array` tag). The
/// deterministic store reads back its own fixed schema, so the section contents
/// only need to be skipped to reach the chunks.
fn skip_tag_section<R: Read>(r: &mut R) -> Result<()> {
    let num_tags = read_u16(r)?;
    for _ in 0..num_tags {
        let name_len = read_u16(r)? as usize;
        let mut name = vec![0u8; name_len];
        r.read_exact(&mut name)?;
        let type_id = read_u8(r)?;
        if type_id == 7 {
            // Array: length-type id + value-type id
            let _ = read_u8(r)?;
            let _ = read_u8(r)?;
        }
    }
    Ok(())
}

/// Read the prelude and return `num_chunks`.
fn read_prelude<R: Read>(r: &mut R) -> Result<u64> {
    let _is_paired = read_u8(r)?;
    let ref_count = read_u64(r)?;
    for _ in 0..ref_count {
        let len = read_u16(r)? as usize;
        let mut name = vec![0u8; len];
        r.read_exact(&mut name)?;
    }
    let num_chunks = read_u64(r)?;
    skip_tag_section(r).context("reading RAD file tag section")?;
    skip_tag_section(r).context("reading RAD read tag section")?;
    skip_tag_section(r).context("reading RAD alignment tag section")?;
    Ok(num_chunks)
}

/// Read one mapped fragment record (the inverse of [`write_record`]).
fn read_record<R: Read>(r: &mut R) -> Result<(u128, Vec<ScoredMapping>)> {
    let na = read_u32(r)? as usize;
    let key = read_u128(r)?;
    let mut maps = Vec::with_capacity(na);
    for _ in 0..na {
        let tid = read_u32(r)?;
        let weight = read_f64(r)?;
        let status = status_from_u8(read_u8(r)?)?;
        let is_fw = read_u8(r)? != 0;
        let fragment_len = read_u32(r)? as i32;
        let read_len = read_u32(r)? as i32;
        let ref_pos = read_u32(r)? as i32;
        let fw_pos = read_u32(r)? as i32;
        let rc_pos = read_u32(r)? as i32;
        let fmt_id = read_u8(r)?;
        let format = (fmt_id != FORMAT_NONE).then(|| LibraryFormat::from_format_id(fmt_id));
        maps.push(ScoredMapping {
            tid,
            is_fw,
            status,
            score: 0,
            fragment_len,
            read_len,
            weight,
            ref_pos,
            fw_pos,
            rc_pos,
            format,
            // SAM-output-only fields; not used by pass-2 inference.
            r1_pos: -1,
            r2_pos: -1,
            r2_fw: false,
            r1_score: 0,
        });
    }
    Ok((key, maps))
}

/// Read a full RAD file written by [`write_rad`] back into `(key, mappings)`
/// fragments (inference fields only).
pub fn read_rad(path: &Path) -> Result<Vec<(u128, Vec<ScoredMapping>)>> {
    let f = File::open(path).with_context(|| format!("opening RAD file {}", path.display()))?;
    let mut r = BufReader::new(f);
    let num_chunks = read_prelude(&mut r)?;
    let mut out = Vec::new();
    let read_one_chunk = |r: &mut BufReader<File>, out: &mut Vec<_>| -> Result<()> {
        let nrec = read_u32(r)?;
        out.reserve(nrec as usize);
        for _ in 0..nrec {
            out.push(read_record(r)?);
        }
        Ok(())
    };
    if num_chunks == 0 {
        // "stream to EOF" sentinel: keep reading chunks until the file ends
        // (the parallel writer does not know the chunk count up front).
        while read_u32_opt(&mut r)?.is_some() {
            read_one_chunk(&mut r, &mut out)?;
        }
    } else {
        for _ in 0..num_chunks {
            let _nbytes = read_u32(&mut r)?;
            read_one_chunk(&mut r, &mut out)?;
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_mapping(tid: u32, status: MateStatus, fmt: Option<LibraryFormat>) -> ScoredMapping {
        ScoredMapping {
            tid,
            is_fw: tid % 2 == 0,
            status,
            score: 0,
            fragment_len: 250 + tid as i32,
            read_len: 100,
            weight: 0.5 + (tid as f64) * 0.125,
            ref_pos: 42 + tid as i32,
            fw_pos: tid as i32,
            rc_pos: -1,
            format: fmt,
            r1_pos: -1,
            r2_pos: -1,
            r2_fw: false,
            r1_score: 0,
        }
    }

    /// Assert that the fields pass-2 inference reads survive a RAD round trip.
    fn assert_inference_fields_eq(a: &ScoredMapping, b: &ScoredMapping) {
        assert_eq!(a.tid, b.tid);
        assert_eq!(a.is_fw, b.is_fw);
        assert_eq!(a.status, b.status);
        assert_eq!(a.fragment_len, b.fragment_len);
        assert_eq!(a.read_len, b.read_len);
        assert_eq!(a.weight.to_bits(), b.weight.to_bits());
        assert_eq!(a.ref_pos, b.ref_pos);
        assert_eq!(a.fw_pos, b.fw_pos);
        assert_eq!(a.rc_pos, b.rc_pos);
        assert_eq!(a.format, b.format);
    }

    #[test]
    fn record_round_trips_through_buffer() {
        let fmt = Some(LibraryFormat::paired_default());
        let maps = vec![
            sample_mapping(0, MateStatus::PairedEndPaired, fmt),
            sample_mapping(7, MateStatus::PairedEndLeft, None),
            sample_mapping(123_456, MateStatus::SingleEnd, fmt),
        ];
        let key: u128 = 0x0123_4567_89AB_CDEF_FEDC_BA98_7654_3210;

        let mut buf = Vec::new();
        write_record(&mut buf, key, &maps);
        let mut cur = std::io::Cursor::new(buf);
        let (rkey, rmaps) = read_record(&mut cur).expect("read record");

        assert_eq!(rkey, key);
        assert_eq!(rmaps.len(), maps.len());
        for (a, b) in maps.iter().zip(&rmaps) {
            assert_inference_fields_eq(a, b);
        }
    }

    #[test]
    fn format_none_sentinel_round_trips() {
        // FORMAT_NONE must not collide with any real format_id.
        let maps = vec![sample_mapping(1, MateStatus::PairedEndPaired, None)];
        let mut buf = Vec::new();
        write_record(&mut buf, 1, &maps);
        let mut cur = std::io::Cursor::new(buf);
        let (_k, rmaps) = read_record(&mut cur).unwrap();
        assert_eq!(rmaps[0].format, None);
    }

    #[test]
    fn tag_sections_skip_is_consistent_with_write() {
        // A prelude's tag sections must be skippable back to byte alignment:
        // write read+aln sections, then skip them, and confirm we land exactly at
        // end of buffer.
        let mut buf = Vec::new();
        read_tag_section().write(&mut buf).unwrap();
        aln_tag_section().write(&mut buf).unwrap();
        let total = buf.len() as u64;
        let mut cur = std::io::Cursor::new(buf);
        skip_tag_section(&mut cur).unwrap();
        skip_tag_section(&mut cur).unwrap();
        assert_eq!(cur.position(), total, "tag section skip misaligned");
    }
}
