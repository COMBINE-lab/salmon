//! `--writeBam` BAM output.
//!
//! The BAM analogue of [`crate::sam`]: the same spoofed-CIGAR alignment records
//! (a full-length `<readLen>M` with transcript-end soft-clips, no real
//! traceback), written as BGZF-compressed BAM via `noodles-bam`. BAM stores read
//! sequences verbatim, so unlike CRAM it needs no reference repository.
//!
//! Records are built per worker thread into a `Vec<RecordBuf>` and flushed to the
//! shared writer under a mutex at batch boundaries, mirroring the SAM path's
//! contention profile. `finish` writes the BGZF EOF block.

use std::io;
use std::num::NonZero;
use std::sync::Mutex;

use bstr::BString;
use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_sam::alignment::io::Write as _;
use noodles_sam::alignment::record::cigar::op::{Kind, Op};
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::record::Flags;
use noodles_sam::alignment::record_buf::data::field::Value;
use noodles_sam::alignment::record_buf::{Cigar, Data, RecordBuf, Sequence};
use noodles_sam::header::record::value::map::ReferenceSequence;
use noodles_sam::header::record::value::Map;
use noodles_sam::{self as sam};

use salmon_core::MateStatus;
use salmon_core::RefProvider as _;
use salmon_index::SalmonIndex;
use salmon_map::ScoredMapping;

// SAM flag bits (same meanings as in `crate::sam`).
const PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
const MATE_UNMAPPED: u16 = 0x8;
const IS_RC: u16 = 0x10;
const MATE_RC: u16 = 0x20;
const READ1: u16 = 0x40;
const READ2: u16 = 0x80;
const SECONDARY: u16 = 0x100;

type Inner = bam::io::Writer<bgzf::io::Writer<std::io::BufWriter<std::fs::File>>>;

/// Thread-safe BAM sink: a `noodles-bam` writer behind a mutex, plus the SAM
/// header it needs for every record.
pub struct BamWriter {
    inner: Mutex<Inner>,
    header: sam::Header,
}

impl BamWriter {
    /// Open `path` for BAM output, build the header (one reference sequence per
    /// index reference, including decoys, plus a `@PG`), and write the BAM
    /// header.
    pub fn create(
        path: &std::path::Path,
        salmon: &SalmonIndex,
        _cmd: &str,
    ) -> anyhow::Result<Self> {
        use anyhow::Context as _;
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("creating BAM output dir {}", parent.display()))?;
            }
        }

        let n = salmon.num_refs();
        let mut ref_seqs = Vec::with_capacity(n);
        for tid in 0..n {
            let name = salmon.ref_name(tid);
            let len = (salmon.ref_len(tid) as usize).max(1);
            ref_seqs.push((
                BString::from(name.as_bytes()),
                Map::<ReferenceSequence>::new(NonZero::<usize>::new(len).unwrap()),
            ));
        }

        let header = sam::Header::builder()
            .set_reference_sequences(ref_seqs.into_iter().collect())
            .build();

        let file = std::fs::File::create(path)
            .with_context(|| format!("creating BAM output {}", path.display()))?;
        let mut writer = bam::io::Writer::new(std::io::BufWriter::new(file));
        writer.write_header(&header).context("writing BAM header")?;

        Ok(Self {
            inner: Mutex::new(writer),
            header,
        })
    }

    /// Write a worker thread's accumulated records under one lock.
    pub fn write_records(&self, records: &[RecordBuf]) -> io::Result<()> {
        if records.is_empty() {
            return Ok(());
        }
        let mut w = self.inner.lock().unwrap();
        for record in records {
            w.write_alignment_record(&self.header, record)?;
        }
        Ok(())
    }

    /// Write the BGZF EOF block (call once at end of run).
    pub fn finish(&self) -> io::Result<()> {
        self.inner.lock().unwrap().try_finish()
    }
}

/// Reverse-complement a read (SEQ written in forward-reference orientation for
/// reverse-strand mappings, matching the SAM path).
fn rc(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Processed read name: up to the first space, with a trailing `/1` or `/2`
/// trimmed (matches `crate::sam`'s record naming).
fn read_name(id: &[u8]) -> &[u8] {
    let end = id
        .iter()
        .position(|&b| b == b' ' || b == b'\t')
        .unwrap_or(id.len());
    let n = &id[..end];
    if n.len() > 2 && n[n.len() - 2] == b'/' && matches!(n[n.len() - 1], b'1' | b'2') {
        &n[..n.len() - 2]
    } else {
        n
    }
}

/// Spoofed CIGAR ops with transcript-end overhang clipping (`adjustOverhang`):
/// `<readLen>M`, becoming `<S>...<M>` where the read hangs off either end.
/// Returns `(ops, adjusted_pos)` (the analogue of `sam::overhang_cigar`).
fn overhang_ops(pos: i32, read_len: i32, txp_len: i32) -> (Vec<Op>, i32) {
    let rl = read_len.max(0) as usize;
    if pos + read_len < 0 {
        (vec![Op::new(Kind::SoftClip, rl)], 0)
    } else if pos < 0 {
        let m = (read_len + pos).max(0) as usize;
        (
            vec![Op::new(Kind::SoftClip, rl - m), Op::new(Kind::Match, m)],
            0,
        )
    } else if pos > txp_len {
        (vec![Op::new(Kind::SoftClip, rl)], pos)
    } else if pos + read_len > txp_len {
        let m = (txp_len - pos).max(0) as usize;
        (
            vec![Op::new(Kind::Match, m), Op::new(Kind::SoftClip, rl - m)],
            pos,
        )
    } else {
        (vec![Op::new(Kind::Match, rl)], pos)
    }
}

/// `NH`/`HI`/`XT`/`AS` tags, matching the SAM record's data.
fn tags(nh: usize, hi: usize, xt: u8, score: i32) -> Data {
    [
        (Tag::ALIGNMENT_HIT_COUNT, Value::Int32(nh as i32)),
        (Tag::HIT_INDEX, Value::Int32(hi as i32)),
        (Tag::new(b'X', b'T'), Value::Character(xt)),
        (Tag::ALIGNMENT_SCORE, Value::Int32(score)),
    ]
    .into_iter()
    .collect()
}

#[inline]
fn position(pos1: i32) -> Option<Position> {
    Position::new(pos1.max(1) as usize)
}

/// Build the BAM records for one fragment's mappings (the eq-class members; the
/// first is primary, the rest secondary). The analogue of `sam::write_fragment`,
/// returning owned [`RecordBuf`]s for the per-thread buffer.
pub fn build_records(
    salmon: &SalmonIndex,
    r1_id: &[u8],
    r1_seq: &[u8],
    r2: Option<(&[u8], &[u8])>,
    maps: &[ScoredMapping],
) -> Vec<RecordBuf> {
    let name1 = read_name(r1_id).to_vec();
    let (name2, r2_seq) = match r2 {
        Some((id, seq)) => (read_name(id).to_vec(), seq),
        None => (name1.clone(), &b""[..]),
    };
    let r1_len = r1_seq.len() as i32;
    let r2_len = r2_seq.len() as i32;
    let nh = maps.len();

    let mut out = Vec::with_capacity(maps.len());
    for (idx, m) in maps.iter().enumerate() {
        let secondary = if idx == 0 { 0 } else { SECONDARY };
        let hi = idx + 1;
        let txp_len = salmon.ref_len(m.tid as usize) as i32;
        let tid = m.tid as usize;
        let xt = if salmon.is_decoy(m.tid) { b'D' } else { b'T' };

        match m.status {
            MateStatus::PairedEndPaired => {
                let (read1_first, min_pos) = (m.r1_pos <= m.r2_pos, m.r1_pos.min(m.r2_pos));
                let mut frag_len = m.fragment_len;
                if min_pos + frag_len > txp_len {
                    frag_len = txp_len - min_pos;
                }
                let mut f1 = PAIRED | PROPER_PAIR | READ1 | secondary;
                let mut f2 = PAIRED | PROPER_PAIR | READ2 | secondary;
                if !m.is_fw {
                    f1 |= IS_RC;
                    f2 |= MATE_RC;
                }
                if !m.r2_fw {
                    f2 |= IS_RC;
                    f1 |= MATE_RC;
                }
                let (c1, p1) = overhang_ops(m.r1_pos, r1_len, txp_len);
                let (c2, p2) = overhang_ops(m.r2_pos, r2_len, txp_len);
                let seq1 = if m.is_fw { r1_seq.to_vec() } else { rc(r1_seq) };
                let seq2 = if m.r2_fw { r2_seq.to_vec() } else { rc(r2_seq) };
                let tlen1 = if read1_first { frag_len } else { -frag_len };

                let mut b1 = RecordBuf::builder()
                    .set_name(name1.clone())
                    .set_flags(Flags::from_bits_retain(f1))
                    .set_reference_sequence_id(tid)
                    .set_cigar(Cigar::from(c1))
                    .set_mate_reference_sequence_id(tid)
                    .set_template_length(tlen1)
                    .set_sequence(Sequence::from(seq1))
                    .set_data(tags(nh, hi, xt, m.r1_score));
                if let Some(p) = position(p1 + 1) {
                    b1 = b1.set_alignment_start(p);
                }
                if let Some(p) = position(p2 + 1) {
                    b1 = b1.set_mate_alignment_start(p);
                }
                out.push(b1.build());

                let mut b2 = RecordBuf::builder()
                    .set_name(name2.clone())
                    .set_flags(Flags::from_bits_retain(f2))
                    .set_reference_sequence_id(tid)
                    .set_cigar(Cigar::from(c2))
                    .set_mate_reference_sequence_id(tid)
                    .set_template_length(-tlen1)
                    .set_sequence(Sequence::from(seq2))
                    .set_data(tags(nh, hi, xt, m.score - m.r1_score));
                if let Some(p) = position(p2 + 1) {
                    b2 = b2.set_alignment_start(p);
                }
                if let Some(p) = position(p1 + 1) {
                    b2 = b2.set_mate_alignment_start(p);
                }
                out.push(b2.build());
            }
            MateStatus::SingleEnd => {
                let mut f = secondary;
                if !m.is_fw {
                    f |= IS_RC;
                }
                let (c, p) = overhang_ops(m.r1_pos, r1_len, txp_len);
                let seq = if m.is_fw { r1_seq.to_vec() } else { rc(r1_seq) };
                let mut b = RecordBuf::builder()
                    .set_name(name1.clone())
                    .set_flags(Flags::from_bits_retain(f))
                    .set_reference_sequence_id(tid)
                    .set_cigar(Cigar::from(c))
                    .set_sequence(Sequence::from(seq))
                    .set_data(tags(nh, hi, xt, m.r1_score));
                if let Some(pos) = position(p + 1) {
                    b = b.set_alignment_start(pos);
                }
                out.push(b.build());
            }
            MateStatus::PairedEndLeft | MateStatus::PairedEndRight => {
                let left_mapped = matches!(m.status, MateStatus::PairedEndLeft);
                let (mapped_pos, mapped_len, mapped_seq, mapped_name, is_r1) = if left_mapped {
                    (m.r1_pos, r1_len, r1_seq, &name1, true)
                } else {
                    (m.r2_pos, r2_len, r2_seq, &name2, false)
                };
                let mut f = PAIRED | secondary | MATE_UNMAPPED | if is_r1 { READ1 } else { READ2 };
                if !m.is_fw {
                    f |= IS_RC;
                }
                let (c, p) = overhang_ops(mapped_pos, mapped_len, txp_len);
                let seq = if m.is_fw {
                    mapped_seq.to_vec()
                } else {
                    rc(mapped_seq)
                };
                let mut b = RecordBuf::builder()
                    .set_name(mapped_name.clone())
                    .set_flags(Flags::from_bits_retain(f))
                    .set_reference_sequence_id(tid)
                    .set_cigar(Cigar::from(c))
                    .set_sequence(Sequence::from(seq))
                    .set_data(tags(nh, hi, xt, m.r1_score));
                if let Some(pos) = position(p + 1) {
                    b = b.set_alignment_start(pos);
                }
                out.push(b.build());
            }
        }
    }
    out
}
