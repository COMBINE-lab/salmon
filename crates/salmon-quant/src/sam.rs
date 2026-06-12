//! `--writeMappings` SAM output.
//!
//! Mirrors salmon/pufferfish's `writeAlignmentsToStream` (`SAMWriter.hpp`): under
//! standard parameters salmon does **not** compute a real CIGAR — it spoofs a
//! full-length `<readLen>M` (soft-clipping only where the read overhangs the
//! transcript end, via `adjustOverhang`), so no traceback is needed and this is
//! compatible with our score-only aligner. Each fragment emits one record per
//! retained mapping (the eq-class members), with `NH`/`HI`/`XT`/`AS` tags.
//!
//! Records are accumulated per worker thread into a string and flushed to the
//! shared output under a mutex at batch boundaries to bound lock contention.

use std::io::Write;
use std::sync::Mutex;

use salmon_core::MateStatus;
use salmon_core::RefProvider as _;
use salmon_index::SalmonIndex;
use salmon_map::ScoredMapping;

const SAM_VERSION: &str = "1.0";

/// Thread-safe SAM sink: a buffered file behind a mutex. Per-thread record
/// strings are appended in bulk via [`SamWriter::write_block`].
pub struct SamWriter {
    inner: Mutex<Box<dyn Write + Send>>,
}

impl SamWriter {
    /// Open `path` for SAM output and write the header (`@HD`, one `@SQ` per
    /// reference — including decoys, matching salmon's full ref table — `@PG`).
    pub fn create(path: &std::path::Path, salmon: &SalmonIndex, cmd: &str) -> anyhow::Result<Self> {
        use anyhow::Context as _;
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("creating SAM output dir {}", parent.display()))?;
            }
        }
        let f = std::fs::File::create(path)
            .with_context(|| format!("creating SAM output {}", path.display()))?;
        let mut w: Box<dyn Write + Send> = Box::new(std::io::BufWriter::with_capacity(1 << 20, f));
        writeln!(w, "@HD\tVN:{SAM_VERSION}\tSO:unknown")?;
        for tid in 0..salmon.num_refs() {
            writeln!(
                w,
                "@SQ\tSN:{}\tLN:{}",
                salmon.ref_name(tid),
                salmon.ref_len(tid)
            )?;
        }
        writeln!(
            w,
            "@PG\tID:salmon\tPN:salmon\tVN:{}\tCL:{}",
            crate::output::SALMON_VERSION,
            cmd
        )?;
        Ok(Self {
            inner: Mutex::new(w),
        })
    }

    /// Append a worker thread's accumulated SAM text in one locked write.
    pub fn write_block(&self, block: &str) -> std::io::Result<()> {
        if block.is_empty() {
            return Ok(());
        }
        let mut g = self.inner.lock().unwrap();
        g.write_all(block.as_bytes())
    }

    /// Flush the underlying writer (call once at end of run).
    pub fn flush(&self) -> std::io::Result<()> {
        self.inner.lock().unwrap().flush()
    }
}

// SAM flag bits.
const PAIRED: u16 = 0x1;
const PROPER_PAIR: u16 = 0x2;
#[allow(dead_code)] // kept for completeness of the SAM flag-bit set
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const IS_RC: u16 = 0x10;
const MATE_RC: u16 = 0x20;
const READ1: u16 = 0x40;
const READ2: u16 = 0x80;
const SECONDARY: u16 = 0x100;

/// Reverse-complement a read; reverse the qualities (so the SAM SEQ/QUAL are in
/// forward-reference orientation, as salmon writes them for reverse-strand mates).
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

/// Spoofed CIGAR with transcript-end overhang clipping (pufferfish `adjustOverhang`):
/// `<readLen>M`, becoming `<S>...<M>` where the read hangs off either end.
/// Returns `(cigar, adjusted_pos)`.
fn overhang_cigar(pos: i32, read_len: i32, txp_len: i32) -> (String, i32) {
    if pos + read_len < 0 {
        (format!("{read_len}S"), 0)
    } else if pos < 0 {
        let match_len = read_len + pos;
        let clip = read_len - match_len;
        (format!("{clip}S{match_len}M"), 0)
    } else if pos > txp_len {
        (format!("{read_len}S"), pos)
    } else if pos + read_len > txp_len {
        let match_len = txp_len - pos;
        let clip = read_len - match_len;
        (format!("{match_len}M{clip}S"), pos)
    } else {
        (format!("{read_len}M"), pos)
    }
}

/// Processed read name: up to the first space, with a trailing `/1` or `/2` trimmed.
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

/// Append SAM records for one fragment's mappings to `buf`.
///
/// `r2` is `None` for single-end. `maps` are the fragment's retained mappings
/// (eq-class members); the first is primary, the rest carry the secondary flag.
pub fn write_fragment(
    buf: &mut String,
    salmon: &SalmonIndex,
    r1_id: &[u8],
    r1_seq: &[u8],
    r2: Option<(&[u8], &[u8])>,
    maps: &[ScoredMapping],
) {
    use std::fmt::Write as _;
    let name1 = String::from_utf8_lossy(read_name(r1_id)).into_owned();
    let (name2, r2_seq) = match r2 {
        Some((id, seq)) => (String::from_utf8_lossy(read_name(id)).into_owned(), seq),
        None => (name1.clone(), &b""[..]),
    };
    let r1_len = r1_seq.len() as i32;
    let r2_len = r2_seq.len() as i32;
    let nh = maps.len();

    for (idx, m) in maps.iter().enumerate() {
        let secondary = if idx == 0 { 0 } else { SECONDARY };
        let hi = idx + 1;
        let txp_len = salmon.ref_len(m.tid as usize) as i32;
        let refname = salmon.ref_name(m.tid as usize).to_string();
        let xt = if salmon.is_decoy(m.tid) { 'D' } else { 'T' };

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
                let (c1, p1) = overhang_cigar(m.r1_pos, r1_len, txp_len);
                let (c2, p2) = overhang_cigar(m.r2_pos, r2_len, txp_len);
                let seq1 = if m.is_fw { r1_seq.to_vec() } else { rc(r1_seq) };
                let seq2 = if m.r2_fw { r2_seq.to_vec() } else { rc(r2_seq) };
                let tlen1 = if read1_first { frag_len } else { -frag_len };
                let as1 = m.r1_score;
                let as2 = m.score - m.r1_score;
                // QUAL is `*` (salmon's default writeMappings omits qualities).
                let _ = writeln!(
                    buf,
                    "{name1}\t{f1}\t{refname}\t{}\t1\t{c1}\t=\t{}\t{tlen1}\t{}\t*\tNH:i:{nh}\tHI:i:{hi}\tXT:A:{xt}\tAS:i:{as1}",
                    p1 + 1, p2 + 1, String::from_utf8_lossy(&seq1),
                );
                let _ = writeln!(
                    buf,
                    "{name2}\t{f2}\t{refname}\t{}\t1\t{c2}\t=\t{}\t{}\t{}\t*\tNH:i:{nh}\tHI:i:{hi}\tXT:A:{xt}\tAS:i:{as2}",
                    p2 + 1, p1 + 1, -tlen1, String::from_utf8_lossy(&seq2),
                );
            }
            MateStatus::SingleEnd => {
                let mut f1 = secondary;
                if !m.is_fw {
                    f1 |= IS_RC;
                }
                let (c1, p1) = overhang_cigar(m.r1_pos, r1_len, txp_len);
                let seq1 = if m.is_fw { r1_seq.to_vec() } else { rc(r1_seq) };
                let _ = writeln!(
                    buf,
                    "{name1}\t{f1}\t{refname}\t{}\t1\t{c1}\t*\t0\t0\t{}\t*\tNH:i:{nh}\tHI:i:{hi}\tXT:A:{xt}\tAS:i:{}",
                    p1 + 1, String::from_utf8_lossy(&seq1), m.r1_score,
                );
            }
            MateStatus::PairedEndLeft | MateStatus::PairedEndRight => {
                // Orphan: one mate mapped. Emit the mapped mate (the unmapped mate
                // is omitted; salmon writes it with the unmapped flag, but it carries
                // no placement and downstream tools tolerate its absence).
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
                let (c, p) = overhang_cigar(mapped_pos, mapped_len, txp_len);
                let seq = if m.is_fw {
                    mapped_seq.to_vec()
                } else {
                    rc(mapped_seq)
                };
                let _ = writeln!(
                    buf,
                    "{mapped_name}\t{f}\t{refname}\t{}\t1\t{c}\t*\t0\t0\t{}\t*\tNH:i:{nh}\tHI:i:{hi}\tXT:A:{xt}\tAS:i:{}",
                    p + 1, String::from_utf8_lossy(&seq), m.r1_score,
                );
            }
        }
    }
}
