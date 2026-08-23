//! Format-neutral borrowed alignment records shared by SAM and BAM output.
//!
//! # Why this module exists
//!
//! salmon can write its per-fragment mappings as text (SAM) or as the equivalent
//! compressed binary (BAM). The two formats encode exactly the same information,
//! so if each writer derived flags, CIGARs, positions and tags for itself the two
//! would inevitably drift apart. Instead this module derives every field once, in
//! [`emit_fragment_records`], and hands each format a borrowed record to
//! serialize. The header is shared the same way, via [`header_text`].
//!
//! "Borrowed" means the record points at the caller's read name and sequence
//! rather than copying them, so emitting a mapping allocates nothing.
//!
//! # Background: SAM records
//!
//! A SAM record describes one read's placement: which reference, at what
//! position, in which orientation, and how the read's bases line up with the
//! reference. The alignment shape is a CIGAR string — here only `M` (bases
//! consuming both read and reference) and `S` (soft-clipped: present in the read,
//! not aligned). salmon does not compute a base-level alignment, so it emits a
//! spoofed `<readLen>M` with soft clips where the read overhangs a transcript end.

use std::io;

use salmon_core::MateStatus;
use salmon_core::RefProvider as _;
use salmon_index::SalmonIndex;
use salmon_map::ScoredMapping;

// SAM FLAG bits (a bitfield in field 2 of every record). Each names one yes/no
// property of the alignment; a record's flags are these OR-ed together.
/// The read is part of a pair.
pub const PAIRED: u16 = 0x1;
/// Both mates placed consistently (right orientation, plausible distance).
pub const PROPER_PAIR: u16 = 0x2;
/// This read's mate did not map (an orphan).
pub const MATE_UNMAPPED: u16 = 0x8;
/// This read aligns to the reverse strand, so its stored sequence is the
/// reverse complement of what the sequencer produced.
pub const IS_RC: u16 = 0x10;
/// The mate aligns to the reverse strand.
pub const MATE_RC: u16 = 0x20;
/// This is read 1 of the pair.
pub const READ1: u16 = 0x40;
/// This is read 2 of the pair.
pub const READ2: u16 = 0x80;
/// An additional, non-primary placement of the same fragment. Exactly one
/// placement per fragment lacks this bit; downstream tools count primaries to
/// avoid counting a multi-mapping fragment several times.
pub const SECONDARY: u16 = 0x100;

#[derive(Clone, Copy)]
pub enum CigarKind {
    Match,
    SoftClip,
}

#[derive(Clone, Copy)]
pub struct CigarOp {
    pub kind: CigarKind,
    pub len: usize,
}

/// A CIGAR of at most two operations, stored inline.
///
/// salmon's spoofed alignments never need more than two (`M`, or a clip plus a
/// match at an overhanging end), so a fixed array avoids allocating one `Vec` per
/// emitted record.
#[derive(Clone, Copy)]
pub struct Cigar {
    ops: [CigarOp; 2],
    len: usize,
}

impl Cigar {
    pub fn as_slice(&self) -> &[CigarOp] {
        &self.ops[..self.len]
    }
}

/// One read's placement, with every field already in the form both writers need.
/// The `'a` lifetime ties the borrowed name and sequence to the caller's buffers.
pub struct AlignmentRecord<'a> {
    pub name: &'a [u8],
    pub flags: u16,
    pub reference_id: usize,
    pub position: i32,
    pub mapping_quality: u8,
    pub cigar: Cigar,
    pub mate_reference_id: Option<usize>,
    pub mate_position: Option<i32>,
    pub template_length: i32,
    pub sequence: &'a [u8],
    pub reverse_complement: bool,
    pub nh: usize,
    pub hi: usize,
    pub xt: u8,
    pub alignment_score: i32,
}

/// SAM version reported in `@HD VN:`. Shared by SAM and BAM so the two formats
/// cannot drift; kept at `1.0` to match C++ salmon's `--writeMappings` header
/// byte-for-byte.
const SAM_VERSION: &str = "1.0";

/// The header text shared by both mapping-output formats: `@HD`, one `@SQ` per
/// reference (including decoys, matching salmon's full ref table), and `@PG`.
///
/// This is the single source of truth for the header, mirroring
/// [`emit_fragment_records`] for record bodies. SAM writes the text verbatim;
/// BAM stores it as the `text` block and then repeats the reference table in
/// binary form (see `bam::encode_header`).
pub fn header_text(salmon: &SalmonIndex, command: &str) -> String {
    use std::fmt::Write as _;
    // ~48 B/@SQ line for typical transcript names, plus the @HD and @PG lines.
    let mut text = String::with_capacity(salmon.num_refs() * 48 + command.len() + 64);
    let _ = writeln!(text, "@HD\tVN:{SAM_VERSION}\tSO:unknown");
    for tid in 0..salmon.num_refs() {
        let _ = writeln!(
            text,
            "@SQ\tSN:{}\tLN:{}",
            salmon.ref_name(tid),
            salmon.ref_len(tid)
        );
    }
    let _ = writeln!(
        text,
        "@PG\tID:salmon\tPN:salmon\tVN:{}\tCL:{}",
        crate::output::SALMON_VERSION,
        command
    );
    text
}

/// The read name to report: the FASTQ header up to the first space or tab, with
/// a trailing `/1` or `/2` mate suffix removed.
///
/// Both mates of a pair must carry the *same* name for a SAM/BAM reader to pair
/// them up, but many FASTQ files distinguish them with that suffix — so it has to
/// go.
pub fn read_name(id: &[u8]) -> &[u8] {
    let end = id
        .iter()
        .position(|&b| b == b' ' || b == b'\t')
        .unwrap_or(id.len());
    let name = &id[..end];
    if name.len() > 2 && name[name.len() - 2] == b'/' && matches!(name[name.len() - 1], b'1' | b'2')
    {
        &name[..name.len() - 2]
    } else {
        name
    }
}

/// Complement one base (unknown bases become `N`), used when emitting a
/// reverse-strand read's sequence.
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

/// The `MAPQ` to report for a fragment with `nh` placements.
///
/// SAM's `MAPQ` is a phred-scaled probability that the placement is wrong, and
/// salmon does not compute one: its whole design is to keep every plausible
/// placement and let the EM apportion them. What it does know is *how many*
/// placements a fragment has, which is the same thing STAR reports, so this uses
/// STAR's mapping: a uniquely placed fragment gets 255 and the value falls off
/// as the fragment becomes more ambiguous. Matching STAR matters because the
/// RNA-seq tools downstream of a salmon BAM are the ones already tuned to it: to
/// all of them `-q 255` means "uniquely placed".
///
/// The previous hardcoded `1` made every record look ambiguous, so any such
/// filter discarded the whole file. It is the one placeholder in this output
/// that costs nothing to replace: the count is already in hand when the records
/// are written, so no alignment work is implied (COMBINE-lab/salmon#1140,
/// #1141).
pub fn mapping_quality(nh: usize) -> u8 {
    match nh {
        // No placements is not a case a mapped record can reach (the writer
        // skips a fragment with none), so this arm exists for the record kinds
        // that have no placement by construction: an unmapped read, once such
        // records are written. 255 there would claim certainty about a
        // placement that does not exist; 0 is what samtools convention expects
        // on an unmapped record.
        0 => 0,
        1 => 255,
        2 => 3,
        3..=4 => 1,
        _ => 0,
    }
}

/// The CIGAR and clamped position for a read placed at `pos` on a transcript of
/// length `txp_len`.
///
/// A read can hang off either end of a transcript — the mapper places it by seed
/// position, and the fragment may genuinely extend past the annotated boundary.
/// SAM cannot express a negative position or an alignment past the reference end,
/// so the overhanging bases are soft-clipped and the position clamped into range.
/// Returns the CIGAR and the position to report.
fn overhang_cigar(pos: i32, read_len: i32, txp_len: i32) -> (Cigar, i32) {
    let read_len = read_len.max(0) as usize;
    let empty = CigarOp {
        kind: CigarKind::Match,
        len: 0,
    };
    let one = |kind, len| Cigar {
        ops: [CigarOp { kind, len }, empty],
        len: 1,
    };
    let two = |a, b| Cigar {
        ops: [a, b],
        len: 2,
    };
    if pos + (read_len as i32) < 0 {
        (one(CigarKind::SoftClip, read_len), 0)
    } else if pos < 0 {
        let matched = (read_len as i32 + pos).max(0) as usize;
        (
            two(
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: read_len - matched,
                },
                CigarOp {
                    kind: CigarKind::Match,
                    len: matched,
                },
            ),
            0,
        )
    } else if pos > txp_len {
        (one(CigarKind::SoftClip, read_len), pos)
    } else if pos + read_len as i32 > txp_len {
        let matched = (txp_len - pos).max(0) as usize;
        (
            two(
                CigarOp {
                    kind: CigarKind::Match,
                    len: matched,
                },
                CigarOp {
                    kind: CigarKind::SoftClip,
                    len: read_len - matched,
                },
            ),
            pos,
        )
    } else {
        (one(CigarKind::Match, read_len), pos)
    }
}

/// Emit borrowed records immediately; no record, name, sequence, CIGAR, or tag
/// allocation is retained between callback invocations.
///
/// This is the single source of truth for record bodies: one call per fragment
/// produces every record that fragment implies (two per placement for a proper
/// pair, one otherwise), and the caller's closure serializes them as SAM text or
/// BAM binary. The callback style is what keeps it allocation-free — nothing is
/// collected into a vector to be handed back.
pub fn emit_fragment_records(
    salmon: &SalmonIndex,
    r1_id: &[u8],
    r1_seq: &[u8],
    r2: Option<(&[u8], &[u8])>,
    maps: &[ScoredMapping],
    mut emit: impl FnMut(&AlignmentRecord<'_>) -> io::Result<()>,
) -> io::Result<()> {
    let name1 = read_name(r1_id);
    let (name2, r2_seq) = r2.map_or((name1, &[][..]), |(id, seq)| (read_name(id), seq));
    let nh = maps.len();
    let mapq = mapping_quality(nh);

    for (index, mapping) in maps.iter().enumerate() {
        // The first placement is primary; the rest are marked secondary so a
        // downstream counter does not count this fragment once per placement.
        let secondary = if index == 0 { 0 } else { SECONDARY };
        // HI is 1-based, unlike the loop index.
        let hi = index + 1;
        let tid = mapping.tid as usize;
        let txp_len = salmon.ref_len(tid) as i32;
        // salmon's XT tag: `T` for a transcript placement, `D` for a decoy, so a
        // reader can tell which placements were sinks rather than real targets.
        let xt = if salmon.is_decoy(mapping.tid) {
            b'D'
        } else {
            b'T'
        };

        match mapping.status {
            MateStatus::PairedEndPaired => {
                let (c1, p1) = overhang_cigar(mapping.r1_pos, r1_seq.len() as i32, txp_len);
                let (c2, p2) = overhang_cigar(mapping.r2_pos, r2_seq.len() as i32, txp_len);
                // TLEN is the observed template (fragment) length, truncated at
                // the transcript end for the same reason positions are clamped.
                let mut fragment_length = mapping.fragment_len;
                let min_pos = mapping.r1_pos.min(mapping.r2_pos);
                if min_pos + fragment_length > txp_len {
                    fragment_length = txp_len - min_pos;
                }
                // SAM signs TLEN: positive on the leftmost mate, negative on the
                // rightmost, so the two records sum to zero.
                let tlen1 = if mapping.r1_pos <= mapping.r2_pos {
                    fragment_length
                } else {
                    -fragment_length
                };
                let mut f1 = PAIRED | PROPER_PAIR | READ1 | secondary;
                let mut f2 = PAIRED | PROPER_PAIR | READ2 | secondary;
                // Each mate's own strand sets its IS_RC bit and the *other*
                // record's MATE_RC bit, so the pair stays mutually consistent.
                if !mapping.is_fw {
                    f1 |= IS_RC;
                    f2 |= MATE_RC;
                }
                if !mapping.r2_fw {
                    f2 |= IS_RC;
                    f1 |= MATE_RC;
                }
                emit(&AlignmentRecord {
                    name: name1,
                    flags: f1,
                    reference_id: tid,
                    position: p1,
                    mapping_quality: mapq,
                    cigar: c1,
                    mate_reference_id: Some(tid),
                    mate_position: Some(p2),
                    template_length: tlen1,
                    sequence: r1_seq,
                    reverse_complement: !mapping.is_fw,
                    nh,
                    hi,
                    xt,
                    alignment_score: mapping.r1_score,
                })?;
                emit(&AlignmentRecord {
                    name: name2,
                    flags: f2,
                    reference_id: tid,
                    position: p2,
                    mapping_quality: mapq,
                    cigar: c2,
                    mate_reference_id: Some(tid),
                    mate_position: Some(p1),
                    template_length: -tlen1,
                    sequence: r2_seq,
                    reverse_complement: !mapping.r2_fw,
                    nh,
                    hi,
                    xt,
                    alignment_score: mapping.score - mapping.r1_score,
                })?;
            }
            MateStatus::SingleEnd => {
                let (cigar, position) =
                    overhang_cigar(mapping.r1_pos, r1_seq.len() as i32, txp_len);
                emit(&AlignmentRecord {
                    name: name1,
                    flags: secondary | if mapping.is_fw { 0 } else { IS_RC },
                    reference_id: tid,
                    position,
                    mapping_quality: mapq,
                    cigar,
                    mate_reference_id: None,
                    mate_position: None,
                    template_length: 0,
                    sequence: r1_seq,
                    reverse_complement: !mapping.is_fw,
                    nh,
                    hi,
                    xt,
                    alignment_score: mapping.r1_score,
                })?;
            }
            // Orphan: one mate placed, the other not. One record is emitted, for
            // whichever mate mapped, flagged so a reader knows the pair is
            // incomplete rather than that this read is single-end.
            MateStatus::PairedEndLeft | MateStatus::PairedEndRight => {
                let left = mapping.status == MateStatus::PairedEndLeft;
                let (name, sequence, position, forward, score, read_flag) = if left {
                    (
                        name1,
                        r1_seq,
                        mapping.r1_pos,
                        mapping.is_fw,
                        mapping.r1_score,
                        READ1,
                    )
                } else {
                    (
                        name2,
                        r2_seq,
                        mapping.r2_pos,
                        mapping.r2_fw,
                        mapping.score - mapping.r1_score,
                        READ2,
                    )
                };
                let (cigar, position) = overhang_cigar(position, sequence.len() as i32, txp_len);
                emit(&AlignmentRecord {
                    name,
                    flags: PAIRED
                        | MATE_UNMAPPED
                        | read_flag
                        | secondary
                        | if forward { 0 } else { IS_RC },
                    reference_id: tid,
                    position,
                    mapping_quality: mapq,
                    cigar,
                    mate_reference_id: None,
                    mate_position: None,
                    template_length: 0,
                    sequence,
                    reverse_complement: !forward,
                    nh,
                    hi,
                    xt,
                    alignment_score: score,
                })?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// MAPQ has to distinguish a uniquely placed fragment from an ambiguous one.
    #[test]
    fn mapping_quality_follows_the_star_convention() {
        assert_eq!(mapping_quality(1), 255);
        assert_eq!(mapping_quality(2), 3);
        assert_eq!(mapping_quality(3), 1);
        assert_eq!(mapping_quality(4), 1);
        assert_eq!(mapping_quality(5), 0);
        // An unmapped record has no placement to be confident about.
        assert_eq!(mapping_quality(0), 0);
    }
}
