//! Format-neutral borrowed alignment records shared by SAM and BAM output.

use std::io;

use salmon_core::MateStatus;
use salmon_core::RefProvider as _;
use salmon_index::SalmonIndex;
use salmon_map::ScoredMapping;

pub const PAIRED: u16 = 0x1;
pub const PROPER_PAIR: u16 = 0x2;
pub const MATE_UNMAPPED: u16 = 0x8;
pub const IS_RC: u16 = 0x10;
pub const MATE_RC: u16 = 0x20;
pub const READ1: u16 = 0x40;
pub const READ2: u16 = 0x80;
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

pub fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

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

    for (index, mapping) in maps.iter().enumerate() {
        let secondary = if index == 0 { 0 } else { SECONDARY };
        let hi = index + 1;
        let tid = mapping.tid as usize;
        let txp_len = salmon.ref_len(tid) as i32;
        let xt = if salmon.is_decoy(mapping.tid) {
            b'D'
        } else {
            b'T'
        };

        match mapping.status {
            MateStatus::PairedEndPaired => {
                let (c1, p1) = overhang_cigar(mapping.r1_pos, r1_seq.len() as i32, txp_len);
                let (c2, p2) = overhang_cigar(mapping.r2_pos, r2_seq.len() as i32, txp_len);
                let mut fragment_length = mapping.fragment_len;
                let min_pos = mapping.r1_pos.min(mapping.r2_pos);
                if min_pos + fragment_length > txp_len {
                    fragment_length = txp_len - min_pos;
                }
                let tlen1 = if mapping.r1_pos <= mapping.r2_pos {
                    fragment_length
                } else {
                    -fragment_length
                };
                let mut f1 = PAIRED | PROPER_PAIR | READ1 | secondary;
                let mut f2 = PAIRED | PROPER_PAIR | READ2 | secondary;
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
                    mapping_quality: 1,
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
                    mapping_quality: 1,
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
                    mapping_quality: 1,
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
                    mapping_quality: 1,
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
