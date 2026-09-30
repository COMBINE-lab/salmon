//! Projection of genome-decoy alignments onto unspliced targets.
//!
//! An index built with `salmon index --unspliced` can hold its unspliced
//! (intron / gene-body) targets in two ways. In the *sequence* layout each
//! target's bases are indexed next to the transcripts; every intronic k-mer
//! then occurs twice, once in its target and once in the genome decoy, and
//! candidate generation pays for both. In the *projection* layout the targets
//! carry no k-mers: a fragment is aligned to the genome decoy as usual, and an
//! alignment that lies wholly inside a target's genomic interval is copied onto
//! that target, with the same score, before decoy filtering. Since a target is
//! a verbatim slice of the genome, the copy is exactly the alignment the
//! sequence layout would have found on the target itself.
//!
//! Coordinates move with the target: shifted for a plus-strand gene, reflected
//! for a minus-strand one (positions, both mates, strands and the observed
//! library format), so bias models, the fragment-length model and library
//! compatibility see the fragment as it lies on the target.

use salmon_core::{LibraryFormat, MateStatus, ReadOrientation, ReadStrandedness, ReadType};

use crate::score::RawMapping;

/// One unspliced target's interval on its genome decoy, 0-based half-open.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProjectionTarget {
    /// the target's reference id in the index
    pub tid: u32,
    pub start: u32,
    pub end: u32,
    /// the gene (and so the target sequence) reads on the minus strand
    pub minus: bool,
}

/// Bin width (log2) of the per-decoy interval lookup: 64 kb.
const BIN_SHIFT: u32 = 16;

/// Targets of one decoy reference, bucketed by 64 kb bins.
#[derive(Debug, Default)]
struct DecoyTargets {
    targets: Vec<ProjectionTarget>,
    /// `bins[b]` lists (indices into `targets`) the targets overlapping bin `b`
    bins: Vec<Vec<u32>>,
}

/// Unspliced target intervals of every genome decoy, for projection.
#[derive(Debug, Default)]
pub struct UnsplicedProjection {
    first_decoy: u32,
    decoys: Vec<DecoyTargets>,
    num_targets: usize,
}

impl UnsplicedProjection {
    /// Build from `(decoy tid, target)` pairs. Decoys occupy the contiguous
    /// block `[first_decoy, first_decoy + num_decoys)`.
    pub fn new(
        first_decoy: u32,
        num_decoys: usize,
        targets: impl IntoIterator<Item = (u32, ProjectionTarget)>,
    ) -> Self {
        let mut decoys: Vec<DecoyTargets> =
            (0..num_decoys).map(|_| DecoyTargets::default()).collect();
        let mut num_targets = 0;
        for (decoy, t) in targets {
            let d = &mut decoys[(decoy - first_decoy) as usize];
            let i = d.targets.len() as u32;
            let (b0, b1) = (t.start >> BIN_SHIFT, (t.end.max(1) - 1) >> BIN_SHIFT);
            if d.bins.len() <= b1 as usize {
                d.bins.resize_with(b1 as usize + 1, Vec::new);
            }
            for b in b0..=b1 {
                d.bins[b as usize].push(i);
            }
            d.targets.push(t);
            num_targets += 1;
        }
        Self {
            first_decoy,
            decoys,
            num_targets,
        }
    }

    /// Number of targets.
    pub fn len(&self) -> usize {
        self.num_targets
    }

    /// Whether there are no targets at all.
    pub fn is_empty(&self) -> bool {
        self.num_targets == 0
    }

    /// Call `f` for every target whose interval contains `[lo, hi)` on decoy
    /// `decoy_tid`.
    fn containing(&self, decoy_tid: u32, lo: i32, hi: i32, mut f: impl FnMut(&ProjectionTarget)) {
        let Some(d) = decoy_tid
            .checked_sub(self.first_decoy)
            .and_then(|i| self.decoys.get(i as usize))
        else {
            return;
        };
        if lo < 0 || hi <= lo {
            return;
        }
        let Some(bin) = d.bins.get((lo as u32 >> BIN_SHIFT) as usize) else {
            return;
        };
        for &i in bin {
            let t = &d.targets[i as usize];
            if t.start as i32 <= lo && hi <= t.end as i32 {
                f(t);
            }
        }
    }
}

/// Genomic extent `[lo, hi)` a decoy mapping covers, from the mates it placed.
fn extent(m: &RawMapping) -> Option<(i32, i32)> {
    match m.status {
        MateStatus::PairedEndPaired => {
            // Both mates are placed (the partner of a recovered pair is
            // estimated from the fragment length); the fragment starts at the
            // leftmost of them.
            let lo = match (m.r1_pos >= 0, m.r2_pos >= 0) {
                (true, true) => m.r1_pos.min(m.r2_pos),
                (true, false) => m.r1_pos,
                (false, true) => m.r2_pos,
                (false, false) => m.ref_pos,
            };
            Some((lo, lo + m.fragment_len))
        }
        MateStatus::PairedEndLeft | MateStatus::PairedEndRight => {
            Some((m.ref_pos, m.ref_pos + m.read_len))
        }
        MateStatus::SingleEnd => Some((m.r1_pos, m.r1_pos + m.read_len)),
    }
}

/// `m` (a decoy mapping inside `t`) expressed on target `t`.
/// `len1` / `len2` are the lengths of read 1 and read 2 (`len2` is unused for
/// single-end reads); `hi` is the end of the mapping's genomic extent.
fn transform(m: &RawMapping, hi: i32, t: &ProjectionTarget, len1: i32, len2: i32) -> RawMapping {
    let mut p = *m;
    p.tid = t.tid;
    p.is_decoy = false;
    let s = t.start as i32;
    if !t.minus {
        let shift = |x: i32| if x >= 0 { x - s } else { x };
        p.ref_pos = shift(m.ref_pos);
        p.fw_pos = shift(m.fw_pos);
        p.rc_pos = shift(m.rc_pos);
        p.r1_pos = shift(m.r1_pos);
        p.r2_pos = shift(m.r2_pos);
        return p;
    }
    // Minus strand: the target is the reverse complement of `[start, end)`.
    // A base at genomic `x` sits at `end - 1 - x` on the target, so a half-open
    // genomic interval `[a, b)` becomes `[end - b, end - a)`.
    let e = t.end as i32;
    let new_lo = e - hi;
    let flip_start = |pos: i32, len: i32| if pos >= 0 { e - (pos + len) } else { pos };
    p.is_fw = !m.is_fw;
    match m.status {
        MateStatus::PairedEndPaired => {
            p.r2_fw = !m.r2_fw;
            p.ref_pos = new_lo;
            p.r1_pos = flip_start(m.r1_pos, len1);
            p.r2_pos = flip_start(m.r2_pos, len2);
            if m.fw_pos >= 0 && m.rc_pos >= 0 {
                // Concordant pair: 5' fragment start and 3' fragment end.
                p.fw_pos = new_lo;
                p.rc_pos = new_lo + m.fragment_len - 1;
            } else if m.fw_pos >= 0 {
                // Recovered pair: the anchor's own start, on its (now reverse)
                // strand.
                let alen = if m.r1_pos == m.fw_pos { len1 } else { len2 };
                p.rc_pos = e - (m.fw_pos + alen);
                p.fw_pos = -1;
            } else if m.rc_pos >= 0 {
                let alen = if m.r1_pos == m.rc_pos { len1 } else { len2 };
                p.fw_pos = e - (m.rc_pos + alen);
                p.rc_pos = -1;
            }
            p.format = m
                .format
                .map(|_| salmon_core::observed_paired_format(p.is_fw, p.r2_fw));
        }
        MateStatus::SingleEnd => {
            // `ref_pos` is the read's 5' boundary: its start forward, its end
            // reverse. Reflection maps one onto the other.
            p.ref_pos = e - m.ref_pos;
            p.r1_pos = new_lo;
            (p.fw_pos, p.rc_pos) = if m.fw_pos >= 0 {
                (-1, new_lo)
            } else {
                (new_lo, -1)
            };
            p.format = m.format.map(|_| {
                LibraryFormat::new(
                    ReadType::SingleEnd,
                    ReadOrientation::None,
                    if p.is_fw {
                        ReadStrandedness::S
                    } else {
                        ReadStrandedness::A
                    },
                )
            });
        }
        _ => {
            // Orphan: one mate placed, leftmost start in `ref_pos`.
            p.ref_pos = new_lo;
            (p.fw_pos, p.rc_pos) = if m.fw_pos >= 0 {
                (-1, new_lo)
            } else {
                (new_lo, -1)
            };
            if m.r1_pos >= 0 {
                p.r1_pos = new_lo;
            }
            if m.r2_pos >= 0 {
                p.r2_pos = new_lo;
            }
        }
    }
    p
}

/// Append to `raw` a copy of every decoy mapping onto each unspliced target
/// whose interval contains it. Run before decoy filtering and per-target
/// collapse, so the copies compete like any transcript mapping: they tie with
/// the decoy mapping they come from, and a tie keeps the target.
pub fn project_decoy_mappings(
    raw: &mut Vec<RawMapping>,
    proj: &UnsplicedProjection,
    len1: i32,
    len2: i32,
) {
    let n = raw.len();
    for i in 0..n {
        let m = raw[i];
        if !m.is_decoy {
            continue;
        }
        let Some((lo, hi)) = extent(&m) else {
            continue;
        };
        proj.containing(m.tid, lo, hi, |t| {
            raw.push(transform(&m, hi, t, len1, len2));
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DECOY: u32 = 10;

    fn proj() -> UnsplicedProjection {
        UnsplicedProjection::new(
            DECOY,
            2,
            [
                // plus-strand target on decoy 10
                (
                    DECOY,
                    ProjectionTarget {
                        tid: 3,
                        start: 1_000,
                        end: 5_000,
                        minus: false,
                    },
                ),
                // minus-strand target overlapping it, spanning two bins
                (
                    DECOY,
                    ProjectionTarget {
                        tid: 4,
                        start: 4_000,
                        end: 70_000,
                        minus: true,
                    },
                ),
                // target on the second decoy
                (
                    DECOY + 1,
                    ProjectionTarget {
                        tid: 5,
                        start: 0,
                        end: 300,
                        minus: false,
                    },
                ),
            ],
        )
    }

    fn pair(pos1: i32, pos2: i32, flen: i32) -> RawMapping {
        let lo = pos1.min(pos2);
        RawMapping {
            tid: DECOY,
            is_fw: true,
            status: MateStatus::PairedEndPaired,
            score: 150,
            fragment_len: flen,
            read_len: 0,
            is_decoy: true,
            ref_pos: lo,
            fw_pos: lo,
            rc_pos: lo + flen - 1,
            format: Some(salmon_core::observed_paired_format(true, false)),
            r1_pos: pos1,
            r2_pos: pos2,
            r2_fw: false,
            r1_score: 75,
        }
    }

    #[test]
    fn projects_only_fully_contained_fragments() {
        let p = proj();
        // [1500, 1750): inside target 3 only
        let mut raw = vec![pair(1_500, 1_675, 250)];
        project_decoy_mappings(&mut raw, &p, 75, 75);
        assert_eq!(raw.iter().map(|m| m.tid).collect::<Vec<_>>(), [DECOY, 3]);
        // [4100, 4350): inside both targets
        let mut raw = vec![pair(4_100, 4_275, 250)];
        project_decoy_mappings(&mut raw, &p, 75, 75);
        assert_eq!(raw.iter().map(|m| m.tid).collect::<Vec<_>>(), [DECOY, 3, 4]);
        // [4900, 5150): runs past target 3's end, inside target 4
        let mut raw = vec![pair(4_900, 5_075, 250)];
        project_decoy_mappings(&mut raw, &p, 75, 75);
        assert_eq!(raw.iter().map(|m| m.tid).collect::<Vec<_>>(), [DECOY, 4]);
        // second bin of target 4, and the second decoy
        let mut raw = vec![
            pair(66_000, 66_175, 250),
            RawMapping {
                tid: DECOY + 1,
                ..pair(10, 185, 250)
            },
        ];
        project_decoy_mappings(&mut raw, &p, 75, 75);
        assert_eq!(
            raw.iter().map(|m| m.tid).collect::<Vec<_>>(),
            [DECOY, DECOY + 1, 4, 5]
        );
        // transcript mappings are never projected
        let mut raw = vec![RawMapping {
            is_decoy: false,
            tid: 0,
            ..pair(1_500, 1_675, 250)
        }];
        project_decoy_mappings(&mut raw, &p, 75, 75);
        assert_eq!(raw.len(), 1);
    }

    #[test]
    fn plus_strand_projection_is_a_shift() {
        let mut raw = vec![pair(1_500, 1_675, 250)];
        project_decoy_mappings(&mut raw, &proj(), 75, 75);
        let m = raw[1];
        assert!(!m.is_decoy);
        assert_eq!((m.ref_pos, m.r1_pos, m.r2_pos), (500, 500, 675));
        assert_eq!((m.fw_pos, m.rc_pos), (500, 749));
        assert_eq!(
            (m.is_fw, m.r2_fw, m.score, m.fragment_len),
            (true, false, 150, 250)
        );
        assert_eq!(m.format, raw[0].format);
    }

    /// A fragment on the genome, reverse-complemented into a minus-strand
    /// target: mates swap sides and strands, and the observed format flips
    /// between sense-first and antisense-first.
    #[test]
    fn minus_strand_projection_reflects_the_fragment() {
        // target 4 is [4000, 70000) on the minus strand
        let mut raw = vec![pair(10_000, 10_175, 250)];
        project_decoy_mappings(&mut raw, &proj(), 75, 75);
        let m = raw[1];
        assert_eq!(m.tid, 4);
        let e = 70_000;
        // genomic fragment [10000, 10250) -> target [e - 10250, e - 10000)
        assert_eq!(m.ref_pos, e - 10_250);
        // read 1 [10000, 10075) -> [e - 10075, e - 10000); read 2 likewise
        assert_eq!(m.r1_pos, e - 10_075);
        assert_eq!(m.r2_pos, e - 10_250);
        assert_eq!((m.is_fw, m.r2_fw), (false, true));
        assert_eq!((m.fw_pos, m.rc_pos), (e - 10_250, e - 10_250 + 249));
        assert_eq!(
            m.format,
            Some(salmon_core::observed_paired_format(false, true))
        );
        assert_ne!(m.format, raw[0].format);
    }

    #[test]
    fn minus_strand_orphans_and_single_end_reads() {
        let e = 70_000;
        let orphan = RawMapping {
            status: MateStatus::PairedEndLeft,
            read_len: 75,
            fragment_len: 0,
            ref_pos: 20_000,
            fw_pos: 20_000,
            rc_pos: -1,
            r1_pos: 20_000,
            r2_pos: -1,
            format: None,
            ..pair(0, 0, 0)
        };
        let mut raw = vec![orphan];
        project_decoy_mappings(&mut raw, &proj(), 75, 75);
        let m = raw[1];
        assert_eq!(
            (m.ref_pos, m.r1_pos, m.r2_pos),
            (e - 20_075, e - 20_075, -1)
        );
        assert_eq!((m.fw_pos, m.rc_pos, m.is_fw), (-1, e - 20_075, false));
        assert_eq!(m.format, None);

        // single end, reverse read: 5' boundary is its end (20075)
        let se = RawMapping {
            status: MateStatus::SingleEnd,
            is_fw: false,
            read_len: 75,
            fragment_len: 0,
            ref_pos: 20_075,
            fw_pos: -1,
            rc_pos: 20_000,
            r1_pos: 20_000,
            r2_pos: -1,
            format: Some(LibraryFormat::new(
                ReadType::SingleEnd,
                ReadOrientation::None,
                ReadStrandedness::A,
            )),
            ..pair(0, 0, 0)
        };
        let mut raw = vec![se];
        project_decoy_mappings(&mut raw, &proj(), 75, 0);
        let m = raw[1];
        // now forward on the target, starting at its 5' end
        assert!(m.is_fw);
        assert_eq!(
            (m.ref_pos, m.r1_pos, m.fw_pos, m.rc_pos),
            (e - 20_075, e - 20_075, e - 20_075, -1)
        );
        assert_eq!(
            m.format,
            Some(LibraryFormat::new(
                ReadType::SingleEnd,
                ReadOrientation::None,
                ReadStrandedness::S
            ))
        );
    }
}
