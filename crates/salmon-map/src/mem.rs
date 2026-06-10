//! Maximal exact match (MEM) anchors.
//!
//! A [`Mem`] is an exact match between a stretch of the read and a stretch of a
//! reference, in a single orientation. In salmon/pufferfish these "uni-MEMs"
//! are produced by extending shared k-mers along a unitig and then projecting
//! the unitig occurrence onto each reference that contains it. Here we keep the
//! anchor purely geometric — read offset, reference offset, and length — so the
//! [chaining](crate::chain) DP can operate on anchors regardless of how they
//! were collected.
//!
//! All coordinates are 0-based and expressed in the *reference's forward*
//! frame: for a reverse-complement mapping the collection step is responsible
//! for transforming read coordinates so that a chain is colinear in increasing
//! `(read_start, ref_start)` order. Positions are `i32` to match piscem's hit
//! coordinates and to allow signed gap arithmetic.

/// An exact-match anchor between the read and a reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mem {
    /// 0-based start offset on the read
    pub read_start: i32,
    /// 0-based start offset on the reference
    pub ref_start: i32,
    /// length of the exact match (bases)
    pub len: i32,
}

impl Mem {
    pub fn new(read_start: i32, ref_start: i32, len: i32) -> Self {
        debug_assert!(len > 0, "MEM length must be positive");
        Self {
            read_start,
            ref_start,
            len,
        }
    }

    /// One past the last matched read base.
    #[inline]
    pub fn read_end(&self) -> i32 {
        self.read_start + self.len
    }

    /// One past the last matched reference base.
    #[inline]
    pub fn ref_end(&self) -> i32 {
        self.ref_start + self.len
    }

    /// The diagonal `ref_start - read_start`. Anchors on the same diagonal are
    /// gap-free continuations of one another.
    #[inline]
    pub fn diagonal(&self) -> i32 {
        self.ref_start - self.read_start
    }
}

/// Union length of read intervals covered by a set of MEMs (given by index into
/// `mems`). Used to compute a chain's read coverage for the consensus-fraction
/// filter. Overlapping anchors are counted once.
pub(crate) fn covered_read_bases(mems: &[Mem], indices: &[usize]) -> i32 {
    if indices.is_empty() {
        return 0;
    }
    // collect [start, end) intervals and merge.
    let mut intervals: Vec<(i32, i32)> = indices
        .iter()
        .map(|&i| (mems[i].read_start, mems[i].read_end()))
        .collect();
    intervals.sort_unstable();
    let mut covered = 0;
    let (mut cur_start, mut cur_end) = intervals[0];
    for &(s, e) in &intervals[1..] {
        if s > cur_end {
            covered += cur_end - cur_start;
            cur_start = s;
            cur_end = e;
        } else if e > cur_end {
            cur_end = e;
        }
    }
    covered += cur_end - cur_start;
    covered
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mem_geometry() {
        let m = Mem::new(5, 100, 20);
        assert_eq!(m.read_end(), 25);
        assert_eq!(m.ref_end(), 120);
        assert_eq!(m.diagonal(), 95);
    }

    #[test]
    fn coverage_disjoint() {
        let mems = [Mem::new(0, 0, 10), Mem::new(20, 20, 10)];
        assert_eq!(covered_read_bases(&mems, &[0, 1]), 20);
    }

    #[test]
    fn coverage_overlapping_counts_once() {
        // [0,15) and [10,25) -> union [0,25) = 25
        let mems = [Mem::new(0, 0, 15), Mem::new(10, 10, 15)];
        assert_eq!(covered_read_bases(&mems, &[0, 1]), 25);
    }
}
