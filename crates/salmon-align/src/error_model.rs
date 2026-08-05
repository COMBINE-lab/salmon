//! Alignment error model — a port of salmon's `AlignmentModel`
//! (`src/alignment/AlignmentModel.cpp`).
//!
//! # What it is for
//!
//! In alignment mode the aligner's own score (`AS`) is the only quality signal,
//! and different aligners compute it differently — so weighting placements by it
//! makes results depend on the aligner's scoring scheme. This model replaces it
//! with salmon's own judgement: given a reference sequence, how *likely* is this
//! particular pattern of matches, mismatches and gaps?
//!
//! The idea is that real sequencing errors have structure. They are not uniformly
//! distributed along the read (later cycles are worse), they depend on which base
//! was expected and which was seen, and they cluster. A model learned from the
//! data itself therefore distinguishes "a read from this transcript with typical
//! errors" from "a read from somewhere else that happens to align" far better
//! than a fixed penalty can.
//!
//! The score handed downstream is a *log-likelihood ratio*: how much more likely
//! the observed alignment is under the error model than under a perfect-match
//! baseline. That is a quantity salmon defines, so it is comparable across
//! aligners and reproducible.
//!
//! A position-binned, first-order Markov model over per-base *alignment states*.
//! Each aligned position is the pair `(refBase, readBase)` encoded as
//! `refBase*9 + readBase` over the 9 symbols {A,C,G,T,DASH,SOFT,HARD,PAD,REF_SKIP}
//! (81 states + a START state = 82). A separate `readBins × {left,right}` set of
//! 82×82 log-space transition matrices captures how error patterns vary along
//! the read and differ between mates. Walking an alignment's CIGAR against the
//! reference yields a state path whose log transition probability is the
//! foreground likelihood `fg`; the background `bg` accumulates the match
//! self-transition `(0,0)`. The per-alignment score used downstream is `fg − bg`.

/// Alignment symbols (salmon's `AlignmentModelChar`).
const ALN_DASH: usize = 4;
const ALN_SOFT_CLIP: usize = 5;
const ALN_HARD_CLIP: usize = 6;
const ALN_PAD: usize = 7;
const ALN_REF_SKIP: usize = 8;

const NUM_STATES: usize = 9;
/// `9*9 + 1` (the START state).
const NUM_ALN_STATES: usize = 82;
const START_STATE: usize = 81;

/// 2-bit encode a reference base (`A=0,C=1,G=2,T=3`; else 0).
#[inline]
fn ref_2bit(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// A CIGAR operation (the subset salmon distinguishes).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlnOp {
    Match,       // M
    SeqMatch,    // =
    SeqMismatch, // X
    Ins,         // I
    Del,         // D
    RefSkip,     // N
    SoftClip,    // S
    HardClip,    // H
    Pad,         // P
}

impl AlnOp {
    #[inline]
    fn consume_seq(self) -> bool {
        matches!(
            self,
            AlnOp::Match | AlnOp::SeqMatch | AlnOp::SeqMismatch | AlnOp::Ins | AlnOp::SoftClip
        )
    }
    #[inline]
    fn consume_ref(self) -> bool {
        matches!(
            self,
            AlnOp::Match | AlnOp::SeqMatch | AlnOp::SeqMismatch | AlnOp::Del | AlnOp::RefSkip
        )
    }
    /// salmon's `setBasesFromCIGAROp_`: adjust the (ref, read) symbols for the op.
    #[inline]
    fn set_bases(self, cur_ref: &mut usize, cur_read: &mut usize) {
        match self {
            AlnOp::Ins => *cur_ref = ALN_DASH,
            AlnOp::Del => *cur_read = ALN_DASH,
            AlnOp::RefSkip => *cur_read = ALN_REF_SKIP,
            AlnOp::SoftClip => *cur_ref = ALN_SOFT_CLIP,
            AlnOp::HardClip => {
                *cur_ref = ALN_HARD_CLIP;
                *cur_read = ALN_HARD_CLIP;
            }
            AlnOp::Pad => {
                *cur_ref = ALN_PAD;
                *cur_read = ALN_PAD;
            }
            _ => {}
        }
    }
}

#[inline]
fn log_add(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (lo - hi).exp().ln_1p()
}

/// An 82×82 log-space transition matrix with maintained per-row sums, so a
/// query returns the (live) normalized log transition probability. Mirrors
/// salmon's `AtomicMatrix` (single-threaded here; training is sequential).
#[derive(Clone)]
struct TransMatrix {
    storage: Vec<f64>, // NUM_ALN_STATES * NUM_ALN_STATES, log space
    rowsums: Vec<f64>,
}

impl TransMatrix {
    #[allow(dead_code)] // retained for unit tests / API symmetry
    fn new(alpha: f64) -> Self {
        Self {
            storage: vec![alpha.ln(); NUM_ALN_STATES * NUM_ALN_STATES],
            rowsums: vec![(NUM_ALN_STATES as f64 * alpha).ln(); NUM_ALN_STATES],
        }
    }
    /// A zero-mass (`-inf` everywhere) delta matrix: a per-thread accumulator
    /// whose increments are later `log_add`-folded into an alpha-seeded global
    /// (so the pseudocount baseline is counted exactly once, in the global).
    fn empty() -> Self {
        Self {
            storage: vec![f64::NEG_INFINITY; NUM_ALN_STATES * NUM_ALN_STATES],
            rowsums: vec![f64::NEG_INFINITY; NUM_ALN_STATES],
        }
    }
    #[inline]
    fn increment(&mut self, prev: usize, cur: usize, amt: f64) {
        let k = prev * NUM_ALN_STATES + cur;
        self.storage[k] = log_add(self.storage[k], amt);
        self.rowsums[prev] = log_add(self.rowsums[prev], amt);
    }
    /// Reset all masses to `-inf` (log-zero) in place, without reallocating, so
    /// a delta accumulator can be reused after flushing.
    fn clear(&mut self) {
        self.storage.iter_mut().for_each(|v| *v = f64::NEG_INFINITY);
        self.rowsums.iter_mut().for_each(|v| *v = f64::NEG_INFINITY);
    }
    /// Fold another matrix's masses into this one element-wise (`log_add`).
    /// Combining an alpha-seeded global with `empty()`-seeded per-thread deltas
    /// reconstructs `alpha + Σ deltas` exactly.
    fn combine(&mut self, other: &TransMatrix) {
        for (s, o) in self.storage.iter_mut().zip(&other.storage) {
            *s = log_add(*s, *o);
        }
        for (s, o) in self.rowsums.iter_mut().zip(&other.rowsums) {
            *s = log_add(*s, *o);
        }
    }
    /// Normalized log transition probability `P(cur | prev)`.
    #[inline]
    #[allow(dead_code)] // retained for unit tests / API symmetry
    fn get(&self, prev: usize, cur: usize) -> f64 {
        self.storage[prev * NUM_ALN_STATES + cur] - self.rowsums[prev]
    }
}

/// The alignment error model: `read_bins` position bins × {left, right} mate.
#[derive(Clone)]
pub struct AlignmentModel {
    left: Vec<TransMatrix>,
    right: Vec<TransMatrix>,
    read_bins: usize,
}

impl AlignmentModel {
    /// salmon's default: `alpha` pseudocount per cell, 4 read-position bins.
    #[allow(dead_code)] // retained for unit tests / API symmetry
    pub fn new(alpha: f64, read_bins: usize) -> Self {
        Self {
            left: (0..read_bins).map(|_| TransMatrix::new(alpha)).collect(),
            right: (0..read_bins).map(|_| TransMatrix::new(alpha)).collect(),
            read_bins,
        }
    }

    /// Reset all per-bin matrices to log-zero in place (reuse after flushing).
    pub fn clear(&mut self) {
        self.left.iter_mut().for_each(|m| m.clear());
        self.right.iter_mut().for_each(|m| m.clear());
    }

    /// A zero-mass per-thread delta accumulator with the same shape. Its updates
    /// are merged into an alpha-seeded global via [`AlignmentModel::combine`].
    pub fn empty(read_bins: usize) -> Self {
        Self {
            left: (0..read_bins).map(|_| TransMatrix::empty()).collect(),
            right: (0..read_bins).map(|_| TransMatrix::empty()).collect(),
            read_bins,
        }
    }

    /// Fold another model's masses into this one (element-wise `log_add` over
    /// every position bin and mate). Used to merge per-thread deltas.
    pub fn combine(&mut self, other: &AlignmentModel) {
        for (s, o) in self.left.iter_mut().zip(&other.left) {
            s.combine(o);
        }
        for (s, o) in self.right.iter_mut().zip(&other.right) {
            s.combine(o);
        }
    }

    /// Walk an alignment's CIGAR producing the `(read_pos_bin, prev, cur)`
    /// transitions, invoking `f(bin, prev, cur)` for each. `read_2bit` is the
    /// read's 2-bit bases (reference-forward orientation, as stored in the BAM);
    /// `ref_bytes` the transcript's ASCII bases; `pos` the 0-based alignment start.
    fn walk<F: FnMut(usize, usize, usize)>(
        read_bins: usize,
        read_2bit: &[u8],
        ref_bytes: &[u8],
        pos: usize,
        ops: &[(AlnOp, usize)],
        mut f: F,
    ) {
        let read_len = read_2bit.len();
        if read_len == 0 || ref_bytes.is_empty() {
            return;
        }
        let inv_len = read_bins as f64 / read_len as f64;
        let mut read_idx = 0usize;
        let mut ref_idx = pos;
        let mut prev = START_STATE;
        for &(op, op_len) in ops {
            for _ in 0..op_len {
                if op.consume_seq() && read_idx >= read_len {
                    return; // inconsistent CIGAR
                }
                if op.consume_ref() && ref_idx >= ref_bytes.len() {
                    return;
                }
                let mut cur_read = if op.consume_seq() {
                    read_2bit[read_idx] as usize
                } else {
                    0
                };
                let mut cur_ref = if op.consume_ref() {
                    ref_2bit(ref_bytes[ref_idx])
                } else {
                    0
                };
                op.set_bases(&mut cur_ref, &mut cur_read);
                let bin = ((read_idx as f64 * inv_len) as usize).min(read_bins - 1);
                let cur = cur_ref * NUM_STATES + cur_read;
                f(bin, prev, cur);
                prev = cur;
                if op.consume_seq() {
                    read_idx += 1;
                }
                if op.consume_ref() {
                    ref_idx += 1;
                }
            }
        }
    }

    /// Accumulate the alignment's transitions into the model, weighted by
    /// `log_weight` (log space). `is_left` selects the mate's matrices.
    pub fn update(
        &mut self,
        read_2bit: &[u8],
        ref_bytes: &[u8],
        pos: usize,
        ops: &[(AlnOp, usize)],
        is_left: bool,
        log_weight: f64,
    ) {
        // Collect transitions first (immutable walk), then apply (mutable).
        let mut trans: Vec<(usize, usize, usize)> = Vec::new();
        Self::walk(
            self.read_bins,
            read_2bit,
            ref_bytes,
            pos,
            ops,
            |bin, prev, cur| {
                trans.push((bin, prev, cur));
            },
        );
        let mats = if is_left {
            &mut self.left
        } else {
            &mut self.right
        };
        for (bin, prev, cur) in trans {
            mats[bin].increment(prev, cur, log_weight);
        }
    }

    /// Foreground/background log-likelihoods `(fg, bg)` of the alignment under
    /// the current model. The per-alignment score is `fg − bg`.
    #[allow(dead_code)] // retained for unit tests; the live path reads SharedAlignmentModel
    pub fn log_likelihood(
        &self,
        read_2bit: &[u8],
        ref_bytes: &[u8],
        pos: usize,
        ops: &[(AlnOp, usize)],
        is_left: bool,
    ) -> (f64, f64) {
        let mats = if is_left { &self.left } else { &self.right };
        let mut fg = 0.0; // LOG_1
        let mut bg = 0.0;
        Self::walk(
            self.read_bins,
            read_2bit,
            ref_bytes,
            pos,
            ops,
            |bin, prev, cur| {
                fg += mats[bin].get(prev, cur);
                bg += mats[bin].get(0, 0);
            },
        );
        (fg, bg)
    }
}

/// An **order-independent** counting error model: a first-order transition
/// *count* per `(read_bin, is_left, prev, cur)` accumulated as `u64` integers.
/// Integer addition is associative and commutative, so merging per-thread
/// counters ([`combine`](Self::combine)) yields the same totals regardless of
/// how fragments were partitioned across threads — the property the online
/// [`SharedAlignmentModel`] cannot offer (its log-space, posterior-weighted,
/// concurrently-flushed updates depend on processing order). Once training is
/// done, [`finalize`](Self::finalize) normalizes the counts once (adding the
/// `alpha` pseudocount exactly, in the denominator too) into a fixed log-space
/// [`AlignmentModel`], so scoring reuses the same `log_likelihood` walk. This is
/// the fidelity model behind deterministic alignment mode (`--deterministic`
/// with `-t`): trained uniformly (one count per placement-mate, no posterior
/// weighting) so the whole pipeline stays byte-reproducible.
#[derive(Clone)]
pub struct CountingAlignmentModel {
    // read_bins vectors, each NUM_ALN_STATES*NUM_ALN_STATES transition counts.
    left: Vec<Vec<u64>>,
    right: Vec<Vec<u64>>,
    read_bins: usize,
}

impl CountingAlignmentModel {
    /// A zeroed counting model with `read_bins` position bins per mate.
    pub fn new(read_bins: usize) -> Self {
        let cells = NUM_ALN_STATES * NUM_ALN_STATES;
        Self {
            left: (0..read_bins).map(|_| vec![0u64; cells]).collect(),
            right: (0..read_bins).map(|_| vec![0u64; cells]).collect(),
            read_bins,
        }
    }

    /// Tally one alignment's transitions (each `+1`). `is_left` selects the
    /// mate's matrices; args match [`AlignmentModel::update`] minus the weight.
    pub fn count(
        &mut self,
        read_2bit: &[u8],
        ref_bytes: &[u8],
        pos: usize,
        ops: &[(AlnOp, usize)],
        is_left: bool,
    ) {
        let mats = if is_left {
            &mut self.left
        } else {
            &mut self.right
        };
        AlignmentModel::walk(
            self.read_bins,
            read_2bit,
            ref_bytes,
            pos,
            ops,
            |bin, prev, cur| {
                mats[bin][prev * NUM_ALN_STATES + cur] += 1;
            },
        );
    }

    /// Add another counting model's transition counts into this one (integer,
    /// order-independent). Both must share the same `read_bins`.
    pub fn combine(&mut self, other: &CountingAlignmentModel) {
        debug_assert_eq!(self.read_bins, other.read_bins);
        for (s, o) in self.left.iter_mut().zip(&other.left) {
            for (a, b) in s.iter_mut().zip(o) {
                *a += *b;
            }
        }
        for (s, o) in self.right.iter_mut().zip(&other.right) {
            for (a, b) in s.iter_mut().zip(o) {
                *a += *b;
            }
        }
    }

    /// Normalize the integer counts into a fixed log-space [`AlignmentModel`],
    /// adding an `alpha` Laplace pseudocount to every cell (and `NUM_ALN_STATES *
    /// alpha` to each row denominator) so unobserved transitions keep a small
    /// mass — matching the online model's alpha-seeded matrices. The result is a
    /// deterministic function of the totals alone.
    pub fn finalize(&self, alpha: f64) -> AlignmentModel {
        let build = |counts: &[Vec<u64>]| -> Vec<TransMatrix> {
            counts
                .iter()
                .map(|cell| {
                    let mut storage = vec![0.0f64; NUM_ALN_STATES * NUM_ALN_STATES];
                    let mut rowsums = vec![0.0f64; NUM_ALN_STATES];
                    for prev in 0..NUM_ALN_STATES {
                        let base = prev * NUM_ALN_STATES;
                        let rowcount: u64 = cell[base..base + NUM_ALN_STATES].iter().sum();
                        rowsums[prev] = (rowcount as f64 + NUM_ALN_STATES as f64 * alpha).ln();
                        for cur in 0..NUM_ALN_STATES {
                            storage[base + cur] = (cell[base + cur] as f64 + alpha).ln();
                        }
                    }
                    TransMatrix { storage, rowsums }
                })
                .collect()
        };
        AlignmentModel {
            left: build(&self.left),
            right: build(&self.right),
            read_bins: self.read_bins,
        }
    }
}

/// An atomic, log-space transition matrix shared across worker threads. Reads
/// (`get`) are lock-free relaxed loads (free on x86); updates arrive as
/// occasional bulk merges of a per-thread plain [`TransMatrix`] delta, so the
/// hot match-state cell is not CAS-contended on every base.
struct SharedTransMatrix {
    storage: Vec<salmon_core::atomic::AtomicF64>,
    rowsums: Vec<salmon_core::atomic::AtomicF64>,
}

impl SharedTransMatrix {
    fn new(alpha: f64) -> Self {
        Self {
            storage: (0..NUM_ALN_STATES * NUM_ALN_STATES)
                .map(|_| salmon_core::atomic::AtomicF64::new(alpha.ln()))
                .collect(),
            rowsums: (0..NUM_ALN_STATES)
                .map(|_| salmon_core::atomic::AtomicF64::new((NUM_ALN_STATES as f64 * alpha).ln()))
                .collect(),
        }
    }
    #[inline]
    fn get(&self, prev: usize, cur: usize) -> f64 {
        self.storage[prev * NUM_ALN_STATES + cur].load() - self.rowsums[prev].load()
    }
    /// Fold a per-thread plain delta into this shared matrix atomically (only the
    /// touched, non-`-inf` cells).
    fn flush_from(&self, delta: &TransMatrix) {
        for (a, &d) in self.storage.iter().zip(&delta.storage) {
            if d != f64::NEG_INFINITY {
                a.log_add_assign(d);
            }
        }
        for (a, &d) in self.rowsums.iter().zip(&delta.rowsums) {
            if d != f64::NEG_INFINITY {
                a.log_add_assign(d);
            }
        }
    }
}

/// Shared, atomic counterpart of [`AlignmentModel`]: read concurrently for the
/// foreground/background likelihood (`basis`) during the online pass while
/// worker threads periodically flush their private deltas into it. Matches
/// salmon's shared `AtomicMatrix` error model, but with batched (per-minibatch)
/// flushing to keep update contention negligible.
pub struct SharedAlignmentModel {
    left: Vec<SharedTransMatrix>,
    right: Vec<SharedTransMatrix>,
    read_bins: usize,
}

impl SharedAlignmentModel {
    pub fn new(alpha: f64, read_bins: usize) -> Self {
        Self {
            left: (0..read_bins)
                .map(|_| SharedTransMatrix::new(alpha))
                .collect(),
            right: (0..read_bins)
                .map(|_| SharedTransMatrix::new(alpha))
                .collect(),
            read_bins,
        }
    }

    /// Foreground/background log-likelihoods `(fg, bg)` under the current shared
    /// model (lock-free atomic reads).
    pub fn log_likelihood(
        &self,
        read_2bit: &[u8],
        ref_bytes: &[u8],
        pos: usize,
        ops: &[(AlnOp, usize)],
        is_left: bool,
    ) -> (f64, f64) {
        let mats = if is_left { &self.left } else { &self.right };
        let mut fg = 0.0;
        let mut bg = 0.0;
        AlignmentModel::walk(
            self.read_bins,
            read_2bit,
            ref_bytes,
            pos,
            ops,
            |bin, prev, cur| {
                fg += mats[bin].get(prev, cur);
                bg += mats[bin].get(0, 0);
            },
        );
        (fg, bg)
    }

    /// Atomically merge a per-thread plain delta model into the shared model.
    pub fn flush_from(&self, delta: &AlignmentModel) {
        for (s, d) in self.left.iter().zip(&delta.left) {
            s.flush_from(d);
        }
        for (s, d) in self.right.iter().zip(&delta.right) {
            s.flush_from(d);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A perfect 5-base match: A C G T A against the same reference.
    fn perfect() -> (Vec<u8>, Vec<u8>, Vec<(AlnOp, usize)>) {
        let read = vec![0u8, 1, 2, 3, 0];
        let refs = b"ACGTA".to_vec();
        (read, refs, vec![(AlnOp::Match, 5)])
    }

    #[test]
    /// Once trained on real alignments, a clean alignment must score above a
    /// mismatched one — the basic property the whole model exists to provide.
    fn matches_score_higher_than_mismatches_after_training() {
        let mut m = AlignmentModel::new(1.0, 4);
        let (read, refs, ops) = perfect();
        // train on many perfect matches
        for _ in 0..200 {
            m.update(&read, &refs, 0, &ops, true, 0.0);
        }
        // a read with a mismatch (read[2] G->A vs ref G)
        let mut bad = read.clone();
        bad[2] = 0; // A where ref is G
        let (fg_good, bg_good) = m.log_likelihood(&read, &refs, 0, &ops, true);
        let (fg_bad, bg_bad) = m.log_likelihood(&bad, &refs, 0, &ops, true);
        // the all-match read has higher foreground likelihood ...
        assert!(fg_good > fg_bad, "fg_good {fg_good} !> fg_bad {fg_bad}");
        // ... and a higher final score (fg - bg); bg is identical (same #positions)
        assert!(
            (fg_good - bg_good) > (fg_bad - bg_bad),
            "perfect score {} !> mismatch score {}",
            fg_good - bg_good,
            fg_bad - bg_bad
        );
    }

    #[test]
    /// Before training the model must not express a preference, so an early
    /// fragment is not scored against an opinion the model has not yet earned.
    fn untrained_model_is_neutral() {
        // With no training every transition is the uniform prior, so fg == bg
        // (score 0): no alignment is preferred until errors are observed.
        let m = AlignmentModel::new(1.0, 4);
        let (read, refs, ops) = perfect();
        let (fg, bg) = m.log_likelihood(&read, &refs, 0, &ops, true);
        assert!((fg - bg).abs() < 1e-9, "untrained score {} not 0", fg - bg);
    }

    #[test]
    /// Same obligation for the deterministic counting variant.
    fn counting_model_prefers_matches_after_training() {
        // The integer counting model, once normalized, ranks a clean match above
        // a mismatch just like the float model.
        let mut c = CountingAlignmentModel::new(4);
        let (read, refs, ops) = perfect();
        for _ in 0..200 {
            c.count(&read, &refs, 0, &ops, true);
        }
        let m = c.finalize(1.0);
        let mut bad = read.clone();
        bad[2] = 0; // A where ref is G
        let (fg_good, bg_good) = m.log_likelihood(&read, &refs, 0, &ops, true);
        let (fg_bad, bg_bad) = m.log_likelihood(&bad, &refs, 0, &ops, true);
        assert!(
            (fg_good - bg_good) > (fg_bad - bg_bad),
            "match score {} !> mismatch score {}",
            fg_good - bg_good,
            fg_bad - bg_bad
        );
    }

    #[test]
    /// Per-thread models must merge to the same result whatever order the threads
    /// finish in; otherwise the deterministic mode would not be deterministic.
    fn counting_model_merge_is_order_independent() {
        // Splitting the same observations across two counters and merging in
        // either order yields byte-identical normalized transition tables — the
        // determinism guarantee for multi-threaded training.
        let (read, refs, ops) = perfect();
        let mut bad = read.clone();
        bad[2] = 0;
        let mut a = CountingAlignmentModel::new(4);
        let mut b = CountingAlignmentModel::new(4);
        for _ in 0..120 {
            a.count(&read, &refs, 0, &ops, true);
        }
        for _ in 0..80 {
            b.count(&bad, &refs, 0, &ops, true);
            b.count(&read, &refs, 0, &ops, false);
        }
        let mut ab = a.clone();
        ab.combine(&b);
        let mut ba = b.clone();
        ba.combine(&a);
        // Same totals ⇒ identical log-space matrices ⇒ identical likelihoods.
        let ll = |m: &AlignmentModel, is_left| m.log_likelihood(&read, &refs, 0, &ops, is_left);
        let (m_ab, m_ba) = (ab.finalize(1.0), ba.finalize(1.0));
        assert_eq!(ll(&m_ab, true), ll(&m_ba, true));
        assert_eq!(ll(&m_ab, false), ll(&m_ba, false));
        // And equal to a single counter fed everything.
        let mut all = CountingAlignmentModel::new(4);
        all.combine(&a);
        all.combine(&b);
        let m_all = all.finalize(1.0);
        assert_eq!(ll(&m_ab, true), ll(&m_all, true));
    }
}
