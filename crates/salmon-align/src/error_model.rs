//! Alignment error model — a port of salmon's `AlignmentModel`
//! (`src/alignment/AlignmentModel.cpp`).
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

/// A CIGAR operation (the subset salmon distinguishes).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlnOp {
    Match,        // M
    SeqMatch,     // =
    SeqMismatch,  // X
    Ins,          // I
    Del,          // D
    RefSkip,      // N
    SoftClip,     // S
    HardClip,     // H
    Pad,          // P
}

impl AlnOp {
    #[inline]
    fn consume_seq(self) -> bool {
        matches!(self, AlnOp::Match | AlnOp::SeqMatch | AlnOp::SeqMismatch | AlnOp::Ins | AlnOp::SoftClip)
    }
    #[inline]
    fn consume_ref(self) -> bool {
        matches!(self, AlnOp::Match | AlnOp::SeqMatch | AlnOp::SeqMismatch | AlnOp::Del | AlnOp::RefSkip)
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
    fn new(alpha: f64) -> Self {
        Self {
            storage: vec![alpha.ln(); NUM_ALN_STATES * NUM_ALN_STATES],
            rowsums: vec![(NUM_ALN_STATES as f64 * alpha).ln(); NUM_ALN_STATES],
        }
    }
    #[inline]
    fn increment(&mut self, prev: usize, cur: usize, amt: f64) {
        let k = prev * NUM_ALN_STATES + cur;
        self.storage[k] = log_add(self.storage[k], amt);
        self.rowsums[prev] = log_add(self.rowsums[prev], amt);
    }
    /// Normalized log transition probability `P(cur | prev)`.
    #[inline]
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
    pub fn new(alpha: f64, read_bins: usize) -> Self {
        Self {
            left: (0..read_bins).map(|_| TransMatrix::new(alpha)).collect(),
            right: (0..read_bins).map(|_| TransMatrix::new(alpha)).collect(),
            read_bins,
        }
    }

    /// Walk an alignment's CIGAR producing the `(read_pos_bin, prev, cur)`
    /// transitions, invoking `f(bin, prev, cur)` for each. `read_2bit` is the
    /// read's 2-bit bases (reference-forward orientation, as stored in the BAM);
    /// `ref_2bit` the transcript's 2-bit bases; `pos` the 0-based alignment start.
    fn walk<F: FnMut(usize, usize, usize)>(
        &self,
        read_2bit: &[u8],
        ref_2bit: &[u8],
        pos: usize,
        ops: &[(AlnOp, usize)],
        mut f: F,
    ) {
        let read_len = read_2bit.len();
        if read_len == 0 || ref_2bit.is_empty() {
            return;
        }
        let inv_len = self.read_bins as f64 / read_len as f64;
        let mut read_idx = 0usize;
        let mut ref_idx = pos;
        let mut prev = START_STATE;
        for &(op, op_len) in ops {
            for _ in 0..op_len {
                if op.consume_seq() && read_idx >= read_len {
                    return; // inconsistent CIGAR
                }
                if op.consume_ref() && ref_idx >= ref_2bit.len() {
                    return;
                }
                let mut cur_read = if op.consume_seq() { read_2bit[read_idx] as usize } else { 0 };
                let mut cur_ref = if op.consume_ref() { ref_2bit[ref_idx] as usize } else { 0 };
                op.set_bases(&mut cur_ref, &mut cur_read);
                let bin = ((read_idx as f64 * inv_len) as usize).min(self.read_bins - 1);
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
        ref_2bit: &[u8],
        pos: usize,
        ops: &[(AlnOp, usize)],
        is_left: bool,
        log_weight: f64,
    ) {
        // Collect transitions first (immutable walk), then apply (mutable).
        let mut trans: Vec<(usize, usize, usize)> = Vec::new();
        self.walk(read_2bit, ref_2bit, pos, ops, |bin, prev, cur| {
            trans.push((bin, prev, cur));
        });
        let mats = if is_left { &mut self.left } else { &mut self.right };
        for (bin, prev, cur) in trans {
            mats[bin].increment(prev, cur, log_weight);
        }
    }

    /// Foreground/background log-likelihoods `(fg, bg)` of the alignment under
    /// the current model. The per-alignment score is `fg − bg`.
    pub fn log_likelihood(
        &self,
        read_2bit: &[u8],
        ref_2bit: &[u8],
        pos: usize,
        ops: &[(AlnOp, usize)],
        is_left: bool,
    ) -> (f64, f64) {
        let mats = if is_left { &self.left } else { &self.right };
        let mut fg = 0.0; // LOG_1
        let mut bg = 0.0;
        self.walk(read_2bit, ref_2bit, pos, ops, |bin, prev, cur| {
            fg += mats[bin].get(prev, cur);
            bg += mats[bin].get(0, 0);
        });
        (fg, bg)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A perfect 5-base match: A C G T A against the same reference.
    fn perfect() -> (Vec<u8>, Vec<u8>, Vec<(AlnOp, usize)>) {
        let read = vec![0u8, 1, 2, 3, 0];
        let refs = vec![0u8, 1, 2, 3, 0];
        (read, refs, vec![(AlnOp::Match, 5)])
    }

    #[test]
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
    fn untrained_model_is_neutral() {
        // With no training every transition is the uniform prior, so fg == bg
        // (score 0): no alignment is preferred until errors are observed.
        let m = AlignmentModel::new(1.0, 4);
        let (read, refs, ops) = perfect();
        let (fg, bg) = m.log_likelihood(&read, &refs, 0, &ops, true);
        assert!((fg - bg).abs() < 1e-9, "untrained score {} not 0", fg - bg);
    }
}
