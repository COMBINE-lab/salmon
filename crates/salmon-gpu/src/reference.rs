//! Pure-Rust reference of the banded affine-gap DP that the GPU kernel computes.
//!
//! salmon's full-length selective alignment scores a read against a reference
//! window with `ksw2`'s banded `extz2` (`KSW_EZ_RIGHT | KSW_EZ_SCORE_ONLY`),
//! reading the score reaching the query end (`mqe`), or the global max when the
//! query end is unreachable within the band. `ksw2` computes this with the
//! Suzuki-Kasahara *difference* recurrence (a `u/v/x/y` encoding tuned for
//! SIMD). That representation is numerically equivalent to a plain Gotoh DP run
//! over the *same set of band cells*: same scoring, same band, same extraction
//! give the same integer scores.
//!
//! This module is that plain Gotoh DP. It exists for two reasons:
//!   1. it is trivially portable to a WGSL compute shader (the GPU kernel is a
//!      near line-by-line transliteration of [`run_dp`]), and
//!   2. it is the bit-exact oracle the GPU backend is differential-tested
//!      against, and (where it matches `ksw2rs`, see the tests) it lets the GPU
//!      path transparently accelerate `--fullLengthAlignment`.
//!
//! The band and scoring mirror `ksw2`'s fast (non-`GENERIC_SC`) path, including
//! the detail that a wildcard base (`N`) scores `-gap_extend` rather than the
//! mismatch penalty.

/// Integer "negative infinity" floor, matching the magnitude `ksw2` uses
/// (`KSW_NEG_INF = -0x40000000`) so saturating adds never wrap.
pub const NEG_INF: i32 = -0x4000_0000;

/// Affine-gap scoring for the banded DP, as positive magnitudes (salmon's
/// `--ma/--mp/--go/--ge`) plus the band half-width (`--dpBandwidth`). A gap of
/// length `L` costs `gap_open + L * gap_extend`.
#[derive(Debug, Clone, Copy)]
pub struct BandedParams {
    pub match_score: i32,
    pub mismatch_pen: i32,
    pub gap_open_pen: i32,
    pub gap_extend_pen: i32,
    pub bandwidth: i32,
}

/// DNA5 code for a base: A=0, C=1, G=2, T=3, anything else (N) = 4.
#[inline]
pub fn dna5(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// Score of aligning query base `qb` against target base `tb` (both DNA5 codes),
/// matching `ksw2`'s fast-row scoring: a wildcard on either side scores
/// `-gap_extend`, otherwise match/mismatch.
#[inline]
fn sub_score(qb: u8, tb: u8, p: &BandedParams) -> i32 {
    const WILDCARD: u8 = 4;
    if qb == WILDCARD || tb == WILDCARD {
        -p.gap_extend_pen
    } else if qb == tb {
        p.match_score
    } else {
        -p.mismatch_pen
    }
}

/// `ksw2`'s per-anti-diagonal band: on anti-diagonal `r = iq + it` (0-based
/// query/target indices), the valid target columns are `[st, en]`. Returns the
/// inclusive `(st, en)` for this anti-diagonal, clamped to the matrix.
#[inline]
fn band_bounds(r: i32, qlen: i32, tlen: i32, w: i32) -> (i32, i32) {
    let mut st = 0;
    let mut en = tlen - 1;
    if st < r - qlen + 1 {
        st = r - qlen + 1;
    }
    if en > r {
        en = r;
    }
    if st < (r - w + 1) >> 1 {
        st = (r - w + 1) >> 1;
    }
    if en > (r + w) >> 1 {
        en = (r + w) >> 1;
    }
    (st, en)
}

#[inline]
fn sat_sub(a: i32, b: i32) -> i32 {
    if a <= NEG_INF / 2 {
        NEG_INF
    } else {
        a - b
    }
}

/// Banded affine-gap alignment score of `query` against `target` (raw ACGTN
/// bytes), reproducing salmon's `ksw2_align_score`: the score reaching the query
/// end within the band, or the global band max if the query end is unreachable.
pub fn banded_extz_score(query: &[u8], target: &[u8], p: &BandedParams) -> i32 {
    let q5: Vec<u8> = query.iter().map(|&b| dna5(b)).collect();
    let t5: Vec<u8> = target.iter().map(|&b| dna5(b)).collect();
    run_dp(&q5, &t5, p)
}

/// As [`banded_extz_score`] but on pre-encoded DNA5 slices (0..=4) — the exact
/// computation the WGSL kernel mirrors, kept separate so the encoding cost is
/// paid once when batching.
pub fn banded_extz_score_dna5(q5: &[u8], t5: &[u8], p: &BandedParams) -> i32 {
    run_dp(q5, t5, p)
}

/// Banded Gotoh DP with explicit running gap state, so each cell is O(1) and the
/// data flow matches the WGSL kernel one-to-one.
///
/// `H` is the alignment score ending at a cell; `E` is the running "gap in the
/// query" (deletion, consuming target) swept left-to-right; `F` is the
/// per-column "gap in the target" (insertion, consuming query) swept
/// top-to-bottom. A virtual column `it = -1` carries the leading-insertion edge.
fn run_dp(q5: &[u8], t5: &[u8], p: &BandedParams) -> i32 {
    let qn = q5.len();
    let tn = t5.len();
    if qn == 0 || tn == 0 {
        return NEG_INF;
    }
    let qlen = qn as i32;
    let tlen = tn as i32;
    let go = p.gap_open_pen;
    let ge = p.gap_extend_pen;
    let w = p.bandwidth;

    // H of the previous query row. Row iq = 0's "previous row" is the virtual
    // query row -1: the leading-deletion edge, H(-1, it) = -(go + (it+1)*ge)
    // (consume it+1 target bases with no query). Seeding this is what lets a read
    // align past a leading deletion of the reference window; a NEG_INF seed would
    // forbid that path and over-charge by a gap open.
    let mut h_prev: Vec<i32> = (0..tn).map(|it| -(go + (it as i32 + 1) * ge)).collect();
    let mut h_cur = vec![NEG_INF; tn]; // H of the current query row
    let mut f_col = vec![NEG_INF; tn]; // F per target column, carried across rows

    let mut global_max = NEG_INF;
    let mut mqe = NEG_INF;

    for iq in 0..qn {
        // Virtual column it = -1: an all-insertion path from the origin.
        // Diagonal source for it = 0 is H(iq-1, -1) (the origin when iq == 0).
        let mut diag_src = if iq == 0 { 0 } else { -(go + iq as i32 * ge) };
        let mut h_left = -(go + (iq as i32 + 1) * ge); // H(iq, it-1), seeded at H(iq,-1)
        let mut e_run = NEG_INF; // E at it = -1 is unreachable

        for it in 0..tn {
            let r = iq as i32 + it as i32;
            let (st, en) = band_bounds(r, qlen, tlen, w);
            let in_band = (it as i32) >= st && (it as i32) <= en;

            let diag = diag_src; // H(iq-1, it-1)
            diag_src = h_prev[it]; // becomes H(iq-1, it) for the next column's diagonal

            if !in_band {
                h_cur[it] = NEG_INF;
                f_col[it] = NEG_INF;
                e_run = NEG_INF;
                h_left = NEG_INF;
                continue;
            }

            let s = sub_score(q5[iq], t5[it], p);
            let m = if diag <= NEG_INF / 2 {
                NEG_INF
            } else {
                diag + s
            };

            // E: gap in query (consume target), from the left neighbour.
            e_run = sat_sub(h_left, go + ge).max(sat_sub(e_run, ge));
            // F: gap in target (consume query), from the cell above.
            let f = sat_sub(h_prev[it], go + ge).max(sat_sub(f_col[it], ge));
            f_col[it] = f;

            let h = m.max(e_run).max(f).max(NEG_INF);
            h_cur[it] = h;
            h_left = h;

            if h > global_max {
                global_max = h;
            }
            if iq + 1 == qn && h > mqe {
                mqe = h;
            }
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
    }

    // Match `ksw2_align_score`: return the score reaching the query end, else the
    // global band max. `ksw2`'s `ez.max` is a `u32` seeded at 0, so the fallback
    // is floored at 0 (an all-negative band reports 0, never a negative max).
    if mqe > NEG_INF {
        mqe
    } else {
        global_max.max(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: BandedParams = BandedParams {
        match_score: 2,
        mismatch_pen: 4,
        gap_open_pen: 6,
        gap_extend_pen: 2,
        bandwidth: 15,
    };

    #[test]
    fn exact_match_is_perfect() {
        let t = b"ACGTACGTTTGGCCAAACGTACAC";
        let q = &t[..16];
        assert_eq!(banded_extz_score(q, t, &P), 16 * 2);
    }

    #[test]
    fn single_mismatch() {
        let t = b"ACGTACGTTTGGCCAAACGTACAC";
        let mut q = t[..16].to_vec();
        q[8] = if q[8] == b'A' { b'C' } else { b'A' };
        // 15 matches, 1 mismatch
        assert_eq!(banded_extz_score(&q, t, &P), 15 * 2 - 4);
    }

    #[test]
    fn internal_deletion_costs_one_affine_gap() {
        // read = reference window with a 3-base internal deletion -> a single gap
        // of length 3 in the read: cost = gap_open + 3 * gap_extend.
        let t = b"ACGTACGTTTGGCCAAACGTACACGGGTTT";
        let mut q = t[..24].to_vec();
        q.drain(10..13); // remove 3 bases
        let perfect = q.len() as i32 * 2;
        let gap = 6 + 3 * 2;
        assert_eq!(banded_extz_score(&q, t, &P), perfect - gap);
    }

    #[test]
    fn leading_deletion_is_reachable() {
        // The read matches the reference only after skipping the first 4 bases of
        // the window: a leading deletion of 4. Regression guard for the virtual
        // top-row (query row -1) boundary seeding.
        let t = b"TTTTACGTACGTTTGGCCAAACGT";
        let q = &t[4..20]; // 16 bp, aligns after a 4-base leading deletion
        let perfect = 16 * 2;
        let gap = 6 + 4 * 2;
        // wide enough band to reach the off-diagonal start
        let p = BandedParams { bandwidth: 31, ..P };
        assert_eq!(banded_extz_score(q, t, &p), perfect - gap);
    }

    #[test]
    fn wildcard_scores_as_gap_extend() {
        // An N in the read scores -gap_extend (ksw2 fast-row rule), not -mismatch.
        let t = b"ACGTACGTTTGGCCAAACGTACAC";
        let mut q = t[..16].to_vec();
        q[8] = b'N';
        assert_eq!(banded_extz_score(&q, t, &P), 15 * 2 - 2);
    }
}
