//! Alignment validation of mapping candidates.
//!
//! Salmon's selective alignment doesn't trust a MEM chain on its own — it
//! validates each candidate by aligning the read against the reference window
//! the chain implies, and keeps the candidate only if the alignment score
//! clears `minScoreFraction * perfectScore`. This mirrors salmon's
//! `PuffAligner::calculateAlignments`, scored with [`ksw2rs`] (a Rust port of
//! `ksw_extz2_sse` with runtime-dispatched SIMD kernels), matching C++ salmon's
//! banded scoring.
//!
//! The wrapper takes the reference sequence bytes as input so it has no
//! dependency on where those bytes are stored; the full mapper supplies them.
//!
//! ## Scoring
//! salmon's affine scoring is given as positive magnitudes (`--ma/--mp/--go/--ge`);
//! ksw2 charges `gapo + gape * len` per gap, which we pass through directly.

use crate::chain::MemChain;

/// Affine-gap scoring and validation thresholds (salmon `--ma/--mp/--go/--ge`,
/// `--minScoreFraction`). Penalties are positive magnitudes.
#[derive(Debug, Clone)]
pub struct AlignConfig {
    pub match_score: i8,
    pub mismatch_pen: i8,
    pub gap_open_pen: i8,
    pub gap_extend_pen: i8,
    /// minimum score as a fraction of the perfect score for a valid alignment
    pub min_score_fraction: f32,
    /// extra reference bases beyond the read length to include in the window,
    /// giving the aligner room for net deletions in the read
    pub indel_margin: usize,
    /// ksw2 DP bandwidth (salmon's `dpBandwidth`)
    pub bandwidth: i32,
    /// score the full read with one DP (`true`) vs. PuffAligner-style scoring
    /// of only the inter-MEM gaps and flanks, treating MEMs as exact matches
    /// (`false`, salmon's default = faster). Salmon's `--fullLengthAlignment`.
    pub full_length_alignment: bool,
}

impl Default for AlignConfig {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_pen: 4,
            gap_open_pen: 6,
            gap_extend_pen: 2,
            min_score_fraction: 0.65,
            indel_margin: 32,
            bandwidth: 15,
            // anchored (inter-MEM-gap) scoring by default, matching salmon.
            full_length_alignment: false,
        }
    }
}

/// Result of validating a mapping candidate by alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Alignment {
    /// optimal alignment score over the window
    pub score: i32,
    /// whether `score >= min_score_fraction * perfect_score`
    pub valid: bool,
    /// absolute reference position where the aligned window began
    pub ref_window_start: i32,
}

/// Reverse-complement a DNA byte slice (ACGT, case-insensitive output upper).
pub(crate) fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
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

/// The perfect (all-match) score for a read of length `read_len`.
#[inline]
pub fn perfect_score(read_len: usize, cfg: &AlignConfig) -> i32 {
    read_len as i32 * cfg.match_score as i32
}

/// Validate a mapping candidate by banded alignment of the read against the
/// reference window implied by the chain.
///
/// `ref_seq` is the full forward sequence of the candidate's reference.
/// `chain.is_fw` selects whether the read or its reverse complement is the
/// query (the chain is always expressed in the reference-forward frame).
/// Returns `None` if a usable window cannot be formed (e.g. read shorter than
/// nothing, or the window is empty).
pub fn align_chain(
    read: &[u8],
    ref_seq: &[u8],
    chain: &MemChain,
    cfg: &AlignConfig,
) -> Option<Alignment> {
    if read.is_empty() || ref_seq.is_empty() {
        return None;
    }

    // Query in the reference-forward frame.
    let query: Vec<u8> = if chain.is_fw {
        read.to_vec()
    } else {
        revcomp(read)
    };
    let qlen = query.len() as i32;
    let diag_origin = chain.ref_start() - chain.read_start();
    let win_start = diag_origin.max(0);

    let score = if cfg.full_length_alignment {
        // Full-read DP over the implied window.
        let win_end = (win_start + qlen + cfg.indel_margin as i32).min(ref_seq.len() as i32);
        if win_end <= win_start {
            return None;
        }
        let rwin = &ref_seq[win_start as usize..win_end as usize];
        ksw2_align_score(&query, rwin, cfg)
    } else {
        // PuffAligner-style: exact-MEM segments + DP of inter-MEM gaps/flanks.
        anchored_align_score(&query, ref_seq, &chain.mems, cfg)
    };

    let perfect = perfect_score(query.len(), cfg);
    let valid = (score as f32) >= cfg.min_score_fraction * (perfect as f32);

    Some(Alignment {
        score,
        valid,
        ref_window_start: win_start,
    })
}

thread_local! {
    /// Per-thread cache of gap/flank DP scores, keyed by the exact (query, ref)
    /// substring pair — salmon's alignment cache, to avoid re-running identical
    /// DPs (common for shared gaps/flanks across the candidate placements of one
    /// read). salmon scopes this cache per read; we instead bound it by size
    /// (see [`GAP_CACHE_CAP`]) so it cannot grow without limit over a whole run
    /// — an unbounded thread-local here cost tens of GB on a 36M-read library.
    static GAP_CACHE: std::cell::RefCell<ahash::AHashMap<(Box<[u8]>, Box<[u8]>), i32>> =
        std::cell::RefCell::new(ahash::AHashMap::new());
}

/// Maximum number of cached (query, ref) DP scores per thread. The cache exists
/// to dedup identical inter-MEM gap / flank alignments across one read's
/// candidate targets, where the working set is tiny; cross-read reuse is rare
/// (flanks are read-specific). When the cap is hit we clear the whole cache,
/// keeping peak memory at a few MB per thread instead of unbounded growth.
const GAP_CACHE_CAP: usize = 32_768;

/// Insert into the gap/flank cache, clearing it first if it has reached the cap.
#[inline]
fn gap_cache_insert(
    c: &std::cell::RefCell<ahash::AHashMap<(Box<[u8]>, Box<[u8]>), i32>>,
    key: (Box<[u8]>, Box<[u8]>),
    score: i32,
) {
    let mut m = c.borrow_mut();
    if m.len() >= GAP_CACHE_CAP {
        m.clear();
    }
    m.insert(key, score);
}

/// PuffAligner-style score: each MEM contributes an exact-match score and only
/// the inter-MEM gaps and the read flanks are DP-aligned (with ksw2), reusing
/// cached results for repeated substrings. Treating MEMs as fixed exact matches
/// is what makes this cheaper than a full-read DP; it can be marginally lower
/// than the global optimum when re-aligning across a MEM boundary would help.
fn anchored_align_score(query: &[u8], ref_seq: &[u8], mems: &[crate::mem::Mem], cfg: &AlignConfig) -> i32 {
    if mems.is_empty() {
        return i32::MIN / 4;
    }
    let qlen = query.len() as i32;
    let reflen = ref_seq.len() as i32;
    let m = cfg.match_score as i32;
    let mut score = 0i32;

    // 5' flank: read before the first MEM, anchored on its right at the MEM.
    let first = &mems[0];
    if first.read_start > 0 {
        let qf = &query[0..first.read_start as usize];
        let r_lo = (first.ref_start - first.read_start - cfg.indel_margin as i32).max(0);
        let r_hi = first.ref_start.max(0);
        let tf = &ref_seq[r_lo as usize..r_hi as usize];
        score += cached_flank_score(qf, tf, cfg, /*anchor_right=*/ true);
    }

    // MEM exact-match segments + inter-MEM gaps.
    for (i, mem) in mems.iter().enumerate() {
        score += mem.len * m;
        if i + 1 < mems.len() {
            let nxt = &mems[i + 1];
            let (rg_s, rg_e) = (mem.read_end(), nxt.read_start);
            let (tg_s, tg_e) = (mem.ref_end(), nxt.ref_start);
            if rg_e >= rg_s && tg_e >= tg_s {
                let qg = &query[rg_s as usize..rg_e as usize];
                let tg = &ref_seq[tg_s as usize..tg_e as usize];
                score += cached_gap_score(qg, tg, cfg);
            } else {
                // Overlapping MEMs: remove the double-counted exact bases.
                let ov = (rg_s - rg_e).max(tg_s - tg_e).max(0);
                score -= ov * m;
            }
        }
    }

    // 3' flank: read after the last MEM, anchored on its left at the MEM.
    let last = mems.last().unwrap();
    if last.read_end() < qlen {
        let qf = &query[last.read_end() as usize..qlen as usize];
        let r_hi = (last.ref_end() + (qlen - last.read_end()) + cfg.indel_margin as i32).min(reflen);
        let tf = &ref_seq[last.ref_end() as usize..r_hi as usize];
        score += cached_flank_score(qf, tf, cfg, /*anchor_right=*/ false);
    }

    score
}

/// 5x5 DNA5 scoring matrix (match diagonal, mismatch off-diagonal, N-N = 0).
fn dna5_mat(cfg: &AlignConfig) -> [i8; 25] {
    let mut mat = [-cfg.mismatch_pen; 25];
    for i in 0..5 {
        mat[i * 5 + i] = cfg.match_score;
    }
    mat[24] = 0;
    mat
}

/// Global DP score of an inter-MEM gap (both substrings fully consumed). The
/// band is widened to guarantee the alignment can reach the corner.
fn ksw2_gap_global(qg: &[u8], tg: &[u8], cfg: &AlignConfig) -> i32 {
    use ksw2rs::{extz2, Extz, Extz2Input, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF};
    if qg.is_empty() {
        return -(cfg.gap_open_pen as i32 + tg.len() as i32 * cfg.gap_extend_pen as i32);
    }
    if tg.is_empty() {
        return -(cfg.gap_open_pen as i32 + qg.len() as i32 * cfg.gap_extend_pen as i32);
    }
    let q = dna5_encode(qg);
    let t = dna5_encode(tg);
    let mat = dna5_mat(cfg);
    let w = cfg
        .bandwidth
        .max((qg.len() as i32 - tg.len() as i32).abs() + 4);
    let input = Extz2Input {
        query: &q,
        target: &t,
        m: 5,
        mat: &mat,
        q: cfg.gap_open_pen,
        e: cfg.gap_extend_pen,
        w,
        zdrop: -1,
        end_bonus: 0,
        // Mapping only needs the score, not the CIGAR — score-only mode skips the
        // O(qlen·tlen) traceback-matrix fill and backtrack pass.
        flag: KSW_EZ_RIGHT | KSW_EZ_SCORE_ONLY,
    };
    let mut ez = Extz::default();
    ez.reset();
    extz2(&input, &mut ez);
    if ez.score > KSW_NEG_INF {
        ez.score
    } else {
        // Corner unreachable: approximate by matching the common prefix length.
        let common = qg.len().min(tg.len()) as i32;
        common * cfg.match_score as i32 - (common * cfg.mismatch_pen as i32)
    }
}

/// Extension score of a read flank, anchored at the MEM-adjacent end: the read
/// flank must be fully consumed; the reference may extend freely outward. For a
/// 5' flank (`anchor_right`) the sequences are reversed so the anchor is on the
/// left, matching ksw2's left-anchored extension.
fn ksw2_flank_extend(qf: &[u8], tf: &[u8], cfg: &AlignConfig, anchor_right: bool) -> i32 {
    use ksw2rs::{extz2, Extz, Extz2Input, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF};
    if qf.is_empty() {
        return 0;
    }
    if tf.is_empty() {
        return -(cfg.gap_open_pen as i32 + qf.len() as i32 * cfg.gap_extend_pen as i32);
    }
    let (q, t): (Vec<u8>, Vec<u8>) = if anchor_right {
        (
            dna5_encode(&qf.iter().rev().copied().collect::<Vec<u8>>()),
            dna5_encode(&tf.iter().rev().copied().collect::<Vec<u8>>()),
        )
    } else {
        (dna5_encode(qf), dna5_encode(tf))
    };
    let mat = dna5_mat(cfg);
    let input = Extz2Input {
        query: &q,
        target: &t,
        m: 5,
        mat: &mat,
        q: cfg.gap_open_pen,
        e: cfg.gap_extend_pen,
        w: cfg.bandwidth.max(qf.len() as i32),
        zdrop: -1,
        end_bonus: 0,
        // Score-only: we read `mqe`/`max`, never the CIGAR (see ksw2_gap_global).
        flag: KSW_EZ_RIGHT | KSW_EZ_SCORE_ONLY,
    };
    let mut ez = Extz::default();
    ez.reset();
    extz2(&input, &mut ez);
    // read flank fully consumed, reference free to overhang -> mqe
    if ez.mqe > KSW_NEG_INF {
        ez.mqe
    } else {
        ez.max as i32
    }
}

fn cached_gap_score(qg: &[u8], tg: &[u8], cfg: &AlignConfig) -> i32 {
    if qg.is_empty() && tg.is_empty() {
        return 0;
    }
    GAP_CACHE.with(|c| {
        let key = (
            qg.to_vec().into_boxed_slice(),
            tg.to_vec().into_boxed_slice(),
        );
        if let Some(&s) = c.borrow().get(&key) {
            return s;
        }
        let s = ksw2_gap_global(qg, tg, cfg);
        gap_cache_insert(c, key, s);
        s
    })
}

fn cached_flank_score(qf: &[u8], tf: &[u8], cfg: &AlignConfig, anchor_right: bool) -> i32 {
    if qf.is_empty() {
        return 0;
    }
    GAP_CACHE.with(|c| {
        // distinguish flank orientation in the key
        let mut qk = qf.to_vec();
        qk.push(if anchor_right { 1 } else { 0 });
        let key = (qk.into_boxed_slice(), tf.to_vec().into_boxed_slice());
        if let Some(&s) = c.borrow().get(&key) {
            return s;
        }
        let s = ksw2_flank_extend(qf, tf, cfg, anchor_right);
        gap_cache_insert(c, key, s);
        s
    })
}

/// 2-bit/DNA5 encode (A=0,C=1,G=2,T=3,N/other=4) for ksw2.
fn dna5_encode(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .map(|&b| match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => 4,
        })
        .collect()
}

/// ksw2 score, matching C++ salmon's banded `ksw_extz2_sse` configuration
/// (`KSW_EZ_RIGHT | KSW_EZ_SCORE_ONLY`, gapo/gape, bandwidth, z-drop off). Uses
/// `mqe` — the max score reaching the query end — i.e. the read aligned in full
/// with the reference free to overhang, within the band around the chain
/// diagonal. This narrow band is what makes salmon reject off-diagonal /
/// large-indel placements that a wide-band aligner would accept.
fn ksw2_align_score(query: &[u8], rwin: &[u8], cfg: &AlignConfig) -> i32 {
    use ksw2rs::{extz2, Extz, Extz2Input, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF};
    let q_enc = dna5_encode(query);
    let t_enc = dna5_encode(rwin);
    // 5x5 DNA5 scoring matrix: match on the diagonal, mismatch off, N-N = 0.
    let mut mat = [-cfg.mismatch_pen; 25];
    for i in 0..5 {
        mat[i * 5 + i] = cfg.match_score;
    }
    mat[24] = 0;
    let input = Extz2Input {
        query: &q_enc,
        target: &t_enc,
        m: 5,
        mat: &mat,
        q: cfg.gap_open_pen,
        e: cfg.gap_extend_pen,
        w: cfg.bandwidth,
        zdrop: -1,
        end_bonus: 0,
        flag: KSW_EZ_RIGHT | KSW_EZ_SCORE_ONLY,
    };
    let mut ez = Extz::default();
    ez.reset(); // initialize mqe/mte/score to KSW_NEG_INF as ksw2 expects
    extz2(&input, &mut ez);
    // `mqe` is the best score with the entire query consumed (the read aligned
    // fully). Fall back to the local max if the query end was never reached
    // within the band.
    if ez.mqe > KSW_NEG_INF {
        ez.mqe
    } else {
        ez.max as i32
    }
}

/// Result of placing a read anywhere within a reference window.
#[derive(Debug, Clone, Copy)]
pub struct WindowAlignment {
    pub score: i32,
    pub valid: bool,
    /// reference column (within the window) where the read alignment ends
    pub end_col: usize,
}

/// Align `query` as a substring of `window` (free reference gaps at both ends),
/// used for orphan recovery where the mate's position within a fragment-length
/// window is unknown. The query must align in full; the reference may overhang
/// on either side for free, and the best reference end column is reported.
///
/// This is a "fitting" affine-gap alignment (global on the query, free ends on
/// the reference). It's a plain scalar DP — orphan recovery is opt-in
/// (`recover_orphans`, off by default) and infrequent, so it is not on the hot
/// path; keeping it scalar avoids any SIMD/platform dependency here. Returns
/// `None` if either sequence is empty.
pub fn align_in_window(query: &[u8], window: &[u8], cfg: &AlignConfig) -> Option<WindowAlignment> {
    if query.is_empty() || window.is_empty() {
        return None;
    }
    let n = query.len();
    let m = window.len();
    let ma = cfg.match_score as i32;
    let mp = cfg.mismatch_pen as i32;
    let go = cfg.gap_open_pen as i32;
    let ge = cfg.gap_extend_pen as i32;
    const NEG: i32 = i32::MIN / 4;

    // Rolling DP rows over reference columns 0..=m.
    //  h[j]  = best score with query[0..i] consumed, alignment ending at ref col j
    //  ix[j] = best such score ending in a gap that consumed a query base (no ref)
    // Free reference start: row i=0 is 0 for every column, so the read may begin
    // opposite any reference position at no cost.
    let mut h_prev = vec![0i32; m + 1];
    let mut ix_prev = vec![NEG; m + 1];
    let mut h_cur = vec![0i32; m + 1];
    let mut ix_cur = vec![NEG; m + 1];

    for i in 1..=n {
        // Column 0: the read cannot be placed entirely left of the window.
        h_cur[0] = NEG;
        ix_cur[0] = NEG;
        let mut iy = NEG; // gap consuming a reference base (deletion in the read)
        let qi = query[i - 1];
        for j in 1..=m {
            let s = if qi == window[j - 1] { ma } else { -mp };
            let diag = h_prev[j - 1].saturating_add(s);
            ix_cur[j] = (h_prev[j].saturating_sub(go))
                .max(ix_prev[j].saturating_sub(ge))
                .max(NEG);
            iy = (h_cur[j - 1].saturating_sub(go))
                .max(iy.saturating_sub(ge))
                .max(NEG);
            h_cur[j] = diag.max(ix_cur[j]).max(iy).max(NEG);
        }
        std::mem::swap(&mut h_prev, &mut h_cur);
        std::mem::swap(&mut ix_prev, &mut ix_cur);
    }

    // Query fully consumed; reference free to end anywhere → best over all columns.
    let mut best = NEG;
    let mut end_col = 0usize;
    for j in 1..=m {
        if h_prev[j] > best {
            best = h_prev[j];
            end_col = j;
        }
    }
    let perfect = perfect_score(n, cfg);
    Some(WindowAlignment {
        score: best,
        valid: (best as f32) >= cfg.min_score_fraction * (perfect as f32),
        end_col,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mem::Mem;

    fn gen_seq(n: usize, seed: u64) -> Vec<u8> {
        const B: [u8; 4] = *b"ACGT";
        let mut x = seed;
        let mut s = Vec::with_capacity(n);
        for _ in 0..n {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            s.push(B[((x >> 33) & 3) as usize]);
        }
        s
    }

    /// A forward chain whose single MEM places read offset 0 at ref `r0`.
    fn fw_chain(r0: i32, len: i32) -> MemChain {
        MemChain::new(vec![Mem::new(0, r0, len)], len as f32, true)
    }

    #[test]
    fn exact_match_scores_perfect() {
        let reference = gen_seq(300, 7);
        let read = &reference[50..130]; // 80 bp exact
        let chain = fw_chain(50, 31);
        let cfg = AlignConfig::default();
        let aln = align_chain(read, &reference, &chain, &cfg).unwrap();
        assert_eq!(aln.score, perfect_score(read.len(), &cfg));
        assert!(aln.valid);
        assert_eq!(aln.ref_window_start, 50);
    }

    #[test]
    fn anchored_and_full_agree_on_clean_read() {
        let reference = gen_seq(300, 71);
        let read = &reference[50..130];
        let chain = fw_chain(50, 31);
        let anchored = AlignConfig::default(); // full_length_alignment = false
        let full = AlignConfig {
            full_length_alignment: true,
            ..AlignConfig::default()
        };
        let a = align_chain(read, &reference, &chain, &anchored).unwrap();
        let f = align_chain(read, &reference, &chain, &full).unwrap();
        assert_eq!(a.score, perfect_score(read.len(), &anchored));
        assert_eq!(a.score, f.score);
    }

    #[test]
    fn anchored_scores_two_mem_chain_with_exact_gap() {
        // read exactly matches ref[100..160]; two MEMs with an exact gap between.
        let reference = gen_seq(400, 73);
        let read = reference[100..160].to_vec();
        let chain = MemChain::new(vec![Mem::new(0, 100, 25), Mem::new(30, 130, 30)], 55.0, true);
        let cfg = AlignConfig::default();
        let a = align_chain(&read, &reference, &chain, &cfg).unwrap();
        assert_eq!(a.score, perfect_score(read.len(), &cfg), "got {}", a.score);
        assert!(a.valid);
    }

    #[test]
    fn single_mismatch_reduces_score_but_stays_valid() {
        let reference = gen_seq(300, 11);
        let mut read = reference[50..130].to_vec();
        // flip one base in the middle to a different base
        let mid = read.len() / 2;
        read[mid] = if read[mid] == b'A' { b'C' } else { b'A' };
        let cfg = AlignConfig::default();
        let aln = align_chain(&read, &reference, &fw_chain(50, 31), &cfg).unwrap();
        let perfect = perfect_score(read.len(), &cfg);
        assert!(aln.score < perfect, "score {} should be < {perfect}", aln.score);
        // one mismatch: (len-1)*match + (-mismatch)
        let expected = (read.len() as i32 - 1) * cfg.match_score as i32 - cfg.mismatch_pen as i32;
        assert_eq!(aln.score, expected);
        assert!(aln.valid);
    }

    #[test]
    fn reverse_complement_read_aligns_in_fw_frame() {
        let reference = gen_seq(300, 13);
        let window = &reference[50..130];
        let rc_read = revcomp(window); // the read is the RC of the ref window
        let chain = MemChain::new(vec![Mem::new(0, 50, 31)], 31.0, false);
        let cfg = AlignConfig::default();
        let aln = align_chain(&rc_read, &reference, &chain, &cfg).unwrap();
        // wrapper revcomps the read back to the forward window -> perfect match
        assert_eq!(aln.score, perfect_score(rc_read.len(), &cfg));
        assert!(aln.valid);
    }

    #[test]
    fn unrelated_sequence_is_invalid() {
        let reference = gen_seq(300, 17);
        let foreign = gen_seq(80, 999); // unrelated
        let cfg = AlignConfig::default();
        let aln = align_chain(&foreign, &reference, &fw_chain(50, 31), &cfg).unwrap();
        assert!(!aln.valid, "unrelated read should not validate (score {})", aln.score);
    }

    #[test]
    fn fitting_window_places_read_in_middle() {
        // Read exactly matches window[40..120]; free ref-start must find it and
        // report the read's end column (120) with a perfect score.
        let window = gen_seq(200, 91);
        let read = window[40..120].to_vec();
        let cfg = AlignConfig::default();
        let a = align_in_window(&read, &window, &cfg).unwrap();
        assert_eq!(a.score, perfect_score(read.len(), &cfg), "score {}", a.score);
        assert!(a.valid);
        assert_eq!(a.end_col, 120, "end_col {}", a.end_col);
    }

    #[test]
    fn fitting_window_one_mismatch() {
        let window = gen_seq(200, 92);
        let mut read = window[40..120].to_vec();
        let mid = read.len() / 2;
        read[mid] = if read[mid] == b'A' { b'C' } else { b'A' };
        let cfg = AlignConfig::default();
        let a = align_in_window(&read, &window, &cfg).unwrap();
        let expected = (read.len() as i32 - 1) * cfg.match_score as i32 - cfg.mismatch_pen as i32;
        assert_eq!(a.score, expected, "score {}", a.score);
        assert!(a.valid);
        assert_eq!(a.end_col, 120);
    }

    #[test]
    fn fitting_window_rejects_unrelated() {
        let window = gen_seq(200, 93);
        let foreign = gen_seq(80, 54321);
        let cfg = AlignConfig::default();
        let a = align_in_window(&foreign, &window, &cfg).unwrap();
        assert!(!a.valid, "unrelated read should not validate (score {})", a.score);
    }

    #[test]
    fn small_deletion_in_read_handled_by_band() {
        // Read equals ref window with a 2bp deletion in the middle.
        let reference = gen_seq(400, 23);
        let win = &reference[100..200];
        let mut read = Vec::new();
        read.extend_from_slice(&win[..50]);
        read.extend_from_slice(&win[52..]); // drop 2 bases
        let cfg = AlignConfig::default();
        let aln = align_chain(&read, &reference, &fw_chain(100, 31), &cfg).unwrap();
        // a 2bp gap: most bases still match, should remain valid
        assert!(aln.valid, "score {} should validate a 2bp deletion", aln.score);
        assert!(aln.score < perfect_score(read.len(), &cfg));
    }
}
