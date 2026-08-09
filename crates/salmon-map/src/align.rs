//! Alignment validation of mapping candidates.
//!
//! # Why candidates are verified
//!
//! Chaining says a read *could* align here; it does not say how well. Transcripts
//! sharing a repeat, a domain or a low-complexity stretch produce plausible
//! chains for reads that did not come from them. Alignment settles it by lining
//! the read's bases up against the reference and scoring the result: matches earn
//! points, mismatches and gaps cost them.
//!
//! A candidate is kept only if it scores close enough to a perfect alignment,
//! which is what makes salmon's mapping *selective* — sensitive enough to keep
//! genuine multi-mappings, strict enough to reject coincidental ones.
//!
//! # Affine gaps, briefly
//!
//! A gap of length `n` costs `gapo + gape * n`, not `n` times a flat penalty.
//! The one-off opening cost makes a single long gap much cheaper than several
//! short ones, which matches how insertions and deletions actually occur.
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
    /// allow soft-clipping read ends in flank extension (salmon's `--softclip`):
    /// a read flank takes `max(mqe, mte, 0)` instead of requiring full
    /// consumption (`mqe`). Implies [`Self::softclip_overhangs`].
    pub softclip: bool,
    /// allow soft-clipping only read bases that overhang a transcript end
    /// (salmon's `--softclipOverhangs`): a read flank takes `max(mqe, mte)`.
    pub softclip_overhangs: bool,
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
            softclip: false,
            softclip_overhangs: false,
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

/// Reverse-complement a DNA base (ACGT, case-insensitive; anything else → N).
#[inline]
fn comp(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

/// Reverse-complement a DNA byte slice (ACGT, case-insensitive output upper).
pub(crate) fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}

/// Reverse-complement `s` into `out` (cleared first), reusing its allocation.
#[inline]
fn revcomp_into(s: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(s.iter().rev().map(|&b| comp(b)));
}

/// The perfect (all-match) score for a read of length `read_len`.
#[inline]
pub fn perfect_score(read_len: usize, cfg: &AlignConfig) -> i32 {
    read_len as i32 * cfg.match_score as i32
}

/// Minimum alignment score for a read of `read_len` to be accepted, matching
/// salmon/pufferfish's `minAcceptedScore = floor(minScoreFraction · matchScore ·
/// readLen)` with acceptance `score > minAcceptedScore` (strict). Using `floor`
/// with strict `>` (rather than `score (f32) >= frac·perfect`) avoids an
/// off-by-one at read lengths where `frac·matchScore·len` is an integer.
#[inline]
pub fn min_accepted_score(read_len: usize, cfg: &AlignConfig) -> i32 {
    (cfg.min_score_fraction as f64 * cfg.match_score as f64 * read_len as f64).floor() as i32
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

    let diag_origin = chain.ref_start() - chain.read_start();
    let win_start = diag_origin.max(0);
    let reflen_i = ref_seq.len() as i32;

    // Build the oriented query into a reused thread-local buffer — `align_chain`
    // runs once per candidate (millions of times), so a fresh
    // `read.to_vec()`/`revcomp` per call was pure allocator traffic.
    ALIGN_QUERY.with(|qbuf| {
        let mut query = qbuf.borrow_mut();
        if chain.is_fw {
            query.clear();
            query.extend_from_slice(read);
        } else {
            revcomp_into(read, &mut query);
        }
        let qlen = query.len() as i32;

        // Invariant boundary: every MEM must lie within the read and the reference.
        // Seeding guarantees this; a violation is an upstream bug. We `debug_assert`
        // so it surfaces loudly in tests/dev (rather than being silently masked —
        // this is how issue #1038, a contig-overshoot in piscem-rs's skip search,
        // was found), and degrade to "skip candidate" in release so a future
        // coordinate bug can never turn into an out-of-bounds slice panic.
        let out_of_bounds = chain.mems.iter().any(|m| {
            m.read_start < 0
                || m.ref_start < 0
                || m.len <= 0
                || m.read_end() > qlen
                || m.ref_end() > reflen_i
        });
        debug_assert!(
            !out_of_bounds,
            "align_chain: chain MEM outside read/reference bounds (qlen={qlen} reflen={reflen_i}) — upstream seeding bug"
        );
        if out_of_bounds {
            return None;
        }

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

        let valid = score > min_accepted_score(query.len(), cfg);

        Some(Alignment {
            score,
            valid,
            ref_window_start: win_start,
        })
    })
}

thread_local! {
    /// Per-thread cache of gap/flank DP scores, keyed by the exact (kind, query,
    /// ref) triple — salmon's alignment cache, to avoid re-running identical DPs
    /// (common for shared gaps/flanks across the candidate placements of one
    /// read). salmon scopes this cache per read; we instead bound it by size
    /// (see [`GAP_CACHE_CAP`]) so it cannot grow without limit over a whole run
    /// — an unbounded thread-local here cost tens of GB on a 36M-read library.
    static GAP_CACHE: std::cell::RefCell<GapCache> =
        std::cell::RefCell::new(GapCache::default());
    /// Per-thread ksw2 scratch. This reuses both the ksw DP workspace/result and
    /// the DNA5-encoded query/target buffers between small gap/flank DPs.
    static KSW_SCRATCH: std::cell::RefCell<KswScratch> =
        std::cell::RefCell::new(KswScratch::default());
    /// Per-thread query buffer for [`align_chain`], reused across the millions of
    /// per-candidate alignments instead of allocating a fresh `Vec` each call.
    static ALIGN_QUERY: std::cell::RefCell<Vec<u8>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

struct KswScratch {
    aligner: ksw2rs::Aligner,
    query: Vec<u8>,
    target: Vec<u8>,
}

impl Default for KswScratch {
    fn default() -> Self {
        Self {
            aligner: ksw2rs::Aligner::new(),
            query: Vec::new(),
            target: Vec::new(),
        }
    }
}

/// Maximum number of cached (kind, query, ref) DP scores per thread. The cache
/// dedups identical inter-MEM gap / flank alignments, both across one read's
/// candidate targets and across reads. When the cap is hit we clear the whole
/// cache, keeping peak memory at a few MB per thread instead of unbounded growth.
///
/// Measured hit rates on 10M paired reads against a GRCh38 cDNA index:
///
/// | cache | probes      | hits       | hit rate |
/// |-------|-------------|------------|----------|
/// | gap   |   5,044,438 |  4,550,348 |    90.2% |
/// | flank | 113,786,498 | 57,614,418 |    50.6% |
///
/// An earlier version of this comment claimed "cross-read reuse is rare (flanks
/// are read-specific)" and on that basis the flank cache looked close to
/// useless. It is not: it avoids ~57.6M DP calls per run, and flank probes
/// outnumber gap probes 22:1, so the flank path is where both the saved work and
/// the allocator traffic actually are.
///
/// Note these caches are on the default (anchored) scoring path only.
/// `--fullLengthAlignment` takes the whole-read `ksw2_align_score` branch and
/// never reaches them (measured: zero probes).
const GAP_CACHE_CAP: usize = 32_768;

/// Which DP a cache entry holds. Gap and flank scores come from *different*
/// functions (`ksw2_gap_global` vs `ksw2_flank_extend`) over the same kind of
/// byte pair, so the discriminant is part of the key. It used to be smuggled in
/// as a trailing `0`/`1` byte appended to the query, which worked only because
/// no ACGT base equals `0x00` or `0x01` — an unasserted assumption about read
/// content guarding a silently-wrong-score failure mode.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
enum GapKind {
    /// inter-MEM gap, both substrings fully consumed
    Gap = 0,
    /// read flank anchored on its left (3' flank)
    FlankLeftAnchored = 1,
    /// read flank anchored on its right (5' flank)
    FlankRightAnchored = 2,
}

/// One cached DP score. `seqs` holds the query bytes followed by the target
/// bytes in a single allocation (rather than two `Box<[u8]>`), so a miss costs
/// one allocation instead of two.
struct GapEntry {
    kind: GapKind,
    seqs: Box<[u8]>,
    qlen: u32,
    score: i32,
}

/// Per-thread exact memo table for gap/flank DP scores.
///
/// Uses `hashbrown::HashTable` rather than a `HashMap` so a lookup can probe
/// with the *borrowed* `(kind, query, target)` triple. A `HashMap` keyed by owned
/// bytes forces the caller to materialize the key before it can ask whether the
/// entry exists, which meant two heap allocations and two frees on **every**
/// probe, cache hits included — pure overhead on the hottest inner loop in the
/// mapper. Only a miss allocates now, and only once.
///
/// The comparison is still byte-exact, so a hit always returns the score that
/// the DP would have produced. (A digest-only key would have been cheaper still,
/// but it trades exactness for a collision probability, and a wrong alignment
/// score is not a trade worth making here.)
struct GapCache {
    table: hashbrown::HashTable<GapEntry>,
    hasher: ahash::RandomState,
}

impl Default for GapCache {
    fn default() -> Self {
        Self {
            table: hashbrown::HashTable::new(),
            // Fixed seeds: the lookup is byte-exact so the seed cannot affect
            // results, but pinning it keeps bucket layout (and therefore
            // profiling runs) reproducible.
            hasher: ahash::RandomState::with_seeds(
                0x51_7c_c1_b7,
                0x27_22_0a_95,
                0x9e_37_79_b9,
                0x85_eb_ca_6b,
            ),
        }
    }
}

/// Hash a borrowed cache key. Must stay in lockstep with [`GapCache::eq`].
///
/// A free function rather than a method so the rehash-on-growth closure can call
/// it while the table itself is mutably borrowed.
#[inline]
fn gap_hash(hasher: &ahash::RandomState, kind: GapKind, q: &[u8], t: &[u8]) -> u64 {
    use std::hash::{BuildHasher, Hash, Hasher};
    let mut h = hasher.build_hasher();
    (kind as u8).hash(&mut h);
    q.hash(&mut h);
    t.hash(&mut h);
    h.finish()
}

impl GapCache {
    /// Exact equality between a stored entry and a borrowed key.
    #[inline]
    fn eq(e: &GapEntry, kind: GapKind, q: &[u8], t: &[u8]) -> bool {
        e.kind == kind
            && e.qlen as usize == q.len()
            && e.seqs.len() == q.len() + t.len()
            && &e.seqs[..q.len()] == q
            && &e.seqs[q.len()..] == t
    }

    /// Cached score for this key, if present. Allocates nothing.
    #[inline]
    fn get(&self, kind: GapKind, q: &[u8], t: &[u8]) -> Option<i32> {
        let h = gap_hash(&self.hasher, kind, q, t);
        self.table
            .find(h, |e| Self::eq(e, kind, q, t))
            .map(|e| e.score)
    }

    /// Record a score, clearing the whole table first if it has reached the cap.
    fn insert(&mut self, kind: GapKind, q: &[u8], t: &[u8], score: i32) {
        let GapCache { table, hasher } = self;
        if table.len() >= GAP_CACHE_CAP {
            table.clear();
        }
        let h = gap_hash(hasher, kind, q, t);
        let mut seqs = Vec::with_capacity(q.len() + t.len());
        seqs.extend_from_slice(q);
        seqs.extend_from_slice(t);
        let entry = GapEntry {
            kind,
            seqs: seqs.into_boxed_slice(),
            qlen: q.len() as u32,
            score,
        };
        // The hash is a pure function of the entry's own key, so rehashing on
        // growth just recomputes it.
        table.insert_unique(h, entry, |e| {
            let (q, t) = e.seqs.split_at(e.qlen as usize);
            gap_hash(hasher, e.kind, q, t)
        });
    }
}

/// PuffAligner-style score: each MEM contributes an exact-match score and only
/// the inter-MEM gaps and the read flanks are DP-aligned (with ksw2), reusing
/// cached results for repeated substrings. Treating MEMs as fixed exact matches
/// is what makes this cheaper than a full-read DP; it can be marginally lower
/// than the global optimum when re-aligning across a MEM boundary would help.
fn anchored_align_score(
    query: &[u8],
    ref_seq: &[u8],
    mems: &[crate::mem::Mem],
    cfg: &AlignConfig,
) -> i32 {
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
            // Per-axis overlap of this MEM with the next (negative ⇒ a true gap).
            let read_ov = mem.read_end() - nxt.read_start;
            let ref_ov = mem.ref_end() - nxt.ref_start;
            // Trim the overlapping exact bases back to a common boundary on BOTH
            // axes, then DP-align whatever residual remains between the MEMs:
            //  * same diagonal (read_ov == ref_ov): residual is empty on both
            //    axes, so this degenerates to "remove the double-counted bases"
            //    (the cached gap DP returns 0 for empty/empty);
            //  * gap-separated (both overlaps <= 0): ov is 0 and we DP the gap;
            //  * different diagonals (read_ov != ref_ov): trimming to the larger
            //    overlap exposes the implied indel as a one-sided residual, which
            //    the gap DP penalizes — instead of the indel being silently
            //    absorbed as a plain overlap (which inflated the score).
            let ov = read_ov.max(ref_ov).max(0);
            score -= ov * m;
            let q_lo = (mem.read_end() - ov)
                .max(mem.read_start)
                .min(nxt.read_start) as usize;
            let t_lo = (mem.ref_end() - ov).max(mem.ref_start).min(nxt.ref_start) as usize;
            let qg = &query[q_lo..nxt.read_start as usize];
            let tg = &ref_seq[t_lo..nxt.ref_start as usize];
            score += cached_gap_score(qg, tg, cfg);
        }
    }

    // 3' flank: read after the last MEM, anchored on its left at the MEM.
    let last = mems.last().unwrap();
    if last.read_end() < qlen {
        let qf = &query[last.read_end() as usize..qlen as usize];
        let r_hi =
            (last.ref_end() + (qlen - last.read_end()) + cfg.indel_margin as i32).min(reflen);
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
    use ksw2rs::{Extz2Input, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF};
    if qg.is_empty() {
        return -(cfg.gap_open_pen as i32 + tg.len() as i32 * cfg.gap_extend_pen as i32);
    }
    if tg.is_empty() {
        return -(cfg.gap_open_pen as i32 + qg.len() as i32 * cfg.gap_extend_pen as i32);
    }
    let mat = dna5_mat(cfg);
    let w = cfg
        .bandwidth
        .max((qg.len() as i32 - tg.len() as i32).abs() + 4);
    let score = KSW_SCRATCH.with(|cell| {
        let KswScratch {
            aligner,
            query,
            target,
        } = &mut *cell.borrow_mut();
        dna5_encode_into(qg, query);
        dna5_encode_into(tg, target);
        let input = Extz2Input {
            query,
            target,
            m: 5,
            mat: &mat,
            q: cfg.gap_open_pen,
            e: cfg.gap_extend_pen,
            w,
            zdrop: -1,
            end_bonus: 0,
            // Mapping only needs the score, not the CIGAR — score-only mode skips
            // the O(qlen·tlen) traceback-matrix fill and backtrack pass.
            flag: KSW_EZ_RIGHT | KSW_EZ_SCORE_ONLY,
        };
        aligner.align(&input).score
    });
    if score > KSW_NEG_INF {
        score
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
    use ksw2rs::{Extz2Input, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF};
    if qf.is_empty() {
        return 0;
    }
    if tf.is_empty() {
        return -(cfg.gap_open_pen as i32 + qf.len() as i32 * cfg.gap_extend_pen as i32);
    }
    let mat = dna5_mat(cfg);
    let (mqe, mte, max_score) = KSW_SCRATCH.with(|cell| {
        let KswScratch {
            aligner,
            query,
            target,
        } = &mut *cell.borrow_mut();
        if anchor_right {
            dna5_encode_rev_into(qf, query);
            dna5_encode_rev_into(tf, target);
        } else {
            dna5_encode_into(qf, query);
            dna5_encode_into(tf, target);
        }
        let input = Extz2Input {
            query,
            target,
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
        let ez = aligner.align(&input);
        (ez.mqe, ez.mte, ez.max as i32)
    });
    // Soft-clip semantics mirror PuffAligner's flank `part_score` (mqe = read
    // flank fully consumed, reference free to overhang; mte = reference fully
    // consumed, read free to overhang the transcript end):
    //   none      -> mqe              (full read flank must align)
    //   overhangs  -> max(mqe, mte)    (clip only bases past a transcript end)
    //   full       -> max(mqe, mte, 0) (clip read ends freely; never negative)
    // `--softclip` implies `--softclipOverhangs`. `ez.max` is the fallback when
    // neither extension is valid (e.g. band-clipped), matching the prior guard.
    if cfg.softclip {
        let s = mqe.max(mte);
        if s > KSW_NEG_INF {
            s.max(0)
        } else {
            max_score.max(0)
        }
    } else if cfg.softclip_overhangs {
        let s = mqe.max(mte);
        if s > KSW_NEG_INF {
            s
        } else {
            max_score
        }
    } else if mqe > KSW_NEG_INF {
        mqe
    } else {
        max_score
    }
}

fn cached_gap_score(qg: &[u8], tg: &[u8], cfg: &AlignConfig) -> i32 {
    if qg.is_empty() && tg.is_empty() {
        return 0;
    }
    GAP_CACHE.with(|c| {
        // Probe with the borrowed slices; nothing is allocated unless we miss.
        if let Some(s) = c.borrow().get(GapKind::Gap, qg, tg) {
            return s;
        }
        let s = ksw2_gap_global(qg, tg, cfg);
        c.borrow_mut().insert(GapKind::Gap, qg, tg, s);
        s
    })
}

fn cached_flank_score(qf: &[u8], tf: &[u8], cfg: &AlignConfig, anchor_right: bool) -> i32 {
    if qf.is_empty() {
        return 0;
    }
    // The anchor side selects a different DP, so it is part of the key.
    let kind = if anchor_right {
        GapKind::FlankRightAnchored
    } else {
        GapKind::FlankLeftAnchored
    };
    GAP_CACHE.with(|c| {
        if let Some(s) = c.borrow().get(kind, qf, tf) {
            return s;
        }
        let s = ksw2_flank_extend(qf, tf, cfg, anchor_right);
        c.borrow_mut().insert(kind, qf, tf, s);
        s
    })
}

/// DNA5 lookup table: `DNA5_LUT[b]` is the 2-bit code for ACGT (upper/lower) and
/// 4 for anything else. A single table load per base beats the per-base `match`
/// (a chain of compares) on the alignment hot path, where encoding the read and
/// reference window is a measurable share of full-length scoring.
static DNA5_LUT: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

fn dna5_encode_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(seq.len());
    out.extend(seq.iter().map(|&b| DNA5_LUT[b as usize]));
}

/// Like [`dna5_encode_into`] but emits the codes in reverse order. Note this
/// only reverses — it does *not* complement — so the per-base LUT applies
/// unchanged (the caller has already oriented the sequence).
fn dna5_encode_rev_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(seq.len());
    out.extend(seq.iter().rev().map(|&b| DNA5_LUT[b as usize]));
}

/// ksw2 score, matching C++ salmon's banded `ksw_extz2_sse` configuration
/// (`KSW_EZ_RIGHT | KSW_EZ_SCORE_ONLY`, gapo/gape, bandwidth, z-drop off). Uses
/// `mqe` — the max score reaching the query end — i.e. the read aligned in full
/// with the reference free to overhang, within the band around the chain
/// diagonal. This narrow band is what makes salmon reject off-diagonal /
/// large-indel placements that a wide-band aligner would accept.
fn ksw2_align_score(query: &[u8], rwin: &[u8], cfg: &AlignConfig) -> i32 {
    use ksw2rs::{Extz2Input, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF};
    // 5x5 DNA5 scoring matrix: match on the diagonal, mismatch off, N-N = 0.
    let mut mat = [-cfg.mismatch_pen; 25];
    for i in 0..5 {
        mat[i * 5 + i] = cfg.match_score;
    }
    mat[24] = 0;
    let (mqe, max_score) = KSW_SCRATCH.with(|cell| {
        let KswScratch {
            aligner,
            query: q_enc,
            target: t_enc,
        } = &mut *cell.borrow_mut();
        dna5_encode_into(query, q_enc);
        dna5_encode_into(rwin, t_enc);
        let input = Extz2Input {
            query: q_enc,
            target: t_enc,
            m: 5,
            mat: &mat,
            q: cfg.gap_open_pen,
            e: cfg.gap_extend_pen,
            w: cfg.bandwidth,
            zdrop: -1,
            end_bonus: 0,
            flag: KSW_EZ_RIGHT | KSW_EZ_SCORE_ONLY,
        };
        let ez = aligner.align(&input);
        (ez.mqe, ez.max as i32)
    });
    // `mqe` is the best score with the entire query consumed (the read aligned
    // fully). Fall back to the local max if the query end was never reached
    // within the band.
    if mqe > KSW_NEG_INF {
        mqe
    } else {
        max_score
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
    Some(WindowAlignment {
        score: best,
        valid: best > min_accepted_score(n, cfg),
        end_col,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mem::Mem;

    #[test]
    /// A read hanging off a transcript end must be soft-clipped rather than scored
    /// as mismatches, or a legitimate placement at a boundary would be rejected.
    fn flank_softclip_recovers_read_overhang() {
        // Read flank (10 b) matches the reference flank (5 b) then overhangs the
        // end. Without soft-clipping the full read must align (mqe); with overhang
        // / full soft-clipping the 5 matched bases score (mte) and the overhang is
        // clipped, so the score is no worse and reflects the matched prefix.
        let qf = b"ACGTACGTAC";
        let tf = b"ACGTA";
        let none = AlignConfig::default();
        let overhang = AlignConfig {
            softclip_overhangs: true,
            ..AlignConfig::default()
        };
        let full = AlignConfig {
            softclip: true,
            ..AlignConfig::default()
        };
        let s_none = ksw2_flank_extend(qf, tf, &none, false);
        let s_oh = ksw2_flank_extend(qf, tf, &overhang, false);
        let s_full = ksw2_flank_extend(qf, tf, &full, false);
        // overhang recovers the 5 matched bases (mte = 5 * match_score = 10)
        assert_eq!(s_oh, 10, "overhang flank should score the matched prefix");
        assert!(s_oh >= s_none, "overhang {s_oh} >= none {s_none}");
        assert!(s_full >= 0, "full soft-clip never negative: {s_full}");
        assert!(s_full >= s_oh, "full {s_full} >= overhang {s_oh}");
    }

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
    /// A read identical to its reference window must reach exactly the perfect
    /// score, which is the denominator every acceptance threshold is relative to.
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
    /// The cheap anchored alignment and the full one must agree on an easy read;
    /// the anchored path is an optimization, not a different answer.
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
    /// A deletion must cost something, or a stitched-together alignment would
    /// score as well as a contiguous one.
    fn read_overlap_with_reference_gap_penalizes_deletion() {
        // Divergent-paralog geometry (read ERR/sim mate vs a partial paralog):
        // two exact MEMs that overlap by 1 base in the read but are separated by a
        // 60 bp REFERENCE gap (different diagonals) — a large deletion. The pre-fix
        // overlap branch absorbed this as a 1-base overlap (a perfect 152) and
        // over-credited the placement past the score threshold; the residual gap DP
        // must charge the deletion so the divergent placement is rejected.
        let cfg = AlignConfig::default();
        let query = gen_seq(76, 3);
        let reference = gen_seq(400, 9);
        let mems = vec![Mem::new(0, 198, 43), Mem::new(42, 301, 34)];
        let score = anchored_align_score(&query, &reference, &mems, &cfg);
        let m = cfg.match_score as i32;
        // 43*2 + 34*2 - 1bp read overlap, minus a 61 bp deletion gap (open + 61*ext).
        let gap = cfg.gap_open_pen as i32 + 61 * cfg.gap_extend_pen as i32;
        let expected = (43 + 34 - 1) * m - gap; // 152 - 128 = 24
        assert_eq!(score, expected);
        assert!(
            score < min_accepted_score(query.len(), &cfg),
            "divergent placement must be rejected (score {score})"
        );
    }

    #[test]
    /// A chain whose anchors are separated by an exactly-matching stretch should
    /// score as one clean alignment.
    fn anchored_scores_two_mem_chain_with_exact_gap() {
        // read exactly matches ref[100..160]; two MEMs with an exact gap between.
        let reference = gen_seq(400, 73);
        let read = reference[100..160].to_vec();
        let chain = MemChain::new(
            vec![Mem::new(0, 100, 25), Mem::new(30, 130, 30)],
            55.0,
            true,
        );
        let cfg = AlignConfig::default();
        let a = align_chain(&read, &reference, &chain, &cfg).unwrap();
        assert_eq!(a.score, perfect_score(read.len(), &cfg), "got {}", a.score);
        assert!(a.valid);
    }

    #[test]
    /// Anchors on different diagonals imply an indel between them, which must be
    /// charged for.
    fn overlapping_mems_on_shifted_diagonals_incur_gap_penalty() {
        // Two MEMs that overlap by 1 base in the read but by 4 in the reference
        // (diagonals 3 apart) — a hidden 3 bp indel. This is the exact geometry
        // observed for ERR188044.600028 vs ENST00000677249.1. Pre-fix the
        // overlap branch scored it as a plain overlap (sum*match - max_overlap),
        // silently absorbing the indel and inflating the score; now the residual
        // is DP-aligned and the diagonal jump costs an affine gap.
        let cfg = AlignConfig::default();
        let query = gen_seq(76, 5);
        let reference = gen_seq(900, 9);
        // mem0: read[0..45]@ref790 (diag 745); mem1: read[44..76]@ref831 (diag 787).
        // read_ov = 45-44 = 1 ; ref_ov = (790+45)-831 = 4  -> implied 3 bp indel.
        let mems = vec![Mem::new(0, 790, 45), Mem::new(44, 831, 32)];
        let m = cfg.match_score as i32;
        let score = anchored_align_score(&query, &reference, &mems, &cfg);
        // Old (buggy) behavior absorbed the indel as a 4-base overlap only:
        let absorbed = (45 + 32) * m - 4 * m; // = 146
                                              // Fixed: 73 matched bases minus a 3 bp affine gap (open + 3*extend).
        let gap = cfg.gap_open_pen as i32 + 3 * cfg.gap_extend_pen as i32;
        let expected = (45 + 32 - 4) * m - gap; // 146 - 12 = 134
        assert_eq!(score, expected, "expected indel-penalized score");
        assert!(
            score < absorbed,
            "indel was not penalized: {score} >= {absorbed}"
        );
    }

    #[test]
    /// One sequencing error must lower the score without disqualifying the read —
    /// the whole point of a fractional threshold rather than requiring perfection.
    fn single_mismatch_reduces_score_but_stays_valid() {
        let reference = gen_seq(300, 11);
        let mut read = reference[50..130].to_vec();
        // flip one base in the middle to a different base
        let mid = read.len() / 2;
        read[mid] = if read[mid] == b'A' { b'C' } else { b'A' };
        let cfg = AlignConfig::default();
        let aln = align_chain(&read, &reference, &fw_chain(50, 31), &cfg).unwrap();
        let perfect = perfect_score(read.len(), &cfg);
        assert!(
            aln.score < perfect,
            "score {} should be < {perfect}",
            aln.score
        );
        // one mismatch: (len-1)*match + (-mismatch)
        let expected = (read.len() as i32 - 1) * cfg.match_score as i32 - cfg.mismatch_pen as i32;
        assert_eq!(aln.score, expected);
        assert!(aln.valid);
    }

    #[test]
    /// A read from the opposite strand must align just as well once transformed
    /// into the reference's frame.
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
    /// The rejection case: a read that did not come from here must fail the
    /// threshold, which is what makes the alignment check worth doing.
    fn unrelated_sequence_is_invalid() {
        let reference = gen_seq(300, 17);
        let foreign = gen_seq(80, 999); // unrelated
        let cfg = AlignConfig::default();
        let aln = align_chain(&foreign, &reference, &fw_chain(50, 31), &cfg).unwrap();
        assert!(
            !aln.valid,
            "unrelated read should not validate (score {})",
            aln.score
        );
    }

    #[test]
    /// Window-fitting must locate the read within the window rather than assuming
    /// it starts at the window's edge.
    fn fitting_window_places_read_in_middle() {
        // Read exactly matches window[40..120]; free ref-start must find it and
        // report the read's end column (120) with a perfect score.
        let window = gen_seq(200, 91);
        let read = window[40..120].to_vec();
        let cfg = AlignConfig::default();
        let a = align_in_window(&read, &window, &cfg).unwrap();
        assert_eq!(
            a.score,
            perfect_score(read.len(), &cfg),
            "score {}",
            a.score
        );
        assert!(a.valid);
        assert_eq!(a.end_col, 120, "end_col {}", a.end_col);
    }

    #[test]
    /// Same, with a sequencing error present.
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
    /// And it must still reject a read that does not belong in the window.
    fn fitting_window_rejects_unrelated() {
        let window = gen_seq(200, 93);
        let foreign = gen_seq(80, 54321);
        let cfg = AlignConfig::default();
        let a = align_in_window(&foreign, &window, &cfg).unwrap();
        assert!(
            !a.valid,
            "unrelated read should not validate (score {})",
            a.score
        );
    }

    #[test]
    /// A short indel must fall inside the alignment band; too narrow a band would
    /// silently miss legitimate alignments.
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
        assert!(
            aln.valid,
            "score {} should validate a 2bp deletion",
            aln.score
        );
        assert!(aln.score < perfect_score(read.len(), &cfg));
    }

    /// All three DP kinds must stay distinct for the same bytes. Previously the
    /// kind was smuggled into the key as a trailing marker byte on the query —
    /// nothing for a gap, `0` or `1` for the two flank orientations — which
    /// separated them only because no ACGT base equals `0x00` or `0x01`. That is
    /// an unasserted assumption about read content standing between two
    /// different DPs and a silently wrong score.
    #[test]
    fn every_dp_kind_is_distinct_for_the_same_bytes() {
        use GapKind::*;
        let kinds = [Gap, FlankLeftAnchored, FlankRightAnchored];
        let (q, t) = (b"ACGTACGT".as_slice(), b"ACGTACGTAC".as_slice());

        let mut c = GapCache::default();
        for (i, &k) in kinds.iter().enumerate() {
            c.insert(k, q, t, -(i as i32) - 1);
        }
        for (i, &k) in kinds.iter().enumerate() {
            assert_eq!(
                c.get(k, q, t),
                Some(-(i as i32) - 1),
                "kind {i} read back another kind's score for identical bytes"
            );
        }

        // And a kind that was never inserted must miss rather than borrow one.
        let mut c = GapCache::default();
        c.insert(Gap, q, t, -3);
        assert_eq!(c.get(FlankLeftAnchored, q, t), None);
        assert_eq!(c.get(FlankRightAnchored, q, t), None);
    }

    /// The entry stores query and target concatenated, so the split point is
    /// part of the key too: the same bytes divided differently are a different
    /// DP. `eq` compares `qlen`, and this is what says so.
    #[test]
    fn the_query_target_split_is_part_of_the_key() {
        let mut c = GapCache::default();
        c.insert(GapKind::Gap, b"AC", b"GTA", -1);
        assert_eq!(c.get(GapKind::Gap, b"AC", b"GTA"), Some(-1));
        assert_eq!(
            c.get(GapKind::Gap, b"ACG", b"TA"),
            None,
            "same concatenation, different split — must not alias"
        );
        assert_eq!(c.get(GapKind::Gap, b"A", b"CGTA"), None);
    }

    /// Entries must survive the table growing. The rehash closure has to derive
    /// each entry's hash from *its own* bytes; deriving it from the key being
    /// inserted would leave older entries stored under a hash no probe computes,
    /// so they would never be found again — a cache that silently stops caching.
    #[test]
    fn entries_survive_rehashing_on_growth() {
        let mut c = GapCache::default();
        let n = 4_096usize;
        for i in 0..n {
            let q = format!("Q{i}");
            let t = format!("T{i}");
            c.insert(GapKind::Gap, q.as_bytes(), t.as_bytes(), i as i32);
        }
        for i in 0..n {
            let q = format!("Q{i}");
            let t = format!("T{i}");
            assert_eq!(
                c.get(GapKind::Gap, q.as_bytes(), t.as_bytes()),
                Some(i as i32),
                "entry {i} was lost across growth"
            );
        }
    }

    /// The cap bounds per-thread memory by clearing wholesale rather than
    /// evicting, so the entry that trips it survives and the rest do not.
    #[test]
    fn reaching_the_cap_clears_the_table() {
        let mut c = GapCache::default();
        for i in 0..GAP_CACHE_CAP {
            let q = format!("q{i}");
            c.insert(GapKind::Gap, q.as_bytes(), b"t", i as i32);
        }
        assert_eq!(c.get(GapKind::Gap, b"q0", b"t"), Some(0));
        c.insert(GapKind::Gap, b"trip", b"t", -7);
        assert_eq!(c.get(GapKind::Gap, b"q0", b"t"), None, "cap should clear");
        assert_eq!(c.get(GapKind::Gap, b"trip", b"t"), Some(-7));
    }
}
