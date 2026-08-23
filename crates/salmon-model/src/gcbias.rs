//! Fragment-GC bias model (`GCFragModel`).
//!
//! # The effect being modelled
//!
//! PCR amplification, which every library prep relies on, does not copy all
//! fragments equally well. Fragments that are very GC-rich or very GC-poor
//! amplify less efficiently than middling ones, so they end up under-represented
//! in the sequencer's output. A transcript with unusual GC content therefore
//! looks less abundant than it is — an artefact of chemistry, not biology.
//!
//! A port of salmon's `GCFragModel` (`include/.../model/GCFragModel.hpp`): a
//! `condBins × numGCBins` table of fragment counts, where a fragment is keyed by
//! its **GC fraction** (0–100, → `numGCBins` bins) and a **conditioning context**
//! (a coarse sequence/context fraction, → `condBins` bins).
//!
//! **Why condition on the context.** GC content immediately *around* a fragment's
//! ends also affects priming and ligation. Splitting the table by that context
//! keeps the GC effect from absorbing what is really a sequence-composition
//! effect at the ends.
//!
//! An *observed* model (from mapped fragments) and an *expected* model (from the
//! transcriptome) are each normalized per conditioning row into a distribution,
//! then divided — [`gc_ratio`] — to give the per-`(context, GC)` bias factor used
//! to weight fragments in the bias-corrected effective length. The expected model
//! answers "what GC contents would we see if there were no GC bias at all, given
//! these transcripts and abundances?", so the ratio isolates the bias itself.
//!
//! salmon accumulates the observed model in log space and the expected in linear
//! space, but both are normalized to linear distributions before the ratio is
//! taken, so the normalized result is independent of the accumulation space.
//! We therefore accumulate in **linear** space (summing weights) throughout —
//! numerically equivalent after [`GcFragModel::normalize`], and simpler.

/// salmon's default number of conditioning (context) bins (`numConditionalGCBins`).
pub const DEFAULT_COND_BINS: usize = 3;
/// salmon's default number of fragment-GC bins for accumulation
/// (`salmon::defaults::numFragGCBins`).
pub const DEFAULT_GC_BINS: usize = 25;
/// salmon's `ratio()` clamp (`gcBias = gcCounts.ratio(transcriptGCDist, 1000.0)`).
///
/// A GC bin observed a handful of times against a near-zero expectation would
/// otherwise produce an enormous factor from pure noise; the clamp bounds how
/// far any single cell can move a transcript's effective length.
pub const GC_MAX_RATIO: f64 = 1000.0;
/// salmon's per-row normalization prior.
///
/// A pseudocount added to every bin before normalizing, so an unobserved GC bin
/// gets a small probability instead of exactly zero (which would make its ratio a
/// division by zero).
const NORM_PRIOR: f64 = 0.1;

/// Bin a 0–100 fraction into `n` bins (salmon's `GCDesc::fragBin(n)`/`contextBin(n)`:
/// `min(n-1, floor(frac / (100/n)))`). With `n == 101` this is the identity.
///
/// Binning trades resolution for statistical strength: 25 bins each see plenty of
/// fragments, whereas 101 would be noisy on a small sample.
#[inline]
pub fn bin_frac(frac: i32, n: usize) -> usize {
    if n == 101 {
        return (frac.clamp(0, 100)) as usize;
    }
    let w = 100.0 / n as f64;
    let b = (frac as f64 / w) as i32;
    b.clamp(0, n as i32 - 1) as usize
}

/// `(condBins × numGCBins)` table of fragment-GC counts (then a normalized
/// distribution / ratio after [`normalize`](Self::normalize) / [`gc_ratio`]).
#[derive(Debug, Clone)]
pub struct GcFragModel {
    cond_bins: usize,
    gc_bins: usize,
    /// row-major `counts[ctx * gc_bins + gc]`. Until [`normalize`], this is the
    /// f64 *view* materialized from `counts_fp`; the live accumulator during
    /// collection is the integer `counts_fp` (for thread-count-independent sums).
    counts: Vec<f64>,
    /// Fixed-point integer accumulator (`mass * BIAS_WEIGHT_SCALE`, truncated);
    /// summed order-independently across worker threads, materialized into
    /// `counts` once at [`normalize`]. See [`crate::BIAS_WEIGHT_SCALE`].
    counts_fp: Vec<u64>,
    normalized: bool,
    /// Precomputed `[0..=100] -> bin` maps for the context and GC fractions
    /// (both always lie in `[0, 100]`). [`bin_frac`] is ~44% of the eff-length
    /// convolution self-time (two `100.0/n` + `frac/w` divides per [`get`]/[`inc`]
    /// call, called per fragment); these LUTs replace both divides with a byte
    /// load. Each entry is exactly `bin_frac(frac, n)`, so the result is
    /// bit-identical to the divide path (verified in tests).
    ctx_lut: Vec<u8>,
    gc_lut: Vec<u8>,
}

impl GcFragModel {
    pub fn new(cond_bins: usize, gc_bins: usize) -> Self {
        let ctx_lut = (0..=100).map(|f| bin_frac(f, cond_bins) as u8).collect();
        let gc_lut = (0..=100).map(|f| bin_frac(f, gc_bins) as u8).collect();
        Self {
            cond_bins,
            gc_bins,
            counts: vec![0.0; cond_bins * gc_bins],
            counts_fp: vec![0u64; cond_bins * gc_bins],
            normalized: false,
            ctx_lut,
            gc_lut,
        }
    }

    #[inline]
    fn idx(&self, ctx_frac: i32, gc_frac: i32) -> usize {
        // Clamp-then-LUT is bit-identical to `bin_frac` for every i32 input:
        // for frac in [0,100] the LUT entry *is* `bin_frac(frac, n)`, and a frac
        // outside [0,100] clamps to the same boundary bin `bin_frac` would
        // produce (0 or n-1). Inputs are always lrint(100·ratio) in [0,100].
        let ctx = if self.cond_bins > 1 {
            self.ctx_lut[ctx_frac.clamp(0, 100) as usize] as usize
        } else {
            0
        };
        let gc = self.gc_lut[gc_frac.clamp(0, 100) as usize] as usize;
        ctx * self.gc_bins + gc
    }

    /// The flattened `cond_bins × gc_bins` table (`counts[ctx * gc_bins + gc]`,
    /// row-major), for dumping to the aux bias files.
    pub fn dump(&self) -> &[f64] {
        &self.counts
    }

    /// Accumulate `weight` for a fragment with the given GC and context fractions.
    pub fn inc(&mut self, gc_frac: i32, ctx_frac: i32, weight: f64) {
        debug_assert!(!self.normalized, "cannot inc a normalized model");
        let i = self.idx(ctx_frac, gc_frac);
        self.counts_fp[i] += crate::bias_mass_to_fp(weight);
    }

    /// Value at a `(context, GC)` cell (a normalized density after `normalize`,
    /// or a clamped ratio for a model produced by [`gc_ratio`]).
    #[inline]
    pub fn get(&self, gc_frac: i32, ctx_frac: i32) -> f64 {
        self.counts[self.idx(ctx_frac, gc_frac)]
    }

    /// Merge another (compatible) model's counts. Both must be pre-normalization.
    pub fn combine_counts(&mut self, other: &GcFragModel) {
        debug_assert!(!self.normalized && !other.normalized);
        // Integer sum — associative, so the merged model is independent of how
        // fragments were partitioned across worker threads.
        for (a, b) in self.counts_fp.iter_mut().zip(&other.counts_fp) {
            *a += *b;
        }
    }

    /// Normalize each conditioning row into a distribution over GC bins, with a
    /// pseudocount `prior` (salmon's default 0.1). Idempotent.
    ///
    /// Per *row*, not over the whole table: the ratio compares like with like
    /// within one context, so each context's GC distribution must sum to 1 on its
    /// own.
    pub fn normalize(&mut self) {
        if self.normalized {
            return;
        }
        // Materialize the integer accumulator into the f64 `counts` the rest of
        // the model (rows, `get`, `gc_ratio`, dumps) reads.
        for (c, &fp) in self.counts.iter_mut().zip(&self.counts_fp) {
            *c = fp as f64 / crate::BIAS_WEIGHT_SCALE;
        }
        for r in 0..self.cond_bins {
            let base = r * self.gc_bins;
            let row = &mut self.counts[base..base + self.gc_bins];
            let row_mass: f64 = row.iter().map(|c| NORM_PRIOR + *c).sum();
            if row_mass > 0.0 {
                let norm = 1.0 / row_mass;
                for c in row.iter_mut() {
                    *c = (NORM_PRIOR + *c) * norm;
                }
            }
        }
        self.normalized = true;
    }
}

/// Per-cell bias ratio `observed / expected`, clamped to `[1/max_ratio, max_ratio]`
/// (salmon's `GCFragModel::ratio`). Both models are normalized first.
///
/// A value above 1 means fragments of that GC/context were seen *more* often than
/// chance predicts, so such fragments are easy to sequence; below 1 means they are
/// suppressed. The effective-length convolution multiplies each candidate fragment
/// by its factor, so a transcript made of suppressed fragments gets a smaller
/// effective length and therefore a higher inferred abundance.
pub fn gc_ratio(
    observed: &mut GcFragModel,
    expected: &mut GcFragModel,
    max_ratio: f64,
) -> GcFragModel {
    observed.normalize();
    expected.normalize();
    let min_ratio = 1.0 / max_ratio;
    let mut out = GcFragModel::new(observed.cond_bins, observed.gc_bins);
    for i in 0..out.counts.len() {
        let e = expected.counts[i];
        let rat = if e != 0.0 {
            observed.counts[i] / e
        } else {
            max_ratio
        };
        out.counts[i] = rat.clamp(min_ratio, max_ratio);
    }
    out.normalized = true; // a ratio model is used directly, not re-normalized
    out
}

// ===========================================================================
// Fragment-GC content (salmon's `Transcript::GCCount_` / `gcFrac` / `gcDesc`)
// ===========================================================================

use crate::seqbias::{conditional_cdf, MIN_ALPHA, MIN_CDF_MASS};

/// salmon's fragment-length sampling stride for the GC convolution
/// (`pdfSampFactor` = `biasSpeedSamp`, default 5).
pub const GC_SAMP_STRIDE: usize = 5;

/// `gcDesc` context-window geometry (salmon `Transcript::gcDesc`):
/// `outsideContext = 3`, `insideContext = 2`.
const OUTSIDE_5P: i32 = 4; // outsideContext + 1
const OUTSIDE_3P: i32 = 3; // outsideContext
const INSIDE_5P: i32 = 1; // insideContext - 1
const INSIDE_3P: i32 = 2; // insideContext

/// Round to nearest (ties to even), matching C++ `std::lrint` under the default
/// rounding mode.
///
/// On x86-64 this is a single `cvtsd2si` instruction (which rounds per the
/// current — default round-to-nearest-even — MXCSR mode), identical to salmon's
/// `std::lrint`. The portable `f64::round_ties_even() as i32` lowers to a libm
/// `roundeven` *function call* on this baseline build; since this is called
/// ~2× per `gc_desc` over tens of billions of fragment evaluations under
/// `--gcBias`, that call was the dominant per-iteration cost (≈4× the GC-bias
/// runtime vs C++). Other architectures keep the portable path.
#[inline]
fn lrint(x: f64) -> i32 {
    #[cfg(target_arch = "x86_64")]
    {
        // SAFETY: SSE2 is part of the x86-64 baseline, so `_mm_set_sd` /
        // `_mm_cvtsd_si32` are always available. `cvtsd2si` matches `std::lrint`
        // (round-to-nearest-even); our arguments are small (0–100), no overflow.
        #[allow(unsafe_code)]
        unsafe {
            core::arch::x86_64::_mm_cvtsd_si32(core::arch::x86_64::_mm_set_sd(x))
        }
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        x.round_ties_even() as i32
    }
}

/// Cumulative G+C counts (salmon's `Transcript::GCCount_`): `prefix[p]` is the
/// number of `G`/`C` bases in `seq[0..=p]`.
///
/// A *prefix sum*: with it, the GC count of any interval is one subtraction
/// (`prefix[b] - prefix[a-1]`) instead of a scan. The convolution asks that
/// question for every (start, length) pair, so the difference is decisive.
pub fn gc_prefix(seq: &[u8]) -> Vec<u32> {
    let mut prefix = Vec::with_capacity(seq.len());
    let mut acc = 0u32;
    for &b in seq {
        if matches!(b, b'G' | b'g' | b'C' | b'c') {
            acc += 1;
        }
        prefix.push(acc);
    }
    prefix
}

/// GC count in the closed interval `[a, b]` from a cumulative-count prefix.
#[inline]
fn gc_in(prefix: &[u32], a: i32, b: i32) -> i64 {
    let lo = if a > 0 {
        prefix[(a - 1) as usize] as i64
    } else {
        0
    };
    prefix[b as usize] as i64 - lo
}

/// Fragment GC fraction (0–100) over the closed interval `[s, e]`
/// (salmon's `Transcript::gcFrac`): `round(100·GC[s,e] / (e−s+1))`.
#[inline]
pub fn gc_frac(prefix: &[u32], s: i32, e: i32) -> i32 {
    lrint(100.0 * gc_in(prefix, s, e) as f64 / (e - s + 1) as f64)
}

/// A rank-enabled bitvector over the *concatenated* reference sequence marking
/// G/C positions (salmon's `--reduceGCMemory` representation). Replaces the
/// per-transcript dense cumulative-GC arrays (`Vec<Vec<u32>>`, 4 bytes/base) with
/// one bitvector (~1 bit/base + 25% rank overhead) — GC in any interval is a
/// difference of two O(1) ranks. Build once over the whole transcriptome.
///
/// A *rank* query answers "how many bits are set before position `i`?" in
/// constant time, using a small precomputed table of block counts. That is
/// exactly a prefix sum, at roughly 1/25th of the memory: on a human
/// transcriptome plus genome decoy, ~300 MB becomes ~12 MB.
pub struct GcRank {
    rank: sux::rank_sel::Rank9,
}

impl GcRank {
    /// Build the rank bitvector from the concatenated reference bytes
    /// (`seq[i]` is G/C ⇒ bit `i` set).
    pub fn new(concat: &[u8]) -> GcRank {
        use sux::traits::bit_vec_ops::BitVecOpsMut;
        let mut bv = sux::bits::BitVec::new(concat.len());
        for (i, &b) in concat.iter().enumerate() {
            if matches!(b, b'G' | b'g' | b'C' | b'c') {
                bv.set(i, true);
            }
        }
        GcRank {
            rank: sux::rank_sel::Rank9::new(bv),
        }
    }

    /// A cumulative-GC view of the transcript occupying `concat[off..off+len]`.
    #[inline]
    pub fn view(&self, off: usize, len: usize) -> GcView<'_> {
        use sux::traits::Rank;
        GcView::Rank {
            rank: &self.rank,
            off,
            base: self.rank.rank(off),
            len,
        }
    }
}

/// A per-transcript cumulative-GC accessor: either the dense prefix array
/// (default) or a window of the shared [`GcRank`] (`--reduceGCMemory`). Both
/// return the *same* cumulative counts, so results are identical; `cum(p)` is
/// the number of G/C bases in transcript-local `[0, p]`.
///
/// An enum rather than a trait object: the choice is made once per run, and this
/// keeps the call a predictable branch instead of an indirect jump on a path
/// taken billions of times.
#[derive(Clone, Copy)]
pub enum GcView<'a> {
    Dense(&'a [u32]),
    Rank {
        rank: &'a sux::rank_sel::Rank9,
        off: usize,
        /// `rank(off)` — precomputed so `cum` is a single rank query.
        base: usize,
        len: usize,
    },
}

impl GcView<'_> {
    /// Cumulative G/C count in transcript-local `[0, p]` (inclusive).
    #[inline]
    pub fn cum(&self, p: i32) -> i64 {
        match self {
            GcView::Dense(prefix) => prefix[p as usize] as i64,
            GcView::Rank {
                rank, off, base, ..
            } => {
                use sux::traits::Rank;
                (rank.rank(off + p as usize + 1) - base) as i64
            }
        }
    }

    /// Transcript length in bases.
    #[inline]
    pub fn ref_len(&self) -> usize {
        match self {
            GcView::Dense(prefix) => prefix.len(),
            GcView::Rank { len, .. } => *len,
        }
    }
}

/// A whole-transcriptome source of per-transcript [`GcView`]s: either dense
/// per-transcript prefixes (default) or one shared [`GcRank`] window per
/// transcript (`--reduceGCMemory`). Lets callers stay representation-agnostic.
#[derive(Clone, Copy)]
pub enum GcStore<'a> {
    Dense(&'a [Vec<u32>]),
    Rank {
        rank: &'a GcRank,
        offsets: &'a [u64],
    },
}

impl<'a> GcStore<'a> {
    /// The cumulative-GC view for transcript `tid`. The view borrows the
    /// underlying data (lifetime `'a`), not `self`, so it outlives a temporary
    /// `GcStore`.
    #[inline]
    pub fn view(&self, tid: usize) -> GcView<'a> {
        match *self {
            GcStore::Dense(prefixes) => GcView::Dense(&prefixes[tid]),
            GcStore::Rank { rank, offsets } => {
                let off = offsets[tid] as usize;
                let len = (offsets[tid + 1] - offsets[tid]) as usize;
                rank.view(off, len)
            }
        }
    }
}

/// Fragment GC descriptor `(fragFrac, contextFrac)` for the closed fragment
/// `[s, e]` (salmon's `Transcript::gcDesc`). The context fraction is the GC
/// content over a 5-base 5' window `[s−3, s+1]` and a 5-base 3' window
/// `[e−1, e+3]` (edge-clamped). Returns `None` when the context window is empty
/// (matching salmon's `valid = false`).
///
/// The two returned numbers are exactly the table's two axes: how GC-rich the
/// fragment is overall, and how GC-rich the short windows straddling its two ends
/// are. The windows deliberately extend a few bases *outside* the fragment,
/// because priming and ligation see the sequence just beyond the cut as well.
#[inline]
pub fn gc_desc(v: &GcView, s: i32, e: i32) -> Option<(i32, i32)> {
    let last = v.ref_len() as i32 - 1;
    let cs = if s > 0 { v.cum(s - 1) } else { 0 };
    let ce = v.cum(e);

    let fs = s - OUTSIDE_5P;
    let fe = s + INSIDE_5P;
    let ts = e - INSIDE_3P;
    let te = e + OUTSIDE_3P;

    let fp_left = fs >= 0;
    let fp_right = fe <= last;
    let tp_left = ts >= 0;
    let tp_right = te <= last;

    let fps = if fp_left { v.cum(fs) } else { 0 };
    let fpe = if fp_right { v.cum(fe) } else { ce };
    let tps = if tp_left { v.cum(ts) } else { 0 };
    let tpe = if tp_right { v.cum(te) } else { ce };

    let fs_c = fs.max(0);
    let fe_c = fe.min(last);
    let ts_c = ts.max(0);
    let te_c = te.min(last);
    let fp_context_size = if !fp_left { fe_c + 1 } else { fe_c - fs_c };
    let tp_context_size = if !tp_left { te_c + 1 } else { te_c - ts_c };
    let context_size = (fp_context_size + tp_context_size) as f64;
    if context_size == 0.0 {
        return None;
    }

    let frag_frac = lrint(100.0 * (ce - cs) as f64 / (e - s + 1) as f64);
    let context_frac = lrint(100.0 * ((fpe - fps) + (tpe - tps)) as f64 / context_size);
    Some((frag_frac, context_frac))
}

/// Per-position context-GC cache for a single transcript — salmon's
/// `populateContextCounts`. [`gc_desc`] recomputes the 5'/3' context window
/// geometry (four conditional prefix loads, four edge branches, a division)
/// for *every* `(fragStart, fragEnd)` pair in the effective-length convolution;
/// since the 5' term depends only on `fragStart` and the 3' term only on
/// `fragEnd`, we hoist both to O(ref_len) per-position arrays and reduce the
/// inner loop to two array loads + one division + one `lrint` — matching the
/// C++ hot path. [`GcContext::desc`] is bit-identical to [`gc_desc`] for every
/// fragment the convolution actually visits (`fragStart + INSIDE_5P <= last`,
/// i.e. `fp_right` always holds — guaranteed since `fragStart <= ref_len-2`).
pub struct GcContext {
    cum: Vec<u32>,      // cumulative G/C count in [0, p], indexed by position
    fp_count: Vec<i64>, // 5' window GC count, indexed by fragStart
    fp_wlen: Vec<i32>,  // 5' window length, indexed by fragStart
    tp_count: Vec<i64>, // 3' window GC count, indexed by fragEnd
    tp_wlen: Vec<i32>,  // 3' window length, indexed by fragEnd
}

impl GcContext {
    /// Precompute the per-position 5'/3' context arrays (and a transient copy of
    /// the cumulative-GC counts) from any [`GcView`]. Each entry replicates
    /// exactly the corresponding sub-expression of [`gc_desc`] (same
    /// edge-clamping, same `fp_right`/`tp_right` handling). Materializing `cum`
    /// here keeps the hot inner loop ([`GcContext::desc`]) free of any rank
    /// queries or dense lookups beyond cached arrays.
    pub fn build(v: &GcView) -> GcContext {
        let n = v.ref_len();
        let last = n as i32 - 1;
        let cum: Vec<u32> = (0..n).map(|p| v.cum(p as i32) as u32).collect();
        let at = |i: i32| cum[i as usize] as i64;
        let mut fp_count = vec![0i64; n];
        let mut fp_wlen = vec![0i32; n];
        let mut tp_count = vec![0i64; n];
        let mut tp_wlen = vec![0i32; n];
        for p in 0..n {
            let pos = p as i32;
            // 3' window for a fragment *ending* at `pos`.
            let ts = pos - INSIDE_3P;
            let te = pos + OUTSIDE_3P;
            let tp_left = ts >= 0;
            let tp_right = te <= last;
            let tps = if tp_left { at(ts) } else { 0 };
            // `tpe = ce = cum[e]` when the 3' window runs off the end.
            let tpe = if tp_right { at(te) } else { at(pos) };
            let ts_c = ts.max(0);
            let te_c = te.min(last);
            tp_count[p] = tpe - tps;
            tp_wlen[p] = if !tp_left { te_c + 1 } else { te_c - ts_c };

            // 5' window for a fragment *starting* at `pos`. Within the
            // convolution `fp_right` always holds; the `cum[p]` fallback is a
            // placeholder for the never-visited `pos == last` slot.
            let fs = pos - OUTSIDE_5P;
            let fe = pos + INSIDE_5P;
            let fp_left = fs >= 0;
            let fp_right = fe <= last;
            let fps = if fp_left { at(fs) } else { 0 };
            let fpe = if fp_right { at(fe) } else { at(pos) };
            let fs_c = fs.max(0);
            let fe_c = fe.min(last);
            fp_count[p] = fpe - fps;
            fp_wlen[p] = if !fp_left { fe_c + 1 } else { fe_c - fs_c };
        }
        GcContext {
            cum,
            fp_count,
            fp_wlen,
            tp_count,
            tp_wlen,
        }
    }

    /// GC descriptor `(fragFrac, contextFrac)` for the closed fragment `[s, e]`,
    /// bit-identical to [`gc_desc`] but using the cached context arrays. `s` and
    /// `e` must lie in `0..ref_len` with `s + INSIDE_5P <= ref_len - 1`.
    #[inline]
    pub fn desc(&self, s: i32, e: i32) -> Option<(i32, i32)> {
        let context_size = (self.fp_wlen[s as usize] + self.tp_wlen[e as usize]) as f64;
        if context_size == 0.0 {
            return None;
        }
        let cs = if s > 0 {
            self.cum[(s - 1) as usize] as i64
        } else {
            0
        };
        let ce = self.cum[e as usize] as i64;
        let count = (self.fp_count[s as usize] + self.tp_count[e as usize]) as f64;
        let frag_frac = lrint(100.0 * (ce - cs) as f64 / (e - s + 1) as f64);
        let context_frac = lrint(100.0 * count / context_size);
        Some((frag_frac, context_frac))
    }
}

/// Build the *expected* fragment-GC model (salmon's `transcriptGCDist`): slide
/// every fragment `[fragStart, fragStart+fl−1]` over each expressed transcript,
/// weighting each by `(alpha/effLen)·(conditionalCDF(fl) − prevMass)` — the same
/// abundance-density × conditional-FLD weighting used for the expected
/// sequence-bias model. `k` is the leading offset salmon excludes from the
/// fragment-start loop (`9` with `--seqBias`, `1` otherwise); `stride` subsamples
/// fragment lengths ([`GC_SAMP_STRIDE`]).
///
/// In words: enumerate every fragment that *could* have been sequenced, weight it
/// by how abundant its transcript is and how likely its length is, and tally its
/// GC content. The result is the GC distribution an unbiased protocol would have
/// produced. Weighting by abundance matters because a highly expressed transcript
/// dominates the observed distribution and must dominate the expected one too.
///
/// `stride` is a pure speed knob: the fragment-length distribution is smooth, so
/// sampling every fifth length changes the model negligibly while cutting the
/// work fivefold.
#[allow(clippy::too_many_arguments)]
pub fn build_expected_gc<'a, FS, FP>(
    num_targets: usize,
    seq_of: FS,
    view_of: FP,
    alphas: &[f64],
    eff_lens: &[f64],
    cdf: &[f64],
    fld_low: usize,
    fld_high: usize,
    cond_bins: usize,
    gc_bins: usize,
    k: usize,
    stride: usize,
) -> GcFragModel
where
    FS: Fn(usize) -> &'a [u8] + Sync,
    FP: Fn(usize) -> GcView<'a> + Sync,
{
    let stride = stride.max(1) as i32;
    // The expected-GC distribution is a sum of independent per-transcript
    // contributions, each an O(refLen · fldRange/stride) double loop — billions
    // of `gc_desc` calls in total. salmon parallelizes this over transcripts;
    // do the same with rayon (per-thread `GcFragModel` partials, reduced via
    // `combine_counts`). `seq_of`/`prefix_of` must be `Sync` to share across
    // threads (they are: closures over `&[Vec<_>]`).
    use rayon::prelude::*;
    let per_tid = |tid: usize| -> Option<GcFragModel> {
        if alphas[tid] < MIN_ALPHA || eff_lens[tid] <= 0.0 {
            return None;
        }
        let seq = seq_of(tid);
        let ref_len = seq.len();
        if ref_len <= k {
            return None;
        }
        let cdf_max_arg = (cdf.len() - 1).min(ref_len);
        let cdf_max_val = cdf[cdf_max_arg];
        if cdf_max_val < MIN_CDF_MASS {
            return None;
        }
        let view = view_of(tid);
        let ctx = GcContext::build(&view);
        let weight = alphas[tid] / eff_lens[tid];
        let cond = |x: i32| conditional_cdf(cdf, cdf_max_arg, cdf_max_val, x);
        let sp = if fld_low > 0 { fld_low as i32 - 1 } else { 0 };
        let mut model = GcFragModel::new(cond_bins, gc_bins);
        for frag_start in 0..(ref_len - k) {
            let mut prev = cond(sp);
            let mut fl = fld_low as i32;
            while fl <= fld_high as i32 {
                let frag_end = frag_start as i32 + fl - 1;
                if (frag_end as usize) < ref_len {
                    if let Some((ff, cf)) = ctx.desc(frag_start as i32, frag_end) {
                        model.inc(ff, cf, weight * (cond(fl) - prev));
                    }
                    prev = cond(fl);
                } else {
                    break;
                }
                fl += stride;
            }
        }
        Some(model)
    };
    (0..num_targets)
        .into_par_iter()
        .fold(
            || GcFragModel::new(cond_bins, gc_bins),
            |mut acc, tid| {
                if let Some(m) = per_tid(tid) {
                    acc.combine_counts(&m);
                }
                acc
            },
        )
        .reduce(
            || GcFragModel::new(cond_bins, gc_bins),
            |mut a, b| {
                a.combine_counts(&b);
                a
            },
        )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 101 bins over a 0-100 fraction is one bin per value, i.e. no binning.
    #[test]
    fn bin_frac_identity_at_101() {
        assert_eq!(bin_frac(0, 101), 0);
        assert_eq!(bin_frac(57, 101), 57);
        assert_eq!(bin_frac(100, 101), 100);
    }

    /// Coarse binning must floor, and must clamp 100 into the last bin rather
    /// than overflowing past it.
    #[test]
    fn bin_frac_coarse() {
        // 3 bins: width 33.33 -> [0,33),[33,66),[66,100]
        assert_eq!(bin_frac(0, 3), 0);
        assert_eq!(bin_frac(30, 3), 0);
        assert_eq!(bin_frac(40, 3), 1);
        assert_eq!(bin_frac(70, 3), 2);
        assert_eq!(bin_frac(100, 3), 2); // clamped to n-1
    }

    /// Each conditioning row is an independent distribution, so each must sum to
    /// 1 on its own.
    #[test]
    fn normalize_makes_rows_sum_to_one() {
        let mut m = GcFragModel::new(2, 5);
        for gc in 0..5 {
            m.inc(gc * 25, 10, (gc + 1) as f64); // row 0 (ctx 10 -> bin 0)
        }
        m.normalize();
        let row0: f64 = (0..5).map(|c| m.counts[c]).sum();
        assert!((row0 - 1.0).abs() < 1e-9, "row0 sum = {row0}");
    }

    /// No bias in the data must mean no correction: observed identical to
    /// expected has to give a factor of exactly 1 everywhere.
    #[test]
    fn ratio_of_identical_models_is_one() {
        let mut obs = GcFragModel::new(1, 101);
        let mut exp = GcFragModel::new(1, 101);
        for gc in 0..=100 {
            obs.inc(gc, 0, (gc + 1) as f64);
            exp.inc(gc, 0, (gc + 1) as f64);
        }
        let r = gc_ratio(&mut obs, &mut exp, GC_MAX_RATIO);
        for gc in 0..=100 {
            assert!(
                (r.get(gc, 0) - 1.0).abs() < 1e-9,
                "gc {gc}: {}",
                r.get(gc, 0)
            );
        }
    }

    #[test]
    fn gc_frac_basic() {
        // 10-base sequence, 5 G/C -> 50%
        let seq = b"ACGTACGTAC"; // A C G T A C G T A C : GC at 1,2,5,6,9 = 5
        let p = gc_prefix(seq);
        assert_eq!(gc_frac(&p, 0, 9), 50);
        // all-GC window
        let seq2 = b"GGGGCCCC";
        let p2 = gc_prefix(seq2);
        assert_eq!(gc_frac(&p2, 0, 7), 100);
        // single base
        assert_eq!(gc_frac(&p, 2, 2), 100); // 'G'
        assert_eq!(gc_frac(&p, 0, 0), 0); // 'A'
    }

    #[test]
    fn gc_desc_matches_frag_and_context() {
        // 40-base sequence; check an interior fragment
        let seq: Vec<u8> = b"ACGT".iter().cycle().take(40).copied().collect();
        let p = gc_prefix(&seq);
        let (ff, cf) = gc_desc(&GcView::Dense(&p), 10, 29).unwrap();
        // ACGT repeat is exactly 50% GC everywhere
        assert_eq!(ff, 50, "fragFrac");
        assert!((0..=100).contains(&cf), "contextFrac in range: {cf}");
        // fragFrac must equal gc_frac over the same interval
        assert_eq!(ff, gc_frac(&p, 10, 29));
    }

    /// The lookup table is a performance substitution for the divide-based
    /// formula, so it has to agree bit-for-bit on every input the model can see,
    /// including out-of-range ones that must clamp identically.
    #[test]
    fn idx_lut_matches_bin_frac() {
        // The LUT-based idx must equal the original divide-based formula for
        // every (ctx, gc) fraction the model can see, including out-of-[0,100].
        for &(cb, gb) in &[(3usize, 25usize), (1, 101), (3, 101), (4, 50)] {
            let m = GcFragModel::new(cb, gb);
            for cf in -5i32..=105 {
                for ff in -5i32..=105 {
                    let want_ctx = if cb > 1 { bin_frac(cf, cb) } else { 0 };
                    let want = want_ctx * gb + bin_frac(ff, gb);
                    assert_eq!(m.idx(cf, ff), want, "cb={cb} gb={gb} cf={cf} ff={ff}");
                }
            }
        }
    }

    /// Same obligation for the hoisted per-position context cache: it exists only
    /// to be faster, so any disagreement with the direct computation is a bug.
    #[test]
    fn gc_context_matches_gc_desc() {
        // Cached context (GcContext::desc) must be bit-identical to the
        // per-fragment gc_desc for every fragment the convolution can visit
        // (fragStart in 0..ref_len-fl, fragEnd = fragStart+fl-1).
        let bases = *b"ACGTCGATGC";
        let seq: Vec<u8> = (0..200).map(|i| bases[(i * 7 + 3) % bases.len()]).collect();
        let p = gc_prefix(&seq);
        let v = GcView::Dense(&p);
        let ctx = GcContext::build(&v);
        let n = seq.len() as i32;
        for fl in 1..=120 {
            let kmax = n - fl;
            for s in 0..kmax {
                let e = s + fl - 1;
                assert_eq!(
                    ctx.desc(s, e),
                    gc_desc(&v, s, e),
                    "mismatch at fl={fl}, s={s}, e={e}"
                );
            }
        }
    }

    /// And for the rank-bitvector representation: --reduceGCMemory must change
    /// memory use and nothing else. The offset window also checks the transcript
    /// is located correctly inside the concatenated buffer.
    #[test]
    fn gc_rank_view_matches_dense() {
        // The rank-bitvector view must return identical cumulative GC (and hence
        // identical gc_desc) to the dense prefix, for an offset window of a
        // concatenated buffer.
        let bases = *b"ACGTGCAT";
        let pad: Vec<u8> = (0..37).map(|i| bases[(i * 3) % bases.len()]).collect();
        let seq: Vec<u8> = (0..150).map(|i| bases[(i * 5 + 1) % bases.len()]).collect();
        // concat = [pad | seq], so the transcript sits at offset pad.len().
        let mut concat = pad.clone();
        concat.extend_from_slice(&seq);
        let p = gc_prefix(&seq);
        let rank = GcRank::new(&concat);
        let dense = GcView::Dense(&p);
        let rview = rank.view(pad.len(), seq.len());
        for fl in 1..=100 {
            let n = seq.len() as i32;
            for s in 0..(n - fl) {
                let e = s + fl - 1;
                assert_eq!(
                    gc_desc(&dense, s, e),
                    gc_desc(&rview, s, e),
                    "fl={fl} s={s}"
                );
            }
        }
    }

    /// Pin the exact window geometry by making only those windows GC-rich: if the
    /// offsets were off by one the fraction would not come out at exactly 100.
    #[test]
    fn gc_desc_context_window_geometry() {
        // Make the 5' context window [s-3, s+1] all GC and the rest AT, so the
        // contextFrac is dominated by the GC island near the fragment ends.
        let mut seq = vec![b'A'; 60];
        for i in 17..=21 {
            seq[i] = b'G'; // 5' window for s=20 is [17,21]
        }
        for i in 39..=43 {
            seq[i] = b'C'; // 3' window for e=40 is [39,43]
        }
        let p = gc_prefix(&seq);
        let (_ff, cf) = gc_desc(&GcView::Dense(&p), 20, 40).unwrap();
        // both 5-base context windows are fully G/C -> 100%
        assert_eq!(cf, 100, "contextFrac");
    }

    /// A cell enriched in one model and absent from the other would give an
    /// unbounded factor; the clamp must contain it.
    #[test]
    fn ratio_clamped() {
        let mut obs = GcFragModel::new(1, 2);
        let mut exp = GcFragModel::new(1, 2);
        obs.inc(0, 0, 1e6); // bin 0 hugely enriched in obs
        exp.inc(75, 0, 1e6); // bin 1 enriched in exp
        let r = gc_ratio(&mut obs, &mut exp, 1000.0);
        for v in [r.get(0, 0), r.get(75, 0)] {
            assert!(
                (1.0 / 1000.0 - 1e-12..=1000.0 + 1e-6).contains(&v),
                "ratio {v} out of clamp"
            );
        }
    }
}
