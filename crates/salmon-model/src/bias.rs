//! Unified bias-corrected effective length: composes sequence-specific, GC, and
//! positional bias exactly as salmon's `updateEffectiveLengths` does — a single
//! conditional-FLD convolution whose per-fragment factor is the product of the
//! enabled bias terms.
//!
//! # What this computes, in words
//!
//! An uncorrected effective length counts how many positions a fragment *could*
//! start at. A bias-corrected one counts how many positions a fragment could
//! start at **weighted by how likely a fragment there is to actually be
//! sequenced**. A transcript whose sequence, GC content and geometry all
//! discourage sequencing has a smaller corrected effective length, so the same
//! number of observed fragments implies a higher abundance.
//!
//! So the calculation sweeps every (start, length) pair, multiplies together the
//! bias factors that apply to that particular fragment, and weights the result by
//! how probable that fragment length is:
//!
//! `fragFactor = seqFW[fragStart]·seqRC[fragEnd] · gcBias({fragFrac, ctxFrac}) ·
//! posFW[fragStart]·posRC[fragEnd]`, summed over fragment starts and convolved
//! with the conditional fragment-length distribution.
//!
//! "Conditional" because a transcript shorter than the longest possible fragment
//! can only host part of the length distribution, so the distribution is
//! renormalized to the lengths that actually fit.

use crate::gcbias::{GcContext, GcFragModel, GcView};
use crate::posbias::{length_class_index, SimplePosBias, NUM_LENGTH_CLASSES, NUM_POS_BINS};
use crate::seqbias::{
    conditional_cdf, revcomp_bytes, LogBiasTable, CONTEXT_LEFT, CONTEXT_LENGTH, MIN_ALPHA,
    MIN_CDF_MASS,
};

/// salmon's `EPSILON` (mass cutoff for adding positional expected mass).
/// Contributions below it cannot affect a normalized 20-bin model, so skipping
/// them saves work without changing the result.
const EPSILON: f64 = 0.375e-10;

/// Additive (Laplace) smoothing fraction for the positional-bias factor: the
/// smoothing constant is `POS_SMOOTH_FRAC · mean(expected density)`.
const POS_SMOOTH_FRAC: f64 = 0.1;

/// Per-position positional-bias factors `obs/exp`, additively smoothed toward 1.
///
/// Instead of dividing two projected densities that the old `0.001` floor drove
/// to a tiny, noisy denominator in the tails (amplifying small model
/// differences), we add a smoothing constant `c = POS_SMOOTH_FRAC·mean(exp)` to
/// both: `factor = (obs + c) / (exp + c)`. Where the expected density is
/// substantial the ratio is preserved; where it vanishes (uninformative tails)
/// the factor shrinks to 1 (no bias), which is the correct default.
///
/// This is the standard fix for a ratio of two noisy estimates: adding the same
/// constant to numerator and denominator leaves well-measured ratios essentially
/// unchanged while pulling ill-measured ones toward the neutral value.
pub fn positional_factor(obs: &[f64], exp: &[f64]) -> Vec<f64> {
    let mut out = Vec::new();
    positional_factor_into(obs, exp, &mut out);
    out
}

/// [`positional_factor`] writing into a caller-provided buffer (cleared first),
/// so a per-worker `Vec` can be reused across transcripts instead of allocating
/// a fresh result for each one.
///
/// Computes exactly what [`positional_factor`] does; that function is now a thin
/// wrapper, so the two cannot drift apart.
pub fn positional_factor_into(obs: &[f64], exp: &[f64], out: &mut Vec<f64>) {
    out.clear();
    let n = exp.len();
    if n == 0 {
        return;
    }
    // Scaling the constant to the mean makes the smoothing strength independent
    // of the arbitrary units the densities happen to be in.
    let mean_exp: f64 = exp.iter().sum::<f64>() / n as f64;
    let c = (POS_SMOOTH_FRAC * mean_exp).max(f64::MIN_POSITIVE);
    out.extend(obs.iter().zip(exp).map(|(&o, &e)| (o + c) / (e + c)));
}

/// The enabled bias terms for one transcript's effective-length correction.
///
/// Each field is `Option`, so "this correction is off" is represented in the type
/// rather than by a neutral all-ones table — which lets the hot loops below skip
/// the work entirely instead of multiplying by 1.
#[derive(Clone, Copy, Default)]
pub struct BiasInputs<'a> {
    /// `(fw, rc)` precomputed `obs − exp` log-bias tables for the 5′/3′
    /// sequence-bias factors (`--seqBias`). Built once per quant run.
    pub seq: Option<(&'a LogBiasTable, &'a LogBiasTable)>,
    /// GC ratio model + this transcript's cumulative-GC view (`--gcBias`)
    pub gc: Option<(&'a GcFragModel, GcView<'a>)>,
    /// per-position 5'/3' positional-bias factors (`--posBias`), transcript-length sized
    pub pos: Option<(&'a [f64], &'a [f64])>,
}

impl BiasInputs<'_> {
    /// Whether any correction is enabled at all; if not, the whole computation is
    /// skipped.
    fn any(&self) -> bool {
        self.seq.is_some() || self.gc.is_some() || self.pos.is_some()
    }
}

/// Bias-corrected effective length composing every enabled bias term, matching
/// salmon's combined `updateEffectiveLengths` convolution. Floors at the lower
/// barrier `min(elen, max(1, unprocessedLen))` (no upper cap), unless
/// `no_length_threshold` is set (salmon's `--noBiasLengthThreshold`), in which
/// case the corrected length is accepted outright (floored only at 1.0) or the
/// uncorrected length is kept.
///
/// The floor exists because effective length is a *divisor*: an unluckily small
/// corrected length would inflate that transcript's abundance without bound, so
/// the correction is allowed to shrink the length only so far.
#[allow(clippy::too_many_arguments)]
pub fn corrected_effective_length_full(
    seq: &[u8],
    cdf: &[f64],
    fld_low: usize,
    fld_high: usize,
    bias: &BiasInputs,
    elen: f64,
    stride: usize,
    no_length_threshold: bool,
) -> f64 {
    corrected_effective_length_strided(
        seq,
        cdf,
        fld_low,
        fld_high,
        bias,
        elen,
        stride,
        no_length_threshold,
        1,
    )
}

/// Position stride of the bias sweeps over an unspliced target (see
/// [`unspliced_position_stride`]).
pub const UNSPLICED_BIAS_STRIDE: usize = 16;

/// Fewest fragment starts an unspliced target keeps in the bias sweeps; the
/// stride shrinks below [`UNSPLICED_BIAS_STRIDE`] for targets too short to
/// provide them.
pub const UNSPLICED_BIAS_MIN_STARTS: usize = 64;

/// Position stride for the bias sweeps of an unspliced target of length `len`:
/// every 16th fragment start, weighted by 16 (fewer for targets under 1 kb).
///
/// A fragment's GC content and end contexts change little when its start moves
/// by a few bases, since it spans ~250 of them, so the per-start bias factor is
/// smooth and systematic sampling at a stride well below the fragment length
/// estimates the full sweep closely (see the tests). Unspliced targets are
/// what makes the sweep expensive: 2.4 Gb of them against 1.5 Gb of
/// transcripts for GENCODE v50 intron targets. Spliced targets always use 1.
pub fn unspliced_position_stride(len: usize) -> usize {
    (len / UNSPLICED_BIAS_MIN_STARTS).clamp(1, UNSPLICED_BIAS_STRIDE)
}

/// [`corrected_effective_length_full`] with the scalar (GC) sweep visiting only
/// every `pos_stride`-th fragment start, each visit weighted by the stride.
/// `pos_stride = 1` is the exact computation, bit for bit. The FFT path (no GC
/// model) is already O(L log L) and ignores the stride.
#[allow(clippy::too_many_arguments)]
pub fn corrected_effective_length_strided(
    seq: &[u8],
    cdf: &[f64],
    fld_low: usize,
    fld_high: usize,
    bias: &BiasInputs,
    elen: f64,
    stride: usize,
    no_length_threshold: bool,
    pos_stride: usize,
) -> f64 {
    let pos_step = pos_stride.max(1) as i32;
    // Multiplying by exactly 1.0 is the identity, so stride 1 is unchanged.
    let pos_scale = pos_step as f64;
    if !bias.any() {
        return elen;
    }
    // Sequence bias needs a full context window at each end; without it a single
    // base suffices as the minimum workable length.
    let k = if bias.seq.is_some() {
        CONTEXT_LENGTH
    } else {
        1
    };
    let ref_len = seq.len();
    // How much of the transcript the uncorrected effective length has "used up";
    // also the floor's basis below.
    let unprocessed = (ref_len as i32 - elen as i32).max(0);
    let cdf_max_arg = (cdf.len() - 1).min(ref_len);
    let cdf_max_val = cdf[cdf_max_arg];
    // Bail out when there is nothing meaningful to correct: too short for a
    // context, no headroom, or almost no fragment-length mass fits at all.
    if ref_len < k || unprocessed <= 0 || cdf_max_val < MIN_CDF_MASS {
        return elen;
    }
    // The fragment-length CDF renormalized to lengths this transcript can host.
    let cond = |x: i32| conditional_cdf(cdf, cdf_max_arg, cdf_max_val, x);

    // Per-position sequence-bias factors. Only built (and applied in the inner
    // loop) when `--seqBias` is on; otherwise the factors are all 1.0, so we
    // skip both the per-transcript allocation and the per-fragment multiply
    // entirely (the common `--gcBias`-only case).
    let have_seq = bias.seq.is_some();
    let mut fw: Vec<f64> = Vec::new();
    let mut rc: Vec<f64> = Vec::new();
    if let Some((tab_fw, tab_rc)) = bias.seq {
        fw = vec![1.0f64; ref_len];
        rc = vec![1.0f64; ref_len];
        let cu = CONTEXT_LEFT;
        // The 3' factor reads the reverse complement, because the second read of
        // a pair is sequenced from the opposite strand.
        let rc_seq = revcomp_bytes(seq);
        for frag_start in 0..(ref_len - CONTEXT_LENGTH) {
            // The factor is indexed by the *read* start, which sits `CONTEXT_LEFT`
            // bases into the context window.
            let read_start = frag_start + cu;
            if read_start < ref_len {
                fw[read_start] = tab_fw
                    .eval(&seq[frag_start..frag_start + CONTEXT_LENGTH], false)
                    .exp();
                rc[read_start] = tab_rc
                    .eval(&rc_seq[frag_start..frag_start + CONTEXT_LENGTH], false)
                    .exp();
            }
        }
        // `rc` was filled in reverse-complement coordinates; flipping it puts it
        // back in forward coordinates so both arrays index the same way.
        rc.reverse();
    }

    let gc_model = bias.gc.map(|(m, _)| m);
    // Precompute the per-position 5'/3' context-GC arrays once per transcript
    // (salmon's `populateContextCounts`) so the inner convolution does cheap
    // array lookups instead of re-deriving the context geometry per fragment.
    let gc_ctx = bias.gc.map(|(_, view)| GcContext::build(&view));
    let (pos_fw, pos_rc) = match bias.pos {
        Some((a, b)) => (Some(a), Some(b)),
        None => (None, None),
    };

    // No-GC fast path: with GC absent, the per-fragment factor is
    // `(seqFW·posFW)[start] · (seqRC·posRC)[end]` — separable into a start array
    // `a` and an end array `b`, so the length sweep is their cross-correlation,
    // computed once via FFT in O(L log L) (vs the O(L·n_len) scalar convolution).
    // GC's windowed-GC binning couples start×length and is not separable, so any
    // GC run stays on the scalar path below. (`fw`/`rc` are length `ref_len` when
    // `have_seq`, else empty; the positional factors are length `ref_len`.)
    //
    // "Separable" is the key word: if a fragment's weight is (something about its
    // start) × (something about its end), then summing over all fragments of a
    // given length is exactly a cross-correlation of the two arrays — and a
    // cross-correlation at every offset at once is what an FFT computes cheaply.
    // GC breaks this because a fragment's GC content depends on start *and*
    // length jointly.
    if gc_model.is_none() && (have_seq || pos_fw.is_some()) {
        let mut a = vec![1.0f64; ref_len];
        let mut b = vec![1.0f64; ref_len];
        if have_seq {
            a.copy_from_slice(&fw);
            b.copy_from_slice(&rc);
        }
        // Fold the positional factors into the same two arrays, so the FFT sees
        // one combined start array and one combined end array.
        if let (Some(pf), Some(pr)) = (pos_fw, pos_rc) {
            for (((ai, bi), &pfi), &pri) in
                a.iter_mut().zip(b.iter_mut()).zip(pf.iter()).zip(pr.iter())
            {
                *ai *= pfi;
                *bi *= pri;
            }
        }
        return crate::seqbias::eff_len_from_xcorr(
            &a,
            &b,
            cond,
            fld_low,
            fld_high,
            elen,
            unprocessed,
            stride.max(1),
            no_length_threshold,
        );
    }

    // Scalar path: sweep fragment lengths, and for each, every start position.
    let stride = stride.max(1) as i32;
    // No fragment can be longer than the transcript, nor than the FLD's support.
    let max_len = (ref_len as i32).min(fld_high as i32 + 1);
    let mut fl = fld_low as i32;
    let mut done = fl >= max_len;
    let sp = if fl > 0 { fl - 1 } else { 0 };
    // Running CDF value, so each length's probability is one subtraction rather
    // than a fresh lookup pair.
    let mut prev_mass = cond(sp);
    let mut eff = 0.0f64;
    while !done {
        // With a stride > 1 the last step can overshoot; clamp and finish.
        if fl >= max_len {
            done = true;
            fl = max_len - 1;
        }
        // Probability of exactly this fragment length = CDF difference.
        let fl_weight = cond(fl) - prev_mass;
        prev_mass = cond(fl);
        let mut mass = 0.0f64;
        // Hoist the bound: for kstart in [0, kmax) we have
        // frag_end = kstart+fl-1 <= ref_len-2 < ref_len, so the old per-iteration
        // `frag_end < ref_len` guard is always true and is dropped (it kept
        // `ref_len` live in the inner loop, forcing a spill/reload — the single
        // hottest source line in `perf annotate`).
        let kmax = ref_len as i32 - fl;
        // Dispatch the bias combination ONCE per fragment length rather than
        // re-testing every bias model's presence per fragment: the common
        // `--gcBias`-only case gets a tight loop with no per-fragment branches.
        match (have_seq, gc_model.zip(gc_ctx.as_ref()), pos_fw.zip(pos_rc)) {
            // Specialized GC-only loop: one lookup per fragment, no branches on
            // absent models.
            (false, Some((gc, ctx)), None) => {
                let mut kstart = 0i32;
                while kstart < kmax {
                    let frag_end = kstart + fl - 1;
                    if let Some((ff, cf)) = ctx.desc(kstart, frag_end) {
                        mass += gc.get(ff, cf);
                    } else {
                        // No usable GC context (too close to an end): neutral.
                        mass += 1.0;
                    }
                    kstart += pos_step;
                }
            }
            // General loop: multiply in whichever factors are present.
            _ => {
                let mut kstart = 0i32;
                while kstart < kmax {
                    let frag_start = kstart;
                    let frag_end = kstart + fl - 1;
                    let mut frag_factor = if have_seq {
                        fw[frag_start as usize] * rc[frag_end as usize]
                    } else {
                        1.0
                    };
                    if let (Some(gc), Some(ctx)) = (gc_model, gc_ctx.as_ref()) {
                        if let Some((ff, cf)) = ctx.desc(frag_start, frag_end) {
                            frag_factor *= gc.get(ff, cf);
                        }
                    }
                    if let (Some(pf), Some(pr)) = (pos_fw, pos_rc) {
                        frag_factor *= pf[frag_start as usize] * pr[frag_end as usize];
                    }
                    mass += frag_factor;
                    kstart += pos_step;
                }
            }
        }
        // Weight this length's total start-position mass by the length's own
        // probability, and accumulate — this is the convolution.
        eff += fl_weight * (mass * pos_scale);
        fl += stride;
    }

    if no_length_threshold {
        // salmon's `noThreshold` path: accept the bias-corrected length outright
        // (floored only at 1.0), else keep the uncorrected length. `unprocessed`
        // is already > 0 here (the early return handled `unprocessed <= 0`).
        if eff > 1.0 {
            eff
        } else {
            elen
        }
    } else {
        // The lower barrier: never shrink below `min(elen, unprocessed)`, so a
        // pathological correction cannot blow up the abundance that divides by it.
        let offset = (unprocessed as f64).max(1.0);
        eff.max(elen.min(offset))
    }
}

/// Build the *expected* positional-bias models (5'/3', one per length class),
/// mirroring salmon's expected-pos accumulation in `updateEffectiveLengths`:
/// for each expressed transcript and fragment start, add `log(weight·density)`
/// to the length-class bin (forward density = fragments that can start here,
/// reverse density = fragments that can end here). Models are finalized.
///
/// This is the "expected" half of observed-vs-expected: what the positional
/// histogram would look like with *no* positional bias at all, given only the
/// transcripts' lengths, their abundances and the fragment-length distribution.
/// Dividing the observed histogram by this isolates the bias itself.
#[allow(clippy::too_many_arguments)]
pub fn build_expected_pos<FL>(
    num_targets: usize,
    ref_len_of: FL,
    alphas: &[f64],
    eff_lens: &[f64],
    cdf: &[f64],
    quantiles: &[u32],
    k: usize,
) -> (Vec<SimplePosBias>, Vec<SimplePosBias>)
where
    FL: Fn(usize) -> usize + Sync,
{
    use rayon::prelude::*;
    type Partials = (Vec<SimplePosBias>, Vec<SimplePosBias>);
    // Per-transcript contributions are independent, each an O(refLen) sweep;
    // salmon parallelizes this expected-pos accumulation over transcripts and so
    // do we. Per-thread partials use `new_empty` (masses at -inf, the `log_add`
    // identity, carrying *no* pseudocount) so the fold/reduce merge via `combine`
    // is associative; the single `log(1)` pseudocount salmon seeds each bin with
    // is injected once at the end (combine into a `default()` model). `ref_len_of`
    // must be `Sync` to share across threads. `num_targets` excludes decoys (the
    // contiguous tail), which are never expressed and never contribute.
    fn empty() -> Partials {
        (
            (0..NUM_LENGTH_CLASSES)
                .map(|_| SimplePosBias::new_empty(NUM_POS_BINS))
                .collect(),
            (0..NUM_LENGTH_CLASSES)
                .map(|_| SimplePosBias::new_empty(NUM_POS_BINS))
                .collect(),
        )
    }
    let (sum5, sum3) = (0..num_targets)
        .into_par_iter()
        .fold(empty, |mut acc, tid| {
            // Skip transcripts that cannot contribute: unexpressed, degenerate,
            // too short for a context, or with no headroom.
            if alphas[tid] < MIN_ALPHA || eff_lens[tid] <= 0.0 {
                return acc;
            }
            let ref_len = ref_len_of(tid) as i32;
            if (ref_len as usize) <= k {
                return acc;
            }
            let unprocessed = ref_len - eff_lens[tid] as i32;
            if unprocessed <= 0 {
                return acc;
            }
            let cdf_max_arg = (cdf.len() - 1).min(ref_len as usize);
            let cdf_max_val = cdf[cdf_max_arg];
            if cdf_max_val < MIN_CDF_MASS {
                return acc;
            }
            // Fragments per position: abundance spread over the positions that can
            // host a fragment. This is what makes an abundant transcript dominate
            // the expected model, exactly as it dominates the observed one.
            let weight = alphas[tid] / eff_lens[tid];
            let lc = length_class_index(quantiles, ref_len as u32);
            let cond = |x: i32| conditional_cdf(cdf, cdf_max_arg, cdf_max_val, x);
            for frag_start in 0..(ref_len - k as i32) {
                // How much room is left to the right (for a fragment starting
                // here) and to the left (for one ending here); the CDF turns that
                // room into a probability.
                let max_fw = ref_len - frag_start + 1;
                let max_rc = frag_start;
                let density_fw = cond(max_fw);
                let density_rc = cond(max_rc);
                if weight * density_fw > EPSILON {
                    acc.0[lc].add_mass(frag_start, ref_len, weight * density_fw);
                }
                if weight * density_rc > EPSILON {
                    acc.1[lc].add_mass(frag_start, ref_len, weight * density_rc);
                }
            }
            acc
        })
        .reduce(empty, |mut a, b| {
            for (x, y) in a.0.iter_mut().zip(&b.0) {
                x.combine(y);
            }
            for (x, y) in a.1.iter_mut().zip(&b.1) {
                x.combine(y);
            }
            a
        });

    // Inject the single per-bin `log(1)` pseudocount salmon seeds each bin with
    // (`SimplePosBias::default`) by merging the raw observed log-sums into it,
    // then finalize. An empty bin (sum at -inf) collapses to exactly `log(1)`,
    // matching the serial accumulation.
    let seed = |sums: Vec<SimplePosBias>| -> Vec<SimplePosBias> {
        sums.into_iter()
            .map(|s| {
                let mut m = SimplePosBias::default();
                m.combine(&s);
                m.finalize();
                m
            })
            .collect()
    };
    (seed(sum5), seed(sum3))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A GC-biased ratio model, an FLD around 250 and a random sequence with
    /// GC-rich and GC-poor stretches: the ingredients of a scalar GC sweep.
    fn gc_fixture(len: usize) -> (Vec<u8>, Vec<u32>, GcFragModel, Vec<f64>, usize, usize) {
        gc_fixture_blocks(len, 3000)
    }

    fn gc_fixture_blocks(
        len: usize,
        block: usize,
    ) -> (Vec<u8>, Vec<u32>, GcFragModel, Vec<f64>, usize, usize) {
        use crate::gcbias::{gc_prefix, gc_ratio, DEFAULT_COND_BINS, DEFAULT_GC_BINS};
        let mut x = 0x1234_5678u64;
        let seq: Vec<u8> = (0..len)
            .map(|i| {
                x = x
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let r = (x >> 33) % 100;
                // alternate GC-rich and GC-poor stretches of `block` bases
                let gc_rich = (i / block).is_multiple_of(2);
                match (gc_rich, r) {
                    (true, 0..=34) | (false, 0..=14) => b'G',
                    (true, 35..=69) | (false, 15..=29) => b'C',
                    (_, r) if r % 2 == 0 => b'A',
                    _ => b'T',
                }
            })
            .collect();
        let prefix = gc_prefix(&seq);
        let mut obs = GcFragModel::new(DEFAULT_COND_BINS, DEFAULT_GC_BINS);
        let mut exp = GcFragModel::new(DEFAULT_COND_BINS, DEFAULT_GC_BINS);
        for ctx in 0..=100 {
            for gc in 0..=100 {
                obs.inc(gc, ctx, 1.0 + 0.8 * (gc as f64 / 100.0));
                exp.inc(gc, ctx, 1.0);
            }
        }
        let model = gc_ratio(&mut obs, &mut exp, 1000.0);
        let pmf: Vec<f64> = (0..=1000usize)
            .map(|l| (-0.5 * ((l as f64 - 250.0) / 40.0).powi(2)).exp())
            .collect();
        let (cdf, lo, hi) = crate::seqbias::fld_cdf_and_bounds(&pmf);
        (seq, prefix, model, cdf, lo, hi)
    }

    /// Stride 1 is the exact computation, bit for bit, so every existing call
    /// (and every default output) is unchanged.
    #[test]
    fn position_stride_one_is_exact() {
        let (seq, prefix, model, cdf, lo, hi) = gc_fixture(20_000);
        let bias = BiasInputs {
            seq: None,
            gc: Some((&model, GcView::Dense(&prefix))),
            pos: None,
        };
        let full = corrected_effective_length_full(&seq, &cdf, lo, hi, &bias, 19_750.0, 5, false);
        let one =
            corrected_effective_length_strided(&seq, &cdf, lo, hi, &bias, 19_750.0, 5, false, 1);
        assert_eq!(full.to_bits(), one.to_bits());
    }

    /// Sampling every 16th start and weighting by 16 estimates the full sweep
    /// closely, for a typical 5 kb intron target and a long one alike, with
    /// GC content switching every 300 bases (the scale of an Alu repeat).
    #[test]
    fn position_stride_estimates_intron_targets() {
        for len in [5_000usize, 100_000] {
            let (seq, prefix, model, cdf, lo, hi) = gc_fixture_blocks(len, 300);
            let bias = BiasInputs {
                seq: None,
                gc: Some((&model, GcView::Dense(&prefix))),
                pos: None,
            };
            let elen = len as f64 - 250.0;
            let full = corrected_effective_length_full(&seq, &cdf, lo, hi, &bias, elen, 5, false);
            let s = unspliced_position_stride(len);
            assert_eq!(s, UNSPLICED_BIAS_STRIDE);
            let est =
                corrected_effective_length_strided(&seq, &cdf, lo, hi, &bias, elen, 5, false, s);
            let rel = (est - full).abs() / full;
            assert!(
                rel < 2e-3,
                "len {len}: strided {est} vs exact {full} ({rel:e})"
            );
            // and the bias did move the length, so the test measures something
            assert!((full - elen).abs() / elen > 0.01, "{full} vs {elen}");
        }
    }

    /// The expected GC and sequence models: stride 1 is exact, a stride on a
    /// long target lands on the same normalized distribution. The abundance is
    /// set so each increment sits far above the accumulators' fixed-point
    /// resolution (`BIAS_WEIGHT_SCALE`, truncating): at a few reads per
    /// megabase the exact sweep's increments truncate towards zero, which a
    /// stride (larger, fewer increments) mitigates rather than worsens.
    #[test]
    fn expected_models_accept_a_position_stride() {
        use crate::gcbias::{build_expected_gc, build_expected_gc_strided, GcFragModel};
        use crate::seqbias::{build_expected, build_expected_strided};
        let (seq, prefix, _, cdf, lo, hi) = gc_fixture(300_000);
        let alphas = [1e7];
        let eff = [seq.len() as f64 - 250.0];
        let gc = |stride: Option<usize>| -> GcFragModel {
            let mut m = match stride {
                None => build_expected_gc(
                    1,
                    |_| &seq[..],
                    |_| GcView::Dense(&prefix),
                    &alphas,
                    &eff,
                    &cdf,
                    lo,
                    hi,
                    3,
                    25,
                    1,
                    5,
                ),
                Some(s) => build_expected_gc_strided(
                    1,
                    |_| &seq[..],
                    |_| GcView::Dense(&prefix),
                    |_| s,
                    &alphas,
                    &eff,
                    &cdf,
                    lo,
                    hi,
                    3,
                    25,
                    1,
                    5,
                ),
            };
            m.normalize();
            m
        };
        let (exact, one, sampled) = (gc(None), gc(Some(1)), gc(Some(6)));
        assert_eq!(exact.dump(), one.dump());
        let tv: f64 = exact
            .dump()
            .iter()
            .zip(sampled.dump())
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>();
        // three conditioning rows, each a distribution: < 1% of mass moved per row
        assert!(tv < 0.03, "total variation {tv}");

        let (efw, _) = build_expected(1, |_| &seq[..], &alphas, &eff, &cdf);
        let (ofw, _) = build_expected_strided(1, |_| &seq[..], |_| 1, &alphas, &eff, &cdf);
        assert_eq!(efw.dump(), ofw.dump());
    }

    #[test]
    fn unspliced_stride_keeps_short_targets_dense() {
        assert_eq!(unspliced_position_stride(0), 1);
        assert_eq!(unspliced_position_stride(100), 1);
        assert_eq!(unspliced_position_stride(640), 10);
        assert_eq!(unspliced_position_stride(5_000), UNSPLICED_BIAS_STRIDE);
        assert_eq!(unspliced_position_stride(2_000_000), UNSPLICED_BIAS_STRIDE);
    }
    use crate::posbias::compute_length_quantiles;
    use crate::seqbias::fld_cdf_and_bounds;

    type PosModels = (Vec<SimplePosBias>, Vec<SimplePosBias>);

    /// Total absolute difference between two sets of positional models, used as a
    /// scalar "did this change anything?" measure.
    fn mass_diff(a: &PosModels, b: &PosModels) -> f64 {
        a.0.iter()
            .chain(a.1.iter())
            .zip(b.0.iter().chain(b.1.iter()))
            .map(|(pa, pb)| {
                pa.masses()
                    .iter()
                    .zip(pb.masses())
                    .map(|(x, y)| (x - y).abs())
                    .sum::<f64>()
            })
            .sum()
    }

    /// Decoys sit past `num_targets` and must never enter the expected model — a
    /// 250 Mb genome decoy swept in as a "transcript" would swamp it. The
    /// zero-alpha half of the test confirms the other guard (unexpressed
    /// transcripts) works independently.
    #[test]
    fn build_expected_pos_respects_num_targets_bound() {
        // Five 200 nt transcripts plus a 400 nt "decoy". The decoy must change the
        // expected positional model only when `num_targets` includes it; a zeroed
        // alpha must skip it either way.
        let lens = [200usize, 200, 200, 200, 200, 400];
        let num_refs = lens.len();
        let alphas = vec![1.0; num_refs];
        let eff_lens = vec![150.0; num_refs];
        // A point-mass fragment-length distribution keeps the arithmetic simple.
        let mut pmf = vec![0.0; 200];
        pmf[100] = 1.0;
        let (cdf, _lo, _hi) = fld_cdf_and_bounds(&pmf);
        let qlens: Vec<u32> = lens.iter().map(|&l| l as u32).collect();
        let quantiles = compute_length_quantiles(&qlens, NUM_LENGTH_CLASSES);
        let k = 1usize;

        let exclude = build_expected_pos(5, |t| lens[t], &alphas, &eff_lens, &cdf, &quantiles, k);
        let include = build_expected_pos(6, |t| lens[t], &alphas, &eff_lens, &cdf, &quantiles, k);
        assert!(exclude
            .0
            .iter()
            .chain(exclude.1.iter())
            .all(|p| p.masses().iter().all(|v| v.is_finite())));
        let diff = mass_diff(&exclude, &include);
        assert!(
            diff > 1e-9,
            "a target beyond num_targets must not contribute (diff={diff})"
        );

        let mut alphas0 = alphas.clone();
        alphas0[5] = 0.0;
        let include0 = build_expected_pos(6, |t| lens[t], &alphas0, &eff_lens, &cdf, &quantiles, k);
        let diff2 = mass_diff(&exclude, &include0);
        assert!(
            diff2 < 1e-9,
            "zero-alpha target must not contribute (diff={diff2})"
        );
    }
}
