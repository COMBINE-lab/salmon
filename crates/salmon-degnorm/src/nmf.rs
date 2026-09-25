//! Rank-one non-negative over-approximation of a transcript's coverage matrix,
//! and the degradation index derived from it.
//!
//! # The model
//!
//! Stack one transcript's coverage curves for `m` samples into a matrix `K`
//! (`m × n`, one column per position bin). If nothing had degraded, every
//! sample would show the *same shape* scaled by *that sample's abundance* —
//! that is, `K` would be exactly rank one: `K = a bᵀ`, with `a` the per-sample
//! abundance and `b` the shape ("envelope") the protocol produces for this
//! transcript.
//!
//! Degradation only ever removes coverage. So the observed curve lies *under*
//! the intact one, and the right fit is not the closest rank-one matrix but the
//! closest one that stays above the data:
//!
//! ```text
//!   minimise ‖K̃ − a bᵀ‖_F   subject to   K̃ ≥ K,  a ≥ 0,  b ≥ 0
//! ```
//!
//! The degradation index of sample `i` is then how much of the envelope's area
//! is missing from that sample's actual curve:
//!
//! ```text
//!   DIᵢ = 1 − (Σ_l K_il) / (Σ_l aᵢ b_l)
//! ```
//!
//! `DI = 0` is an intact transcript; `DI = 0.5` means half of the fragments the
//! transcript should have produced are gone, so its raw count understates it by
//! a factor of two. That is exactly the correction applied downstream:
//! `adjusted = raw / (1 − DI)`.
//!
//! This is the model behind DegNorm (Xiong *et al.*, *Genome Biology* 2019,
//! <https://nustatbioinfo.github.io/DegNorm/>). What follows implements that
//! model over salmon's own binned transcript coverage; it is not a port of the
//! DegNorm code, and the two differ in the ways `docs/degnorm.md` sets out.
//!
//! # Why it needs several samples
//!
//! With one sample the problem is degenerate: any curve is "rank one" on its
//! own, so `a bᵀ` fits it exactly and every DI is zero. Degradation is only
//! identifiable because *other* samples show the parts of the transcript this
//! one lost. This is why coverage is dumped during `salmon quant` but the fit
//! happens later, across a cohort.
//!
//! # How it is solved
//!
//! Block coordinate descent on `‖K̃ − a bᵀ‖_F`, which is what makes the
//! constrained problem cheap:
//!
//! 1. Hold `K̃` and fit the best rank-one approximation by alternating least
//!    squares (`b ← K̃ᵀa / aᵀa`, `a ← K̃b / bᵀb`). Because `K̃` is non-negative
//!    so are the iterates, which is why no non-negativity projection appears.
//! 2. Hold `a, b` and project back onto the constraint set: `K̃ ← max(K, a bᵀ)`,
//!    the nearest matrix that dominates the data.
//!
//! Step 2 is what lets the fit see through degradation: the entries damage
//! removed are replaced by the model's own prediction, so step 1 is not dragged
//! down by the holes.
//!
//! # Excluding damaged positions from the fit
//!
//! Step 2 alone is not enough. Take three samples of a transcript, one of which
//! lost its first half. Over the second half three curves agree; over the first
//! half only two carry evidence and the third contributes whatever the current
//! model predicts. That is self-consistent at any envelope height, so the
//! iteration settles on an envelope that sags over the damaged half — and a
//! sagging envelope cannot explain the *intact* samples' flat curves without
//! slack, which surfaces as a non-zero index for samples that lost nothing
//! (~0.10 in that example, against a true 0.5 for the truncated one).
//!
//! The fix is to fit the envelope only where there is undamaged evidence for
//! it: an entry whose observed coverage has fallen below
//! [`FitOptions::mask_frac`] of the current fit is treated as damaged and left
//! out of the least squares — never out of the index, which is computed from
//! the original data. A column or row that would lose every entry keeps them
//! all, so the fit never runs out of data. On the example above this recovers
//! `0, 0.5, 0` exactly. DegNorm reaches for the same idea with its
//! baseline-selection step; this is a per-entry version of it, not a
//! reproduction of that procedure.
//!
//! # What the index cannot tell you
//!
//! Two limits are worth stating, because both are properties of the model
//! rather than of this implementation:
//!
//! * **Uniform loss is invisible.** A sample that lost the same fraction of
//!   coverage at *every* position is indistinguishable from a sample with less
//!   of that transcript. Only *uneven* loss is identifiable — which is what
//!   degradation actually looks like.
//! * **The index has a positive floor.** Coverage curves are noisy, and an
//!   envelope constrained to sit above the data sits above the noise too. On
//!   simulated intact libraries the floor runs ~0.1 at 20× and ~0.25 at 5×
//!   coverage, so DI is a comparative measure — between samples, for a given
//!   transcript — and not an absolute percentage of degraded molecules. This is
//!   why [`crate::cohort::CohortOptions::min_mean_coverage`] exists.

/// Tuning for [`fit`].
#[derive(Debug, Clone, Copy)]
pub struct FitOptions {
    /// Maximum outer (fit / project) iterations.
    pub max_iter: usize,
    /// Stop when the relative Frobenius change of `K̃` falls below this.
    pub tol: f64,
    /// Inner alternating-least-squares sweeps per outer iteration. The factors
    /// are warm-started from the previous outer iteration, so a handful is
    /// enough after the first one.
    pub als_iter: usize,
    /// Entries whose observed coverage is below this fraction of the fitted
    /// envelope are held out of the envelope fit as damaged; see the module
    /// docs. `0.0` fits every entry, which biases the index on exactly the
    /// transcripts the model exists to catch.
    pub mask_frac: f64,
}

impl Default for FitOptions {
    fn default() -> Self {
        Self {
            max_iter: 100,
            tol: 1e-4,
            als_iter: 8,
            mask_frac: 0.5,
        }
    }
}

/// Ceiling on a reported DI.
///
/// A DI of exactly 1 says "every fragment this transcript should have produced
/// is missing", which the count adjustment `1 / (1 − DI)` cannot represent and
/// which no finite amount of evidence supports anyway. Capping keeps a
/// near-empty row from becoming an astronomical adjusted count.
pub const MAX_DI: f64 = 0.99;

/// The fitted rank-one model for one transcript.
#[derive(Debug, Clone, PartialEq)]
pub struct Fit {
    /// Per-sample abundance scale `a` (length `m`).
    pub scales: Vec<f64>,
    /// Shared coverage envelope `b` (length `n`).
    pub envelope: Vec<f64>,
    /// Per-sample degradation index, in `[0, MAX_DI]`.
    pub di: Vec<f64>,
    /// Outer iterations actually run (`0` for a degenerate all-zero matrix).
    pub iters: usize,
}

/// Fit the rank-one over-approximation of `k` (`m × n`, row-major) and derive
/// the per-sample degradation index.
///
/// A degenerate matrix — all zero, or a single sample — has no identifiable
/// degradation and comes back with `DI = 0` rather than an error: callers fit
/// hundreds of thousands of transcripts and some of them are always degenerate.
pub fn fit(k: &[f64], m: usize, n: usize, opts: &FitOptions) -> Fit {
    assert_eq!(k.len(), m * n, "coverage matrix shape mismatch");
    if m == 0 || n == 0 || k.iter().sum::<f64>() <= 0.0 {
        return Fit {
            scales: vec![0.0; m],
            envelope: vec![0.0; n],
            di: vec![0.0; m],
            iters: 0,
        };
    }

    // `K̃` starts at the data itself, which trivially satisfies `K̃ ≥ K`, and
    // `a` at the row sums so the first sweep does not start from an arbitrary
    // vector.
    let mut kt = k.to_vec();
    let mut a: Vec<f64> = (0..m).map(|i| k[i * n..(i + 1) * n].iter().sum()).collect();
    let mut b = vec![0.0f64; n];
    let mut mask = vec![true; m * n];
    let mut iters = 0;

    for _ in 0..opts.max_iter {
        iters += 1;
        alternating_least_squares(k, &kt, &mut mask, m, n, &mut a, &mut b, opts);
        // Project onto `K̃ ≥ K`, tracking how far the projection moved.
        let mut delta = 0.0;
        let mut norm = 0.0;
        for i in 0..m {
            for j in 0..n {
                let idx = i * n + j;
                let next = (a[i] * b[j]).max(k[idx]);
                let d = next - kt[idx];
                delta += d * d;
                norm += next * next;
                kt[idx] = next;
            }
        }
        if norm <= 0.0 || (delta / norm).sqrt() < opts.tol {
            break;
        }
    }

    // The index compares each sample's *observed* area with the envelope's, so
    // the masked-out entries are counted here — they are the missing coverage
    // it is measuring.
    let envelope_area: f64 = b.iter().sum();
    let di = (0..m)
        .map(|i| {
            let observed: f64 = k[i * n..(i + 1) * n].iter().sum();
            let expected = a[i] * envelope_area;
            if expected <= 0.0 {
                0.0
            } else {
                (1.0 - observed / expected).clamp(0.0, MAX_DI)
            }
        })
        .collect();

    Fit {
        scales: a,
        envelope: b,
        di,
        iters,
    }
}

/// Flag the entries that have lost more than `1 − frac` of their fitted
/// coverage, so the envelope is fitted to what survived rather than to the
/// damage.
///
/// A column with no surviving entry (every sample lost that position) or a row
/// with none (a sample that lost everything) keeps all of its entries: an empty
/// least-squares system has no solution, and "no evidence anywhere" is not a
/// reason to discard what little there is.
fn update_mask(k: &[f64], m: usize, n: usize, a: &[f64], b: &[f64], frac: f64, mask: &mut [bool]) {
    if frac <= 0.0 || b.iter().all(|&x| x <= 0.0) {
        mask.fill(true);
        return;
    }
    for i in 0..m {
        for j in 0..n {
            mask[i * n + j] = k[i * n + j] >= frac * a[i] * b[j];
        }
    }
    for j in 0..n {
        if (0..m).all(|i| !mask[i * n + j]) {
            for i in 0..m {
                mask[i * n + j] = true;
            }
        }
    }
    for i in 0..m {
        let row = &mut mask[i * n..(i + 1) * n];
        if row.iter().all(|&x| !x) {
            row.fill(true);
        }
    }
}

/// A few alternating-least-squares sweeps towards the best rank-one
/// approximation of the unmasked entries of `kt`, updating `a`, `b` and the
/// damaged-entry `mask` in place. `k` is the original data the mask is derived
/// from; `kt` the over-approximating completion the fit runs on.
fn alternating_least_squares(
    k: &[f64],
    kt: &[f64],
    mask: &mut [bool],
    m: usize,
    n: usize,
    a: &mut [f64],
    b: &mut [f64],
    opts: &FitOptions,
) {
    for _ in 0..opts.als_iter {
        // Re-derived every sweep, not once per outer iteration: the mask is a
        // function of the current fit, and letting it lag leaves damaged
        // entries in the least squares long enough to bias the envelope.
        update_mask(k, m, n, a, b, opts.mask_frac, mask);
        // b_j ← Σ_i a_i K_ij / Σ_i a_i², over the unmasked entries of column j.
        for j in 0..n {
            let mut num = 0.0;
            let mut den = 0.0;
            for i in 0..m {
                if mask[i * n + j] {
                    num += a[i] * kt[i * n + j];
                    den += a[i] * a[i];
                }
            }
            if den > 0.0 {
                b[j] = (num / den).max(0.0);
            }
        }
        // a_i ← Σ_j b_j K_ij / Σ_j b_j², over the unmasked entries of row i.
        for i in 0..m {
            let mut num = 0.0;
            let mut den = 0.0;
            for j in 0..n {
                if mask[i * n + j] {
                    num += b[j] * kt[i * n + j];
                    den += b[j] * b[j];
                }
            }
            if den > 0.0 {
                a[i] = (num / den).max(0.0);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build an `m × n` matrix from a closure over (sample, position).
    fn matrix(m: usize, n: usize, f: impl Fn(usize, usize) -> f64) -> Vec<f64> {
        (0..m)
            .flat_map(|i| (0..n).map(move |j| (i, j)))
            .map(|(i, j)| f(i, j))
            .collect()
    }

    #[test]
    fn an_exactly_rank_one_matrix_has_no_degradation() {
        // Three samples at different abundances, identical shape.
        let shape = [1.0, 3.0, 5.0, 4.0, 2.0, 1.0, 2.0, 3.0];
        let scales = [1.0, 2.5, 0.4];
        let k = matrix(3, 8, |i, j| scales[i] * shape[j]);
        let fit = fit(&k, 3, 8, &FitOptions::default());
        for (i, di) in fit.di.iter().enumerate() {
            assert!(*di < 1e-3, "sample {i} should be intact, got DI={di}");
        }
    }

    #[test]
    fn a_truncated_sample_loses_exactly_the_area_it_is_missing() {
        // Sample 1 lost its first half: half the flat envelope's area is gone.
        let n = 20;
        let k = matrix(3, n, |i, j| if i == 1 && j < n / 2 { 0.0 } else { 10.0 });
        let fit = fit(&k, 3, n, &FitOptions::default());
        assert!(fit.di[0] < 0.01, "intact sample: {}", fit.di[0]);
        assert!(fit.di[2] < 0.01, "intact sample: {}", fit.di[2]);
        assert!(
            (fit.di[1] - 0.5).abs() < 0.01,
            "truncated sample should read 0.5, got {}",
            fit.di[1]
        );
    }

    #[test]
    fn samples_damaged_in_different_places_are_all_seen() {
        // Each sample lost a different 40% stretch; every position survives in
        // two of the three, so the envelope is fully determined.
        let n = 40;
        let holes = [(0, 16), (12, 28), (24, 40)];
        let k = matrix(3, n, |i, j| {
            let (lo, hi) = holes[i];
            if j >= lo && j < hi {
                0.0
            } else {
                100.0
            }
        });
        let fit = fit(&k, 3, n, &FitOptions::default());
        for (i, di) in fit.di.iter().enumerate() {
            assert!(
                (di - 0.4).abs() < 0.02,
                "sample {i}: expected ~0.4, got {di}"
            );
        }
    }

    #[test]
    fn more_degradation_reads_as_a_higher_index() {
        let n = 40;
        // A 3'-biased decay of increasing severity: sample i keeps
        // `exp(-rate_i * (1 - x))` of its coverage.
        let rates = [0.0, 1.0, 3.0];
        let k = matrix(3, n, |i, j| {
            let x = j as f64 / (n - 1) as f64;
            100.0 * (-rates[i] * (1.0 - x)).exp()
        });
        let fit = fit(&k, 3, n, &FitOptions::default());
        assert!(
            fit.di[0] < fit.di[1] && fit.di[1] < fit.di[2],
            "DI must increase with degradation, got {:?}",
            fit.di
        );
        assert!(fit.di[0] < 0.1, "intact sample drifted to {}", fit.di[0]);
        assert!(fit.di[2] > 0.4, "heavy decay under-read as {}", fit.di[2]);
    }

    #[test]
    fn masking_is_what_keeps_intact_samples_at_zero() {
        // The same truncation, fitted without the damaged-entry mask: the
        // envelope sags over the missing half and the intact samples pick up an
        // index they have not earned. This is the failure the mask exists for,
        // pinned here so it cannot come back unnoticed.
        let n = 20;
        let k = matrix(3, n, |i, j| if i == 1 && j < n / 2 { 0.0 } else { 10.0 });
        let unmasked = fit(
            &k,
            3,
            n,
            &FitOptions {
                mask_frac: 0.0,
                ..Default::default()
            },
        );
        assert!(
            unmasked.di[0] > 0.05,
            "expected the unmasked fit to be biased, got {}",
            unmasked.di[0]
        );
        let masked = fit(&k, 3, n, &FitOptions::default());
        assert!(masked.di[0] < unmasked.di[0]);
    }

    #[test]
    fn uniform_loss_is_reported_as_lower_abundance_not_degradation() {
        // Sample 1 lost half of its coverage *everywhere*. Nothing in the data
        // distinguishes that from simply having less of the transcript, and the
        // model must not pretend otherwise.
        let n = 20;
        let k = matrix(3, n, |i, j| {
            let base = 10.0 + (j % 4) as f64;
            if i == 1 {
                base * 0.5
            } else {
                base
            }
        });
        let fit = fit(&k, 3, n, &FitOptions::default());
        assert!(
            fit.di.iter().all(|&d| d < 0.01),
            "uniform loss must not read as degradation, got {:?}",
            fit.di
        );
    }

    #[test]
    fn the_index_stays_inside_its_range() {
        let n = 16;
        let k = matrix(4, n, |i, j| ((i + 1) * (j % 5 + 1)) as f64);
        let fit = fit(&k, 4, n, &FitOptions::default());
        assert!(fit.di.iter().all(|&d| (0.0..=MAX_DI).contains(&d)));
    }

    #[test]
    fn a_degenerate_matrix_returns_zeros_rather_than_nan() {
        let fit = fit(&[0.0; 12], 3, 4, &FitOptions::default());
        assert_eq!(fit.di, vec![0.0; 3]);
        assert_eq!(fit.iters, 0);
        assert!(fit.envelope.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn a_single_sample_cannot_identify_degradation() {
        let k = [5.0, 4.0, 0.0, 0.0, 1.0];
        let fit = fit(&k, 1, 5, &FitOptions::default());
        assert!(fit.di[0] < 1e-6, "got {}", fit.di[0]);
    }

    #[test]
    fn the_index_is_capped_below_one() {
        // One surviving bin against two intact samples.
        let n = 50;
        let k = matrix(3, n, |i, j| if i == 0 && j > 0 { 0.0 } else { 50.0 });
        let fit = fit(&k, 3, n, &FitOptions::default());
        assert!(fit.di[0] <= MAX_DI && fit.di[0] > 0.9, "got {}", fit.di[0]);
    }

    #[test]
    fn a_sample_that_lost_everything_does_not_break_the_fit() {
        let n = 12;
        let k = matrix(3, n, |i, _| if i == 2 { 0.0 } else { 7.0 });
        let fit = fit(&k, 3, n, &FitOptions::default());
        assert!(fit.di.iter().all(|d| d.is_finite()));
        assert!(fit.di[0] < 0.01 && fit.di[1] < 0.01);
    }
}
