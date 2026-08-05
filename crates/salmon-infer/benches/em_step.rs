//! Criterion benchmark for the EM/VBEM step over equivalence classes.
//!
//! Builds a synthetic power-law equivalence-class set (many singletons, a heavy
//! tail of multi-transcript classes drawn with a Zipf-like bias so a handful of
//! "hot" transcripts recur — the scatter-contention shape the real M-step guards
//! against) and times the convergence loop through the public `optimize_packed`
//! entry point:
//!   - a single M-step (`min_iter = max_iter = 1`, no truncation): the surgical
//!     per-iteration throughput measurement,
//!   - a bounded convergence run (fixed `max_iter`): EM vs VBEM, parallel (rayon)
//!     vs sequential.
//!
//! This is the guard that will quantify the SQUAREM win in Phase 2: re-run with an
//! acceleration flag and compare wall-clock here, and iterations-to-convergence via
//! `EmResult::iters` in the `salmon-infer` tests, on this same fixture.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use salmon_eqclass::{EquivalenceClassBuilder, TranscriptGroup};
use salmon_infer::{optimize_packed, EmAccel, EmOptions, PackedEqClasses};

/// Deterministic splitmix64 PRNG — no external rng dependency, fully reproducible.
///
/// A benchmark must measure the same work every run, so the fixture cannot come
/// from a real random source. splitmix64 is a few lines of arithmetic with good
/// enough statistical properties for generating a plausible workload.
struct Rng(u64);
impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: u32) -> u32 {
        (self.next_u64() % n as u64) as u32
    }
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Zipf-ish transcript sampler: square a uniform to skew draws toward low ids, so a
/// small set of transcripts recurs across many classes.
///
/// This shape is the point of the benchmark. Real expression is highly skewed —
/// a few transcripts dominate — so a handful of "hot" transcripts appear in a
/// large share of classes. That is exactly the pattern that makes a shared
/// accumulator contend, and what the sharded M-step exists to handle. A uniform
/// fixture would measure a workload salmon never sees.
fn pick(rng: &mut Rng, num_txps: u32) -> u32 {
    let u = rng.unit();
    (((u * u) * num_txps as f64) as u32) % num_txps
}

/// Synthetic packed eq-class set: `num_txps` transcripts, `num_classes` classes.
///
/// Built through the ordinary public builder rather than by constructing the
/// packed form directly, so the fixture is guaranteed to be one the real pipeline
/// could produce.
///
/// ~70% singletons; the rest span 2..=8 transcripts drawn with the hot-transcript
/// bias above. Effective lengths spread over a realistic range so the
/// multiplicative update does real work.
fn synthetic(num_txps: u32, num_classes: usize, seed: u64) -> PackedEqClasses {
    let mut rng = Rng(seed);
    let b = EquivalenceClassBuilder::new();
    for _ in 0..num_classes {
        let k = if rng.unit() < 0.70 {
            1
        } else {
            2 + rng.below(7) as usize // 2..=8
        };
        let mut txps = Vec::with_capacity(k);
        while txps.len() < k {
            let t = pick(&mut rng, num_txps);
            if !txps.contains(&t) {
                txps.push(t);
            }
        }
        let weights = vec![1.0; txps.len()];
        let count = 1 + rng.below(50) as u64;
        b.add_group(TranscriptGroup::new(txps), weights, count);
    }
    let mut eq = b.finish();
    let eff: Vec<f64> = (0..num_txps).map(|_| 200.0 + rng.unit() * 3000.0).collect();
    eq.update_eff_lengths(&eff);
    PackedEqClasses::from_collapsed(&eq, num_txps as usize)
}

/// One M-step, no post-convergence truncation — isolates the hot kernel cost.
fn single_step(use_vbem: bool) -> EmOptions {
    EmOptions {
        min_iter: 1,
        max_iter: 1,
        min_alpha: 0.0,
        use_vbem,
        ..Default::default()
    }
}

/// A fixed, bounded number of iterations so wall-clock is comparable across
/// EM/VBEM/par/seq (and, in Phase 2, plain vs SQUAREM) rather than dominated by a
/// data-dependent convergence count.
fn bounded(use_vbem: bool, iters: u32) -> EmOptions {
    EmOptions {
        min_iter: iters,
        max_iter: iters,
        use_vbem,
        ..Default::default()
    }
}

fn bench_em(c: &mut Criterion) {
    let sizes = [(20_000u32, 50_000usize), (100_000u32, 300_000usize)];
    let mut group = c.benchmark_group("em");
    // Convergence loops are the expensive case; keep the sample count modest so the
    // whole bench finishes in well under a minute.
    group.sample_size(20);

    for &(num_txps, num_classes) in &sizes {
        let p = synthetic(num_txps, num_classes, 0xDEAD_BEEF);
        group.throughput(Throughput::Elements(num_classes as u64));
        let id = format!("{num_txps}txp_{num_classes}cls");

        // --- single M-step (per-iteration throughput) ---
        let step_em = single_step(false);
        group.bench_with_input(BenchmarkId::new("step_em_par", &id), &p, |b, p| {
            b.iter(|| optimize_packed(p, &step_em, true))
        });
        group.bench_with_input(BenchmarkId::new("step_em_seq", &id), &p, |b, p| {
            b.iter(|| optimize_packed(p, &step_em, false))
        });
        let step_vbem = single_step(true);
        group.bench_with_input(BenchmarkId::new("step_vbem_par", &id), &p, |b, p| {
            b.iter(|| optimize_packed(p, &step_vbem, true))
        });

        // --- bounded convergence loop (EM vs VBEM, par) ---
        let loop_em = bounded(false, 50);
        group.bench_with_input(BenchmarkId::new("loop50_em_par", &id), &p, |b, p| {
            b.iter(|| optimize_packed(p, &loop_em, true))
        });
        let loop_vbem = bounded(true, 50);
        group.bench_with_input(BenchmarkId::new("loop50_vbem_par", &id), &p, |b, p| {
            b.iter(|| optimize_packed(p, &loop_vbem, true))
        });

        // --- full convergence: plain EM vs SQUAREM (the headline comparison) ---
        let conv_plain = EmOptions::default();
        group.bench_with_input(BenchmarkId::new("converge_em_plain", &id), &p, |b, p| {
            b.iter(|| optimize_packed(p, &conv_plain, true))
        });
        let conv_sq = EmOptions {
            accel: EmAccel::Squarem,
            ..Default::default()
        };
        group.bench_with_input(BenchmarkId::new("converge_em_squarem", &id), &p, |b, p| {
            b.iter(|| optimize_packed(p, &conv_sq, true))
        });
        let conv_da = EmOptions {
            accel: EmAccel::Daarem,
            ..Default::default()
        };
        group.bench_with_input(BenchmarkId::new("converge_em_daarem", &id), &p, |b, p| {
            b.iter(|| optimize_packed(p, &conv_da, true))
        });
    }
    group.finish();
}

criterion_group!(benches, bench_em);
criterion_main!(benches);
