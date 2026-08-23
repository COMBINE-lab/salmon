//! Criterion benchmark for the EM/VBEM step over equivalence classes.
//!
//! Builds a synthetic power-law equivalence-class set (many singletons, a heavy
//! tail of multi-transcript classes drawn with a Zipf-like bias so a handful of
//! "hot" transcripts recur — the scatter-contention shape the real M-step guards
//! against) and times the convergence loop through the public `optimize_packed`
//! entry point. Every parallel case owns an explicit Rayon pool, so its benchmark
//! id records the worker count and Criterion compares like with like across
//! revisions instead of accidentally comparing two different `RAYON_NUM_THREADS`
//! environments under the same id.
//!
//! The one-step cases include `ShardPlan` construction and reusable workspace
//! allocation. The 50-step cases amortize that setup and are the primary M-step
//! throughput/scaling measurement. Full-convergence cases additionally capture
//! convergence-check and accelerator costs.

use std::hint::black_box;

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
    const THREADS: [usize; 7] = [1, 2, 4, 8, 16, 32, 64];
    const NUM_TXPS: u32 = 100_000;
    const NUM_CLASSES: usize = 300_000;

    let p = synthetic(NUM_TXPS, NUM_CLASSES, 0xDEAD_BEEF);
    let shape = format!("{NUM_TXPS}txp_{NUM_CLASSES}cls");
    let step_em = single_step(false);
    let step_vbem = single_step(true);
    let loop_em = bounded(false, 50);
    let loop_vbem = bounded(true, 50);

    let mut group = c.benchmark_group("em_scaling");
    // Seven pool sizes make a full run intentionally substantial; ten samples is
    // Criterion's minimum and enough for the paired baseline/candidate comparison.
    group.sample_size(10);

    // The sequential path is what each outer-parallel bootstrap replicate uses.
    group.throughput(Throughput::Elements(NUM_CLASSES as u64 * 50));
    group.bench_with_input(BenchmarkId::new("loop50_vbem_seq", &shape), &p, |b, p| {
        b.iter(|| black_box(optimize_packed(p, &loop_vbem, false)))
    });

    for threads in THREADS {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("benchmark Rayon pool");
        let suffix = format!("t{threads}_{shape}");

        group.throughput(Throughput::Elements(NUM_CLASSES as u64));
        group.bench_with_input(BenchmarkId::new("step_em_par", &suffix), &p, |b, p| {
            b.iter(|| pool.install(|| black_box(optimize_packed(p, &step_em, true))))
        });
        group.bench_with_input(BenchmarkId::new("step_vbem_par", &suffix), &p, |b, p| {
            b.iter(|| pool.install(|| black_box(optimize_packed(p, &step_vbem, true))))
        });

        group.throughput(Throughput::Elements(NUM_CLASSES as u64 * 50));
        group.bench_with_input(BenchmarkId::new("loop50_em_par", &suffix), &p, |b, p| {
            b.iter(|| pool.install(|| black_box(optimize_packed(p, &loop_em, true))))
        });
        group.bench_with_input(BenchmarkId::new("loop50_vbem_par", &suffix), &p, |b, p| {
            b.iter(|| pool.install(|| black_box(optimize_packed(p, &loop_vbem, true))))
        });
    }
    group.finish();

    // Full convergence is intentionally a separate, smaller fixture so it can be
    // filtered independently from the scaling sweep.
    let conv_p = synthetic(20_000, 50_000, 0xC0FF_EE11);
    let conv_plain = EmOptions {
        use_vbem: true,
        ..Default::default()
    };
    let conv_sq = EmOptions {
        use_vbem: true,
        accel: EmAccel::Squarem,
        ..Default::default()
    };
    let conv_da = EmOptions {
        use_vbem: true,
        accel: EmAccel::Daarem,
        ..Default::default()
    };
    let mut convergence = c.benchmark_group("em_convergence");
    convergence.sample_size(10);
    for threads in THREADS {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("benchmark Rayon pool");
        for (name, opts) in [
            ("vbem_none", &conv_plain),
            ("vbem_squarem", &conv_sq),
            ("vbem_daarem", &conv_da),
        ] {
            convergence.bench_with_input(
                BenchmarkId::new(name, format!("t{threads}_20000txp_50000cls")),
                &conv_p,
                |b, p| b.iter(|| pool.install(|| black_box(optimize_packed(p, opts, true)))),
            );
        }
    }
    convergence.finish();
}

criterion_group!(benches, bench_em);
criterion_main!(benches);
