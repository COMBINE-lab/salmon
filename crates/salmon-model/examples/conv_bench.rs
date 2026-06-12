//! Isolated microbenchmark for the GC bias-corrected effective-length
//! convolution (`corrected_effective_length_full`), so it can be profiled and
//! A/B'd without the surrounding 130s quant run (where it's only ~5% of wall).
//!
//! Generates a realistic set of transcripts (human-cDNA-like length spread,
//! random ACGT), a 3×25 GC ratio model, and a fragment-length CDF, then times
//! the convolution over every transcript. Run:
//!   cargo run --release --example conv_bench -- [num_transcripts] [passes]

use std::hint::black_box;
use std::time::Instant;

use salmon_model::gcbias::{DEFAULT_COND_BINS, DEFAULT_GC_BINS};
use salmon_model::{
    corrected_effective_length_full, gc_prefix, gc_ratio, BiasInputs, GcFragModel, GC_SAMP_STRIDE,
};

/// Tiny deterministic LCG so the bench is reproducible without a rand dep.
struct Lcg(u64);
impl Lcg {
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 33) as u32
    }
    fn next_f64(&mut self) -> f64 {
        self.next_u32() as f64 / u32::MAX as f64
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(100_000);
    let passes: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(3);

    let mut rng = Lcg(0x9E3779B97F4A7C15);

    // Build a realistic transcript-length distribution: log-normal-ish around
    // ~1.5 kb with a long tail, clamped to [80, 15000] — close to GRCh38 cDNA.
    let bases = [b'A', b'C', b'G', b'T'];
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(n);
    let mut prefixes: Vec<Vec<u32>> = Vec::with_capacity(n);
    for _ in 0..n {
        let u = rng.next_f64();
        // exp of a normal-ish variate -> heavy right tail
        let len = (7.3 + 0.9 * (u - 0.5) * 4.0).exp() as usize;
        let len = len.clamp(80, 15_000);
        let seq: Vec<u8> = (0..len)
            .map(|_| bases[(rng.next_u32() & 3) as usize])
            .collect();
        prefixes.push(gc_prefix(&seq));
        seqs.push(seq);
    }

    // A normalized 3×25 GC ratio model (the shape doesn't matter for timing; we
    // exercise the same binning + lookup the real run does).
    let mut obs = GcFragModel::new(DEFAULT_COND_BINS, DEFAULT_GC_BINS);
    let mut exp = GcFragModel::new(DEFAULT_COND_BINS, DEFAULT_GC_BINS);
    for ctx in 0..=100 {
        for gc in 0..=100 {
            obs.inc(gc, ctx, 1.0 + 0.3 * ((gc + ctx) as f64).sin());
            exp.inc(gc, ctx, 1.0);
        }
    }
    let gc_model = gc_ratio(&mut obs, &mut exp, 1000.0);

    // Fragment-length CDF: gaussian-ish pmf around mean 250, sd 40, over 0..1000.
    let fld_max = 1000usize;
    let mean = 250.0f64;
    let sd = 40.0f64;
    let pmf: Vec<f64> = (0..=fld_max)
        .map(|l| {
            let z = (l as f64 - mean) / sd;
            (-0.5 * z * z).exp()
        })
        .collect();
    let (cdf, fld_low, fld_high) = salmon_model::seqbias::fld_cdf_and_bounds(&pmf);

    let mut acc = 0.0f64;
    let mut best = f64::INFINITY;
    for p in 0..passes {
        let t = Instant::now();
        for (seq, prefix) in seqs.iter().zip(&prefixes) {
            let ref_len = seq.len() as f64;
            let elen = (ref_len - 200.0).max(1.0); // ensure unprocessed > 0
            let bias = BiasInputs {
                seq: None,
                gc: Some((&gc_model, salmon_model::GcView::Dense(prefix.as_slice()))),
                pos: None,
            };
            acc += corrected_effective_length_full(
                seq,
                &cdf,
                fld_low,
                fld_high,
                &bias,
                elen,
                GC_SAMP_STRIDE,
                false,
            );
        }
        let dt = t.elapsed().as_secs_f64();
        best = best.min(dt);
        eprintln!(
            "pass {p}: {:.3}s  ({:.2} µs/transcript)  acc={}",
            dt,
            dt * 1e6 / n as f64,
            black_box(acc)
        );
    }
    eprintln!(
        "BEST: {:.3}s over {} transcripts ({:.2} µs/transcript)",
        best,
        n,
        best * 1e6 / n as f64
    );
}
