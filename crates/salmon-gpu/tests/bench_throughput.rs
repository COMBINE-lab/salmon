//! Rough throughput probe (not a determinism test): time a realistic mini-batch
//! of full-length alignments on the real CPU backend (ksw2rs, SIMD) vs. the GPU.
//! Run with `cargo test --release -p salmon-gpu --features gpu --test
//! bench_throughput -- --ignored --nocapture`. The CPU baseline is single-
//! threaded, so it understates salmon's real CPU path (all cores); read it as
//! raw alignment throughput, GPU (one device) vs one ksw2rs thread.
#![cfg(feature = "gpu")]

use std::time::Instant;

use ksw2rs::{Aligner, Extz2Input, KSW_EZ_RIGHT, KSW_EZ_SCORE_ONLY, KSW_NEG_INF};
use salmon_gpu::gpu::{GpuContext, GpuTask};
use salmon_gpu::reference::BandedParams;

const MATCH: i8 = 2;
const MISMATCH: i8 = 4;
const GAPO: i8 = 6;
const GAPE: i8 = 2;
const W: i32 = 15;

fn ksw_score(aln: &mut Aligner, q: &[u8], t: &[u8], mat: &[i8; 25]) -> i32 {
    let input = Extz2Input {
        query: q,
        target: t,
        m: 5,
        mat,
        q: GAPO,
        e: GAPE,
        w: W,
        zdrop: -1,
        end_bonus: 0,
        flag: KSW_EZ_RIGHT | KSW_EZ_SCORE_ONLY,
    };
    let ez = aln.align(&input);
    if ez.mqe > KSW_NEG_INF {
        ez.mqe
    } else {
        ez.max as i32
    }
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
        self.0 >> 33
    }
    fn b(&mut self) -> u8 {
        (self.next() % 4) as u8
    }
}

#[test]
#[ignore]
fn throughput_cpu_vs_gpu() {
    let Some(ctx) = GpuContext::new() else {
        eprintln!("no GPU adapter; skipping");
        return;
    };
    let p = BandedParams {
        match_score: MATCH as i32,
        mismatch_pen: MISMATCH as i32,
        gap_open_pen: GAPO as i32,
        gap_extend_pen: GAPE as i32,
        bandwidth: W,
    };
    let n = 50_000usize;
    let read_len = 150usize;
    let win = read_len + 32;
    let mut rng = Rng(0x1234_5678);
    let mut qs: Vec<Vec<u8>> = Vec::with_capacity(n);
    let mut ts: Vec<Vec<u8>> = Vec::with_capacity(n);
    for _ in 0..n {
        let t: Vec<u8> = (0..win).map(|_| rng.b()).collect();
        let mut q: Vec<u8> = t[..read_len].to_vec();
        for _ in 0..3 {
            let pos = (rng.next() as usize) % read_len;
            q[pos] = rng.b();
        }
        qs.push(q);
        ts.push(t);
    }
    let tasks: Vec<GpuTask> = qs
        .iter()
        .zip(&ts)
        .map(|(q, t)| GpuTask {
            query5: q,
            target5: t,
        })
        .collect();

    // CPU: ksw2rs (the real salmon backend), single thread, reusing one aligner.
    let mut mat = [-MISMATCH; 25];
    for i in 0..5 {
        mat[i * 5 + i] = MATCH;
    }
    mat[24] = 0;
    let mut aln = Aligner::new();
    let t0 = Instant::now();
    let mut cpu_sum = 0i64;
    for (q, t) in qs.iter().zip(&ts) {
        cpu_sum += ksw_score(&mut aln, q, t, &mat) as i64;
    }
    let cpu = t0.elapsed();

    // GPU warm-up then timed dispatch.
    let _ = ctx.batch_align(&tasks[..2000.min(n)], &p);
    let t1 = Instant::now();
    let g = ctx.batch_align(&tasks, &p);
    let gpu = t1.elapsed();
    // Keep both result sets live so neither loop is optimized away.
    std::hint::black_box((cpu_sum, &g));

    let mref = n as f64 / cpu.as_secs_f64() / 1e6;
    let mgpu = n as f64 / gpu.as_secs_f64() / 1e6;
    eprintln!("n={n} read_len={read_len} band={W}");
    eprintln!("CPU ksw2rs (1 thread): {:?}  ({mref:.2} M aln/s)", cpu);
    eprintln!("GPU (1 dispatch):      {:?}  ({mgpu:.2} M aln/s)", gpu);
    eprintln!(
        "GPU vs 1 ksw2rs thread: {:.2}x",
        cpu.as_secs_f64() / gpu.as_secs_f64()
    );
}
