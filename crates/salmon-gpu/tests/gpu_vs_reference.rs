//! GPU kernel vs. CPU reference: the wgpu/WGSL banded-DP must return scores
//! bit-identical to [`salmon_gpu::reference`] over a randomized corpus. Run with
//! `--features gpu`; the test no-ops with a note if no GPU adapter is present
//! (e.g. headless CI without a software Vulkan device).

#![cfg(feature = "gpu")]

use salmon_gpu::gpu::{GpuContext, GpuTask};
use salmon_gpu::reference::{banded_extz_score_dna5, dna5, BandedParams};

fn params(w: i32) -> BandedParams {
    BandedParams {
        match_score: 2,
        mismatch_pen: 4,
        gap_open_pen: 6,
        gap_extend_pen: 2,
        bandwidth: w,
    }
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
    fn base5(&mut self) -> u8 {
        self.below(4) as u8
    }
}

#[test]
fn gpu_matches_reference_over_corpus() {
    let Some(ctx) = GpuContext::new() else {
        eprintln!("no GPU adapter available; skipping GPU differential test");
        return;
    };

    let mut rng = Rng(0xA11_6A704);
    let margin = 32usize;

    // Build a mixed batch of DNA5 query/target pairs across scenarios.
    let mut queries: Vec<Vec<u8>> = Vec::new();
    let mut targets: Vec<Vec<u8>> = Vec::new();
    for _ in 0..5000 {
        let qlen = 40 + rng.below(120);
        let target: Vec<u8> = (0..qlen + margin).map(|_| rng.base5()).collect();
        let mut query: Vec<u8> = target[..qlen].to_vec();
        match rng.below(5) {
            0 => {}
            1 => {
                for _ in 0..1 + rng.below(6) {
                    let pos = rng.below(qlen);
                    query[pos] = rng.base5();
                }
            }
            2 => {
                let d = 1 + rng.below(4);
                let at = rng.below(qlen.saturating_sub(d).max(1));
                query.drain(at..(at + d).min(query.len()));
            }
            3 => {
                let ins = 1 + rng.below(4);
                let at = rng.below(qlen);
                let extra: Vec<u8> = (0..ins).map(|_| rng.base5()).collect();
                query.splice(at..at, extra);
            }
            _ => {
                query = (0..qlen).map(|_| rng.base5()).collect();
            }
        }
        if rng.below(5) == 0 && !query.is_empty() {
            let pos = rng.below(query.len());
            query[pos] = dna5(b'N');
        }
        queries.push(query);
        targets.push(target);
    }

    for &w in &[15i32, 31] {
        let p = params(w);
        let tasks: Vec<GpuTask> = queries
            .iter()
            .zip(&targets)
            .map(|(q, t)| GpuTask {
                query5: q,
                target5: t,
            })
            .collect();
        let gpu = ctx.batch_align(&tasks, &p);
        assert_eq!(gpu.len(), tasks.len());
        let mut mismatches = 0;
        for (i, t) in tasks.iter().enumerate() {
            let want = banded_extz_score_dna5(t.query5, t.target5, &p);
            if gpu[i] != want {
                if mismatches < 10 {
                    eprintln!(
                        "w={} task={} qlen={} tlen={} gpu={} reference={}",
                        w,
                        i,
                        t.query5.len(),
                        t.target5.len(),
                        gpu[i],
                        want
                    );
                }
                mismatches += 1;
            }
        }
        assert_eq!(
            mismatches, 0,
            "{mismatches} GPU/reference mismatches at w={w}"
        );
    }
}
