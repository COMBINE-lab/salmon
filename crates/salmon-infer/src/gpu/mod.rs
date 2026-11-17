//! Experimental GPU warm start for the EM/VBEM point estimate (on by default
//! in builds with a backend; `--no-gpu` turns it off).
//!
//! The GPU does not replace the optimizer. It runs the same fixed-point
//! iteration in `f32` until the global relative-difference criterion is met,
//! and hands the result to the unchanged CPU engine as a warm start, which then
//! finishes in `f64` under its own stopping rule. Outputs therefore keep their
//! usual guarantees of convergence and mass, are deterministic for a given GPU
//! (gather-only kernels, no atomics), and are *not* bit-identical to a CPU-only
//! run. See `docs/em-architecture-and-optimization.md`.
//!
//! All backends (Metal, CUDA, ROCm) run the same kernels over [`Problem`]:
//! 1. `et` (VBEM only): `et[t] = exp(ψ(alpha[t] + prior[t]))`, the global ψ
//!    normalizer dropped as in the component engine (it cancels per class);
//! 2. `chunk_dot` over the class CSR, then `class_finish`:
//!    `scale[c] = count[c] / Σ basis[label]·w` per class;
//! 3. `chunk_dot` over the transcript CSR, then `txp_finish`:
//!    `out[t] = single[t] + basis[t] · Σ w·scale[class]` per transcript, where
//!    `single[t]` is the count of `t`'s single-transcript classes, credited
//!    whatever `basis[t]` is (as the CPU kernels do).
//!
//! `chunk_dot` has one thread per chunk of at most [`CHUNK_LEN`] incidences
//! (rows are cut into chunks), so threads stay about equally busy although
//! row lengths span three orders of magnitude; the finish kernels sum each
//! row's chunk partials.
//!
//! shortcut: `exp(ψ(x))` underflows `f32` below `x ≈ 0.0115`, so under VBEM a
//! class whose members *all* sit below that is dropped on the GPU, where `f64`
//! still splits it. A class supplies its own count to its members every step,
//! so this needs a start with such zeros; the `f64` CPU finish restores it.
//! Upgrade to per-class log-space scaling if warm starts like that matter.

use crate::packed::PackedEqClasses;
#[cfg(gpu)]
use crate::{max_rel_diff, EmOptions, GPU_BACKEND};

// `build.rs` sets `cfg(gpu = ...)` to at most one backend.
#[cfg(gpu = "cuda")]
use cuda as backend;
#[cfg(gpu = "metal")]
use metal as backend;
#[cfg(gpu = "rocm")]
use rocm as backend;
#[cfg(gpu = "cuda")]
mod cuda;
#[cfg(any(gpu = "cuda", gpu = "rocm"))]
mod dl;
#[cfg(gpu = "metal")]
mod metal;
#[cfg(gpu = "rocm")]
mod rocm;
#[cfg(gpu)]
pub(crate) use backend::{available, KmerDevice};

/// M-steps submitted between two convergence checks (each check reads back
/// two iterates).
#[cfg(gpu)]
const CHUNK: u32 = 10;

/// Below this many class–transcript incidences the CPU is faster: the GPU's
/// fixed setup (25–50 ms on Metal) outweighs its faster M-steps.
// shortcut: crossover measured once (Apple M4 Max, 16 CPU threads, synthetic
// classes; EM gains from ~150k, VBEM from ~500k); recalibrate per backend.
#[cfg(gpu)]
const MIN_INCIDENCES: usize = if cfg!(test) { 1 } else { 250_000 };

/// The packed classes re-laid out for gather-only GPU kernels, in `f32`.
///
/// Single-transcript classes are folded into `single`; the others keep their
/// packed order, whose neighbouring classes share transcripts, so the gathers
/// stay local (sorting rows by size measured 1.3-1.5x slower on Metal for that
/// loss of locality; the chunks balance the load instead). Every array is
/// indexed with `i32`-sized offsets (cuTile pointer arithmetic), hence the size
/// check in [`Problem::new`].
pub(crate) struct Problem {
    /// transcript id of each incidence
    pub labels: Vec<u32>,
    /// `combined` weight (`weight / effLen`) of each incidence
    pub w: Vec<f32>,
    /// fragment count of each class
    pub count: Vec<f32>,
    /// class index (into `count`) of each incidence, ascending within a transcript
    pub t_cls: Vec<u32>,
    /// the same weight as `w`, in transcript order
    pub t_w: Vec<f32>,
    /// summed count of each transcript's single-transcript classes
    pub single: Vec<f32>,
    /// the class CSR (rows into `labels`/`w`) and its transpose (rows into
    /// `t_cls`/`t_w`), cut into chunks of at most [`CHUNK_LEN`] incidences:
    /// chunk `k` spans `c_chunk[k]..c_chunk[k + 1]` of `labels`, and class `c`
    /// owns chunks `c_first[c]..c_first[c + 1]` (likewise `t_chunk`/`t_first`).
    /// One GPU thread per chunk keeps a SIMD group's threads equally busy:
    /// real transcript degrees run from ~50 to ~45,000.
    pub c_chunk: Vec<u32>,
    pub c_first: Vec<u32>,
    pub t_chunk: Vec<u32>,
    pub t_first: Vec<u32>,
}

/// Incidences per chunk (see [`Problem::c_chunk`]).
pub(crate) const CHUNK_LEN: u32 = 32;

/// Cut the rows of a CSR (`start`, `len = rows + 1`) into chunks of at most
/// [`CHUNK_LEN`]: returns the chunk boundaries and each row's first chunk.
fn chunks(start: &[u32]) -> (Vec<u32>, Vec<u32>) {
    let mut bounds = Vec::new();
    let mut first = Vec::with_capacity(start.len());
    for r in start.windows(2) {
        first.push(bounds.len() as u32);
        bounds.extend((r[0]..r[1]).step_by(CHUNK_LEN as usize));
    }
    first.push(bounds.len() as u32);
    bounds.push(*start.last().unwrap_or(&0));
    (bounds, first)
}

impl Problem {
    pub(crate) fn new(p: &PackedEqClasses, counts: &[u64]) -> Option<Self> {
        let n = p.labels.len();
        let nt = p.num_txps;
        if n > i32::MAX as usize || nt >= i32::MAX as usize {
            return None;
        }
        let mut single = vec![0f32; nt];
        let mut c_start = vec![0u32];
        let mut labels = Vec::with_capacity(n);
        let mut w = Vec::with_capacity(n);
        let mut count = Vec::new();
        for c in 0..p.num_classes() {
            let (tids, ws) = p.class(c);
            match tids.len() {
                0 => {}
                1 => single[tids[0] as usize] += counts[c] as f32,
                _ => {
                    labels.extend_from_slice(tids);
                    w.extend(ws.iter().map(|&x| x as f32));
                    count.push(counts[c] as f32);
                    c_start.push(labels.len() as u32);
                }
            }
        }
        // transpose by counting sort; classes are visited in ascending order
        let mut t_start = vec![0u32; nt + 1];
        for &t in &labels {
            t_start[t as usize + 1] += 1;
        }
        for i in 0..nt {
            t_start[i + 1] += t_start[i];
        }
        let mut fill: Vec<u32> = t_start[..nt].to_vec();
        let mut t_cls = vec![0u32; labels.len()];
        let mut t_w = vec![0f32; labels.len()];
        for c in 0..count.len() {
            for e in c_start[c] as usize..c_start[c + 1] as usize {
                let t = labels[e] as usize;
                t_cls[fill[t] as usize] = c as u32;
                t_w[fill[t] as usize] = w[e];
                fill[t] += 1;
            }
        }
        let (c_chunk, c_first) = chunks(&c_start);
        let (t_chunk, t_first) = chunks(&t_start);
        Some(Self {
            labels,
            w,
            count,
            t_cls,
            t_w,
            single,
            c_chunk,
            c_first,
            t_chunk,
            t_first,
        })
    }

    pub(crate) fn num_classes(&self) -> usize {
        self.count.len()
    }

    pub(crate) fn num_txps(&self) -> usize {
        self.single.len()
    }
}

/// Iterate on the GPU from `init` until the global criterion holds (after
/// `min_iter` steps) or `max_iter`; returns the last iterate and the number of
/// M-steps run. `None` — a problem too small to gain or too large for the
/// kernels' index width, or no usable device — means "continue on the CPU from
/// `init`".
#[cfg(gpu)]
pub(crate) fn converge(
    p: &PackedEqClasses,
    counts: &[u64],
    init: &[f64],
    prior: &[f64],
    opts: &EmOptions,
    min_iter: u32,
) -> Option<(Vec<f64>, u32)> {
    let name = GPU_BACKEND.unwrap_or("gpu");
    if p.labels.len() < MIN_INCIDENCES {
        return None;
    }
    if !backend::available() {
        tracing::debug!("no {name} GPU device; EM on the CPU");
        return None;
    }
    let Some(prob) = Problem::new(p, counts) else {
        tracing::warn!("GPU: problem exceeds the {name} kernels' 32-bit indexing; EM on the CPU");
        return None;
    };
    let run = || -> Result<(Vec<f64>, u32), String> {
        let mut dev = backend::Device::new(&prob, init, prior, opts.use_vbem)?;
        let mut it = 0u32;
        loop {
            let n = CHUNK.min(opts.max_iter - it);
            dev.run(n)?;
            it += n;
            let (prev, cur) = dev.read()?;
            let cur: Vec<f64> = cur.iter().map(|&x| f64::from(x)).collect();
            if it >= opts.max_iter {
                return Ok((cur, it));
            }
            let prev: Vec<f64> = prev.iter().map(|&x| f64::from(x)).collect();
            if it >= min_iter
                && max_rel_diff(&prev, &cur, opts.alpha_check_cutoff) < opts.rel_diff_tol
            {
                return Ok((cur, it));
            }
        }
    };
    match run() {
        Ok((alphas, it)) => {
            tracing::info!("GPU ({name}) warm start ran {it} M-steps; finishing on the CPU");
            Some((alphas, it))
        }
        Err(e) => {
            tracing::warn!("GPU ({name}) failed ({e}); EM on the CPU (--no-gpu skips the attempt)");
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small packed problem with uneven class sizes and one isolated transcript.
    pub(super) fn fixture() -> PackedEqClasses {
        let classes: &[(&[u32], &[f64], u64)] = &[
            (&[0, 1, 2], &[0.3, 0.5, 0.2], 7),
            (&[1], &[1.0], 3),
            (&[2, 3], &[0.25, 0.75], 11),
            (&[0, 1, 2, 3, 5], &[0.1, 0.2, 0.3, 0.15, 0.25], 5),
            (&[5], &[1.0], 2),
        ];
        let mut p = PackedEqClasses {
            labels: Vec::new(),
            starts: vec![0],
            combined: Vec::new(),
            weights: Vec::new(),
            counts: Vec::new(),
            num_txps: 6,
            total_count: 0,
        };
        for &(t, w, c) in classes {
            p.labels.extend_from_slice(t);
            p.combined.extend_from_slice(w);
            p.starts.push(p.labels.len() as u64);
            p.counts.push(c);
            p.total_count += c;
        }
        p
    }

    #[test]
    fn layout_covers_every_incidence_once_in_both_directions() {
        let p = fixture();
        let g = Problem::new(&p, &p.counts).unwrap();
        let row_starts = |bounds: &[u32], first: &[u32]| -> Vec<u32> {
            first.iter().map(|&k| bounds[k as usize]).collect()
        };
        let c_start = row_starts(&g.c_chunk, &g.c_first);
        let t_start = row_starts(&g.t_chunk, &g.t_first);
        assert_eq!(g.num_txps(), p.num_txps);
        // single-transcript classes are folded into `single`
        assert_eq!(g.single, [0.0, 3.0, 0.0, 0.0, 0.0, 2.0]);
        // the others keep their order, members and count
        let seen: Vec<(Vec<u32>, u64)> = (0..g.num_classes())
            .map(|c| {
                let r = c_start[c] as usize..c_start[c + 1] as usize;
                (g.labels[r].to_vec(), g.count[c] as u64)
            })
            .collect();
        let want: Vec<(Vec<u32>, u64)> = (0..p.num_classes())
            .map(|c| (p.class(c).0.to_vec(), p.counts[c]))
            .filter(|(t, _)| t.len() > 1)
            .collect();
        assert_eq!(seen, want);
        // the transcript side is the transpose of the class side
        assert_eq!(*t_start.last().unwrap() as usize, g.labels.len());
        for t in 0..p.num_txps as u32 {
            let r = t_start[t as usize] as usize..t_start[t as usize + 1] as usize;
            let cls = &g.t_cls[r.clone()];
            assert!(cls.windows(2).all(|s| s[0] < s[1]), "ascending class order");
            for (&c, &w) in cls.iter().zip(&g.t_w[r]) {
                let e = (c_start[c as usize]..c_start[c as usize + 1])
                    .find(|&e| g.labels[e as usize] == t)
                    .expect("transcript belongs to the class it lists");
                assert_eq!(g.w[e as usize], w);
            }
        }
        // transcript 4 is in no class
        assert_eq!(t_start[4], t_start[5]);
        check_chunks(&c_start, &g.c_chunk, &g.c_first);
        check_chunks(&t_start, &g.t_chunk, &g.t_first);
    }

    /// Each row's chunks tile exactly its range, each 1..=CHUNK_LEN long.
    fn check_chunks(start: &[u32], bounds: &[u32], first: &[u32]) {
        assert_eq!(first.len(), start.len());
        assert_eq!(bounds.last(), start.last());
        for r in 0..start.len() - 1 {
            let mut e = start[r];
            for k in first[r] as usize..first[r + 1] as usize {
                let (a, b) = (bounds[k], bounds[k + 1]);
                assert!(a == e && a < b && b - a <= CHUNK_LEN, "row {r} chunk {k}");
                e = b;
            }
            assert_eq!(e, start[r + 1], "row {r} fully covered");
        }
    }

    #[test]
    fn long_and_empty_rows_are_chunked() {
        let start = [0, 0, 5, 5 + 2 * CHUNK_LEN + 1, 5 + 2 * CHUNK_LEN + 1];
        let (bounds, first) = chunks(&start);
        check_chunks(&start, &bounds, &first);
        assert_eq!(first, [0, 0, 1, 4, 4]);
    }

    /// A pseudo-random problem with uneven class sizes, large enough to span
    /// many threadgroups / tiles.
    #[cfg(gpu)]
    fn random_fixture(num_txps: usize, num_classes: usize) -> PackedEqClasses {
        let mut seed = 0x9E37_79B9_7F4A_7C15u64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let mut p = PackedEqClasses {
            labels: Vec::new(),
            starts: vec![0],
            combined: Vec::new(),
            weights: Vec::new(),
            counts: Vec::new(),
            num_txps,
            total_count: 0,
        };
        for c in 0..num_classes {
            let n = 1 + (rnd() % [1, 2, 4, 16, 64][c % 5]) as usize;
            let mut tids: Vec<u32> = (0..n).map(|_| (rnd() % num_txps as u64) as u32).collect();
            tids.sort_unstable();
            tids.dedup();
            for _ in &tids {
                p.combined
                    .push(1e-3 * (1.0 + (rnd() % 1000) as f64 / 997.0));
            }
            p.labels.extend_from_slice(&tids);
            p.starts.push(p.labels.len() as u64);
            let count = 1 + rnd() % 50;
            p.counts.push(count);
            p.total_count += count;
        }
        p
    }

    #[cfg(gpu)]
    #[test]
    fn one_gpu_step_matches_the_cpu_kernels() {
        if !backend::available() {
            eprintln!("skipping: no {} device", GPU_BACKEND.unwrap());
            return;
        }
        for p in [fixture(), random_fixture(3000, 8000)] {
            let g = Problem::new(&p, &p.counts).unwrap();
            let prior = vec![1e-2; p.num_txps];
            for vb in [false, true] {
                // EM starts with zero abundances, which single-transcript classes
                // still credit; VBEM stays clear of the f32 `exp(ψ)` underflow.
                let floor = if vb { 0.5 } else { 0.0 };
                let init: Vec<f64> = (0..p.num_txps)
                    .map(|i| floor + (i % 7) as f64 * 1.5)
                    .collect();
                let mut want = vec![0.0; p.num_txps];
                let mut scratch = Vec::new();
                if vb {
                    let mut et = vec![0.0; p.num_txps];
                    crate::packed::vbem_step_seq::<false>(
                        &p,
                        &p.counts,
                        &prior,
                        &init,
                        &mut want,
                        &mut et,
                        &mut scratch,
                    );
                } else {
                    crate::packed::em_step_seq::<false>(
                        &p,
                        &p.counts,
                        &init,
                        &mut want,
                        &mut scratch,
                    );
                }
                let mut dev = backend::Device::new(&g, &init, &prior, vb).unwrap();
                dev.run(1).unwrap();
                let (prev, got) = dev.read().unwrap();
                assert!(prev.iter().zip(&init).all(|(&a, &b)| a == b as f32));
                for (t, (&w, &g)) in want.iter().zip(&got).enumerate() {
                    let err = (f64::from(g) - w).abs() / w.abs().max(1e-3);
                    assert!(err < 1e-4, "vb={vb} t={t}: gpu {g} cpu {w}");
                }
            }
        }
    }

    #[cfg(gpu)]
    #[test]
    fn gpu_warm_start_converges_to_the_cpu_answer() {
        use crate::{optimize_packed_with_init, EffLens, InitAlphas};
        if !backend::available() {
            eprintln!("skipping: no {} device", GPU_BACKEND.unwrap());
            return;
        }
        let p = random_fixture(3000, 8000);
        for vb in [false, true] {
            let cpu_opts = EmOptions {
                use_vbem: vb,
                rel_diff_tol: 1e-6,
                max_iter: 100_000,
                ..Default::default()
            };
            let gpu_opts = EmOptions {
                gpu: true,
                ..cpu_opts.clone()
            };
            let uniform = vec![p.total_count as f64 / p.num_txps as f64; p.num_txps];
            let prior = vec![cpu_opts.vb_prior; p.num_txps];
            assert!(
                converge(&p, &p.counts, &uniform, &prior, &gpu_opts, 50).is_some(),
                "the GPU path must actually run"
            );
            let cpu =
                optimize_packed_with_init(&p, &cpu_opts, true, InitAlphas::NONE, EffLens::NONE);
            let gpu =
                optimize_packed_with_init(&p, &gpu_opts, true, InitAlphas::NONE, EffLens::NONE);
            assert!(cpu.converged && gpu.converged, "vb={vb}");
            assert_eq!((cpu.gpu_backend, gpu.gpu_backend), (None, GPU_BACKEND));
            let mass = |a: &[f64]| a.iter().sum::<f64>();
            let (mc, mg) = (mass(&cpu.alphas), mass(&gpu.alphas));
            assert!((mc - mg).abs() < 1e-6 * mc, "vb={vb}: mass {mc} vs {mg}");
            for (t, (&c, &g)) in cpu.alphas.iter().zip(&gpu.alphas).enumerate() {
                if c > 1.0 {
                    assert!((g - c).abs() / c < 1e-3, "vb={vb} t={t}: gpu {g} cpu {c}");
                }
            }
        }
    }
}
