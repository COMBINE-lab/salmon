//! [`Aligner`] backends that score candidates with the banded full-length DP.
//!
//! These implement salmon-map's [`Aligner`] trait, so they are drop-in for the
//! [`CpuAligner`](salmon_map::CpuAligner) the mapper uses in full-length mode:
//!   * [`RefAligner`] scores on the CPU via [`crate::reference`], and
//!   * [`GpuAligner`] (feature `gpu`) scores a whole batch on the GPU.
//!
//! Both reproduce `ksw2_align_score`'s full-length scoring on every alignment
//! that affects quant output (see the crate's differential tests), so selecting
//! either is equivalent to `--fullLengthAlignment` and preserves determinism.
//! They always score full-length regardless of [`AlignConfig::full_length_alignment`];
//! the GPU path is only wired in when full-length mode is in effect.

use salmon_map::{
    full_length_window, min_accepted_score, AlignConfig, AlignTask, Aligner, Alignment,
};

use crate::reference::{banded_extz_score_dna5, dna5, BandedParams};

fn banded_params(cfg: &AlignConfig) -> BandedParams {
    BandedParams {
        match_score: cfg.match_score as i32,
        mismatch_pen: cfg.mismatch_pen as i32,
        gap_open_pen: cfg.gap_open_pen as i32,
        gap_extend_pen: cfg.gap_extend_pen as i32,
        bandwidth: cfg.bandwidth,
    }
}

/// The DNA5 query/target window for a task, plus the bookkeeping needed to build
/// its [`Alignment`]. `None` when the chain implies no usable window (matching
/// [`align_chain`](salmon_map::align_chain) returning `None`).
struct Prepared {
    q5: Vec<u8>,
    t5: Vec<u8>,
    win_start: i32,
    read_len: usize,
}

fn prepare(task: &AlignTask, cfg: &AlignConfig) -> Option<Prepared> {
    let (query, rwin, win_start) = full_length_window(task.read, task.ref_seq, task.chain, cfg)?;
    Some(Prepared {
        q5: query.iter().map(|&b| dna5(b)).collect(),
        t5: rwin.iter().map(|&b| dna5(b)).collect(),
        win_start,
        read_len: query.len(),
    })
}

fn make_alignment(score: i32, prep: &Prepared, cfg: &AlignConfig) -> Alignment {
    Alignment {
        score,
        valid: score > min_accepted_score(prep.read_len, cfg),
        ref_window_start: prep.win_start,
    }
}

/// CPU full-length aligner using the banded-DP reference. Deterministic and, on
/// every score-affecting alignment, equal to `ksw2_align_score`. Also the
/// fallback when no GPU is present and the oracle the GPU backend is tested
/// against.
pub struct RefAligner;

impl Aligner for RefAligner {
    fn align_batch(&self, tasks: &[AlignTask], cfg: &AlignConfig) -> Vec<Option<Alignment>> {
        let p = banded_params(cfg);
        tasks
            .iter()
            .map(|t| {
                prepare(t, cfg).map(|prep| {
                    let score = banded_extz_score_dna5(&prep.q5, &prep.t5, &p);
                    make_alignment(score, &prep, cfg)
                })
            })
            .collect()
    }
}

#[cfg(feature = "gpu")]
pub use gpu_backend::GpuAligner;

#[cfg(feature = "gpu")]
mod gpu_backend {
    use super::*;
    use crate::gpu::{GpuContext, GpuTask};

    /// GPU full-length aligner: scores an entire batch of candidate alignments in
    /// one dispatch. Build once (it acquires the device); reuse across batches.
    pub struct GpuAligner {
        ctx: GpuContext,
    }

    impl GpuAligner {
        /// Acquire a GPU and build the pipeline, or `None` if no adapter exists.
        pub fn new() -> Option<Self> {
            GpuContext::new().map(|ctx| Self { ctx })
        }
    }

    impl Aligner for GpuAligner {
        fn align_batch(&self, tasks: &[AlignTask], cfg: &AlignConfig) -> Vec<Option<Alignment>> {
            let p = banded_params(cfg);
            // Prepare every task's window; only those with a usable window go to
            // the GPU, the rest stay `None` (as align_chain would return).
            let preps: Vec<Option<Prepared>> = tasks.iter().map(|t| prepare(t, cfg)).collect();
            let mut gpu_tasks = Vec::new();
            let mut origin = Vec::new(); // gpu-task index -> original task index
            for (i, pr) in preps.iter().enumerate() {
                if let Some(prep) = pr {
                    gpu_tasks.push(GpuTask {
                        query5: &prep.q5,
                        target5: &prep.t5,
                    });
                    origin.push(i);
                }
            }
            let scores = self.ctx.batch_align(&gpu_tasks, &p);

            let mut out: Vec<Option<Alignment>> = vec![None; tasks.len()];
            for (gi, &oi) in origin.iter().enumerate() {
                let prep = preps[oi].as_ref().unwrap();
                out[oi] = Some(make_alignment(scores[gi], prep, cfg));
            }
            out
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use salmon_map::{AlignConfig, CpuAligner, Mem, MemChain};

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

    fn full_length_cfg() -> AlignConfig {
        AlignConfig {
            full_length_alignment: true,
            ..AlignConfig::default()
        }
    }

    #[test]
    fn ref_aligner_matches_cpu_full_length() {
        // RefAligner must reproduce the CpuAligner (ksw2) in full-length mode on
        // every accepted alignment and agree on validity. Build a few tasks: an
        // exact read, a read with mismatches, and an unrelated read.
        let reference = gen_seq(400, 7);
        let exact = reference[50..130].to_vec();
        let mut mutated = reference[140..220].to_vec();
        for i in [10, 25, 40] {
            mutated[i] = if mutated[i] == b'A' { b'C' } else { b'A' };
        }
        let foreign = gen_seq(80, 9999);
        let reads = [exact, mutated, foreign];
        let chains = [
            MemChain::new(vec![Mem::new(0, 50, 31)], 31.0, true),
            MemChain::new(vec![Mem::new(0, 140, 31)], 31.0, true),
            MemChain::new(vec![Mem::new(0, 50, 31)], 31.0, true),
        ];
        let cfg = full_length_cfg();
        let tasks: Vec<AlignTask> = reads
            .iter()
            .zip(&chains)
            .map(|(r, c)| AlignTask {
                read: r,
                ref_seq: &reference,
                chain: c,
            })
            .collect();

        let cpu = CpuAligner.align_batch(&tasks, &cfg);
        let refa = RefAligner.align_batch(&tasks, &cfg);
        assert_eq!(cpu.len(), refa.len());
        for (i, (c, r)) in cpu.iter().zip(&refa).enumerate() {
            match (c, r) {
                (Some(c), Some(r)) => {
                    // validity always agrees; score agrees on accepted alignments.
                    assert_eq!(c.valid, r.valid, "validity differs at task {i}");
                    if c.valid {
                        assert_eq!(c.score, r.score, "score differs at task {i}");
                    }
                }
                (None, None) => {}
                _ => panic!("Some/None mismatch at task {i}"),
            }
        }
    }

    #[cfg(feature = "gpu")]
    #[test]
    fn gpu_aligner_matches_ref_aligner() {
        // Through the Aligner trait, the GPU backend produces Alignments
        // bit-identical to the CPU reference backend (same banded DP).
        let Some(gpu) = GpuAligner::new() else {
            eprintln!("no GPU adapter; skipping");
            return;
        };
        let reference = gen_seq(600, 13);
        let mut reads = Vec::new();
        let mut chains = Vec::new();
        for k in 0..40usize {
            let start = 20 + k * 12;
            let mut r = reference[start..start + 70].to_vec();
            if k % 3 == 0 {
                r[35] = if r[35] == b'A' { b'C' } else { b'A' };
            }
            if k % 5 == 0 && r.len() > 40 {
                r.drain(30..33); // small deletion
            }
            chains.push(MemChain::new(
                vec![Mem::new(0, start as i32, 31)],
                31.0,
                true,
            ));
            reads.push(r);
        }
        let cfg = full_length_cfg();
        let tasks: Vec<AlignTask> = reads
            .iter()
            .zip(&chains)
            .map(|(r, c)| AlignTask {
                read: r,
                ref_seq: &reference,
                chain: c,
            })
            .collect();
        let g = gpu.align_batch(&tasks, &cfg);
        let r = RefAligner.align_batch(&tasks, &cfg);
        assert_eq!(g, r, "GPU and reference Alignments must be identical");
    }
}
