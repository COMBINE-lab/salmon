//! `--sketch` mapping with the GPU k-mer lookup (`salmon_infer::kmers`): one
//! GPU call per batch finds every read's raw hits, then each fragment is mapped
//! from them (`salmon_map::sketch_hits`) exactly as the CPU path would; a read past `HITS_CAP` hits, or one the replay declines, is mapped on
//! the CPU instead.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::OnceLock;

use salmon_infer::kmers::{Hit, KmerGpu, HITS_CAP};

/// The device lookup once its build thread has set it, the switch that sends
/// every later batch to the CPU after a GPU failure, and whether a batch ran
/// on the GPU.
#[derive(Clone, Copy)]
pub(crate) struct GpuMap<'a> {
    pub gpu: &'a OnceLock<KmerGpu>,
    pub off: &'a AtomicBool,
    pub used: &'a AtomicBool,
}

/// Raise the calling thread's scheduling class for the GPU lookup's build:
/// the device waits on it while the mapping threads (and on a busy machine
/// other processes) compete for the cores, and finishing it sooner moves
/// more batches to the GPU. macOS QoS only; elsewhere a no-op.
#[allow(unsafe_code)]
pub(crate) fn urgent() {
    #[cfg(target_os = "macos")]
    {
        extern "C" {
            fn pthread_set_qos_class_self_np(qos: u32, relpri: i32) -> i32;
        }
        const QOS_CLASS_USER_INTERACTIVE: u32 = 0x21;
        // SAFETY: changes only the calling thread's QoS class; a failure
        // leaves it at the default.
        unsafe { pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0) };
    }
}

/// One batch's reads and their hits (reused across a worker's batches).
#[derive(Default)]
pub(crate) struct BatchHits {
    seq: Vec<u8>,
    off: Vec<u32>,
    hits: Vec<Hit>,
    counts: Vec<u32>,
}

impl BatchHits {
    /// Look `reads` up on the GPU; `false` means map this batch on the CPU.
    pub fn lookup<R: AsRef<[u8]>>(&mut self, g: GpuMap, reads: impl Iterator<Item = R>) -> bool {
        let Some(gpu) = g.gpu.get().filter(|_| !g.off.load(Ordering::Relaxed)) else {
            return false;
        };
        self.seq.clear();
        self.off.clear();
        self.off.push(0);
        for r in reads {
            self.seq.extend_from_slice(r.as_ref());
            // a batch is a few thousand reads, far below 4 GiB
            self.off.push(self.seq.len() as u32);
        }
        match gpu.hits(&self.seq, &self.off, &mut self.hits, &mut self.counts) {
            Ok(()) => {
                g.used.store(true, Ordering::Relaxed);
                true
            }
            Err(e) => {
                if !g.off.swap(true, Ordering::Relaxed) {
                    tracing::warn!("GPU k-mer lookup failed ({e}); mapping the rest on the CPU");
                }
                false
            }
        }
    }

    /// Read `i`'s hits, `None` when it had more than the GPU kept or met one
    /// of sshash's heavy buckets.
    pub fn read(&self, i: usize) -> Option<&[Hit]> {
        let n = self.counts[i] as usize;
        (n <= HITS_CAP).then(|| &self.hits[i * HITS_CAP..i * HITS_CAP + n])
    }
}
