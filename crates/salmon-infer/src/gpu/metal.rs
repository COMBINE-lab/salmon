//! Metal backend of the GPU warm start (macOS).
//!
//! One thread per chunk of incidences, then one per class or transcript to sum
//! the chunks (see `gpu/mod.rs`). Buffers use shared
//! storage, so on Apple silicon's unified memory a read-back is a slice copy.
#![allow(unsafe_code)]

use std::ptr::NonNull;

use objc2::rc::{autoreleasepool, Retained};
use objc2::runtime::ProtocolObject;
use objc2_foundation::{NSError, NSRange, NSString};
use objc2_metal::{
    MTLBlitCommandEncoder, MTLBuffer, MTLCommandBuffer, MTLCommandBufferStatus, MTLCommandEncoder,
    MTLCommandQueue, MTLCompileOptions, MTLComputeCommandEncoder, MTLComputePipelineState,
    MTLCreateSystemDefaultDevice, MTLDevice, MTLLibrary, MTLResourceOptions, MTLSize,
};

use super::Problem;

/// The kernels of `gpu/mod.rs`. Buffer indices are the argument order
/// [`Device::run`] binds; the element count is always the last argument.
const SRC: &str = r#"
#include <metal_stdlib>
using namespace metal;

// exp(psi(x)) for x > 0: f32 port of components.rs `exp_digamma`.
static float exp_digamma(float x) {
    float y = x, s = 0.0f;
    if (x < 10.0f) {
        float num = 0.0f, den = 1.0f;
        for (int j = 0; j < 10; j += 2) {
            float a = x + float(j), b = a + 1.0f;
            num = num * (a * b) + (a + b) * den;
            den *= a * b;
        }
        y = x + 10.0f;
        s = num / den;
    }
    float r = 1.0f / y, r2 = r * r;
    float u = -0.5f * r - r2 * (1.0f / 12.0f - r2 * (1.0f / 120.0f - r2 * (1.0f / 252.0f
            - r2 * (1.0f / 240.0f - r2 * (1.0f / 132.0f - r2 * (691.0f / 32760.0f))))));
    return y * exp(u - s);
}

kernel void et_kernel(device const float* alpha [[buffer(0)]],
                      device const float* prior [[buffer(1)]],
                      device float* et [[buffer(2)]],
                      constant uint& n [[buffer(3)]],
                      uint t [[thread_position_in_grid]]) {
    if (t >= n) return;
    float ap = alpha[t] + prior[t];
    et[t] = ap > 1e-10f ? exp_digamma(ap) : 0.0f;  // DIGAMMA_MIN
}

// partial[k] = sum of x[idx[e]] * w[e] over chunk k (both CSR directions)
kernel void chunk_dot(device const uint* chunk [[buffer(0)]],
                      device const uint* idx [[buffer(1)]],
                      device const float* w [[buffer(2)]],
                      device const float* x [[buffer(3)]],
                      device float* partial [[buffer(4)]],
                      constant uint& n [[buffer(5)]],
                      uint k [[thread_position_in_grid]]) {
    if (k >= n) return;
    float acc = 0.0f;
    for (uint e = chunk[k]; e < chunk[k + 1]; ++e) {
        float b = x[idx[e]];
        acc += b > 0.0f ? b * w[e] : 0.0f;
    }
    partial[k] = acc;
}

kernel void class_finish(device const uint* first [[buffer(0)]],
                         device const float* partial [[buffer(1)]],
                         device const float* count [[buffer(2)]],
                         device float* scale [[buffer(3)]],
                         constant uint& n [[buffer(4)]],
                         uint c [[thread_position_in_grid]]) {
    if (c >= n) return;
    float denom = 0.0f;
    for (uint k = first[c]; k < first[c + 1]; ++k) denom += partial[k];
    scale[c] = denom > 1.17549435e-38f ? count[c] / denom : 0.0f;  // FLT_MIN
}

kernel void txp_finish(device const uint* first [[buffer(0)]],
                       device const float* partial [[buffer(1)]],
                       device const float* basis [[buffer(2)]],
                       device const float* single [[buffer(3)]],
                       device float* out [[buffer(4)]],
                       constant uint& n [[buffer(5)]],
                       uint t [[thread_position_in_grid]]) {
    if (t >= n) return;
    float acc = 0.0f;
    for (uint k = first[t]; k < first[t + 1]; ++k) acc += partial[k];
    out[t] = single[t] + basis[t] * acc;
}
"#;

type Buf = Retained<ProtocolObject<dyn MTLBuffer>>;
type Pso = Retained<ProtocolObject<dyn MTLComputePipelineState>>;

/// Whether this machine exposes a Metal device at all (else the CPU runs; tests skip).
pub(crate) fn available() -> bool {
    MTLCreateSystemDefaultDevice().is_some()
}

pub(super) struct Device {
    queue: Retained<ProtocolObject<dyn MTLCommandQueue>>,
    /// `Some` under VBEM
    et: Option<Pso>,
    dot: Pso,
    class_finish: Pso,
    txp_finish: Pso,
    c_chunk: Buf,
    c_first: Buf,
    labels: Buf,
    w: Buf,
    count: Buf,
    t_chunk: Buf,
    t_first: Buf,
    t_cls: Buf,
    t_w: Buf,
    single: Buf,
    prior: Buf,
    basis: Buf,
    class_scale: Buf,
    /// per-chunk sums, reused by the class pass and then the transcript pass
    partial: Buf,
    /// ping-pong iterates; `alpha[cur]` is the latest
    alpha: [Buf; 2],
    cur: usize,
    nc: u32,
    nt: u32,
    /// chunks on the class and on the transcript side
    ncc: u32,
    ntc: u32,
}

fn describe(e: Retained<NSError>) -> String {
    e.localizedDescription().to_string()
}

fn buffer<T: Copy>(dev: &ProtocolObject<dyn MTLDevice>, data: &[T]) -> Result<Buf, String> {
    let opts = MTLResourceOptions::StorageModeShared;
    let made = if data.is_empty() {
        // Metal refuses zero-length buffers; nothing ever reads this one.
        dev.newBufferWithLength_options(4, opts)
    } else {
        // SAFETY: `data` is a live slice of exactly this many bytes; Metal copies it.
        unsafe {
            dev.newBufferWithBytes_length_options(
                NonNull::from(data).cast(),
                std::mem::size_of_val(data),
                opts,
            )
        }
    };
    made.ok_or_else(|| "Metal buffer allocation failed".to_string())
}

fn pipeline(
    dev: &ProtocolObject<dyn MTLDevice>,
    lib: &ProtocolObject<dyn MTLLibrary>,
    name: &str,
) -> Result<Pso, String> {
    let f = lib
        .newFunctionWithName(&NSString::from_str(name))
        .ok_or_else(|| format!("Metal kernel {name} missing"))?;
    dev.newComputePipelineStateWithFunction_error(&f)
        .map_err(describe)
}

/// Bind `bufs` to indices `0..`, the element count `n` after them, and launch
/// one thread per element.
fn dispatch(enc: &ProtocolObject<dyn MTLComputeCommandEncoder>, pso: &Pso, bufs: &[&Buf], n: u32) {
    if n == 0 {
        return;
    }
    enc.setComputePipelineState(pso);
    for (i, b) in bufs.iter().enumerate() {
        // SAFETY: index `i` is the kernel's `[[buffer(i)]]` argument of this
        // element type (see SRC), and the buffer outlives the command buffer.
        unsafe { enc.setBuffer_offset_atIndex(Some(b), 0, i) };
    }
    let count = n;
    // SAFETY: 4 bytes of a live `u32`, copied by Metal at encode time.
    unsafe {
        enc.setBytes_length_atIndex(NonNull::from(&count).cast(), 4, bufs.len());
    }
    let width = pso.maxTotalThreadsPerThreadgroup().min(256);
    enc.dispatchThreads_threadsPerThreadgroup(
        MTLSize {
            width: n as usize,
            height: 1,
            depth: 1,
        },
        MTLSize {
            width,
            height: 1,
            depth: 1,
        },
    );
}

impl Device {
    pub(super) fn new(
        prob: &Problem,
        init: &[f64],
        prior: &[f64],
        vb: bool,
    ) -> Result<Self, String> {
        autoreleasepool(|_| {
            let dev = MTLCreateSystemDefaultDevice().ok_or("no Metal device")?;
            let opts = MTLCompileOptions::new();
            // IEEE math rather than fast math: the kernels are memory-bound, so
            // fast math would only cost accuracy.
            #[allow(deprecated)]
            opts.setFastMathEnabled(false);
            let lib = dev
                .newLibraryWithSource_options_error(&NSString::from_str(SRC), Some(&opts))
                .map_err(describe)?;
            let f32s = |v: &[f64]| v.iter().map(|&x| x as f32).collect::<Vec<f32>>();
            let nt = prob.num_txps();
            let zeros_t = vec![0f32; nt];
            let (ncc, ntc) = (prob.c_chunk.len() - 1, prob.t_chunk.len() - 1);
            Ok(Self {
                queue: dev.newCommandQueue().ok_or("no Metal command queue")?,
                et: if vb {
                    Some(pipeline(&dev, &lib, "et_kernel")?)
                } else {
                    None
                },
                dot: pipeline(&dev, &lib, "chunk_dot")?,
                class_finish: pipeline(&dev, &lib, "class_finish")?,
                txp_finish: pipeline(&dev, &lib, "txp_finish")?,
                c_chunk: buffer(&dev, &prob.c_chunk)?,
                c_first: buffer(&dev, &prob.c_first)?,
                labels: buffer(&dev, &prob.labels)?,
                w: buffer(&dev, &prob.w)?,
                count: buffer(&dev, &prob.count)?,
                t_chunk: buffer(&dev, &prob.t_chunk)?,
                t_first: buffer(&dev, &prob.t_first)?,
                t_cls: buffer(&dev, &prob.t_cls)?,
                t_w: buffer(&dev, &prob.t_w)?,
                single: buffer(&dev, &prob.single)?,
                prior: buffer(&dev, &f32s(prior))?,
                basis: buffer(&dev, &zeros_t)?,
                class_scale: buffer(&dev, &vec![0f32; prob.num_classes()])?,
                partial: buffer(&dev, &vec![0f32; ncc.max(ntc)])?,
                alpha: [buffer(&dev, &f32s(init))?, buffer(&dev, &zeros_t)?],
                cur: 0,
                nc: prob.num_classes() as u32,
                nt: nt as u32,
                ncc: ncc as u32,
                ntc: ntc as u32,
            })
        })
    }

    /// Run `n` M-steps in one command buffer and wait for them.
    pub(super) fn run(&mut self, n: u32) -> Result<(), String> {
        autoreleasepool(|_| {
            let cb = self
                .queue
                .commandBuffer()
                .ok_or("no Metal command buffer")?;
            let enc = cb
                .computeCommandEncoder()
                .ok_or("no Metal compute encoder")?;
            // A serial compute encoder runs dispatches in order, each seeing
            // the previous one's writes.
            for _ in 0..n {
                let (src, dst) = (&self.alpha[self.cur], &self.alpha[1 - self.cur]);
                let basis = match &self.et {
                    Some(et) => {
                        dispatch(&enc, et, &[src, &self.prior, &self.basis], self.nt);
                        &self.basis
                    }
                    None => src,
                };
                let p = &self.partial;
                let class_dot = [&self.c_chunk, &self.labels, &self.w, basis, p];
                dispatch(&enc, &self.dot, &class_dot, self.ncc);
                let class_fin = [&self.c_first, p, &self.count, &self.class_scale];
                dispatch(&enc, &self.class_finish, &class_fin, self.nc);
                let txp_dot = [&self.t_chunk, &self.t_cls, &self.t_w, &self.class_scale, p];
                dispatch(&enc, &self.dot, &txp_dot, self.ntc);
                let txp_fin = [&self.t_first, p, basis, &self.single, dst];
                dispatch(&enc, &self.txp_finish, &txp_fin, self.nt);
                self.cur = 1 - self.cur;
            }
            enc.endEncoding();
            cb.commit();
            cb.waitUntilCompleted();
            if cb.status() == MTLCommandBufferStatus::Completed {
                Ok(())
            } else {
                Err(cb
                    .error()
                    .map(describe)
                    .unwrap_or_else(|| "Metal command buffer failed".to_string()))
            }
        })
    }

    /// `(previous, latest)` iterates.
    pub(super) fn read(&self) -> Result<(Vec<f32>, Vec<f32>), String> {
        let read = |b: &Buf| {
            // SAFETY: a shared-storage buffer of `nt` f32s (or a 4-byte stub
            // when `nt == 0`), and `run` waited for all GPU work on it.
            unsafe {
                std::slice::from_raw_parts(b.contents().as_ptr().cast::<f32>(), self.nt as usize)
            }
            .to_vec()
        };
        Ok((read(&self.alpha[1 - self.cur]), read(&self.alpha[self.cur])))
    }
}

/// The k-mer lookup kernel (`kmers.rs`): one thread per read.
fn kmer_src() -> String {
    let prelude = "#include <metal_stdlib>\nusing namespace metal;\n\
        typedef uchar u8; typedef uint u32; typedef ulong u64; typedef long i64;\n\
        #define G device\n#define PRIV thread\n#define FN static inline\n";
    let entry = r#"
kernel void kmer_hits(device const u8* seq [[buffer(0)]],
                      device const u32* off [[buffer(1)]],
                      device const u64* keys [[buffer(2)]],
                      device const u32* vals [[buffer(3)]],
                      device const u32* kstart [[buffer(4)]],
                      device const u64* useq [[buffer(5)]],
                      device const u32* prm [[buffer(6)]],
                      device u32* hits [[buffer(7)]],
                      device u32* counts [[buffer(8)]],
                      constant uint& n [[buffer(9)]],
                      uint i [[thread_position_in_grid]]) {
    if (i >= n) return;
    map_mate(i, seq, off, keys, vals, kstart, useq, prm, hits, counts);
}
// the sshash lookup's minimizer table (`KmerTable::build_sshash`): bucket `b`
// claims a slot with a CAS on its fingerprint word (never all ones), then
// writes its id; `bad` when another minimizer's equal fingerprint is ahead
kernel void kmer_table(device const u32* aux [[buffer(0)]],
                       device const u64* words [[buffer(1)]],
                       device const u32* prm [[buffer(2)]],
                       device u32* table [[buffer(3)]],
                       device atomic_uint* bad [[buffer(4)]],
                       constant uint& n [[buffer(5)]],
                       uint b [[thread_position_in_grid]]) {
    if (b >= n || aux[b] == aux[b + 1]) return;
    int m = (int)prm[5];
    u32 cap = prm[1];
    u64 v = word_at(words, (u64)aux[prm[8] + aux[b]], m), r = 0, x = v;
    for (int j = 0; j < m; ++j) { r = (r << 2) | ((x & 3) ^ 2); x >>= 2; }
    u64 h = fmix64(v < r ? v : r);
    u32 s = (u32)(((h >> 32) * (u64)cap) >> 32), fp = (u32)h < 0xFFFFFFFFu ? (u32)h : 0xFFFFFFFEu;
    device atomic_uint* t = (device atomic_uint*)table;
    for (;;) {
        uint seen = 0xFFFFFFFFu;
        if (atomic_compare_exchange_weak_explicit(&t[2 * s + 1], &seen, fp,
                                                  memory_order_relaxed, memory_order_relaxed)) {
            table[2 * s] = b;
            return;
        }
        if (seen == 0xFFFFFFFFu) continue;
        if (seen == fp) { atomic_store_explicit(bad, 1u, memory_order_relaxed); return; }
        s = s + 1 == cap ? 0 : s + 1;
    }
}
"#;
    format!("{prelude}{}{entry}", crate::kmers::KMER_BODY)
}

/// The sshash lookup's minimizer table, built on the device from the bucket
/// positions (`aux`) and the text (`words`): `prm[1]` slots for the
/// `prm[8] - 1` buckets.
fn build_table(
    dev: &ProtocolObject<dyn MTLDevice>,
    queue: &ProtocolObject<dyn MTLCommandQueue>,
    lib: &ProtocolObject<dyn MTLLibrary>,
    aux_b: &Buf,
    words_b: &Buf,
    prm_b: &Buf,
    prm: &[u32],
) -> Result<Buf, String> {
    let (cap, nb) = (prm[1] as usize, prm[8] - 1);
    let table = dev
        .newBufferWithLength_options(cap * 8, MTLResourceOptions::StorageModeShared)
        .ok_or("Metal buffer allocation failed")?;
    let bad = buffer(dev, &[0u32])?;
    let cb = queue.commandBuffer().ok_or("no Metal command buffer")?;
    let blit = cb.blitCommandEncoder().ok_or("no Metal blit encoder")?;
    blit.fillBuffer_range_value(&table, NSRange::new(0, cap * 8), 0xFF);
    blit.endEncoding();
    let enc = cb
        .computeCommandEncoder()
        .ok_or("no Metal compute encoder")?;
    dispatch(
        &enc,
        &pipeline(dev, lib, "kmer_table")?,
        &[aux_b, words_b, prm_b, &table, &bad],
        nb,
    );
    enc.endEncoding();
    cb.commit();
    cb.waitUntilCompleted();
    if cb.status() != MTLCommandBufferStatus::Completed {
        return Err(cb
            .error()
            .map(describe)
            .unwrap_or_else(|| "Metal command buffer failed".to_string()));
    }
    // SAFETY: a 4-byte shared-storage buffer, and the GPU work on it has completed.
    if unsafe { *bad.contents().as_ptr().cast::<u32>() } != 0 {
        return Err("two minimizers share a table fingerprint".into());
    }
    Ok(table)
}

/// The k-mer table on the device, shared by every mapping thread.
pub(crate) struct KmerDevice {
    dev: Retained<ProtocolObject<dyn MTLDevice>>,
    queue: Retained<ProtocolObject<dyn MTLCommandQueue>>,
    pso: Pso,
    keys: Buf,
    vals: Buf,
    kstart: Buf,
    useq: Buf,
    prm: Buf,
}

// SAFETY: MTLDevice, MTLCommandQueue, MTLComputePipelineState and MTLBuffer
// are thread-safe in Metal; each `hits` call encodes into its own command
// buffer and only reads the shared buffers.
unsafe impl Send for KmerDevice {}
unsafe impl Sync for KmerDevice {}

impl KmerDevice {
    pub(crate) fn new(
        keys: &[u64],
        vals: &[u32],
        kstart: &[u32],
        useq: &[u64],
        prm: &[u32],
    ) -> Result<Self, String> {
        autoreleasepool(|_| {
            let dev = MTLCreateSystemDefaultDevice().ok_or("no Metal device")?;
            let lib = dev
                .newLibraryWithSource_options_error(&NSString::from_str(&kmer_src()), None)
                .map_err(describe)?;
            let queue = dev.newCommandQueue().ok_or("no Metal command queue")?;
            let (kstart, useq, prm_b) = (
                buffer(&dev, kstart)?,
                buffer(&dev, useq)?,
                buffer(&dev, prm)?,
            );
            let keys = if keys.is_empty() {
                build_table(&dev, &queue, &lib, &kstart, &useq, &prm_b, prm)?
            } else {
                buffer(&dev, keys)?
            };
            Ok(Self {
                queue,
                pso: pipeline(&dev, &lib, "kmer_hits")?,
                keys,
                vals: buffer(&dev, vals)?,
                kstart,
                useq,
                prm: prm_b,
                dev,
            })
        })
    }

    /// See `KmerGpu::hits`; `hits` and `counts` are already sized.
    pub(crate) fn hits(
        &self,
        seq: &[u8],
        off: &[u32],
        hits: &mut [u32],
        counts: &mut [u32],
    ) -> Result<(), String> {
        autoreleasepool(|_| {
            let n = counts.len() as u32;
            let seq_b = buffer(&self.dev, seq)?;
            let off_b = buffer(&self.dev, off)?;
            let hits_b = buffer(&self.dev, &*hits)?;
            let counts_b = buffer(&self.dev, &*counts)?;
            let cb = self
                .queue
                .commandBuffer()
                .ok_or("no Metal command buffer")?;
            let enc = cb
                .computeCommandEncoder()
                .ok_or("no Metal compute encoder")?;
            let bufs = [
                &seq_b,
                &off_b,
                &self.keys,
                &self.vals,
                &self.kstart,
                &self.useq,
                &self.prm,
                &hits_b,
                &counts_b,
            ];
            dispatch(&enc, &self.pso, &bufs, n);
            enc.endEncoding();
            cb.commit();
            cb.waitUntilCompleted();
            if cb.status() != MTLCommandBufferStatus::Completed {
                return Err(cb
                    .error()
                    .map(describe)
                    .unwrap_or_else(|| "Metal command buffer failed".to_string()));
            }
            let copy = |b: &Buf, out: &mut [u32]| {
                // SAFETY: a shared-storage buffer made from `out`, so of its
                // length, and the GPU work on it has completed.
                let src = unsafe {
                    std::slice::from_raw_parts(b.contents().as_ptr().cast::<u32>(), out.len())
                };
                out.copy_from_slice(src);
            };
            copy(&hits_b, hits);
            copy(&counts_b, counts);
            Ok(())
        })
    }
}
