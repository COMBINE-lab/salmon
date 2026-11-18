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
use crate::kmers::{HITS_CAP, SEL_CAP};

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

/// A shared buffer of `len` `u32`s the GPU fills (nothing copied in).
fn out_buffer(dev: &ProtocolObject<dyn MTLDevice>, len: usize) -> Result<Buf, String> {
    dev.newBufferWithLength_options((len * 4).max(4), MTLResourceOptions::StorageModeShared)
        .ok_or_else(|| "Metal buffer allocation failed".to_string())
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
    dispatch_lanes(enc, pso, bufs, n, 1);
}

/// [`dispatch`] with `lanes` threads per item, in threadgroups of a multiple
/// of `lanes` (a power of two up to 32).
fn dispatch_lanes(
    enc: &ProtocolObject<dyn MTLComputeCommandEncoder>,
    pso: &Pso,
    bufs: &[&Buf],
    n: u32,
    lanes: usize,
) {
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
    let width = pso.maxTotalThreadsPerThreadgroup().min(256) & !(lanes - 1);
    enc.dispatchThreads_threadsPerThreadgroup(
        MTLSize {
            width: n as usize * lanes,
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
// `select_ref` (`kmers.rs`), after `kmer_hits`: one SIMD group per read. Each
// lane takes a target of the read's shortest occurrence list, binary-searches
// every hit's list for it (a hit on the previous hit's unitig reuses its
// entries) and applies them in hit order; a prefix sum over the lanes writes
// the survivors in target order. `occ_off`/`occ_ent` are `Occurrences`, `sp`:
// ref_len_bits, max_occ, max_read_occ.
#define SEL_CAP 128u
#define SEL_CPU 0xFFFFFFFFu
#define SEL_OVER 0xFFFFFFFEu
struct Tgt { int lf; int lr; int af; int ar; u32 fh; u32 rh; };
static inline void sel_apply(thread Tgt& t, device const u32* h, u64 e, u64 pmask, u32 ulen, int rl, int k) {
    bool ori = (h[1] & 1u) != 0u, cfw = (e & 1) != 0;
    u32 cpos = h[2], bp = (u32)((e >> 1) & pmask);
    int rp = (int)h[0];
    int rpos = (int)(cfw ? bp + cpos : bp + ulen - (cpos + (u32)k));
    if (cfw == ori) {
        if (rp > t.lf) { if (t.lf == -1) t.af = rpos - rp; t.lf = rp; t.fh += 1; }
    } else if (rp > t.lr) {
        t.ar = rpos - (rl - (rp + k)); t.rh += 1; t.lr = rp;
    }
}
kernel void kmer_select(device const u32* off [[buffer(0)]],
                        device const u32* hits [[buffer(1)]],
                        device const u32* counts [[buffer(2)]],
                        device const u64* keys [[buffer(3)]],
                        device const u32* vals [[buffer(4)]],
                        device const u32* kstart [[buffer(5)]],
                        device const u64* useq [[buffer(6)]],
                        device const u32* prm [[buffer(7)]],
                        device const u32* occ_off [[buffer(8)]],
                        device const u64* occ_ent [[buffer(9)]],
                        device const u32* sp [[buffer(10)]],
                        device u32* sel [[buffer(11)]],
                        device u32* nsel [[buffer(12)]],
                        constant uint& nr [[buffer(13)]],
                        uint gid [[thread_position_in_grid]]) {
    // the dispatch is 32 threads per read in threadgroups of a multiple of 32
    u32 i = gid >> 5, lane = gid & 31u;
    if (i >= nr) return;
    u32 n = counts[i], hcap = prm[3];
    if (n > hcap || n == 0) { if (lane == 0u) nsel[i] = n == 0u ? 0u : SEL_CPU; return; }
    device const u32* h = hits + (u64)i * hcap * 3;
    u32 shift = sp[0] + 1, max_occ = sp[1], max_read_occ = sp[2];
    u64 pmask = (((u64)1) << sp[0]) - 1;
    int k = (int)prm[0], rl = (int)(off[i + 1] - off[i]);
    u32 d = 0u, dlen = 0xFFFFFFFFu;
    for (u32 j = 0; j < n; ++j) {
        u32 u = h[3 * j + 1] >> 1, len = occ_off[u + 1] - occ_off[u];
        if (len > max_occ) { if (lane == 0u) nsel[i] = SEL_CPU; return; }
        if (len < dlen) { dlen = len; d = j; }
    }
    u32 ud = h[3 * d + 1] >> 1, s0 = occ_off[ud], e0 = occ_off[ud + 1];
    device u32* out = sel + (u64)i * SEL_CAP * 3;
    u32 cnt = 0u;
    for (u32 base = s0; base < e0 && cnt <= max_read_occ; base += 32u) {
        u32 c = base + lane, nout = 0u, t = 0u;
        Tgt st; st.lf = -1; st.lr = -1; st.af = -1; st.ar = -1; st.fh = 0u; st.rh = 0u;
        if (c < e0) t = (u32)(occ_ent[c] >> shift);
        if (c < e0 && (c == s0 || (u32)(occ_ent[c - 1] >> shift) != t)) {
            u32 pu2 = 0xFFFFFFFFu, a = 0u, b = 0u, len = 0u;
            bool alive = true;
            for (u32 j = 0; j < n && alive; ++j) {
                u32 u = h[3 * j + 1] >> 1;
                if (u != pu2) {
                    u32 lo = occ_off[u], hi = occ_off[u + 1], end = hi;
                    while (lo < hi) {
                        u32 mid = (lo + hi) >> 1;
                        if ((u32)(occ_ent[mid] >> shift) < t) lo = mid + 1; else hi = mid;
                    }
                    a = lo; b = lo;
                    while (b < end && (u32)(occ_ent[b] >> shift) == t) b += 1;
                    len = (u32)(unitig_base(u + 1, keys, vals, kstart, useq, prm) - unitig_base(u, keys, vals, kstart, useq, prm));
                    pu2 = u;
                }
                for (u32 e = a; e < b; ++e)
                    if ((st.fh > st.rh ? st.fh : st.rh) >= j) sel_apply(st, h + 3 * j, occ_ent[e], pmask, len, rl, k);
                alive = (st.fh > st.rh ? st.fh : st.rh) > j;
            }
            // piscem `push_if_accepted`, with every hit valid
            if (alive) nout = st.fh > st.rh ? (st.fh >= n ? 1u : 0u)
                            : st.fh < st.rh ? (st.rh >= n ? 1u : 0u)
                            : (st.fh >= n ? 2u : 0u);
        }
        u32 pos = cnt + simd_prefix_exclusive_sum(nout);
        bool fw_first = st.fh >= st.rh;
        for (u32 q = 0; q < nout; ++q, ++pos) {
            if (pos >= SEL_CAP) break;
            bool fw = (q == 0u) == fw_first;
            out[3 * pos] = t << 1 | (fw ? 1u : 0u);
            out[3 * pos + 1] = (u32)(fw ? st.af : st.ar);
            out[3 * pos + 2] = fw ? st.fh : st.rh;
        }
        cnt += simd_sum(nout);
    }
    if (lane == 0u) nsel[i] = cnt > max_read_occ ? SEL_OVER : (cnt > SEL_CAP ? SEL_CPU : cnt);
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
    sel_pso: Pso,
    /// `set_occurrences`: the occurrence offsets, entries and parameters
    occ: Option<[Buf; 3]>,
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
                sel_pso: pipeline(&dev, &lib, "kmer_select")?,
                occ: None,
                keys,
                vals: buffer(&dev, vals)?,
                kstart,
                useq,
                prm: prm_b,
                dev,
            })
        })
    }

    /// See `KmerGpu::device_bytes`.
    pub(crate) fn device_bytes(&self) -> Option<u64> {
        Some(self.dev.currentAllocatedSize() as u64)
    }

    /// See `KmerGpu::set_occurrences`.
    pub(crate) fn set_occurrences(
        &mut self,
        off: &[u32],
        ent: &[u64],
        prm: &[u32; 3],
    ) -> Result<(), String> {
        // `kmer_select` takes each SIMD group as one read's 32 lanes
        if self.sel_pso.threadExecutionWidth() != 32 {
            return Err("the selection kernel needs 32-wide SIMD groups".into());
        }
        autoreleasepool(|_| {
            self.occ = Some([
                buffer(&self.dev, off)?,
                buffer(&self.dev, ent)?,
                buffer(&self.dev, prm)?,
            ]);
            Ok(())
        })
    }

    /// See `KmerGpu::hits` and `KmerGpu::hits_select` (`sel`); every output
    /// is already sized.
    pub(crate) fn run(
        &self,
        seq: &[u8],
        off: &[u32],
        hits: &mut [u32],
        counts: &mut [u32],
        mut sel: Option<(&mut [u32], &mut [u32])>,
    ) -> Result<(), String> {
        autoreleasepool(|_| {
            let n = counts.len() as u32;
            let seq_b = buffer(&self.dev, seq)?;
            let off_b = buffer(&self.dev, off)?;
            let hits_b = out_buffer(&self.dev, hits.len())?;
            let counts_b = out_buffer(&self.dev, counts.len())?;
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
            // the default (serial) encoder runs this after the lookup
            let sel_bufs = match (&sel, &self.occ) {
                (Some((s, ns)), Some([occ_off, occ_ent, sp])) => {
                    let (s_b, ns_b) = (
                        out_buffer(&self.dev, s.len())?,
                        out_buffer(&self.dev, ns.len())?,
                    );
                    let bufs = [
                        &off_b,
                        &hits_b,
                        &counts_b,
                        &self.keys,
                        &self.vals,
                        &self.kstart,
                        &self.useq,
                        &self.prm,
                        occ_off,
                        occ_ent,
                        sp,
                        &s_b,
                        &ns_b,
                    ];
                    dispatch_lanes(&enc, &self.sel_pso, &bufs, n, 32);
                    Some((s_b, ns_b))
                }
                (Some(_), None) => return Err("no occurrence lists on the device".into()),
                _ => None,
            };
            enc.endEncoding();
            cb.commit();
            cb.waitUntilCompleted();
            if cb.status() != MTLCommandBufferStatus::Completed {
                return Err(cb
                    .error()
                    .map(describe)
                    .unwrap_or_else(|| "Metal command buffer failed".to_string()));
            }
            // SAFETY: a shared-storage buffer of `out`'s length (`out_buffer`),
            // and the GPU work on it has completed.
            let src = |b: &Buf, out: &[u32]| unsafe {
                std::slice::from_raw_parts(b.contents().as_ptr().cast::<u32>(), out.len())
            };
            // each read's written prefix only (`w` u32s per entry, `cap`
            // entries per read): the rest of a row is never read
            let rows = |b: &Buf, out: &mut [u32], cap: usize, w: usize, len: &[u32]| {
                let src = src(b, out);
                for (i, &c) in len.iter().enumerate() {
                    let r = i * cap * w..i * cap * w + (c as usize).min(cap) * w;
                    out[r.clone()].copy_from_slice(&src[r]);
                }
            };
            counts.copy_from_slice(src(&counts_b, counts));
            rows(&hits_b, hits, HITS_CAP, 3, counts);
            if let (Some((s, ns)), Some((s_b, ns_b))) = (sel.as_mut(), sel_bufs) {
                ns.copy_from_slice(src(&ns_b, ns));
                // `SEL_CPU`/`SEL_OVER` exceed the cap: nothing to copy
                let ns: &[u32] = ns;
                let keep: Vec<u32> = ns
                    .iter()
                    .map(|&c| if c as usize > SEL_CAP { 0 } else { c })
                    .collect();
                rows(&s_b, s, SEL_CAP, 3, &keep);
            }
            Ok(())
        })
    }
}
