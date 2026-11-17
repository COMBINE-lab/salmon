//! Metal backend of the GPU warm start (macOS).
//!
//! One thread per class or transcript. Buffers use shared
//! storage, so on Apple silicon's unified memory a read-back is a slice copy.
#![allow(unsafe_code)]

use std::ptr::NonNull;

use objc2::rc::{autoreleasepool, Retained};
use objc2::runtime::ProtocolObject;
use objc2_foundation::{NSError, NSString};
use objc2_metal::{
    MTLBuffer, MTLCommandBuffer, MTLCommandBufferStatus, MTLCommandEncoder, MTLCommandQueue,
    MTLCompileOptions, MTLComputeCommandEncoder, MTLComputePipelineState,
    MTLCreateSystemDefaultDevice, MTLDevice, MTLLibrary, MTLResourceOptions, MTLSize,
};

use super::Problem;

/// The three kernels of `gpu/mod.rs`. Buffer indices are the argument order
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

kernel void class_scale(device const uint* c_start [[buffer(0)]],
                        device const uint* labels [[buffer(1)]],
                        device const float* w [[buffer(2)]],
                        device const float* count [[buffer(3)]],
                        device const float* basis [[buffer(4)]],
                        device float* scale [[buffer(5)]],
                        constant uint& n [[buffer(6)]],
                        uint c [[thread_position_in_grid]]) {
    if (c >= n) return;
    float denom = 0.0f;
    for (uint e = c_start[c]; e < c_start[c + 1]; ++e) {
        float b = basis[labels[e]];
        denom += b > 0.0f ? b * w[e] : 0.0f;
    }
    scale[c] = denom > 1.17549435e-38f ? count[c] / denom : 0.0f;  // FLT_MIN
}

kernel void txp_gather(device const uint* t_start [[buffer(0)]],
                       device const uint* t_cls [[buffer(1)]],
                       device const float* t_w [[buffer(2)]],
                       device const float* scale [[buffer(3)]],
                       device const float* basis [[buffer(4)]],
                       device const float* single [[buffer(5)]],
                       device float* out [[buffer(6)]],
                       constant uint& n [[buffer(7)]],
                       uint t [[thread_position_in_grid]]) {
    if (t >= n) return;
    float acc = 0.0f;
    for (uint e = t_start[t]; e < t_start[t + 1]; ++e) {
        acc += t_w[e] * scale[t_cls[e]];
    }
    out[t] = single[t] + basis[t] * acc;
}
"#;

type Buf = Retained<ProtocolObject<dyn MTLBuffer>>;
type Pso = Retained<ProtocolObject<dyn MTLComputePipelineState>>;

/// Whether this machine exposes a Metal device at all (tests skip otherwise).
#[cfg(test)]
pub(super) fn available() -> bool {
    MTLCreateSystemDefaultDevice().is_some()
}

pub(super) struct Device {
    queue: Retained<ProtocolObject<dyn MTLCommandQueue>>,
    /// `Some` under VBEM
    et: Option<Pso>,
    scale: Pso,
    gather: Pso,
    c_start: Buf,
    labels: Buf,
    w: Buf,
    count: Buf,
    t_start: Buf,
    t_cls: Buf,
    t_w: Buf,
    single: Buf,
    prior: Buf,
    basis: Buf,
    class_scale: Buf,
    /// ping-pong iterates; `alpha[cur]` is the latest
    alpha: [Buf; 2],
    cur: usize,
    nc: u32,
    nt: u32,
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
            Ok(Self {
                queue: dev.newCommandQueue().ok_or("no Metal command queue")?,
                et: if vb {
                    Some(pipeline(&dev, &lib, "et_kernel")?)
                } else {
                    None
                },
                scale: pipeline(&dev, &lib, "class_scale")?,
                gather: pipeline(&dev, &lib, "txp_gather")?,
                c_start: buffer(&dev, &prob.c_start)?,
                labels: buffer(&dev, &prob.labels)?,
                w: buffer(&dev, &prob.w)?,
                count: buffer(&dev, &prob.count)?,
                t_start: buffer(&dev, &prob.t_start)?,
                t_cls: buffer(&dev, &prob.t_cls)?,
                t_w: buffer(&dev, &prob.t_w)?,
                single: buffer(&dev, &prob.single)?,
                prior: buffer(&dev, &f32s(prior))?,
                basis: buffer(&dev, &zeros_t)?,
                class_scale: buffer(&dev, &vec![0f32; prob.num_classes()])?,
                alpha: [buffer(&dev, &f32s(init))?, buffer(&dev, &zeros_t)?],
                cur: 0,
                nc: prob.num_classes() as u32,
                nt: nt as u32,
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
                let class_bufs = [
                    &self.c_start,
                    &self.labels,
                    &self.w,
                    &self.count,
                    basis,
                    &self.class_scale,
                ];
                dispatch(&enc, &self.scale, &class_bufs, self.nc);
                let txp_bufs = [
                    &self.t_start,
                    &self.t_cls,
                    &self.t_w,
                    &self.class_scale,
                    basis,
                    &self.single,
                    dst,
                ];
                dispatch(&enc, &self.gather, &txp_bufs, self.nt);
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
