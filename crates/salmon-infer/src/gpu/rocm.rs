//! ROCm/HIP backend of the GPU warm start (AMD GPUs, Linux).
//!
//! The Metal kernels in HIP C++, compiled at run time for the current device
//! with hiprtc and launched with the HIP module API, one thread per chunk of
//! incidences, then one per class or transcript to sum the chunks. `libamdhip64` and `libhiprtc` are opened at run time, so the
//! feature builds without ROCm installed, and a machine without it logs a
//! warning and stays on the CPU.
//!
//! Untested: written without AMD hardware; CI only compiles it.
#![allow(unsafe_code)]

use std::ffi::{c_char, c_int, c_uint, c_void, CStr};
use std::ptr::null_mut;

use libloading::Library;

use super::dl::{open, rtc_compile, sym};
use super::Problem;

/// Same kernels, math and argument order as `gpu/metal.rs`; the count is last.
const SRC: &CStr = cr#"
// exp(psi(x)) for x > 0: f32 port of components.rs `exp_digamma`.
static __device__ float exp_digamma(float x) {
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
    return y * expf(u - s);
}

extern "C" __global__ void et_kernel(const float* alpha, const float* prior, float* et,
                                     unsigned n) {
    unsigned t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    float ap = alpha[t] + prior[t];
    et[t] = ap > 1e-10f ? exp_digamma(ap) : 0.0f;  // DIGAMMA_MIN
}

// partial[k] = sum of x[idx[e]] * w[e] over chunk k (both CSR directions)
extern "C" __global__ void chunk_dot(const unsigned* chunk, const unsigned* idx,
                                     const float* w, const float* x, float* partial,
                                     unsigned n) {
    unsigned k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= n) return;
    float acc = 0.0f;
    for (unsigned e = chunk[k]; e < chunk[k + 1]; ++e) {
        float b = x[idx[e]];
        acc += b > 0.0f ? b * w[e] : 0.0f;
    }
    partial[k] = acc;
}

extern "C" __global__ void class_finish(const unsigned* first, const float* partial,
                                        const float* count, float* scale, unsigned n) {
    unsigned c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= n) return;
    float denom = 0.0f;
    for (unsigned k = first[c]; k < first[c + 1]; ++k) denom += partial[k];
    scale[c] = denom > 1.17549435e-38f ? count[c] / denom : 0.0f;  // FLT_MIN
}

extern "C" __global__ void txp_finish(const unsigned* first, const float* partial,
                                      const float* basis, const float* single, float* out,
                                      unsigned n) {
    unsigned t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n) return;
    float acc = 0.0f;
    for (unsigned k = first[t]; k < first[t + 1]; ++k) acc += partial[k];
    out[t] = single[t] + basis[t] * acc;
}
"#;

/// Device pointer, `hipModule_t`, `hipFunction_t`, `hiprtcProgram`.
type Ptr = *mut c_void;
const THREADS: u32 = 256;
const HOST_TO_DEVICE: c_int = 1;
const DEVICE_TO_HOST: c_int = 2;

/// The HIP runtime entry points, the module and every device allocation;
/// dropping it frees them.
struct Hip {
    malloc: unsafe extern "C" fn(*mut Ptr, usize) -> c_int,
    free: unsafe extern "C" fn(Ptr) -> c_int,
    memcpy: unsafe extern "C" fn(Ptr, *const c_void, usize, c_int) -> c_int,
    launch: unsafe extern "C" fn(
        Ptr,
        c_uint,
        c_uint,
        c_uint,
        c_uint,
        c_uint,
        c_uint,
        c_uint,
        Ptr,
        *mut Ptr,
        *mut Ptr,
    ) -> c_int,
    sync: unsafe extern "C" fn() -> c_int,
    unload: unsafe extern "C" fn(Ptr) -> c_int,
    error_string: unsafe extern "C" fn(c_int) -> *const c_char,
    module: Ptr,
    allocs: Vec<Ptr>,
    devices: c_int,
    lib: Library,
}

impl Hip {
    fn load() -> Result<Self, String> {
        let lib = open(&[
            "libamdhip64.so",
            "libamdhip64.so.7",
            "libamdhip64.so.6",
            "/opt/rocm/lib/libamdhip64.so",
        ])?;
        let init: unsafe extern "C" fn(c_uint) -> c_int = sym!(lib, "hipInit");
        let count: unsafe extern "C" fn(*mut c_int) -> c_int = sym!(lib, "hipGetDeviceCount");
        let mut hip = Self {
            malloc: sym!(lib, "hipMalloc"),
            free: sym!(lib, "hipFree"),
            memcpy: sym!(lib, "hipMemcpy"),
            launch: sym!(lib, "hipModuleLaunchKernel"),
            sync: sym!(lib, "hipDeviceSynchronize"),
            unload: sym!(lib, "hipModuleUnload"),
            error_string: sym!(lib, "hipGetErrorString"),
            module: null_mut(),
            allocs: Vec::new(),
            devices: 0,
            lib,
        };
        let mut devices = 0;
        // SAFETY: plain C calls with valid out-pointers.
        unsafe {
            hip.check(init(0), "hipInit")?;
            hip.check(count(&mut devices), "hipGetDeviceCount")?;
        }
        hip.devices = devices;
        Ok(hip)
    }

    fn check(&self, code: c_int, what: &str) -> Result<(), String> {
        if code == 0 {
            return Ok(());
        }
        // SAFETY: hipGetErrorString returns a static C string (or null).
        let msg = unsafe { (self.error_string)(code) };
        let msg = if msg.is_null() {
            format!("error {code}")
        } else {
            // SAFETY: non-null static C string from the runtime.
            unsafe { CStr::from_ptr(msg) }
                .to_string_lossy()
                .into_owned()
        };
        Err(format!("{what}: {msg}"))
    }

    fn upload<T: Copy>(&mut self, data: &[T]) -> Result<Ptr, String> {
        let bytes = std::mem::size_of_val(data);
        let mut p = null_mut();
        // SAFETY: valid out-pointer; at least 4 bytes so an empty slice still
        // gets a real (never read) allocation.
        self.check(unsafe { (self.malloc)(&mut p, bytes.max(4)) }, "hipMalloc")?;
        self.allocs.push(p);
        // SAFETY: `p` holds at least `bytes` bytes and `data` is a live slice of them.
        let copied = unsafe { (self.memcpy)(p, data.as_ptr().cast(), bytes, HOST_TO_DEVICE) };
        self.check(copied, "hipMemcpy")?;
        Ok(p)
    }

    /// Compile `src` for the current device and load it as the module.
    fn load_module(&mut self, src: &CStr) -> Result<(), String> {
        // ROCm 7 moved hiprtc to its own library; earlier ones export it from libamdhip64.
        let own = open(&[
            "libhiprtc.so",
            "libhiprtc.so.7",
            "libhiprtc.so.6",
            "/opt/rocm/lib/libhiprtc.so",
        ])
        .ok();
        let rtc = own.as_ref().unwrap_or(&self.lib);
        // No --offload-arch: hiprtc targets the current device. Default (IEEE)
        // math, not -ffast-math, as on Metal.
        let code = rtc_compile(rtc, "hiprtc", "Code", src, &[c"-O3"])?;
        let load: unsafe extern "C" fn(*mut Ptr, *const c_void) -> c_int =
            sym!(self.lib, "hipModuleLoadData");
        let mut module = null_mut();
        // SAFETY: `code` is the code object hiprtc just produced.
        self.check(
            unsafe { load(&mut module, code.as_ptr().cast()) },
            "hipModuleLoadData",
        )?;
        self.module = module;
        Ok(())
    }

    fn function(&self, name: &CStr) -> Result<Ptr, String> {
        let get: unsafe extern "C" fn(*mut Ptr, Ptr, *const c_char) -> c_int =
            sym!(self.lib, "hipModuleGetFunction");
        let mut f = null_mut();
        // SAFETY: valid out-pointer, loaded module, NUL-terminated name.
        let found = unsafe { get(&mut f, self.module, name.as_ptr()) };
        self.check(found, "hipModuleGetFunction")?;
        Ok(f)
    }

    /// Launch `f` with the device pointers `bufs` then the count `n`, one thread
    /// per element, on the default stream (so launches run in order).
    fn dispatch(&self, f: Ptr, bufs: &[Ptr], n: u32) -> Result<(), String> {
        if n == 0 {
            return Ok(());
        }
        let mut args: Vec<Ptr> = bufs.iter().map(|b| b as *const Ptr as Ptr).collect();
        args.push(&n as *const u32 as Ptr);
        // SAFETY: `args` points at one live value per kernel parameter, in the
        // kernel's order and types (see SRC); HIP copies them at launch.
        let launched = unsafe {
            (self.launch)(
                f,
                n.div_ceil(THREADS),
                1,
                1,
                THREADS,
                1,
                1,
                0,
                null_mut(),
                args.as_mut_ptr(),
                null_mut(),
            )
        };
        self.check(launched, "hipModuleLaunchKernel")
    }

    fn download(&self, p: Ptr, n: usize) -> Result<Vec<f32>, String> {
        let mut out = vec![0f32; n];
        self.download_into(p, &mut out)?;
        Ok(out)
    }

    fn download_into<T: Copy>(&self, p: Ptr, out: &mut [T]) -> Result<(), String> {
        // SAFETY: `p` holds at least `out.len()` `T`s and `out` has room for them.
        let copied = unsafe {
            (self.memcpy)(
                out.as_mut_ptr().cast(),
                p,
                std::mem::size_of_val(out),
                DEVICE_TO_HOST,
            )
        };
        self.check(copied, "hipMemcpy")
    }

    /// Free the allocations made since there were `from`.
    fn release(&mut self, from: usize) {
        for p in self.allocs.drain(from..) {
            // SAFETY: from hipMalloc, freed once (drained from the list).
            unsafe { (self.free)(p) };
        }
    }
}

impl Drop for Hip {
    fn drop(&mut self) {
        // SAFETY: each pointer came from hipMalloc / hipModuleLoadData once;
        // hipFree waits for pending work on the device.
        unsafe {
            for &p in &self.allocs {
                (self.free)(p);
            }
            if !self.module.is_null() {
                (self.unload)(self.module);
            }
        }
    }
}

/// Whether this machine has the HIP runtime and an AMD GPU (else the CPU runs; tests skip).
pub(crate) fn available() -> bool {
    Hip::load().is_ok_and(|h| h.devices > 0)
}

pub(super) struct Device {
    /// `Some` under VBEM
    et: Option<Ptr>,
    dot: Ptr,
    class_finish: Ptr,
    txp_finish: Ptr,
    c_chunk: Ptr,
    c_first: Ptr,
    labels: Ptr,
    w: Ptr,
    count: Ptr,
    t_chunk: Ptr,
    t_first: Ptr,
    t_cls: Ptr,
    t_w: Ptr,
    single: Ptr,
    prior: Ptr,
    basis: Ptr,
    class_scale: Ptr,
    /// per-chunk sums, reused by the class pass and then the transcript pass
    partial: Ptr,
    /// ping-pong iterates; `alpha[cur]` is the latest
    alpha: [Ptr; 2],
    cur: usize,
    nc: u32,
    nt: u32,
    /// chunks on the class and on the transcript side
    ncc: u32,
    ntc: u32,
    hip: Hip,
}

impl Device {
    pub(super) fn new(
        prob: &Problem,
        init: &[f64],
        prior: &[f64],
        vb: bool,
    ) -> Result<Self, String> {
        let mut hip = Hip::load()?;
        if hip.devices == 0 {
            return Err("no ROCm device".into());
        }
        hip.load_module(SRC)?;
        let f32s = |v: &[f64]| v.iter().map(|&x| x as f32).collect::<Vec<f32>>();
        let nt = prob.num_txps();
        let zeros_t = vec![0f32; nt];
        let (ncc, ntc) = (prob.c_chunk.len() - 1, prob.t_chunk.len() - 1);
        Ok(Self {
            et: if vb {
                Some(hip.function(c"et_kernel")?)
            } else {
                None
            },
            dot: hip.function(c"chunk_dot")?,
            class_finish: hip.function(c"class_finish")?,
            txp_finish: hip.function(c"txp_finish")?,
            c_chunk: hip.upload(&prob.c_chunk)?,
            c_first: hip.upload(&prob.c_first)?,
            labels: hip.upload(&prob.labels)?,
            w: hip.upload(&prob.w)?,
            count: hip.upload(&prob.count)?,
            t_chunk: hip.upload(&prob.t_chunk)?,
            t_first: hip.upload(&prob.t_first)?,
            t_cls: hip.upload(&prob.t_cls)?,
            t_w: hip.upload(&prob.t_w)?,
            single: hip.upload(&prob.single)?,
            prior: hip.upload(&f32s(prior))?,
            basis: hip.upload(&zeros_t)?,
            class_scale: hip.upload(&vec![0f32; prob.num_classes()])?,
            partial: hip.upload(&vec![0f32; ncc.max(ntc)])?,
            alpha: [hip.upload(&f32s(init))?, hip.upload(&zeros_t)?],
            cur: 0,
            nc: prob.num_classes() as u32,
            nt: nt as u32,
            ncc: ncc as u32,
            ntc: ntc as u32,
            hip,
        })
    }

    /// Run `n` M-steps and wait for them.
    pub(super) fn run(&mut self, n: u32) -> Result<(), String> {
        for _ in 0..n {
            let (src, dst) = (self.alpha[self.cur], self.alpha[1 - self.cur]);
            let basis = match self.et {
                Some(et) => {
                    self.hip
                        .dispatch(et, &[src, self.prior, self.basis], self.nt)?;
                    self.basis
                }
                None => src,
            };
            let (hip, p) = (&self.hip, self.partial);
            let class_dot = [self.c_chunk, self.labels, self.w, basis, p];
            hip.dispatch(self.dot, &class_dot, self.ncc)?;
            let class_fin = [self.c_first, p, self.count, self.class_scale];
            hip.dispatch(self.class_finish, &class_fin, self.nc)?;
            let txp_dot = [self.t_chunk, self.t_cls, self.t_w, self.class_scale, p];
            hip.dispatch(self.dot, &txp_dot, self.ntc)?;
            let txp_fin = [self.t_first, p, basis, self.single, dst];
            hip.dispatch(self.txp_finish, &txp_fin, self.nt)?;
            self.cur = 1 - self.cur;
        }
        // SAFETY: plain C call.
        let synced = unsafe { (self.hip.sync)() };
        self.hip.check(synced, "hipDeviceSynchronize")
    }

    /// `(previous, latest)` iterates.
    pub(super) fn read(&self) -> Result<(Vec<f32>, Vec<f32>), String> {
        let nt = self.nt as usize;
        Ok((
            self.hip.download(self.alpha[1 - self.cur], nt)?,
            self.hip.download(self.alpha[self.cur], nt)?,
        ))
    }
}

/// The k-mer table on the device (`kmers.rs`).
pub(crate) struct KmerDevice {
    /// one caller at a time: per-call buffers live in `Hip::allocs`
    // shortcut: calls from the mapping threads serialize here (and on the null
    // stream); give each a stream and its own buffers if the GPU idles.
    hip: std::sync::Mutex<Hip>,
    f: Ptr,
    keys: Ptr,
    vals: Ptr,
    kstart: Ptr,
    prm: Ptr,
}

// SAFETY: the HIP runtime is thread-safe and its device pointers are valid
// from any host thread; `Hip` is only used under the mutex.
unsafe impl Send for KmerDevice {}
unsafe impl Sync for KmerDevice {}

impl KmerDevice {
    pub(crate) fn new(k: u32, keys: &[u64], vals: &[u32], kstart: &[u32]) -> Result<Self, String> {
        let mut hip = Hip::load()?;
        if hip.devices == 0 {
            return Err("no ROCm device".into());
        }
        hip.load_module(&super::dl::kmer_src())?;
        let prm = [
            k,
            keys.len() as u32,
            (kstart.len() - 1) as u32,
            crate::kmers::HITS_CAP as u32,
        ];
        Ok(Self {
            f: hip.function(c"kmer_hits")?,
            keys: hip.upload(keys)?,
            vals: hip.upload(vals)?,
            kstart: hip.upload(kstart)?,
            prm: hip.upload(&prm)?,
            hip: std::sync::Mutex::new(hip),
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
        let mut hip = self.hip.lock().map_err(|_| "ROCm device poisoned")?;
        let base = hip.allocs.len();
        let mut run = |hip: &mut Hip| -> Result<(), String> {
            let seq_d = hip.upload(seq)?;
            let off_d = hip.upload(off)?;
            let hits_d = hip.upload(&*hits)?;
            let counts_d = hip.upload(&*counts)?;
            let bufs = [
                seq_d,
                off_d,
                self.keys,
                self.vals,
                self.kstart,
                self.prm,
                hits_d,
                counts_d,
            ];
            hip.dispatch(self.f, &bufs, counts.len() as u32)?;
            // SAFETY: plain C call.
            let synced = unsafe { (hip.sync)() };
            hip.check(synced, "hipDeviceSynchronize")?;
            hip.download_into(hits_d, hits)?;
            hip.download_into(counts_d, counts)
        };
        let res = run(&mut hip);
        hip.release(base);
        res
    }
}
