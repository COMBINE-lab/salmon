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
use std::ptr::{null, null_mut};

use libloading::Library;

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

fn open(names: &[&str]) -> Result<Library, String> {
    let mut last = String::new();
    for name in names {
        // SAFETY: a ROCm system library, trusted like one linked at build time.
        match unsafe { Library::new(name) } {
            Ok(lib) => return Ok(lib),
            Err(e) => last = e.to_string(),
        }
    }
    Err(format!("cannot load {} ({last})", names[0]))
}

/// `*lib.get(name)`: the symbol must have the C signature `T` in the ROCm headers.
macro_rules! sym {
    ($lib:expr, $name:literal) => {
        // SAFETY: `$name` has this signature in ROCm's hip_runtime_api.h / hiprtc.h,
        // and the library outlives every use (it is kept beside the pointer).
        *unsafe { $lib.get(concat!($name, "\0").as_bytes()) }
            .map_err(|e| format!("{}: {e}", $name))?
    };
}

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

    /// Compile `SRC` for the current device and load it as the module.
    fn load_module(&mut self) -> Result<(), String> {
        // ROCm 7 moved hiprtc to its own library; earlier ones export it from libamdhip64.
        let own = open(&[
            "libhiprtc.so",
            "libhiprtc.so.7",
            "libhiprtc.so.6",
            "/opt/rocm/lib/libhiprtc.so",
        ])
        .ok();
        let rtc = own.as_ref().unwrap_or(&self.lib);
        let create: unsafe extern "C" fn(
            *mut Ptr,
            *const c_char,
            *const c_char,
            c_int,
            *const *const c_char,
            *const *const c_char,
        ) -> c_int = sym!(rtc, "hiprtcCreateProgram");
        let compile: unsafe extern "C" fn(Ptr, c_int, *const *const c_char) -> c_int =
            sym!(rtc, "hiprtcCompileProgram");
        let size_of: unsafe extern "C" fn(Ptr, *mut usize) -> c_int =
            sym!(rtc, "hiprtcGetCodeSize");
        let code_of: unsafe extern "C" fn(Ptr, *mut c_char) -> c_int = sym!(rtc, "hiprtcGetCode");
        let log_size: unsafe extern "C" fn(Ptr, *mut usize) -> c_int =
            sym!(rtc, "hiprtcGetProgramLogSize");
        let log_of: unsafe extern "C" fn(Ptr, *mut c_char) -> c_int =
            sym!(rtc, "hiprtcGetProgramLog");
        let destroy: unsafe extern "C" fn(*mut Ptr) -> c_int = sym!(rtc, "hiprtcDestroyProgram");
        let load: unsafe extern "C" fn(*mut Ptr, *const c_void) -> c_int =
            sym!(self.lib, "hipModuleLoadData");
        fn rtc_ok(code: c_int) -> Result<(), String> {
            match code {
                0 => Ok(()),
                _ => Err(format!("hiprtc error {code}")),
            }
        }

        let mut prog = null_mut();
        // SAFETY: NUL-terminated source and name, no headers; `prog` is
        // destroyed on every path below.
        rtc_ok(unsafe {
            create(
                &mut prog,
                SRC.as_ptr(),
                c"salmon_em.hip".as_ptr(),
                0,
                null(),
                null(),
            )
        })?;
        // No --offload-arch: hiprtc targets the current device. Default (IEEE)
        // math, not -ffast-math, as on Metal.
        let opts = [c"-O3".as_ptr()];
        // SAFETY: `prog` is live; the buffers are sized from the size queries.
        let code = unsafe {
            let built = if compile(prog, 1, opts.as_ptr()) != 0 {
                let mut n = 0;
                log_size(prog, &mut n);
                let mut log = vec![0u8; n.max(1)];
                log_of(prog, log.as_mut_ptr().cast());
                Err(format!(
                    "hiprtc compile failed: {}",
                    String::from_utf8_lossy(&log).trim_end_matches('\0')
                ))
            } else {
                let mut n = 0;
                rtc_ok(size_of(prog, &mut n)).and_then(|()| {
                    let mut code = vec![0u8; n];
                    rtc_ok(code_of(prog, code.as_mut_ptr().cast())).map(|()| code)
                })
            };
            destroy(&mut prog);
            built?
        };
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
        // SAFETY: `p` holds at least `n` f32s and `out` has room for them.
        let copied = unsafe {
            (self.memcpy)(
                out.as_mut_ptr().cast(),
                p,
                std::mem::size_of_val(out.as_slice()),
                DEVICE_TO_HOST,
            )
        };
        self.check(copied, "hipMemcpy")?;
        Ok(out)
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
pub(super) fn available() -> bool {
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
        hip.load_module()?;
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
