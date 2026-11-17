//! Run-time loading shared by the ROCm and CUDA backends: the system GPU
//! libraries are opened when first used, so the features build without them.
#![allow(unsafe_code)]

use std::ffi::{c_char, c_int, c_void, CStr};
use std::ptr::{null, null_mut};

use libloading::Library;

pub(super) fn open(names: &[&str]) -> Result<Library, String> {
    let mut last = String::new();
    for name in names {
        // SAFETY: a ROCm/CUDA system library, trusted like one linked at build time.
        match unsafe { Library::new(name) } {
            Ok(lib) => return Ok(lib),
            Err(e) => last = e.to_string(),
        }
    }
    Err(format!("cannot load {} ({last})", names[0]))
}

/// `*lib.get(name)`: the symbol must have the C signature `T` in the vendor headers.
#[cfg(gpu = "rocm")]
macro_rules! sym {
    ($lib:expr, $name:literal) => {
        // SAFETY: `$name` has this signature in the vendor's runtime / rtc header,
        // and the library outlives every use (it is kept beside the pointer).
        *unsafe { $lib.get(concat!($name, "\0").as_bytes()) }
            .map_err(|e| format!("{}: {e}", $name))?
    };
}
#[cfg(gpu = "rocm")]
pub(super) use sym;

/// `prefix + name` from `lib`, typed `T` by the caller.
fn get<T: Copy>(lib: &Library, prefix: &str, name: &str) -> Result<T, String> {
    let full = format!("{prefix}{name}\0");
    // SAFETY: every caller below types `T` as this hiprtc/NVRTC entry point's C
    // signature (both APIs declare the same ones); `lib` outlives the call.
    unsafe { lib.get::<T>(full.as_bytes()) }
        .map(|s| *s)
        .map_err(|e| format!("{}: {e}", full.trim_end_matches('\0')))
}

type Prog = *mut c_void;

/// Compile `src` with hiprtc (`prefix` "hiprtc", `code` "Code": a code object)
/// or NVRTC ("nvrtc", "PTX": PTX text, NUL-terminated), which share one C API.
pub(super) fn rtc_compile(
    lib: &Library,
    prefix: &str,
    code: &str,
    src: &CStr,
    opts: &[&CStr],
) -> Result<Vec<u8>, String> {
    let create: unsafe extern "C" fn(
        *mut Prog,
        *const c_char,
        *const c_char,
        c_int,
        *const *const c_char,
        *const *const c_char,
    ) -> c_int = get(lib, prefix, "CreateProgram")?;
    let compile: unsafe extern "C" fn(Prog, c_int, *const *const c_char) -> c_int =
        get(lib, prefix, "CompileProgram")?;
    let size_of: unsafe extern "C" fn(Prog, *mut usize) -> c_int =
        get(lib, prefix, &format!("Get{code}Size"))?;
    let code_of: unsafe extern "C" fn(Prog, *mut c_char) -> c_int =
        get(lib, prefix, &format!("Get{code}"))?;
    let log_size: unsafe extern "C" fn(Prog, *mut usize) -> c_int =
        get(lib, prefix, "GetProgramLogSize")?;
    let log_of: unsafe extern "C" fn(Prog, *mut c_char) -> c_int =
        get(lib, prefix, "GetProgramLog")?;
    let destroy: unsafe extern "C" fn(*mut Prog) -> c_int = get(lib, prefix, "DestroyProgram")?;
    let ok = |c: c_int| match c {
        0 => Ok(()),
        _ => Err(format!("{prefix} error {c}")),
    };

    let mut prog = null_mut();
    // SAFETY: NUL-terminated source and name, no headers; `prog` is
    // destroyed on every path below.
    ok(unsafe {
        create(
            &mut prog,
            src.as_ptr(),
            c"salmon.cu".as_ptr(),
            0,
            null(),
            null(),
        )
    })?;
    let opts: Vec<*const c_char> = opts.iter().map(|o| o.as_ptr()).collect();
    // SAFETY: `prog` is live; the buffers are sized from the size queries.
    unsafe {
        let built = if compile(prog, opts.len() as c_int, opts.as_ptr()) != 0 {
            let mut n = 0;
            log_size(prog, &mut n);
            let mut log = vec![0u8; n.max(1)];
            log_of(prog, log.as_mut_ptr().cast());
            Err(format!(
                "{prefix} compile failed: {}",
                String::from_utf8_lossy(&log).trim_end_matches('\0')
            ))
        } else {
            let mut n = 0;
            ok(size_of(prog, &mut n)).and_then(|()| {
                let mut out = vec![0u8; n];
                ok(code_of(prog, out.as_mut_ptr().cast())).map(|()| out)
            })
        };
        destroy(&mut prog);
        built
    }
}

/// [`crate::kmers::KMER_BODY`] as CUDA/HIP C++, with its `kmer_hits` entry.
pub(super) fn kmer_src() -> std::ffi::CString {
    let prelude = "typedef unsigned char u8; typedef unsigned int u32;\n\
        typedef unsigned long long u64; typedef long long i64;\n\
        #define G\n#define PRIV\n#define FN static __device__ inline\n";
    let entry = r#"
extern "C" __global__ void kmer_hits(const u8* seq, const u32* off, const u64* keys,
                                     const u32* vals, const u32* kstart, const u32* prm,
                                     u32* hits, u32* counts, unsigned n) {
    unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    map_mate(i, seq, off, keys, vals, kstart, prm, hits, counts);
}
"#;
    std::ffi::CString::new(format!("{prelude}{}{entry}", crate::kmers::KMER_BODY))
        .expect("kernel source has no NUL")
}
