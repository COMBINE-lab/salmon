//! Picks the GPU backend of the EM warm start from the enabled features and the
//! target, as `cfg(gpu)` plus `cfg(gpu = "metal" | "cuda" | "rocm")`: Metal only
//! on macOS, ROCm only on Linux, CUDA wherever it is asked for. Metal wins over
//! CUDA, CUDA over ROCm.
fn main() {
    println!("cargo::rerun-if-changed=build.rs");
    println!("cargo::rustc-check-cfg=cfg(gpu, values(none(), \"metal\", \"cuda\", \"rocm\"))");
    let feature = |f: &str| std::env::var_os(format!("CARGO_FEATURE_{f}")).is_some();
    let os = std::env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let backend = if feature("METAL") && os == "macos" {
        "metal"
    } else if feature("CUDA") {
        "cuda"
    } else if feature("ROCM") && os == "linux" {
        "rocm"
    } else {
        return;
    };
    println!("cargo::rustc-cfg=gpu");
    println!("cargo::rustc-cfg=gpu=\"{backend}\"");
}
