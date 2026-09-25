//! Which global allocator the `salmon` binary is built with.
//!
//! The quant hot path is highly multithreaded and allocation-heavy (per-read
//! scratch buffers allocated on a parser thread and freed on a worker thread),
//! so the allocator is a real, measurable choice rather than a detail. mimalloc
//! is the default, and `docs/allocator-choice.md` records both the measurement
//! behind that default and what a candidate replacement would have to show
//! before it is worth carrying:
//!
//! | build flags           | allocator     |
//! |-----------------------|---------------|
//! | (none)                | mimalloc      |
//! | `--features sysalloc` | system malloc |
//!
//! [`NAME`] reports which one a build actually got, so a benchmark labels its
//! own numbers instead of trusting the flags it was invoked with.
//!
//! This file is also pulled into `examples/cache_bench.rs` via `#[path]`, so
//! the microbenchmark and the real binary can never disagree about which
//! allocator is in play.

/// Name of the allocator this binary was actually built with.
pub const NAME: &str = if cfg!(feature = "sysalloc") {
    "system"
} else {
    "mimalloc"
};

#[cfg(not(feature = "sysalloc"))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;
