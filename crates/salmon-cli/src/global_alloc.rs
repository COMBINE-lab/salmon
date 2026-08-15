//! Which global allocator the `salmon` binary is built with.
//!
//! The quant hot path is highly multithreaded and allocation-heavy (per-read
//! scratch buffers allocated on a parser thread and freed on a worker thread),
//! so the allocator is a real, measurable choice rather than a detail. mimalloc
//! is the default; the other three exist so that choice can be re-measured
//! rather than assumed:
//!
//! | build flags                       | allocator      |
//! |-----------------------------------|----------------|
//! | (none)                            | mimalloc       |
//! | `--features sysalloc`             | system malloc  |
//! | `--features jemalloc`             | jemalloc       |
//! | `--features snmalloc`             | snmalloc       |
//!
//! Cargo features are additive, so a build that somehow enables several of
//! these must still pick exactly one `#[global_allocator]`. The `cfg`s below
//! encode a fixed priority (sysalloc > snmalloc > jemalloc > mimalloc), so any
//! combination resolves to one static. [`NAME`] reports which one won, so a
//! benchmark labels its own numbers instead of trusting the flags it was
//! invoked with.
//!
//! This file is also pulled into `examples/cache_bench.rs` via `#[path]`, so
//! the microbenchmark and the real binary can never disagree about which
//! allocator is in play.

/// Name of the allocator this binary was actually built with.
pub const NAME: &str = if cfg!(feature = "sysalloc") {
    "system"
} else if cfg!(feature = "snmalloc") {
    "snmalloc"
} else if cfg!(feature = "jemalloc") {
    "jemalloc"
} else {
    "mimalloc"
};

#[cfg(all(
    not(feature = "sysalloc"),
    not(feature = "snmalloc"),
    not(feature = "jemalloc")
))]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[cfg(all(
    feature = "jemalloc",
    not(feature = "sysalloc"),
    not(feature = "snmalloc")
))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

#[cfg(all(feature = "snmalloc", not(feature = "sysalloc")))]
#[global_allocator]
static GLOBAL: snmalloc_rs::SnMalloc = snmalloc_rs::SnMalloc;
