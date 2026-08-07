//! Shared live-progress counters for the quantification drivers.

use std::sync::atomic::AtomicU64;

/// Lock-free counters a quant driver updates as it runs, so a UI (e.g. a CLI
/// progress bar) can poll mapping progress while quantification is in flight.
/// Shared by both the reads-mode (`salmon-quant`) and alignment-mode
/// (`salmon-align`) drivers.
///
/// **Why "atomic".** Dozens of mapping threads bump these counters constantly
/// while a reporting thread reads them. An `AtomicU64` can be incremented from
/// several threads at once without a lock and without the increments being lost;
/// a plain `u64` shared this way would be a data race (undefined behaviour in
/// Rust, and rejected by the compiler).
///
/// Readers may observe a slightly stale value — the counters are relaxed,
/// best-effort telemetry, not accounting — which is exactly what a progress
/// display needs and costs nothing on the hot path.
#[derive(Debug, Default)]
pub struct ProgressCounters {
    /// fragments observed so far
    pub processed: AtomicU64,
    /// fragments mapped so far
    pub mapped: AtomicU64,
}
