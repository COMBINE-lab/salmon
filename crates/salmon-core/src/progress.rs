//! Shared live-progress counters for the quantification drivers.

use std::sync::atomic::AtomicU64;

/// Lock-free counters a quant driver updates as it runs, so a UI (e.g. a CLI
/// progress bar) can poll mapping progress while quantification is in flight.
/// Shared by both the reads-mode (`salmon-quant`) and alignment-mode
/// (`salmon-align`) drivers.
#[derive(Debug, Default)]
pub struct ProgressCounters {
    /// fragments observed so far
    pub processed: AtomicU64,
    /// fragments mapped so far
    pub mapped: AtomicU64,
}
