//! Coarse per-phase wall-clock timing for profiling.
//!
//! A [`PhaseTimer`] logs the elapsed time of each phase on the dedicated
//! `salmon::timing` tracing target, so a breakdown (mapping/RAD-read vs.
//! inference vs. posterior vs. output) is available inline at the default
//! `info` level, or in isolation via `RUST_LOG=salmon::timing=info`.
//!
//! This is pure instrumentation — a handful of `Instant::now()` calls around
//! second-scale phases, emitting logs — it never touches quant output or
//! determinism. It is shared by every quantification driver (the reads path in
//! `salmon-quant`, and the RAD / alignment paths in `salmon-align`) so the same
//! phase breakdown is available regardless of input mode.

use std::time::Instant;

/// Emits one `salmon::timing` log line per [`mark`](PhaseTimer::mark), reporting
/// the wall-clock elapsed since the previous mark (or since construction).
pub struct PhaseTimer {
    last: Instant,
}

impl Default for PhaseTimer {
    fn default() -> Self {
        Self::new()
    }
}

impl PhaseTimer {
    /// Start timing; the first [`mark`](PhaseTimer::mark) reports time since now.
    pub fn new() -> Self {
        Self {
            last: Instant::now(),
        }
    }

    /// Log the elapsed time since the previous mark under `phase`, and reset the
    /// clock for the next phase.
    pub fn mark(&mut self, phase: &str) {
        let now = Instant::now();
        let elapsed_s = now.duration_since(self.last).as_secs_f64();
        tracing::info!(target: "salmon::timing", phase, elapsed_s, "phase complete");
        self.last = now;
    }
}
