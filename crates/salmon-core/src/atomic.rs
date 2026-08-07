//! Lock-free `f64` accumulation, used for concurrent mass/count updates.

use crate::math::log_add;
use std::sync::atomic::{AtomicU64, Ordering};

/// An `f64` stored as its bit pattern in an `AtomicU64`, supporting lock-free
/// updates via a CAS loop. Salmon updates per-transcript mass and the fragment
/// length distribution concurrently from worker threads; this is the building
/// block for that.
///
/// **Why the bit pattern.** CPUs provide atomic instructions for integers, not
/// for floating-point values, so there is no `AtomicF64` in the standard
/// library. `f64::to_bits` reinterprets the 64-bit float as a 64-bit integer
/// with no loss, `from_bits` reverses it, and all the atomicity happens on the
/// integer.
///
/// **What a CAS loop is.** "Compare-and-swap": read the current value, compute
/// the new one, then ask the CPU to store it *only if* the slot still holds what
/// we read. If another thread got there first the store fails, we take its value
/// and retry. This gives correct concurrent updates without a mutex, and threads
/// never block each other — a losing thread just loops once more.
#[derive(Debug)]
pub struct AtomicF64(AtomicU64);

impl AtomicF64 {
    /// Wrap an initial value.
    pub fn new(v: f64) -> Self {
        Self(AtomicU64::new(v.to_bits()))
    }

    /// Read the current value.
    ///
    /// `Ordering::Relaxed` is the cheapest memory ordering: it guarantees the
    /// read/write itself is atomic but imposes no ordering relative to other
    /// memory operations. That is what we want here, because these are
    /// independent accumulators, not flags that publish other data.
    pub fn load(&self) -> f64 {
        f64::from_bits(self.0.load(Ordering::Relaxed))
    }

    /// Overwrite the current value.
    pub fn store(&self, v: f64) {
        self.0.store(v.to_bits(), Ordering::Relaxed);
    }

    /// Atomically replace the stored value with `f(current)`, returning the new value.
    ///
    /// Note `f` may be called more than once when threads contend, so it must be
    /// a pure function of its argument (no side effects).
    pub fn fetch_update<F: Fn(f64) -> f64>(&self, f: F) -> f64 {
        let mut cur = self.0.load(Ordering::Relaxed);
        loop {
            let next = f(f64::from_bits(cur)).to_bits();
            // `_weak` is allowed to fail spuriously even when the comparison
            // would have succeeded; that is harmless inside a retry loop and
            // compiles to cheaper code on some architectures.
            match self
                .0
                .compare_exchange_weak(cur, next, Ordering::Relaxed, Ordering::Relaxed)
            {
                // Our value was stored: done.
                Ok(_) => return f64::from_bits(next),
                // Someone else won the race; `actual` is their value, so retry
                // the update on top of it rather than clobbering it.
                Err(actual) => cur = actual,
            }
        }
    }

    /// Atomically accumulate `log_inc` into the stored log-space value:
    /// `self = log(exp(self) + exp(log_inc))`. Mirrors salmon's `incLoopLog`.
    ///
    /// In other words: the slot holds the log of a total, and this adds one more
    /// term to that total without ever leaving log space (see [`crate::math`]).
    pub fn log_add_assign(&self, log_inc: f64) {
        self.fetch_update(|cur| log_add(cur, log_inc));
    }

    /// Atomically add `inc` in linear space. Mirrors salmon's `incLoop`.
    pub fn add_assign(&self, inc: f64) {
        self.fetch_update(|cur| cur + inc);
    }
}

impl Default for AtomicF64 {
    fn default() -> Self {
        Self::new(0.0)
    }
}
