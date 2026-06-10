//! Lock-free `f64` accumulation, used for concurrent mass/count updates.

use crate::math::log_add;
use std::sync::atomic::{AtomicU64, Ordering};

/// An `f64` stored as its bit pattern in an `AtomicU64`, supporting lock-free
/// updates via a CAS loop. Salmon updates per-transcript mass and the fragment
/// length distribution concurrently from worker threads; this is the building
/// block for that.
#[derive(Debug)]
pub struct AtomicF64(AtomicU64);

impl AtomicF64 {
    pub fn new(v: f64) -> Self {
        Self(AtomicU64::new(v.to_bits()))
    }

    pub fn load(&self) -> f64 {
        f64::from_bits(self.0.load(Ordering::Relaxed))
    }

    pub fn store(&self, v: f64) {
        self.0.store(v.to_bits(), Ordering::Relaxed);
    }

    /// Atomically replace the stored value with `f(current)`, returning the new value.
    pub fn fetch_update<F: Fn(f64) -> f64>(&self, f: F) -> f64 {
        let mut cur = self.0.load(Ordering::Relaxed);
        loop {
            let next = f(f64::from_bits(cur)).to_bits();
            match self
                .0
                .compare_exchange_weak(cur, next, Ordering::Relaxed, Ordering::Relaxed)
            {
                Ok(_) => return f64::from_bits(next),
                Err(actual) => cur = actual,
            }
        }
    }

    /// Atomically accumulate `log_inc` into the stored log-space value:
    /// `self = log(exp(self) + exp(log_inc))`. Mirrors salmon's `incLoopLog`.
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
