//! Why does reusing the sketch MappingCache (vs allocating one per read) not
//! change wall time? Measure the cost of `MappingCache::new` (alloc+free of an
//! AHashMap cap-256 + Vec cap-64) vs `clear()` (reuse), under whichever global
//! allocator the binary is built with. Run with the default (mimalloc) and with
//! `--features sysalloc` to contrast.
//!
//! Note what this does *not* measure: it allocates and frees on one thread,
//! which is every allocator's fast path and not what the quant hot path does.
//! Its ranking does not survive contact with a real run; see
//! `docs/allocator-choice.md`.

// Shared with the real binary so the two can never disagree about which
// allocator a given set of feature flags selects.
#[path = "../src/global_alloc.rs"]
mod global_alloc;

use piscem_rs::mapping::cache::MappingCache;
use piscem_rs::mapping::sketch_hit_simple::SketchHitInfoSimple;
use std::hint::black_box;
use std::time::Instant;

fn main() {
    let alloc_name = global_alloc::NAME;
    // Enough iterations that per-op timing is not dominated by clock resolution.
    let n: u64 = 30_000_000;

    // alloc + free, every iteration (what the current sketch path does per read)
    let t = Instant::now();
    for _ in 0..n {
        let c = MappingCache::<SketchHitInfoSimple>::new(31);
        // `black_box` stops the optimizer deleting the allocation as dead code,
        // which is the entire thing being measured.
        black_box(&c);
    }
    let new_ns = t.elapsed().as_nanos() as f64 / n as f64;

    // reuse: clear an existing cache (what SketchScratch would do)
    let mut c = MappingCache::<SketchHitInfoSimple>::new(31);
    let t = Instant::now();
    for _ in 0..n {
        c.clear();
        black_box(&c);
    }
    let clear_ns = t.elapsed().as_nanos() as f64 / n as f64;

    println!("[{alloc_name}] MappingCache::new = {new_ns:.1} ns/op   clear(reuse) = {clear_ns:.1} ns/op   alloc overhead = {:.1} ns/op", new_ns - clear_ns);
    println!(
        "  per read pair (left+right+out = 3 caches): {:.1} ns of alloc avoided by reuse",
        3.0 * (new_ns - clear_ns)
    );
}
