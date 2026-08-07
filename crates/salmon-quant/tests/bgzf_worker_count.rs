//! Guards that BGZF compression actually runs on the number of workers we ask
//! for.
//!
//! This has silently broken twice on a noodles bump, and neither break was
//! catchable by output comparison, because the bytes are identical no matter how
//! many threads produced them:
//!
//! - noodles 0.48 deprecated `with_worker_count` and made it a no-op, dispatching
//!   every block onto rayon's *global* pool.
//! - noodles 0.50 restored it, but defined `MultithreadedWriter::new` as
//!   `with_worker_count(NonZero::MIN, ..)` — so the previously-correct `new`
//!   became a single compression thread.
//!
//! Both are silent throughput regressions. This test observes the thread count
//! directly so that a third variation fails loudly instead.
//!
//! **This must remain the only test in this file.** It measures process-wide
//! thread counts, so a concurrently running test in the same binary would make
//! the baseline drift.

#![cfg(target_os = "linux")]

use noodles_bgzf as bgzf;
use std::io::Write;
use std::num::NonZeroUsize;

fn thread_count() -> usize {
    std::fs::read_dir("/proc/self/task")
        .expect("reading /proc/self/task")
        .count()
}

#[test]
fn compression_runs_on_the_requested_worker_count() {
    const WORKERS: usize = 8;

    // Mirror salmon's real shape: a global rayon pool is already up and sized to
    // `-p` before any BAM writing starts. `build_global` is deliberately allowed
    // to fail — if something else in the process got there first, the baseline
    // below still accounts for it.
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build_global();
    std::thread::sleep(std::time::Duration::from_millis(300));

    let baseline = thread_count();
    let mut writer = bgzf::io::MultithreadedWriter::with_worker_count(
        NonZeroUsize::new(WORKERS).unwrap(),
        Vec::new(),
    );
    // Enough data that every worker has real work to pick up.
    writer
        .write_all(&vec![b'A'; 8 << 20])
        .expect("writing payload");
    let peak = thread_count();
    writer.finish().expect("finishing the BGZF stream");

    let spawned = peak.saturating_sub(baseline);
    assert!(
        spawned >= WORKERS,
        "asked for {WORKERS} BGZF compression workers but only {spawned} threads appeared \
         (baseline {baseline}, peak {peak}); the worker count is being ignored"
    );
}
