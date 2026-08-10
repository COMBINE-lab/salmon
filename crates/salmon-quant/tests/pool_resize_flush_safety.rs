//! A worker retired mid-run must still commit everything observable outside it.
//!
//! # Why this test exists
//!
//! Before the thread broker, salmon's mapping pool was fixed: every worker ran
//! from the first record to EOF, so "commit my thread-local state at the end"
//! and "commit it when the input runs out" were the same moment.
//!
//! With a resizable pool they are not. The broker shrinks the mapping side to
//! give slots to the decoder, and a worker can retire at any batch boundary
//! while the run continues. Everything [`QuantProcessor`] holds per thread and
//! publishes elsewhere has to survive that:
//!
//!   * **counters** merged into shared atomics (`merge_counters`, `merge_bias`,
//!     `merge_unmapped`) — accumulated locally and committed once, at thread
//!     end. Note `merge_into` does *not* zero the local counters afterwards; it
//!     is safe only because the worker exits immediately after. A future
//!     mid-run merge would double-count, which is exactly why this is pinned by
//!     a test rather than a comment.
//!   * **output buffers** (`flush_sam`, `flush_bam`, `flush_rad`) — a retiring
//!     worker holds up to one batch of records that exist nowhere else. Dropping
//!     it silently truncates the output.
//!
//! paraseq commits both: its worker loop calls `on_thread_complete` after the
//! retirement break, not only at EOF. This test pins that guarantee from
//! salmon's side, so that a paraseq change which stopped flushing retired
//! workers fails *here* — loudly, in salmon's own suite — rather than showing up
//! as quietly missing reads in someone's quantification.
//!
//! The processor below deliberately mirrors `QuantProcessor`'s contract rather
//! than using it directly: same `Clone`-builds-fresh-state shape, same
//! accumulate-locally-commit-at-thread-end counters, same batch-flushed output
//! buffer — with none of the index and model setup a real one needs.

use std::collections::HashSet;
use std::io::Cursor;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use paraseq::fastx::{Collection, CollectionType, Reader, RefRecord};
use paraseq::parallel::ThreadPool;
use paraseq::prelude::*;

const RECORDS: usize = 250_000;
const MAX_THREADS: usize = 16;

fn input() -> Cursor<Vec<u8>> {
    let mut bytes = Vec::new();
    for i in 0..RECORDS {
        bytes.extend_from_slice(format!("@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n").as_bytes());
    }
    Cursor::new(bytes)
}

#[derive(Default)]
struct Shared {
    /// Where per-thread counters land, as salmon's `merge_into` does.
    counted: AtomicU64,
    /// Where per-thread output buffers land, as `flush_sam`/`flush_rad` do.
    flushed: Mutex<Vec<u32>>,
    /// Worker *incarnations* that ran to completion. A fixed pool produces at
    /// most one per thread; anything more means workers really did retire and
    /// respawn, which is what makes this test non-vacuous.
    completions: AtomicU64,
}

/// Same shape as `QuantProcessor`: shared refs plus per-thread state that is
/// published elsewhere only by an explicit commit.
struct Probe {
    shared: Arc<Shared>,
    /// Committed once, at thread end. Never reset — like salmon's counters.
    local_count: u64,
    /// Committed at every batch boundary *and* at thread end, like the output
    /// buffers. Whatever is still here when a worker retires must not be lost.
    local_out: Vec<u32>,
}

impl Probe {
    fn new(shared: Arc<Shared>) -> Self {
        Self {
            shared,
            local_count: 0,
            local_out: Vec::new(),
        }
    }

    fn flush_out(&mut self) {
        if !self.local_out.is_empty() {
            self.shared
                .flushed
                .lock()
                .unwrap()
                .extend(self.local_out.drain(..));
        }
    }
}

// A fresh worker must start from zeroed local state, or a retire/respawn cycle
// would re-commit whatever the template happened to be carrying. This mirrors
// `impl Clone for QuantProcessor`, which rebuilds from `shared`.
impl Clone for Probe {
    fn clone(&self) -> Self {
        Probe::new(Arc::clone(&self.shared))
    }
}

impl<'r> ParallelProcessor<RefRecord<'r>> for Probe {
    fn process_record(&mut self, record: RefRecord<'r>) -> paraseq::Result<()> {
        let id: u32 = std::str::from_utf8(record.id())
            .unwrap()
            .trim_start_matches('r')
            .parse()
            .unwrap();
        self.local_count += 1;
        self.local_out.push(id);
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        self.flush_out();
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        // Both halves, exactly as salmon does: counters into the shared atomic,
        // and any output the worker is still holding.
        self.shared
            .counted
            .fetch_add(self.local_count, Ordering::Relaxed);
        self.shared.completions.fetch_add(1, Ordering::Relaxed);
        self.flush_out();
        Ok(())
    }
}

fn run(churn: bool) -> (u64, Vec<u32>, u64) {
    let shared = Arc::new(Shared::default());
    let pool = ThreadPool::with_max(MAX_THREADS, MAX_THREADS);
    let done = Arc::new(AtomicBool::new(false));

    // Drive the pool up and down throughout the run, so workers retire and
    // respawn repeatedly while records are still being read.
    let resizer = churn.then(|| {
        let pool = pool.clone();
        let done = Arc::clone(&done);
        std::thread::spawn(move || {
            let sizes = [MAX_THREADS, 3, 11, 1, MAX_THREADS, 2, 7];
            let mut i = 0usize;
            while !done.load(Ordering::Acquire) {
                pool.set_threads(sizes[i % sizes.len()]);
                i += 1;
                std::thread::sleep(std::time::Duration::from_micros(300));
            }
        })
    });

    let mut proc = Probe::new(Arc::clone(&shared));
    let reader = Reader::new(input()).unwrap();
    let collection = Collection::new(vec![reader], CollectionType::Single).unwrap();
    collection
        .process_parallel_pool(&mut proc, &pool, None)
        .expect("mapping pass failed");

    done.store(true, Ordering::Release);
    if let Some(h) = resizer {
        h.join().unwrap();
    }

    let counted = shared.counted.load(Ordering::Relaxed);
    let flushed = shared.flushed.lock().unwrap().clone();
    let completions = shared.completions.load(Ordering::Relaxed);
    (counted, flushed, completions)
}

/// Control: with a fixed pool nothing retires early, so this is the baseline
/// that makes a failure under churn attributable to the resizing.
#[test]
fn fixed_pool_counts_and_flushes_everything() {
    let (counted, flushed, _) = run(false);
    assert_eq!(
        counted, RECORDS as u64,
        "counter lost or duplicated records"
    );
    assert_eq!(flushed.len(), RECORDS, "output buffer lost records");
}

/// The real check: workers retiring mid-run must commit both their counters and
/// their pending output.
#[test]
fn retired_workers_commit_counters_and_output() {
    let (counted, flushed, completions) = run(true);

    assert_eq!(
        counted, RECORDS as u64,
        "thread-local counters were lost or double-counted across a resize: a \
         retiring worker either skipped on_thread_complete or committed twice"
    );
    assert_eq!(
        flushed.len(),
        RECORDS,
        "a retiring worker dropped its pending output buffer: {} of {} records \
         reached the shared sink",
        flushed.len(),
        RECORDS
    );

    // Exactly-once, not merely the right total: loss and duplication in equal
    // measure would cancel out in a count but not in the set of ids.
    let unique: HashSet<u32> = flushed.iter().copied().collect();
    assert_eq!(
        unique.len(),
        RECORDS,
        "records were duplicated across a retire/respawn cycle"
    );
    for id in 0..RECORDS as u32 {
        assert!(unique.contains(&id), "record {id} never reached the sink");
    }

    // Guard against a vacuous pass. If the churn never actually retired anyone,
    // the assertions above hold trivially and this test proves nothing -- the
    // same shape of mistake that let a decoder verdict ship three times without
    // ever changing behaviour.
    assert!(
        completions > MAX_THREADS as u64,
        "only {completions} worker incarnations completed with a ceiling of \
         {MAX_THREADS}: no worker retired mid-run, so this test exercised \
         nothing and its green result is meaningless"
    );
}
