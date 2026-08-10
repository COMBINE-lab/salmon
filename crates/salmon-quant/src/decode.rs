//! Adaptive division of `-p` between gzip decoding and mapping.
//!
//! `-p` names execution slots, not mapping threads. Some of them have to inflate
//! the input before any of them can map it, and the right split depends on the
//! compression ratio, on how much work each fragment costs, and on the machine —
//! none of which are known before the run starts, and the second of which
//! differs by almost an order of magnitude between salmon's two mapping modes.
//!
//! Rather than guess, this wires up piscem's [`thread_broker`]: both sides
//! report the thread-nanoseconds they spend, and the controller solves
//!
//! ```text
//! d* = N · busy_producer / (busy_producer + busy_consumer)
//! ```
//!
//! for the number of decode slots. See `piscem-rs/crates/thread-broker/README.md`
//! for why it solves rather than hill-climbs (the starvation signal is one-sided,
//! so a climber walks to the end of the budget).
//!
//! # Engagement is a separate, earlier question
//!
//! Whether to run the parallel decoder *at all* is decided before any of that,
//! by [`EngagementPolicy`], and it is not free to get wrong. The serial decoder
//! inflates inline on the mapping threads, so it is work-conserving — a thread
//! that has just inflated a batch goes straight on to map it and is never idle.
//! A dedicated decode slot can only decode. Serial therefore gives `F` inflate
//! streams for free (`F` = gzip inputs) while parallel gives `d` but *takes* `d`
//! mapping threads, so parallel pays exactly when `d > F`.
//!
//! Substituting the solved `d*` gives the rule this module implements:
//!
//! ```text
//! d* > F   <=>   N > F · (1 + C/P)
//! ```
//!
//! so `min_threads_per_stream = 1 + C/P` — it scales linearly with consumer cost
//! per fragment. **This is why salmon cannot simply inherit piscem's 8.** That 8
//! encodes `C/P = 7` for piscem's sketch-like mapping; selective alignment does
//! substantially more work per fragment, which raises `C`, lowers the producer's
//! cost share, and pushes the threshold up. See [`ENGAGEMENT_SELECTIVE_ALIGNMENT`]
//! for what was measured.

use std::sync::Arc;

use anyhow::{Context, Result};
use thread_broker::{BusyMeter, Consumer, EngagementPolicy, ResizeError, Work};

/// How the engagement threshold varies across salmon's mapping configurations.
///
/// `min_threads_per_stream = 1 + C/P`, so every configuration that changes the
/// per-fragment consumer cost `C` moves the threshold, while `P` — inflating the
/// same bytes — stays put. Two knobs do that, giving four cells:
///
/// * **selective alignment vs sketch.** SA extends and scores alignments per
///   candidate; sketch does not. SA has the larger `C`, so it engages later.
/// * **`--deterministic`.** Its mapping pass does, per `record_discrete`, "no
///   online inference, eq-class assembly, bias collection, or prefix
///   library-type detection" — four categories of per-fragment work removed. It
///   lowers `C` in *both* modes, so it engages *earlier* than its counterpart.
///
/// The ordering is therefore
/// `sa < sa+deterministic` is **false**; it is
/// `sa+deterministic < sa` and `sketch+deterministic < sketch`, with both
/// sketch cells below their SA counterparts.
///
/// Caveat: `--deterministic` only takes the cheap path when the library is
/// paired (`deterministic_fld && is_paired()`), so single-end gets the ordinary
/// cost and the ordinary threshold.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Thresholds {
    pub sketch: usize,
    pub sketch_deterministic: usize,
    pub selective_alignment: usize,
    pub selective_alignment_deterministic: usize,
}

/// Measured on 26.1 M paired human RNA-seq fragments against a 238 k-transcript
/// index, by reading the broker's own solved `producer_cost_share` (the
/// threshold is `1/share`) rather than inferring it from a wall-time crossover.
///
/// | cell | share | threshold | status |
/// |---|---|---|---|
/// | sketch | 0.1004 ± 0.0056 | 10 | **measured** |
/// | sketch + deterministic | — | 8 | provisional |
/// | selective alignment | — | 29 | provisional |
/// | selective alignment + deterministic | — | 20 | provisional |
///
/// Only the first is measured so far. The other three were blocked on a pool
/// deadlock (`notes/pool-lost-wakeup-deadlock.md`, resolved in rapidgzip-core
/// 0.3.1) and can now be taken the same way: run each cell with `--decoder
/// parallel` and read the `thread broker cost model` log line. The provisional
/// values come from the `1 + C/P` relation under a 4x consumer-cost prior;
/// they are placeholders, and the tests below assert only their *ordering*,
/// which follows from the mechanism, not their magnitudes, which do not.
pub const THRESHOLDS: Thresholds = Thresholds {
    sketch: 10,
    sketch_deterministic: 8,
    selective_alignment: 29,
    selective_alignment_deterministic: 20,
};

/// The engagement policy salmon defaults to for a given configuration.
///
/// Public so the tests and any harness assert against the same values the
/// binary uses, rather than restating them.
pub fn default_engagement_policy(sketch: bool, deterministic: bool) -> EngagementPolicy {
    let t = THRESHOLDS;
    EngagementPolicy {
        min_threads_per_stream: match (sketch, deterministic) {
            (true, false) => t.sketch,
            (true, true) => t.sketch_deterministic,
            (false, false) => t.selective_alignment,
            (false, true) => t.selective_alignment_deterministic,
        },
    }
}

/// Resolve the engagement policy for a run.
///
/// A policy file, when given, wins outright — including over the per-mode
/// default, since the whole reason to write one is that the measured default
/// does not suit your data. Absent a file, the mode picks the default.
///
/// Reuses piscem's `ThreadPolicy` rather than defining a parallel format, so a
/// file written for one tool means the same thing to the other, and an unknown
/// field is a hard error in both rather than a silent no-op.
pub fn engagement_policy(
    sketch: bool,
    deterministic: bool,
    policy_file: Option<&std::path::Path>,
) -> Result<EngagementPolicy> {
    match policy_file {
        Some(path) => Ok(piscem_rs::io::policy::ThreadPolicy::load(Some(path))?.parallel_decode),
        None => Ok(default_engagement_policy(sketch, deterministic)),
    }
}

/// The mapping side of the broker: a resizable paraseq pool plus a busy meter.
///
/// Deliberately not piscem's `MappingConsumer`. That one is tied to piscem's
/// `MappingStats`, and all the broker actually needs is "how do I resize you"
/// and "how much thread-time have you spent working" — both of which salmon can
/// answer with its own types.
pub struct MappingConsumer {
    pool: paraseq::parallel::ThreadPool,
    busy: Arc<BusyMeter>,
}

impl MappingConsumer {
    pub fn new(pool: paraseq::parallel::ThreadPool, busy: Arc<BusyMeter>) -> Self {
        Self { pool, busy }
    }
}

impl Consumer for MappingConsumer {
    fn set_threads(&self, n: usize) -> Result<(), ResizeError> {
        self.pool.set_threads(n);
        Ok(())
    }

    /// Workers really running, which lags `set_threads` by up to a batch.
    ///
    /// `total_live()` rather than this pool's own share: `Collection` splits the
    /// pool across readers and the handle held here is the parent, whose share
    /// never runs anything. Reading its own `live()` would report zero for the
    /// whole run and make the consumer look 100% idle however hard it is working.
    fn live_threads(&self) -> usize {
        self.pool.total_live()
    }

    fn work(&self) -> Work {
        self.busy.work()
    }
}


/// A reader over one opened input, whatever decoder it ended up with.
pub type Input = paraseq::fastx::Reader<Box<dyn std::io::Read + Send>>;

/// Everything the mapping pass needs once the decode/map split is decided.
///
/// Owns the decoder pool for the run: dropping it would retire the decode
/// workers out from under readers still pulling from them.
pub struct MappingPlan {
    /// Opened inputs, in the order paraseq expects (paired mates adjacent).
    pub readers: Vec<Input>,
    /// The mapping pool — resizable, and what the broker actually acts on.
    pub map_pool: paraseq::parallel::ThreadPool,
    /// Where mapping threads report the time they spend mapping.
    pub busy: Arc<BusyMeter>,
    /// Advisory by construction: a broker failure records itself and is never
    /// allowed to fail the run or truncate its output. piscem learned this the
    /// hard way — an error propagating out after mapping finished unwound past
    /// the output finalisation and left a file that parsed but was empty.
    pub broker: piscem_rs::io::broker::AdvisoryBroker,
    pub map_threads: usize,
    pub decode_slots: usize,
    pub parallel: bool,
    _decode_pool: Option<rapidgzip_core::DecoderPool>,
}

/// Decide the decoder, size the budget, open every input, and arm the broker.
///
/// `groups` is one entry per fragment source: `[r1]` single-end, `[r1, r2]`
/// paired. Groups rather than a flat path list because paired `fill` reads both
/// mates together, so a group is the unit that can or cannot be decoded in
/// parallel — and because a non-seekable input among them must not cost the
/// regular files their parallel decoder.
pub fn plan(
    groups: &[Vec<std::path::PathBuf>],
    budget: usize,
    pref: piscem_rs::io::calibrate::DecoderPreference,
    policy: EngagementPolicy,
) -> Result<MappingPlan> {
    let decision = piscem_rs::io::calibrate::choose_decoder(groups, budget, pref, policy);
    let num_files: usize = groups.iter().map(|g| g.len()).sum();
    let mut exec = piscem_rs::io::fastx::plan_thread_budget(
        budget,
        num_files,
        decision.parallel,
        pref,
    )?;

    // Sized to the *entire* budget, not the expected split: `workers` is an
    // immutable maximum and a later grant above it would simply be refused,
    // leaving the broker's accounting silently diverged from reality.
    let decode_pool = if exec.parallel_gzip() {
        Some(
            rapidgzip_core::DecoderPool::builder()
                .workers(exec.effective_budget)
                .initial_worker_limit(exec.decode_slots)
                .build()
                .map_err(|e| anyhow::anyhow!("could not create decoder pool: {e}"))?,
        )
    } else {
        None
    };

    let mut readers = Vec::with_capacity(num_files);
    let mut handles = Vec::new();
    for path in groups.iter().flatten() {
        let opened = piscem_rs::io::fastx::open_input_pooled(
            path,
            decode_pool.as_ref(),
            exec.effective_budget,
        )
        .with_context(|| format!("opening {}", path.display()))?;
        handles.extend(opened.handle.clone());
        readers.push(
            piscem_rs::io::fastx::reader_with_batch_size(opened.reader)
                .map_err(|e| anyhow::anyhow!("opening {}: {e}", path.display()))?,
        );
    }

    // Inputs that turned out not to be parallel-decodable (not gzip, or not
    // seekable) produce no handle, so the plan is reconciled against what was
    // really opened rather than what was hoped for.
    exec.reconcile_parallel_decoders(handles.len());
    if exec.parallel_gzip() {
        if let Some(pool) = &decode_pool {
            pool.set_worker_limit(exec.decode_slots)
                .map_err(|e| anyhow::anyhow!("could not apply decoder execution plan: {e}"))?;
        }
    }

    let arity = groups.first().map(|g| g.len()).unwrap_or(1);
    let busy = Arc::new(BusyMeter::new());
    let map_pool =
        paraseq::parallel::ThreadPool::with_max(exec.map_threads, exec.effective_budget);
    let consumer_floor = piscem_rs::io::fastx::collection_share_floor(
        readers.len(),
        arity,
        exec.map_threads,
    );

    let broker = match (&decode_pool, exec.adaptive()) {
        (Some(pool), true) => {
            match piscem_rs::io::broker::DecodeProducer::new(pool.clone(), handles.clone()) {
                Ok(producer) => {
                    let built = thread_broker::ThreadBroker::builder_with(
                        MappingConsumer::new(map_pool.clone(), Arc::clone(&busy)),
                        producer,
                        piscem_rs::io::fastx::broker_config_for_budget(exec.effective_budget),
                    )
                    .budget(exec.effective_budget)
                    .initial_producer_slots(exec.decode_slots)
                    .min_consumer_threads(consumer_floor)
                    .build()
                    .map_err(|e| anyhow::anyhow!("invalid thread broker configuration: {e}"))?;
                    piscem_rs::io::broker::AdvisoryBroker::start(
                        built,
                        exec.map_threads,
                        exec.decode_slots,
                    )
                }
                Err(error) => piscem_rs::io::broker::AdvisoryBroker::failed(
                    piscem_rs::io::broker::BrokerFailureStage::ProducerMeasurementStartup,
                    error,
                    exec.map_threads,
                    exec.decode_slots,
                ),
            }
        }
        _ => piscem_rs::io::broker::AdvisoryBroker::disabled(),
    };

    tracing::info!(
        requested_budget = exec.requested_budget,
        effective_budget = exec.effective_budget,
        mapping_threads = exec.map_threads,
        decode_slots = exec.decode_slots,
        parallel_decode = exec.parallel_gzip(),
        reason = ?decision.reason,
        "thread execution plan"
    );

    Ok(MappingPlan {
        readers,
        map_pool,
        busy,
        broker,
        map_threads: exec.map_threads,
        decode_slots: exec.decode_slots,
        parallel: exec.parallel_gzip(),
        _decode_pool: decode_pool,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn t(sketch: bool, det: bool) -> usize {
        default_engagement_policy(sketch, det).min_threads_per_stream
    }

    /// Only the *ordering* is asserted, because only the ordering follows from
    /// the mechanism. Three of the four magnitudes are still provisional, so
    /// pinning them here would be pinning a guess.
    #[test]
    fn thresholds_order_by_per_fragment_consumer_cost() {
        // More work per fragment => smaller producer cost share => engage later.
        assert!(
            t(false, false) > t(true, false),
            "selective alignment extends and scores alignments per candidate, so \
             it must engage later than sketch"
        );
        assert!(
            t(false, true) > t(true, true),
            "the SA-vs-sketch gap must survive --deterministic"
        );
        // --deterministic removes online inference, eq-class assembly, bias
        // collection and prefix detection, so it must engage *earlier*.
        assert!(
            t(true, true) < t(true, false),
            "--deterministic lightens sketch, so its threshold must fall"
        );
        assert!(
            t(false, true) < t(false, false),
            "--deterministic lightens selective alignment, so its threshold must fall"
        );
    }

    /// Sketch inherits piscem's measured constant; drift here means the two
    /// projects have silently disagreed about the same measurement.
    /// Sketch is the one cell that has been measured on salmon's own workload:
    /// producer_cost_share = 0.1004 +/- 0.0056, so 1/share = 10. It deliberately
    /// does *not* track piscem's default of 8 -- inheriting that constant was an
    /// assumption, and measurement did not bear it out.
    #[test]
    fn sketch_uses_the_value_measured_for_salmon_not_piscem() {
        assert_eq!(t(true, false), 10);
        assert_ne!(
            t(true, false),
            EngagementPolicy::default().min_threads_per_stream,
            "measured on salmon's own workload; piscem's 8 was slightly optimistic"
        );
    }

    /// A budget below the threshold must not engage, and one at it must.
    #[test]
    fn the_threshold_is_per_gzip_input() {
        let p = default_engagement_policy(true, false);
        let n = p.min_threads_per_stream;
        assert!(!p.engages(n * 2 - 1, 2), "just below the threshold per stream");
        assert!(p.engages(n * 2, 2), "exactly at the threshold per stream");
        assert!(!p.engages(1024, 0), "no gzip input: only ever costs threads");
    }
}

