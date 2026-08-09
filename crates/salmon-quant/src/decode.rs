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

/// Engagement threshold for sketch (pseudoalignment) mode.
///
/// Sketch does per-fragment work close to what piscem does, so it inherits
/// piscem's measured default. Recorded there across 16 cells: the serial decoder
/// won every cell in which the parallel decoder could not add concurrency, by
/// 8–28%, and lost by up to 4.2× in every cell where it could.
pub const ENGAGEMENT_SKETCH: usize = 8;

/// Engagement threshold for selective-alignment mode (salmon's default).
///
/// Selective alignment extends and scores alignments per candidate, so `C` is
/// much larger than in sketch mode while `P` — inflating the same bytes — is
/// unchanged. By the `1 + C/P` relation above the threshold must rise, and the
/// question is only by how much.
///
/// **Provisional.** This constant is a placeholder until the crossover sweep in
/// `scripts/decode_crossover.sh` has been run on both modes; it must not ship
/// unmeasured. The prior from a 4× consumer-cost ratio is `1 + 4·7 = 29`, which
/// is well above the 12–16 one would guess by eye — precisely the reason to
/// measure rather than assume. Override at runtime with `--thread-policy`.
pub const ENGAGEMENT_SELECTIVE_ALIGNMENT: usize = 29;

/// The engagement policy salmon defaults to for a given mapping mode.
///
/// Split out and public so the crossover harness can assert against the same
/// values the binary uses, rather than restating them.
pub fn default_engagement_policy(sketch: bool) -> EngagementPolicy {
    EngagementPolicy {
        min_threads_per_stream: if sketch {
            ENGAGEMENT_SKETCH
        } else {
            ENGAGEMENT_SELECTIVE_ALIGNMENT
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
    policy_file: Option<&std::path::Path>,
) -> Result<EngagementPolicy> {
    match policy_file {
        Some(path) => Ok(piscem_rs::io::policy::ThreadPolicy::load(Some(path))?.parallel_decode),
        None => Ok(default_engagement_policy(sketch)),
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

    /// The two modes must not share a threshold — that is the entire point of
    /// this module. A refactor that collapses them should fail here.
    #[test]
    fn selective_alignment_engages_later_than_sketch() {
        let sketch = default_engagement_policy(true);
        let sa = default_engagement_policy(false);
        assert!(
            sa.min_threads_per_stream > sketch.min_threads_per_stream,
            "selective alignment does more work per fragment, so its producer \
             cost share is lower and the parallel decoder must engage later"
        );
    }

    /// Sketch inherits piscem's measured constant; drift here means the two
    /// projects have silently disagreed about the same measurement.
    #[test]
    fn sketch_matches_piscems_measured_default() {
        assert_eq!(
            default_engagement_policy(true),
            EngagementPolicy::default(),
            "sketch should track piscem's default rather than restate it"
        );
    }

    /// A budget below the threshold must not engage, and one at it must.
    #[test]
    fn the_threshold_is_per_gzip_input() {
        let p = default_engagement_policy(true);
        assert!(!p.engages(8 * 2 - 1, 2), "just below 8 per stream");
        assert!(p.engages(8 * 2, 2), "exactly 8 per stream");
        assert!(!p.engages(1024, 0), "no gzip input: only ever costs threads");
    }
}

