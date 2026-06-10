//! paraseq parallel processor that maps reads and accumulates equivalence
//! classes, fragment lengths, observed library formats, and (for `--seqBias`)
//! observed sequence-bias contexts.
//!
//! Shared `&'a` references hold the index and thread-safe accumulators; the
//! per-thread scratch (a [`HitSearcher`] and, when bias-correcting, observed
//! sequence-bias models) is rebuilt on `Clone` (paraseq clones the processor
//! per worker) and merged into shared state in `on_thread_complete`.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

use paraseq::fastx::RefRecord;
use paraseq::parallel::{PairedParallelProcessor, ParallelProcessor};
use paraseq::Record;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};

use salmon_core::{is_compatible, LibraryFormat, MateStatus};
use salmon_eqclass::{range_factorize_bins, EquivalenceClassBuilder, TranscriptGroup};
use salmon_index::SalmonIndex;
use salmon_map::{
    map_read_pair, map_read_pair_sketch, map_single_read, map_single_read_sketch, MapConfig,
    ScoredMapping,
};
use salmon_model::seqbias::{SBModel, CONTEXT_LEFT, CONTEXT_LENGTH, CONTEXT_RIGHT};
use salmon_model::{FragmentLengthDistribution, LibraryTypeDetector};

/// Shared, thread-safe state (all `Copy` references).
#[derive(Clone, Copy)]
pub(crate) struct Shared<'a> {
    pub salmon: &'a SalmonIndex,
    pub eq: &'a EquivalenceClassBuilder,
    pub fld: &'a FragmentLengthDistribution,
    /// library-type detector when auto-detecting (`-l A`); else `None`
    pub detector: Option<&'a LibraryTypeDetector>,
    pub map_cfg: &'a MapConfig,
    pub sketch: bool,
    pub skip: SkippingStrategy,
    /// range-factorization bins (0 = disabled)
    pub range_factorization_bins: u32,
    /// expected library format for strand-compatibility filtering; `None` when
    /// auto-detecting (filtering only applies with an explicit `-l`)
    pub expected_format: Option<LibraryFormat>,
    /// linear-space weight multiplier for incompatible mappings
    pub incompat_prior: f64,
    /// drop incompatible mappings entirely (set when `incompat_prior <= 0`)
    pub ignore_incompat: bool,
    /// collect observed sequence-bias contexts (`--seqBias`)
    pub collect_seqbias: bool,
    /// shared observed (fw, rc) sequence-bias models to merge per-thread results into
    pub seqbias_obs: Option<&'a Mutex<(SBModel, SBModel)>>,
    pub num_processed: &'a AtomicU64,
    pub num_mapped: &'a AtomicU64,
}

/// Per-thread mapping processor.
pub(crate) struct QuantProcessor<'a> {
    pub shared: Shared<'a>,
    pub hs: Option<HitSearcher<'a>>,
    /// per-thread observed (fw, rc) sequence-bias models (when collecting)
    pub seqbias: Option<(SBModel, SBModel)>,
}

impl<'a> QuantProcessor<'a> {
    pub fn new(shared: Shared<'a>) -> Self {
        let seqbias = shared
            .collect_seqbias
            .then(|| (SBModel::new(), SBModel::new()));
        Self {
            shared,
            hs: None,
            seqbias,
        }
    }
}

impl Clone for QuantProcessor<'_> {
    fn clone(&self) -> Self {
        // Per-thread scratch is rebuilt fresh in each worker.
        QuantProcessor::new(self.shared)
    }
}

/// Collect the orientation-aware 5'/3' sequence-bias contexts for one
/// representative mapping. A forward read's 5' context feeds the forward model;
/// a reverse read's 5' context (reverse-complemented) feeds the RC model. For a
/// paired fragment the opposite end feeds the other model.
fn collect_context(salmon: &SalmonIndex, m: &ScoredMapping, obs: &mut (SBModel, SBModel)) {
    let seq = salmon.ref_seq(m.tid);
    let n = seq.len() as i32;
    let (cl, cr, k) = (CONTEXT_LEFT as i32, CONTEXT_RIGHT as i32, CONTEXT_LENGTH as i32);

    // Forward-strand 5' context window starts CONTEXT_LEFT before the position;
    // a reverse-strand context window is the reverse-complement around it.
    let mut add_fwd = |obs: &mut SBModel, five_prime: i32| {
        let s = five_prime - cl;
        if s >= 0 && s + k <= n {
            obs.add_context(&seq[s as usize..(s + k) as usize], false, 1.0);
        }
    };
    let mut add_rev = |obs: &mut SBModel, five_prime: i32| {
        // window centered so the 5' (rightmost) base sits CONTEXT_LEFT from the
        // RC window's start; reverse-complement on insertion.
        let s = five_prime - cr;
        if s >= 0 && s + k <= n {
            obs.add_context(&seq[s as usize..(s + k) as usize], true, 1.0);
        }
    };

    if m.is_fw {
        add_fwd(&mut obs.0, m.ref_pos); // 5' -> forward model
        if m.fragment_len > 0 {
            add_rev(&mut obs.1, m.ref_pos + m.fragment_len - 1); // 3' -> RC model
        }
    } else {
        add_rev(&mut obs.1, m.ref_pos); // reverse read's 5' -> RC model
        if m.fragment_len > 0 {
            add_fwd(&mut obs.0, m.ref_pos - m.fragment_len + 1); // 3' -> forward model
        }
    }
}

/// Record one fragment's weighted mappings into the shared accumulators.
fn record(sh: &Shared, maps: &[ScoredMapping], seqbias: Option<&mut (SBModel, SBModel)>) {
    sh.num_processed.fetch_add(1, Ordering::Relaxed);
    if maps.is_empty() {
        return;
    }
    sh.num_mapped.fetch_add(1, Ordering::Relaxed);

    // Sample the observed library format for auto-detection (auto mode only).
    if let Some(det) = sh.detector {
        if det.is_active() {
            if let Some(m) = maps
                .iter()
                .filter(|m| m.format.is_some())
                .max_by(|a, b| a.weight.total_cmp(&b.weight))
            {
                det.add_sample(m.format.unwrap());
            }
        }
    }

    // Strand-compatibility filtering against an explicit expected format.
    let mut pairs: Vec<(u32, f64)> = maps
        .iter()
        .filter_map(|m| match sh.expected_format {
            Some(exp) => {
                if is_compatible(exp, m.format, m.is_fw, m.status) {
                    Some((m.tid, m.weight))
                } else if sh.ignore_incompat {
                    None
                } else {
                    Some((m.tid, m.weight * sh.incompat_prior))
                }
            }
            None => Some((m.tid, m.weight)),
        })
        .collect();
    if pairs.is_empty() {
        return; // no compatible mapping -> fragment is unassigned
    }

    // Observe a fragment length from the best concordant compatible pair.
    let best_concordant = maps
        .iter()
        .filter(|m| {
            m.status == MateStatus::PairedEndPaired
                && m.fragment_len > 0
                && sh
                    .expected_format
                    .is_none_or(|exp| is_compatible(exp, m.format, m.is_fw, m.status))
        })
        .max_by(|a, b| a.weight.total_cmp(&b.weight));
    if let Some(best) = best_concordant {
        sh.fld.add_val(best.fragment_len as usize, 0.0);
    }

    // Collect the observed sequence-bias context from a representative mapping.
    if let Some(obs) = seqbias {
        let rep = best_concordant.or_else(|| {
            maps.iter().max_by(|a, b| a.weight.total_cmp(&b.weight))
        });
        if let Some(m) = rep {
            collect_context(sh.salmon, m, obs);
        }
    }

    // Build the equivalence class: sorted, de-duplicated transcript ids + weights.
    pairs.sort_by_key(|p| p.0);
    pairs.dedup_by_key(|p| p.0);
    let tids: Vec<u32> = pairs.iter().map(|p| p.0).collect();
    let mut weights: Vec<f64> = pairs.iter().map(|p| p.1).collect();

    // Normalize to per-fragment conditional probabilities (sum to 1).
    let wsum: f64 = weights.iter().sum();
    if wsum > 0.0 {
        for w in &mut weights {
            *w /= wsum;
        }
    }

    let group = if sh.range_factorization_bins > 0 {
        let bins = range_factorize_bins(&weights, sh.range_factorization_bins);
        TranscriptGroup::with_bins(tids, bins)
    } else {
        TranscriptGroup::from_sorted(tids)
    };
    sh.eq.add_group(group, weights, 1);
}

/// Merge this thread's observed sequence-bias models into the shared models.
fn merge_seqbias(proc: &QuantProcessor) {
    if let (Some(local), Some(shared_mtx)) = (proc.seqbias.as_ref(), proc.shared.seqbias_obs) {
        let mut g = shared_mtx.lock().unwrap();
        g.0.combine_counts(&local.0);
        g.1.combine_counts(&local.1);
    }
}

impl<'a, 'r> PairedParallelProcessor<RefRecord<'r>> for QuantProcessor<'a> {
    fn process_record_pair_batch(
        &mut self,
        pairs: impl Iterator<Item = (RefRecord<'r>, RefRecord<'r>)>,
    ) -> paraseq::Result<()> {
        let QuantProcessor { shared, hs, seqbias } = self;
        let sh = *shared;
        let idx = sh.salmon.inner();
        if hs.is_none() {
            *hs = Some(HitSearcher::new(idx));
        }
        let hs = hs.as_mut().unwrap();
        for (r1, r2) in pairs {
            let s1 = r1.seq();
            let s2 = r2.seq();
            let maps = if sh.sketch {
                map_read_pair_sketch(idx, hs, s1.as_ref(), s2.as_ref(), sh.skip)
            } else {
                map_read_pair(idx, hs, sh.salmon, s1.as_ref(), s2.as_ref(), sh.map_cfg)
            };
            record(&sh, &maps, seqbias.as_mut());
        }
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        merge_seqbias(self);
        Ok(())
    }
}

impl<'a, 'r> ParallelProcessor<RefRecord<'r>> for QuantProcessor<'a> {
    fn process_record_batch(
        &mut self,
        records: impl Iterator<Item = RefRecord<'r>>,
    ) -> paraseq::Result<()> {
        let QuantProcessor { shared, hs, seqbias } = self;
        let sh = *shared;
        let idx = sh.salmon.inner();
        if hs.is_none() {
            *hs = Some(HitSearcher::new(idx));
        }
        let hs = hs.as_mut().unwrap();
        for rec in records {
            let s = rec.seq();
            let maps = if sh.sketch {
                map_single_read_sketch(idx, hs, s.as_ref(), sh.skip)
            } else {
                map_single_read(idx, hs, sh.salmon, s.as_ref(), sh.map_cfg)
            };
            record(&sh, &maps, seqbias.as_mut());
        }
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        merge_seqbias(self);
        Ok(())
    }
}
