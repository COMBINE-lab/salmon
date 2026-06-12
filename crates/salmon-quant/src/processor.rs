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
use salmon_infer::OnlineInference;
use salmon_model::seqbias::{SBModel, CONTEXT_LEFT, CONTEXT_LENGTH, CONTEXT_RIGHT};
use salmon_model::{
    gc_desc, FragmentLengthDistribution, GcFragModel, LibraryTypeDetector, SimplePosBias,
    NUM_LENGTH_CLASSES,
};

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
    /// sketch mode: only orphan a pair when the other mate had no matching
    /// k-mers (strict), instead of the default empty-accepted-target rule
    /// (`--sketchStrictOrphans`)
    pub sketch_strict_orphan: bool,
    /// fragments processed without the auxiliary (FLD) model before it is folded
    /// into the online posterior (salmon's `--numPreAuxModelSamples`)
    pub pre_burnin: u64,
    /// discard a fragment mapping to more than this many places (`--maxReadOcc`)
    pub max_read_occ: usize,
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
    /// collect observed fragment-GC contexts (`--gcBias`)
    pub collect_gcbias: bool,
    /// GC bias model bin counts (`--conditionalGCBins` × `--numGCBins`)
    pub cond_gc_bins: usize,
    pub gc_bins: usize,
    /// per-transcript cumulative-GC source (dense or rank), when GC-correcting
    pub gc_store: Option<salmon_model::GcStore<'a>>,
    /// shared observed fragment-GC model to merge per-thread results into
    pub gcbias_obs: Option<&'a Mutex<GcFragModel>>,
    /// collect observed positional bias (`--posBias`)
    pub collect_posbias: bool,
    /// per-transcript length-class index (for positional bias), when pos-correcting
    pub length_class: Option<&'a [usize]>,
    /// shared observed (5', 3') positional-bias models (per length class) to merge into
    pub posbias_obs: Option<&'a Mutex<(Vec<SimplePosBias>, Vec<SimplePosBias>)>>,
    /// online (dual-phase) inference state; when present, observed bias models
    /// are collected with abundance-aware posteriors instead of score-only weights
    pub online: Option<&'a OnlineInference>,
    /// whether the library is paired (for the online orphan penalty)
    pub paired_lib: bool,
    pub num_processed: &'a AtomicU64,
    pub num_mapped: &'a AtomicU64,
    /// mapped fragments whose representative mapping is an orphan (only one mate
    /// of a paired-end fragment mapped); feeds `lib_format_counts.json`'s
    /// `num_frags_with_inconsistent_or_orphan_mappings`.
    pub num_orphan: &'a AtomicU64,
    /// selective-alignment meta counters (salmon's `aux_info/meta_info.json`):
    /// decoy-dominated fragments, dovetail fragments, fragments that had
    /// candidates but no surviving mapping, and (for mapped fragments) candidate
    /// alignments dropped below the score threshold.
    pub num_decoy: &'a AtomicU64,
    pub num_dovetail: &'a AtomicU64,
    pub num_frags_filtered_vm: &'a AtomicU64,
    pub num_below_threshold_vm: &'a AtomicU64,
    /// when set, collect the names of unmapped fragments (`--writeUnmappedNames`)
    pub unmapped_names: Option<&'a Mutex<Vec<String>>>,
    /// when set, write per-mapping SAM records (`--writeMappings`)
    pub sam: Option<&'a crate::sam::SamWriter>,
}

/// Per-thread mapping processor.
pub(crate) struct QuantProcessor<'a> {
    pub shared: Shared<'a>,
    pub hs: Option<HitSearcher<'a>>,
    /// per-thread observed (fw, rc) sequence-bias models (when collecting)
    pub seqbias: Option<(SBModel, SBModel)>,
    /// per-thread observed fragment-GC model (when collecting)
    pub gcbias: Option<GcFragModel>,
    /// per-thread observed (5', 3') positional-bias models, per length class
    pub posbias: Option<(Vec<SimplePosBias>, Vec<SimplePosBias>)>,
    /// per-thread collected unmapped-fragment names (merged in on_thread_complete)
    pub unmapped: Vec<String>,
    /// per-thread SAM record buffer (flushed to the shared writer per batch)
    pub sam_buf: String,
}

impl<'a> QuantProcessor<'a> {
    pub fn new(shared: Shared<'a>) -> Self {
        let seqbias = shared
            .collect_seqbias
            .then(|| (SBModel::new(), SBModel::new()));
        let gcbias = shared
            .collect_gcbias
            .then(|| GcFragModel::new(shared.cond_gc_bins, shared.gc_bins));
        let posbias = shared.collect_posbias.then(|| {
            (
                (0..NUM_LENGTH_CLASSES).map(|_| SimplePosBias::default()).collect(),
                (0..NUM_LENGTH_CLASSES).map(|_| SimplePosBias::default()).collect(),
            )
        });
        Self {
            shared,
            hs: None,
            seqbias,
            gcbias,
            posbias,
            unmapped: Vec::new(),
            sam_buf: String::new(),
        }
    }
}

impl Clone for QuantProcessor<'_> {
    fn clone(&self) -> Self {
        // Per-thread scratch is rebuilt fresh in each worker.
        QuantProcessor::new(self.shared)
    }
}

/// salmon's `LOG_EPSILON = log(0.375e-10)` — the small log-probability assigned
/// to implausible placements (out-of-range fragment lengths, unexpected orphans).
const LOG_EPSILON: f64 = -23.998_158_637_57; // (0.375e-10f64).ln()
/// salmon's `numPreBurninFrags`: fold the fragment-length term into the online
/// posterior only after this many fragments have been assigned.
pub(crate) const NUM_PRE_BURNIN: u64 = 5000;

/// Collect the orientation-aware 5'/3' sequence-bias contexts for one mapping,
/// weighted by `weight` (the fragment-transcript posterior). A forward read's 5'
/// context feeds the forward model; a reverse read's 5' context
/// (reverse-complemented) feeds the RC model. For a paired fragment the opposite
/// end feeds the other model.
fn collect_context(salmon: &SalmonIndex, m: &ScoredMapping, weight: f64, obs: &mut (SBModel, SBModel)) {
    let seq = salmon.ref_seq(m.tid);
    let n = seq.len() as i32;
    let (cl, cr, k) = (CONTEXT_LEFT as i32, CONTEXT_RIGHT as i32, CONTEXT_LENGTH as i32);

    let mut add_fwd = |obs: &mut SBModel, five_prime: i32| {
        let s = five_prime - cl;
        if s >= 0 && s + k <= n {
            obs.add_context(&seq[s as usize..(s + k) as usize], false, weight);
        }
    };
    let mut add_rev = |obs: &mut SBModel, five_prime: i32| {
        let s = five_prime - cr;
        if s >= 0 && s + k <= n {
            obs.add_context(&seq[s as usize..(s + k) as usize], true, weight);
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

/// Collect observed fragment-GC contexts for one fragment's compatible paired
/// mappings, each weighted by `bias_w[i]` (the abundance-aware posterior under
/// online inference, else the normalized aux weight). Mirrors salmon's
/// per-alignment `observedGCMass.inc(gcDesc(start, stop), logProb)`.
fn collect_gc(sh: &Shared, compat: &[(&ScoredMapping, f64)], bias_w: &[f64], gc: &mut GcFragModel) {
    let Some(store) = sh.gc_store else { return };
    for (i, (m, _)) in compat.iter().enumerate() {
        if m.status != MateStatus::PairedEndPaired || m.fragment_len <= 0 {
            continue;
        }
        let start = m.ref_pos;
        let stop = m.ref_pos + m.fragment_len - 1;
        let view = store.view(m.tid as usize);
        let ref_len = view.ref_len() as i32;
        if start < 0 || stop >= ref_len {
            continue;
        }
        if let Some((ff, cf)) = gc_desc(&view, start, stop) {
            gc.inc(ff, cf, bias_w[i]);
        }
    }
}

/// Collect observed positional-bias mass for one fragment's compatible mappings,
/// each weighted by `bias_w[i]` (log space). Mirrors salmon's
/// `observedPosBias{Fwd,RC}[lengthClass].addMass(pos, refLen, logProb)`.
fn collect_pos(
    sh: &Shared,
    compat: &[(&ScoredMapping, f64)],
    bias_w: &[f64],
    pos: &mut (Vec<SimplePosBias>, Vec<SimplePosBias>),
) {
    let Some(length_class) = sh.length_class else { return };
    for (i, (m, _)) in compat.iter().enumerate() {
        if bias_w[i] <= 0.0 {
            continue;
        }
        let lc = length_class[m.tid as usize];
        let ref_len = sh.salmon.ref_len(m.tid as usize) as i32;
        let mass = bias_w[i].ln();
        if m.fw_pos >= 0 {
            pos.0[lc].add_mass(m.fw_pos, ref_len, mass);
        }
        if m.rc_pos >= 0 {
            pos.1[lc].add_mass(m.rc_pos, ref_len, mass);
        }
    }
}

/// Record one fragment's weighted mappings into the shared accumulators.
/// `log_fm` is the current minibatch's forgetting mass (used only when online
/// inference is active).
fn record(
    sh: &Shared,
    maps: &[ScoredMapping],
    log_fm: f64,
    seqbias: Option<&mut (SBModel, SBModel)>,
    gcbias: Option<&mut GcFragModel>,
    posbias: Option<&mut (Vec<SimplePosBias>, Vec<SimplePosBias>)>,
) {
    sh.num_processed.fetch_add(1, Ordering::Relaxed);
    if maps.is_empty() {
        return;
    }
    sh.num_mapped.fetch_add(1, Ordering::Relaxed);

    // Classify the fragment as orphan when its representative mapping has only
    // one mate of a pair placed (PairedEndLeft / PairedEndRight). Single-end
    // libraries never count as orphan here.
    if matches!(
        maps[0].status,
        MateStatus::PairedEndLeft | MateStatus::PairedEndRight
    ) {
        sh.num_orphan.fetch_add(1, Ordering::Relaxed);
    }

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
    let compat: Vec<(&ScoredMapping, f64)> = maps
        .iter()
        .filter_map(|m| match sh.expected_format {
            Some(exp) => {
                if is_compatible(exp, m.format, m.is_fw, m.status) {
                    Some((m, m.weight))
                } else if sh.ignore_incompat {
                    None
                } else {
                    Some((m, m.weight * sh.incompat_prior))
                }
            }
            None => Some((m, m.weight)),
        })
        .collect();
    if compat.is_empty() {
        return; // no compatible mapping -> fragment is unassigned
    }

    // Per-fragment bias-collection weights. With online (dual-phase) inference
    // these are abundance-aware posteriors `softmax(mass_t + log w_t)`, which
    // also advances the online masses by `logForgettingMass + log(posterior)`;
    // without it they fall back to the normalized aux (score) weights.
    let collecting = seqbias.is_some() || gcbias.is_some() || posbias.is_some();
    let bias_w: Vec<f64> = if collecting {
        if let Some(online) = sh.online {
            // Per-mapping log auxiliary probability (salmon's
            // `auxProb + startPosProb`, abundance-independent):
            //   logFragCov  = ln(score weight)
            //   startPosProb = proper pair -> -ln(refLen - flen + 1) (flen<=refLen,
            //                  else LOG_EPSILON); otherwise -ln(refLen)
            //   logFragProb  = proper pair (after pre-burn-in) -> live FLD pmf(flen);
            //                  unexpected orphan in a paired library -> LOG_EPSILON.
            // The FLD term discriminates isoforms/paralogs whose implied insert
            // size differs (e.g. alternative splicing), which the length norm alone
            // cannot capture.
            let use_aux = online.num_assigned() >= sh.pre_burnin;
            let mm: Vec<(u32, f64)> = compat
                .iter()
                .map(|(m, w)| {
                    let rl = sh.salmon.ref_len(m.tid as usize).max(1) as f64;
                    let proper = m.status == MateStatus::PairedEndPaired && m.fragment_len > 0;
                    let flen = m.fragment_len as f64;
                    let start_pos_prob = if proper {
                        if flen <= rl {
                            -((rl - flen + 1.0).ln())
                        } else {
                            LOG_EPSILON
                        }
                    } else {
                        -(rl.ln())
                    };
                    let log_frag_prob = if proper {
                        if use_aux {
                            sh.fld.pmf(m.fragment_len as usize)
                        } else {
                            0.0
                        }
                    } else if sh.paired_lib {
                        LOG_EPSILON // unexpected orphan
                    } else {
                        0.0
                    };
                    let log_cov = if *w > 0.0 { w.ln() } else { f64::NEG_INFINITY };
                    (m.tid, log_cov + start_pos_prob + log_frag_prob)
                })
                .collect();
            online.assign_fragment(&mm, log_fm)
        } else {
            let wsum: f64 = compat.iter().map(|(_, w)| *w).sum();
            compat
                .iter()
                .map(|(_, w)| if wsum > 0.0 { *w / wsum } else { 0.0 })
                .collect()
        }
    } else {
        Vec::new()
    };
    // After burn-in salmon freezes model collection (still advances masses).
    let collect_now = collecting && sh.online.is_none_or(|o| o.collecting());

    if collect_now {
        if let Some(obs) = seqbias {
            for (i, (m, _)) in compat.iter().enumerate() {
                collect_context(sh.salmon, m, bias_w[i], obs);
            }
        }
        if let Some(gc) = gcbias {
            collect_gc(sh, &compat, &bias_w, gc);
        }
        if let Some(pos) = posbias {
            collect_pos(sh, &compat, &bias_w, pos);
        }
    }

    let mut pairs: Vec<(u32, f64)> = compat.iter().map(|(m, w)| (m.tid, *w)).collect();

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

/// Fold the most recent fragment's selective-alignment [`MapStats`] into the
/// shared meta counters. Call once per mapped/attempted fragment on the SA path.
fn accumulate_vm_stats(sh: &Shared, maps_empty: bool) {
    let s = salmon_map::take_last_map_stats();
    if maps_empty {
        if s.decoy_dominated {
            sh.num_decoy.fetch_add(1, Ordering::Relaxed);
        } else if s.had_candidates {
            sh.num_frags_filtered_vm.fetch_add(1, Ordering::Relaxed);
        }
    } else if s.alns_below_threshold > 0 {
        sh.num_below_threshold_vm
            .fetch_add(s.alns_below_threshold as u64, Ordering::Relaxed);
    }
    if s.dovetail {
        sh.num_dovetail.fetch_add(1, Ordering::Relaxed);
    }
}

/// Merge this thread's observed bias models (sequence and GC) into the shared
/// accumulators.
/// The read name written to unmapped_names.txt: the id up to the first
/// whitespace (the conventional fragment name), lossily decoded.
fn read_name(id: &[u8]) -> String {
    let end = id.iter().position(|b| b.is_ascii_whitespace()).unwrap_or(id.len());
    String::from_utf8_lossy(&id[..end]).into_owned()
}

/// Merge this thread's collected unmapped-fragment names into the shared list.
fn merge_unmapped(proc: &mut QuantProcessor) {
    if let Some(shared) = proc.shared.unmapped_names {
        if !proc.unmapped.is_empty() {
            shared.lock().unwrap().append(&mut proc.unmapped);
        }
    }
}

fn merge_bias(proc: &QuantProcessor) {
    if let (Some(local), Some(shared_mtx)) = (proc.seqbias.as_ref(), proc.shared.seqbias_obs) {
        let mut g = shared_mtx.lock().unwrap();
        g.0.combine_counts(&local.0);
        g.1.combine_counts(&local.1);
    }
    if let (Some(local), Some(shared_mtx)) = (proc.gcbias.as_ref(), proc.shared.gcbias_obs) {
        let mut g = shared_mtx.lock().unwrap();
        g.combine_counts(local);
    }
    if let (Some(local), Some(shared_mtx)) = (proc.posbias.as_ref(), proc.shared.posbias_obs) {
        let mut g = shared_mtx.lock().unwrap();
        for (a, b) in g.0.iter_mut().zip(&local.0) {
            a.combine(b);
        }
        for (a, b) in g.1.iter_mut().zip(&local.1) {
            a.combine(b);
        }
    }
}

impl<'a, 'r> PairedParallelProcessor<RefRecord<'r>> for QuantProcessor<'a> {
    fn process_record_pair_batch(
        &mut self,
        pairs: impl Iterator<Item = (RefRecord<'r>, RefRecord<'r>)>,
    ) -> paraseq::Result<()> {
        let QuantProcessor { shared, hs, seqbias, gcbias, posbias, unmapped, sam_buf } = self;
        let sh = *shared;
        let idx = sh.salmon.inner();
        if hs.is_none() {
            *hs = Some(HitSearcher::new(idx));
        }
        let hs = hs.as_mut().unwrap();
        // one forgetting-mass timestep per minibatch (online inference)
        let log_fm = sh.online.map_or(0.0, |o| o.next_log_fm());
        for (r1, r2) in pairs {
            let s1 = r1.seq();
            let s2 = r2.seq();
            let mut maps = if sh.sketch {
                map_read_pair_sketch(idx, hs, s1.as_ref(), s2.as_ref(), sh.sketch_strict_orphan, sh.skip, sh.map_cfg.collect.max_hit_occ, sh.max_read_occ)
            } else {
                map_read_pair(idx, hs, sh.salmon, s1.as_ref(), s2.as_ref(), sh.map_cfg)
            };
            // A fragment mapping to too many places is discarded (salmon's
            // tooManyHits / maxReadOccs): treat as unmapped everywhere below.
            if maps.len() > sh.max_read_occ {
                maps.clear();
            }
            if maps.is_empty() && sh.unmapped_names.is_some() {
                unmapped.push(format!("{} u", read_name(r1.id())));
            }
            if !sh.sketch {
                accumulate_vm_stats(&sh, maps.is_empty());
            }
            if sh.sam.is_some() && !maps.is_empty() {
                crate::sam::write_fragment(
                    sam_buf, sh.salmon, r1.id(), s1.as_ref(),
                    Some((r2.id(), s2.as_ref())), &maps,
                );
            }
            record(&sh, &maps, log_fm, seqbias.as_mut(), gcbias.as_mut(), posbias.as_mut());
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        flush_sam(self);
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        merge_bias(self);
        merge_unmapped(self);
        flush_sam(self);
        Ok(())
    }
}

impl<'a, 'r> ParallelProcessor<RefRecord<'r>> for QuantProcessor<'a> {
    fn process_record_batch(
        &mut self,
        records: impl Iterator<Item = RefRecord<'r>>,
    ) -> paraseq::Result<()> {
        let QuantProcessor { shared, hs, seqbias, gcbias, posbias, unmapped, sam_buf } = self;
        let sh = *shared;
        let idx = sh.salmon.inner();
        if hs.is_none() {
            *hs = Some(HitSearcher::new(idx));
        }
        let hs = hs.as_mut().unwrap();
        let log_fm = sh.online.map_or(0.0, |o| o.next_log_fm());
        for rec in records {
            let s = rec.seq();
            let mut maps = if sh.sketch {
                map_single_read_sketch(idx, hs, s.as_ref(), sh.skip, sh.map_cfg.collect.max_hit_occ, sh.max_read_occ)
            } else {
                map_single_read(idx, hs, sh.salmon, s.as_ref(), sh.map_cfg)
            };
            if maps.len() > sh.max_read_occ {
                maps.clear();
            }
            if maps.is_empty() && sh.unmapped_names.is_some() {
                unmapped.push(format!("{} u", read_name(rec.id())));
            }
            if !sh.sketch {
                accumulate_vm_stats(&sh, maps.is_empty());
            }
            if sh.sam.is_some() && !maps.is_empty() {
                crate::sam::write_fragment(
                    sam_buf, sh.salmon, rec.id(), s.as_ref(), None, &maps,
                );
            }
            record(&sh, &maps, log_fm, seqbias.as_mut(), gcbias.as_mut(), posbias.as_mut());
        }
        Ok(())
    }

    fn on_batch_complete(&mut self) -> paraseq::Result<()> {
        flush_sam(self);
        Ok(())
    }

    fn on_thread_complete(&mut self) -> paraseq::Result<()> {
        merge_bias(self);
        merge_unmapped(self);
        flush_sam(self);
        Ok(())
    }
}

/// Flush this thread's accumulated SAM buffer to the shared writer (under lock).
fn flush_sam(proc: &mut QuantProcessor) {
    if let Some(sw) = proc.shared.sam {
        if !proc.sam_buf.is_empty() {
            let _ = sw.write_block(&proc.sam_buf);
            proc.sam_buf.clear();
        }
    }
}
