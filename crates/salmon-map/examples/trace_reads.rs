//! Diagnostic: trace why specific read pairs map in kallisto but not in our
//! sketch path. For each pair we run piscem's `map_read` on each mate
//! independently and report the per-mate accepted-target sets, the
//! `has_matching_kmers` flag, and the target-set intersection — the exact
//! signal the sketch pairing decision is built from.
//!
//! Usage: trace_reads <index_dir> <reads.json>
//! reads.json: { "<ord>": {"r1": "ACGT...", "r2": "ACGT..."}, ... }

use piscem_rs::mapping::cache::MappingCache;
use piscem_rs::mapping::engine::map_read;
use piscem_rs::mapping::filters::PoisonState;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use piscem_rs::mapping::hits::SketchHitInfo;
use piscem_rs::mapping::sketch_hit_simple::SketchHitInfoSimple;
use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;
use salmon_index::SalmonIndex;
use sshash_lib::dispatch_on_k;
use std::collections::BTreeSet;

/// Per-mate trace: accepted target ids, `has_matching_kmers`, the number of raw
/// k-mer projected-hit positions, the number of candidate targets in `hit_map`,
/// and the best single-target support (max fw/rc hits over candidates). When
/// `accepted` is empty but `cand > 0`, no target was supported by *all* k-mers
/// (consensus/intersection conflict); `cand == 0` means k-mers projected to
/// nothing (e.g. only ambiguous/over-occ k-mers or poison).
struct MateTrace {
    tids: BTreeSet<u32>,
    has_kmers: bool,
    raw_hits: usize,
    cand: usize,
    best_support: u32,
}

fn trace_mate(
    idx: &salmon_index::SalmonIndex,
    hs: &mut HitSearcher,
    read: &[u8],
    max_hit_occ: usize,
    max_read_occ: usize,
) -> MateTrace {
    let ridx = idx.inner();
    let k = ridx.k();
    if read.len() < k {
        return MateTrace {
            tids: BTreeSet::new(),
            has_kmers: false,
            raw_hits: 0,
            cand: 0,
            best_support: 0,
        };
    }
    dispatch_on_k!(k, K => {
        let mut q = PiscemStreamingQuery::<K>::new(ridx.dict());
        let mut cache = MappingCache::<SketchHitInfoSimple>::new(k);
        cache.max_hit_occ = max_hit_occ;
        cache.max_read_occ = max_read_occ;
        let mut poison = PoisonState::new(ridx.poison_table());
        map_read::<K, SketchHitInfoSimple, _, _>(
            read, &mut cache, hs, &mut q, ridx, &mut poison, SkippingStrategy::Permissive,
        );
        let tids: BTreeSet<u32> = cache.accepted_hits.iter().map(|h| h.tid).collect();
        let raw_hits = hs.left_hits().len();
        let cand = cache.hit_map.len();
        let best_support = cache.hit_map.values().map(|t| t.max_hits_for_target()).max().unwrap_or(0);
        MateTrace { tids, has_kmers: cache.has_matching_kmers, raw_hits, cand, best_support }
    })
}

/// Per-k-mer-position dump for one mate: for each raw projected hit, print the
/// read position, contig id, occurrence count, and the target set. Reveals
/// whether the consensus target is hit by a contiguous run of positions and
/// which positions disagree (candidates for what kallisto skips).
fn dump_kmers(idx: &salmon_index::SalmonIndex, hs: &mut HitSearcher, read: &[u8], tag: &str) {
    let ridx = idx.inner();
    let k = ridx.k();
    if read.len() < k {
        return;
    }
    dispatch_on_k!(k, K => {
        let mut q = PiscemStreamingQuery::<K>::new(ridx.dict());
        let mut cache = MappingCache::<SketchHitInfoSimple>::new(k);
        cache.max_hit_occ = 1000;
        cache.max_read_occ = 200;
        let mut poison = PoisonState::new(ridx.poison_table());
        map_read::<K, SketchHitInfoSimple, _, _>(
            read, &mut cache, hs, &mut q, ridx, &mut poison, SkippingStrategy::Permissive,
        );
        let enc = ridx.encoding();
        println!("  [{tag}] {} raw-hit positions:", hs.left_hits().len());
        for (read_pos, ph) in hs.left_hits() {
            let mut tids: Vec<u32> = ph.ref_range().iter().map(|e| enc.transcript_id(e)).collect();
            tids.sort_unstable(); tids.dedup();
            let shown: Vec<&str> = tids.iter().take(8).map(|&t| idx.ref_name(t as usize)).collect();
            println!("    read_pos={read_pos:>3} contig={:>8} occ={:>4} ntargets={} {:?}{}",
                ph.contig_id(), ph.ref_range().len(), tids.len(), shown,
                if tids.len() > 8 { ".." } else { "" });
        }
    })
}

// Reads are traced one at a time through the same stages the mapper runs, with
// the intermediate state printed after each, so a read that fails to map can be
// attributed to a specific stage rather than to the pipeline as a whole.
fn main() {
    let mut args = std::env::args().skip(1);
    let index_dir = args.next().expect("index dir");
    let reads_json = args.next().expect("reads json");
    let max_hit_occ: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(1000);
    let max_read_occ: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(200);

    let idx = SalmonIndex::load(index_dir).expect("load index");
    let mut hs = HitSearcher::new(idx.inner());

    let text = std::fs::read_to_string(&reads_json).expect("read json");
    // minimal JSON parse: entries "ord": {"r1":"..","r2":".."}
    let v: serde_json::Value = serde_json::from_str(&text).expect("parse json");
    let obj = v.as_object().unwrap();

    for (ord, rec) in obj {
        let r1 = rec["r1"].as_str().unwrap().as_bytes();
        let r2 = rec["r2"].as_str().unwrap().as_bytes();
        let m1 = trace_mate(&idx, &mut hs, r1, max_hit_occ, max_read_occ);
        let m2 = trace_mate(&idx, &mut hs, r2, max_hit_occ, max_read_occ);
        let (t1, t2) = (&m1.tids, &m2.tids);
        let (km1, km2) = (m1.has_kmers, m2.has_kmers);
        let inter: Vec<u32> = t1.intersection(t2).copied().collect();
        let name = |s: &BTreeSet<u32>| -> Vec<String> {
            s.iter()
                .take(6)
                .map(|&t| idx.ref_name(t as usize).to_string())
                .collect()
        };
        // category: which divergence class this read falls into
        let cat = if !inter.is_empty() {
            "shared_tid"
        } else if (!t1.is_empty() && t2.is_empty() && km2)
            || (!t2.is_empty() && t1.is_empty() && km1)
        {
            // one mate has accepted targets, the other had k-mers but no
            // accepted target -> kallisto orphans, we drop (orphan-rule gap)
            "orphan_rule"
        } else if (!t1.is_empty() && t2.is_empty() && !km2)
            || (!t2.is_empty() && t1.is_empty() && !km1)
        {
            "true_orphan"
        } else if t1.is_empty() && t2.is_empty() {
            // sub-split the seeding set: did either mate have candidate targets
            // in hit_map (consensus/intersection conflict) or none at all?
            if m1.cand > 0 || m2.cand > 0 {
                "seed_consensus_conflict"
            } else if km1 || km2 {
                "seed_kmers_no_projection"
            } else {
                "seed_no_kmers"
            }
        } else {
            "both_hit_disjoint"
        };
        if std::env::var_os("TRACE_TSV").is_some() {
            println!(
                "{ord}\t{cat}\tR1[acc={} km={km1} raw={} cand={} best={}]\tR2[acc={} km={km2} raw={} cand={} best={}]",
                t1.len(), m1.raw_hits, m1.cand, m1.best_support,
                t2.len(), m2.raw_hits, m2.cand, m2.best_support,
            );
        } else {
            println!(
                "--- ord {ord} (ERR188044.{}) [{cat}] ---",
                ord.parse::<u64>().unwrap() + 1
            );
            println!("  R1: {} accepted, has_kmers={km1}, raw_hits={}, cand_targets={}, best_support={}, targets={:?}",
                t1.len(), m1.raw_hits, m1.cand, m1.best_support, name(t1));
            println!("  R2: {} accepted, has_kmers={km2}, raw_hits={}, cand_targets={}, best_support={}, targets={:?}",
                t2.len(), m2.raw_hits, m2.cand, m2.best_support, name(t2));
            println!(
                "  intersection: {} -> {:?}",
                inter.len(),
                inter
                    .iter()
                    .take(6)
                    .map(|&t| idx.ref_name(t as usize))
                    .collect::<Vec<_>>()
            );
            if std::env::var_os("VERBOSE_KMERS").is_some() {
                dump_kmers(&idx, &mut hs, r1, "R1");
                dump_kmers(&idx, &mut hs, r2, "R2");
            }
        }
    }
}
