//! Diagnostic: for read pairs the sketch leaves unmapped, show each mate's
//! accepted hits (tid, orientation, pos) and the merge_se_mappings outcome, to
//! see why merge_lists rejects a pair that aligns concordantly under SW
//! (orientation mismatch? fragment-length? shared tid not co-accepted?).
//! Usage: trace_pair <index_dir> <reads.json>

use piscem_rs::mapping::cache::MappingCache;
use piscem_rs::mapping::engine::map_read;
use piscem_rs::mapping::filters::PoisonState;
use piscem_rs::mapping::hit_searcher::{HitSearcher, SkippingStrategy};
use piscem_rs::mapping::merge_pairs::merge_se_mappings;
use piscem_rs::mapping::sketch_hit_simple::SketchHitInfoSimple;
use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;
use salmon_index::SalmonIndex;
use sshash_lib::dispatch_on_k;
use std::collections::BTreeSet;

fn main() {
    let mut a = std::env::args().skip(1);
    let idxdir = a.next().unwrap();
    let json = a.next().unwrap();
    let idx = SalmonIndex::load(idxdir).expect("idx");
    let ridx = idx.inner();
    let k = ridx.k();
    let mut hs = HitSearcher::new(ridx);
    let v: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(json).unwrap()).unwrap();
    for (ord, rec) in v.as_object().unwrap() {
        let r1 = rec["r1"].as_str().unwrap().as_bytes();
        let r2 = rec["r2"].as_str().unwrap().as_bytes();
        // The index's k-mer length is only known at run time, but the mapping
        // engine is generic over it for speed. `dispatch_on_k!` expands to a match
        // that picks the right monomorphized instantiation, binding it as `K`.
        dispatch_on_k!(k, K => {
            let mut left = MappingCache::<SketchHitInfoSimple>::new(k);
            let mut right = MappingCache::<SketchHitInfoSimple>::new(k);
            let mut out = MappingCache::<SketchHitInfoSimple>::new(k);
            let mut pois = PoisonState::new(ridx.poison_table());
            if r1.len() >= k {
                let mut q = PiscemStreamingQuery::<K>::new(ridx.dict());
                map_read::<K, SketchHitInfoSimple, _, _>(r1, &mut left, &mut hs, &mut q, ridx, &mut pois, SkippingStrategy::Permissive);
            }
            if r2.len() >= k {
                let mut q = PiscemStreamingQuery::<K>::new(ridx.dict());
                map_read::<K, SketchHitInfoSimple, _, _>(r2, &mut right, &mut hs, &mut q, ridx, &mut pois, SkippingStrategy::Permissive);
            }
            // The transcripts each mate accepted. A pair can only be concordant on
            // a transcript both mates hit, so the intersection is the candidate set
            // and an empty one already explains the failure.
            let lt: BTreeSet<u32> = left.accepted_hits.iter().map(|h| h.tid).collect();
            let rt: BTreeSet<u32> = right.accepted_hits.iter().map(|h| h.tid).collect();
            let shared: Vec<u32> = lt.intersection(&rt).copied().collect();
            if std::env::var_os("TRACE_TSV").is_some() {
                // classify the best (lowest-tid) shared placement
                let mut bucket = "no_shared";
                // classify only the best (lowest-tid) shared placement
                if let Some(&t) = shared.first() {
                    let l = left.accepted_hits.iter().find(|h| h.tid==t).unwrap();
                    let r = right.accepted_hits.iter().find(|h| h.tid==t).unwrap();
                    if l.is_fw == r.is_fw { bucket = "same_orientation"; }
                    else {
                        let (fw,rc) = if l.is_fw {(l.pos,r.pos)} else {(r.pos,l.pos)};
                        let fl = rc - fw;
                        bucket = if fl > -32 && fl < 2000 { "concordant_window" }
                                 else if fl <= -32 { "dovetail_neg" }
                                 else { "too_far_pos" };
                    }
                }
                println!("{ord}\t{}\t{}\t{bucket}", left.accepted_hits.len(), right.accepted_hits.len());
            } else {
                println!("--- ord {ord} | R1 acc={} R2 acc={} shared_tids={} ---",
                    left.accepted_hits.len(), right.accepted_hits.len(), shared.len());
                for &t in shared.iter().take(4) {
                    let l: Vec<_> = left.accepted_hits.iter().filter(|h| h.tid==t).map(|h| (h.is_fw, h.pos)).collect();
                    let r: Vec<_> = right.accepted_hits.iter().filter(|h| h.tid==t).map(|h| (h.is_fw, h.pos)).collect();
                    println!("    {} R1(fw,pos)={:?}  R2(fw,pos)={:?}", idx.ref_name(t as usize), l, r);
                }
                merge_se_mappings(&mut left, &mut right, r1.len() as i32, r2.len() as i32, &mut out);
                println!("    merge -> {:?}  ({} hits)", out.map_type, out.accepted_hits.len());
            }
        });
    }
}
