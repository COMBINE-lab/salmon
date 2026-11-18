//! Pseudoalignment from raw hits computed elsewhere: the GPU k-mer lookup
//! (`salmon_infer::kmers`), which finds piscem's raw hits for either
//! `SkippingStrategy`. This module turns them into exactly what
//! [`crate::sketch`] gets from `map_read` with that strategy (which only
//! changes how raw hits are collected): piscem's target selection is replayed
//! on them, then the same pairing and conversion run.
//!
//! The selection below is piscem-rs 0.12.3 `mapping::engine::map_read` from
//! its raw hits on (BSD-3-Clause, COMBINE-lab): `collect_live_targets`, the
//! occurrence-recovery retry, `visit_ambiguous_ecs`, `push_if_accepted` and
//! the `max_read_occ` cut. Its hash-map fallback, taken only for an occurrence
//! list not ascending by transcript (never in an index piscem builds), is not
//! ported: such a read returns `false` and the caller maps it on the CPU.

use piscem_rs::index::eq_classes::ec_entry_transcript_id;
use piscem_rs::index::reference_index::ReferenceIndex;
use piscem_rs::mapping::cache::MappingCache;
use piscem_rs::mapping::hits::{HitDirection, MappingType, SimpleHit, SketchHitInfo};
use piscem_rs::mapping::projected_hits::ProjectedHits;
use piscem_rs::mapping::sketch_hit_simple::SketchHitInfoSimple;
use salmon_core::RefProvider;

use crate::score::ScoredMapping;
use crate::sketch::{finish_pair, prepare_pair, push_single, SketchScratch};

/// `[read_pos, unitig << 1 | fw_on_unitig, kmer_pos_in_unitig]`, as
/// `salmon_infer::kmers::Hit`.
pub type RawHit = [u32; 3];

type Cache = MappingCache<SketchHitInfoSimple>;

/// The inputs of `salmon_infer::kmers::KmerTable::build` for `index`: each
/// unitig's first k-mer id (and the total), whether the unitig occurs in a
/// reference (piscem skips hits on one that does not), and its `j`-th k-mer.
#[allow(clippy::type_complexity)]
pub fn unitig_kmers(
    index: &ReferenceIndex,
) -> (
    Vec<u32>,
    impl Fn(usize) -> bool + Sync + '_,
    impl Fn(usize, u32) -> u64 + Sync + '_,
) {
    let spss = index.dict().spss();
    let (nu, k1) = (spss.num_strings(), index.k() as u64 - 1);
    // shortcut: ids past u32 make `build` refuse the table (CPU path)
    let kstart = (0..=nu)
        .map(|u| u32::try_from(spss.string_offset(u) - u * k1).unwrap_or(u32::MAX))
        .collect();
    let keep = kept(index);
    let k = index.k();
    let kmer = move |u: usize, j: u32| {
        spss.decode_mmer_at((spss.string_offset(u as u64) + u64::from(j)) as usize, k)
    };
    (kstart, keep, kmer)
}

/// Whether unitig `u` occurs in a reference (piscem drops hits on the others).
pub fn kept(index: &ReferenceIndex) -> impl Fn(usize) -> bool + Sync + '_ {
    move |u: usize| !index.contig_table().contig_entries(u as u64).is_empty()
}

/// piscem's sshash dictionary, flattened for `KmerTable::build_sshash`.
pub struct SshashParts<'i> {
    pub m: usize,
    /// the dictionary hasher's magic
    pub magic: u64,
    /// the SPSS text, 2 bits per base
    pub spss: &'i [u8],
    /// each unitig's first base, plus the total
    pub soff: Vec<u64>,
    /// bucket `b`'s super-k-mer positions are `offs[bstart[b]..bstart[b + 1]]`
    pub bstart: Vec<u32>,
    pub offs: Vec<u32>,
}

/// `None` when the offsets do not fit `u32`.
pub fn sshash_parts(index: &ReferenceIndex) -> Option<SshashParts<'_>> {
    use value_traits::slices::SliceByValue;
    let dict = index.dict();
    let (spss, sk) = (dict.spss(), dict.index_ref());
    use rayon::prelude::*;
    let n = sk.offsets.len();
    // random access into sshash's compressed arrays: parallel, as this runs
    // while the CPU maps the first batches
    let ef = &sk.num_super_kmers_before_bucket;
    let mut bstart = (0..sk.num_buckets())
        .into_par_iter()
        .map(|b| u32::try_from(ef.access(b) as usize + b).ok())
        .collect::<Option<Vec<_>>>()?;
    bstart.push(u32::try_from(n).ok()?);
    let offs = (0..n)
        .into_par_iter()
        .map(|i| u32::try_from(sk.offsets.index_value(i)).ok())
        .collect::<Option<_>>()?;
    Some(SshashParts {
        m: dict.m(),
        magic: sshash_lib::hasher::DeterministicHasher::new(dict.seed()).magic(),
        spss: spss.strings_raw(),
        soff: (0..spss.offsets_len())
            .into_par_iter()
            .map(|i| spss.offset_at(i))
            .collect(),
        bstart,
        offs,
    })
}

/// The hits as piscem's `(read_pos, ProjectedHits)` pairs. `unitig_base`,
/// when not empty, is each unitig's first base in the SPSS plus the total
/// ([`SshashParts::soff`] as `u32`), read instead of sshash's Elias-Fano
/// select.
fn project<'i>(
    index: &'i ReferenceIndex,
    hits: &[RawHit],
    unitig_base: &[u32],
) -> Vec<(i32, ProjectedHits<'i>)> {
    let spss = index.dict().spss();
    let k = index.k() as u32;
    hits.iter()
        .map(|&[read_pos, uf, cpos]| {
            let u = uf >> 1;
            let (begin, end) = match unitig_base.get(u as usize..u as usize + 2) {
                Some(&[b, e]) => (u64::from(b), u64::from(e)),
                _ => spss.string_offsets(u64::from(u)),
            };
            let hit = ProjectedHits::new(
                u,
                cpos,
                uf & 1 == 1,
                (end - begin) as u32,
                begin + u64::from(cpos),
                k,
                index.contig_table().contig_entries(u64::from(u)),
            );
            (read_pos as i32, hit)
        })
        .collect()
}

/// piscem `map_read` after its hit collection, on `hits` of a read of length
/// `rl`; `false` for a read piscem would send down its hash-map fallback.
fn select(
    cache: &mut Cache,
    index: &ReferenceIndex,
    hits: &[RawHit],
    rl: usize,
    unitig_base: &[u32],
) -> bool {
    cache.clear();
    // salmon never loads the poison table; `map_read` would scan with one
    if index.poison_table().is_some() {
        return false;
    }
    let raw = project(index, hits, unitig_base);
    cache.has_matching_kmers = !raw.is_empty();
    if cache.has_matching_kmers {
        let ambig = index.has_ec_table();
        let (rl, k) = (rl as i32, cache.k as i32);
        let mut nvh = 0;
        let mut min_occ = u64::MAX;
        let first_pass = (cache.max_hit_occ - 1) as u64;
        let mut lean = live_targets(
            cache,
            &raw,
            &mut nvh,
            &mut min_occ,
            first_pass,
            rl,
            k,
            ambig,
            index,
        );
        if lean.is_some()
            && cache.attempt_occ_recover
            && min_occ >= cache.max_hit_occ as u64
            && min_occ < cache.max_hit_occ_recover as u64
        {
            nvh = 0;
            cache.ambiguous_hit_indices.clear();
            let allowed = min_occ;
            lean = live_targets(
                cache,
                &raw,
                &mut nvh,
                &mut min_occ,
                allowed,
                rl,
                k,
                ambig,
                index,
            );
        }
        if lean.is_none() {
            return false;
        }
        if ambig && !cache.live.is_empty() && !cache.ambiguous_hit_indices.is_empty() {
            nvh += ambiguous_ecs(cache, &raw, index);
        }
        for (tid, target) in &cache.live {
            accept(*tid, target, nvh, &mut cache.accepted_hits);
        }
    }
    if cache.accepted_hits.len() > cache.max_read_occ {
        cache.accepted_hits.clear();
        cache.map_type = MappingType::Unmapped;
    } else if !cache.accepted_hits.is_empty() {
        cache.map_type = MappingType::SingleMapped;
    }
    true
}

/// piscem `collect_live_targets` (engine.rs:720): the targets every valid hit
/// so far has reached, ascending by transcript id; `None` on an unsorted
/// occurrence list.
#[allow(clippy::too_many_arguments)]
fn live_targets(
    cache: &mut Cache,
    raw: &[(i32, ProjectedHits<'_>)],
    nvh_out: &mut u32,
    min_occ: &mut u64,
    max_allowed_occ: u64,
    rl: i32,
    k: i32,
    ambig: bool,
    index: &ReferenceIndex,
) -> Option<()> {
    let live = &mut cache.live;
    live.clear();
    let encoding = index.encoding();
    let apply = |target: &mut SketchHitInfoSimple,
                 read_pos: i32,
                 hit: &ProjectedHits<'_>,
                 entry: u64,
                 nvh: u32| {
        if target.max_hits_for_target() >= nvh {
            let at = hit.decode_hit(entry, &encoding);
            if at.is_fw {
                target.add_fw(at.pos as i32, read_pos, rl, k, rl, 1.0);
            } else {
                target.add_rc(at.pos as i32, read_pos, rl, k, rl, 1.0);
            }
        }
    };
    for (hit_idx, (read_pos, hit)) in raw.iter().enumerate() {
        let occs = hit.ref_range();
        *min_occ = (*min_occ).min(occs.len() as u64);
        if occs.len() as u64 > max_allowed_occ {
            if ambig {
                cache.ambiguous_hit_indices.push(hit_idx as u32);
            }
            continue;
        }
        let nvh = *nvh_out;
        let mut prev: Option<u32> = None;
        if nvh == 0 {
            // every occurrence enters; a transcript listed twice shares its accumulator
            for entry in occs.iter() {
                let tid = encoding.transcript_id(entry);
                match prev {
                    Some(p) if tid < p => return None,
                    Some(p) if tid == p => {}
                    _ => live.push((tid, SketchHitInfoSimple::default())),
                }
                prev = Some(tid);
                apply(&mut live.last_mut().unwrap().1, *read_pos, hit, entry, nvh);
            }
        } else {
            let mut i = 0;
            for entry in occs.iter() {
                let tid = encoding.transcript_id(entry);
                if prev.is_some_and(|p| tid < p) {
                    return None;
                }
                prev = Some(tid);
                while i < live.len() && live[i].0 < tid {
                    i += 1;
                }
                if i < live.len() && live[i].0 == tid {
                    apply(&mut live[i].1, *read_pos, hit, entry, nvh);
                }
            }
        }
        live.retain(|(_, t)| t.max_hits_for_target() > nvh);
        *nvh_out += 1;
        if live.is_empty() {
            return Some(());
        }
    }
    Some(())
}

/// piscem `visit_ambiguous_ecs` + `apply_ec_entry` (engine.rs:612, :580):
/// credit the live targets in the class of each too-repetitive hit; the number
/// of classes visited, each one a valid hit.
fn ambiguous_ecs(
    cache: &mut Cache,
    raw: &[(i32, ProjectedHits<'_>)],
    index: &ReferenceIndex,
) -> u32 {
    let ec_table = index.ec_table().unwrap();
    let max_ec = cache.max_ec_card as usize;
    let Cache {
        live,
        ambiguous_hit_indices,
        observed_ecs,
        ..
    } = cache;
    let mut visit = |ent: u64, fw_on_contig: bool| {
        let tid = ec_entry_transcript_id(ent);
        if let Ok(i) = live.binary_search_by_key(&tid, |&(t, _)| t) {
            let t = &mut live[i].1;
            match (ent & 0x3, fw_on_contig) {
                (0, true) | (1, false) => t.inc_fw_hits(),
                (0, false) | (1, true) => t.inc_rc_hits(),
                _ => {
                    t.inc_fw_hits();
                    t.inc_rc_hits();
                }
            }
        }
    };
    observed_ecs.clear();
    let (mut min_size, mut min_ec, mut min_idx) = (usize::MAX, u64::MAX, 0u32);
    let mut visited = 0;
    for &hit_idx in ambiguous_hit_indices.iter() {
        let hit = &raw[hit_idx as usize].1;
        let fw_on_contig = hit.hit_fw_on_contig();
        let ec = ec_table.ec_for_tile(u64::from(hit.contig_id()));
        if !observed_ecs.insert(ec | if fw_on_contig { 0 } else { 1 << 63 }) {
            continue;
        }
        let entries = ec_table.entries_for_ec(ec);
        if entries.len() < min_size {
            (min_size, min_ec, min_idx) = (entries.len(), ec, hit_idx);
        }
        if entries.len() > max_ec {
            continue;
        }
        visited += 1;
        for ent in entries.iter() {
            visit(ent, fw_on_contig);
        }
    }
    // last resort: the smallest class when none was small enough
    if visited == 0 && min_ec != u64::MAX {
        let fw_on_contig = raw[min_idx as usize].1.hit_fw_on_contig();
        for ent in ec_table.entries_for_ec(min_ec).iter() {
            visit(ent, fw_on_contig);
        }
        visited = 1;
    }
    visited
}

/// piscem `push_if_accepted` (engine.rs:669).
fn accept(tid: u32, target: &SketchHitInfoSimple, nvh: u32, out: &mut Vec<SimpleHit>) {
    let best = target.best_hit_direction();
    let mut h = if best != HitDirection::Rc {
        target.get_fw_hit()
    } else {
        target.get_rc_hit()
    };
    if h.num_hits >= nvh {
        h.tid = tid;
        out.push(h);
        if best == HitDirection::Both {
            let mut rc = target.get_rc_hit();
            rc.tid = tid;
            out.push(rc);
        }
    }
}

/// [`crate::map_single_read_sketch_into`] from the read's raw hits (found with
/// the same skipping strategy); `false` (and `result` untouched) when the CPU
/// must map it. `unitig_base` as in [`project`] (`&[]` reads sshash's).
#[allow(clippy::too_many_arguments)]
pub fn map_single_read_hits_into(
    result: &mut Vec<ScoredMapping>,
    index: &ReferenceIndex,
    scratch: &mut SketchScratch,
    read_len: usize,
    hits: &[RawHit],
    max_hit_occ: usize,
    max_read_occ: usize,
    unitig_base: &[u32],
) -> bool {
    if read_len < index.k() {
        result.clear();
        return true;
    }
    let cache = &mut scratch.left;
    cache.max_hit_occ = max_hit_occ;
    cache.max_read_occ = max_read_occ;
    if !select(cache, index, hits, read_len, unitig_base) {
        return false;
    }
    result.clear();
    push_single(result, &cache.accepted_hits, read_len as i32);
    true
}

/// [`crate::map_read_pair_sketch_into`] from each mate's raw hits (found with
/// the same skipping strategy); `false` (and `result` untouched) when the CPU
/// must map it. `unitig_base` as in [`project`] (`&[]` reads sshash's).
#[allow(clippy::too_many_arguments)]
pub fn map_read_pair_hits_into<R: RefProvider>(
    result: &mut Vec<ScoredMapping>,
    index: &ReferenceIndex,
    scratch: &mut SketchScratch,
    (l1, h1): (usize, &[RawHit]),
    (l2, h2): (usize, &[RawHit]),
    strict_orphan: bool,
    allow_dovetail: bool,
    max_hit_occ: usize,
    max_read_occ: usize,
    refs: &R,
    allow_decoy_orphans: bool,
    unitig_base: &[u32],
) -> bool {
    let k = index.k();
    if l1 < k && l2 < k {
        result.clear();
        return true;
    }
    prepare_pair(
        scratch,
        l1,
        l2,
        strict_orphan,
        allow_dovetail,
        max_hit_occ,
        max_read_occ,
    );
    for (len, hits, cache) in [(l1, h1, &mut scratch.left), (l2, h2, &mut scratch.right)] {
        if len < k {
            cache.clear();
        } else if !select(cache, index, hits, len, unitig_base) {
            return false;
        }
    }
    result.clear();
    finish_pair(
        result,
        scratch,
        l1 as i32,
        l2 as i32,
        refs,
        allow_decoy_orphans,
    );
    true
}
