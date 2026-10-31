//! MEM extraction: extend k-mer anchors into longer exact matches.
//!
//! # Why extend at all
//!
//! A k-mer hit says 31 bases matched. But a read that genuinely came from this
//! transcript usually matches far more than 31 bases in a row — the hit is just
//! the part the index happened to be queried at. Extending each hit outward until
//! it stops matching replaces a scatter of short, fixed-length anchors with a few
//! long ones that describe the same agreement.
//!
//! That helps twice: the chainer sees fewer, more informative anchors, and a long
//! exact match is much stronger evidence than several short ones that could each
//! be coincidence.
//!
//! The bare [`collect`](crate::collect) path chains the *sparse, fixed-length*
//! k-mer anchors that piscem's skipping streaming query emits — one length-`k`
//! anchor per unitig transition (`--sparseSeeds`). This module offers two
//! alternatives that extend each seed before chaining (selected by
//! [`SeedMode`](crate::mapper::SeedMode)):
//!
//! - **Reference MEMs** ([`candidates_from_raw_hits_unimems`]): each seed is
//!   extended against the **reference transcript** until a base mismatch or a
//!   read/reference boundary. These cross unitig junctions freely — a
//!   junction-straddling match becomes a *single* anchor (`--refMEMs`).
//! - **True uni-MEMs** ([`candidates_from_raw_hits_true_unimems`]): extension is
//!   clamped to the seed's **unitig** (`[base_pos, base_pos + contig_len)`),
//!   faithfully reproducing pufferfish's `expandHitEfficient` (which stops at
//!   `CONTIG_END`). A junction-straddling match becomes *one uni-MEM per unitig*,
//!   which the chainer then stitches together. **This is the default** — it is
//!   ~11% faster than sparse seeding and at least as accurate, and the extension
//!   is computed once per orientation and projected to every occurrence (a
//!   uni-MEM is identical wherever its unitig occurs).

use std::collections::HashSet;

use piscem_rs::index::contig_table::EntryEncoding;
use piscem_rs::index::reference_index::ReferenceIndex;
use piscem_rs::mapping::hit_searcher::HitSearcher;
use piscem_rs::mapping::projected_hits::ProjectedHits;
use piscem_rs::mapping::streaming_query::PiscemStreamingQuery;
use sshash_lib::dispatch_on_k;

use salmon_core::RefProvider;

use crate::align::revcomp;
use crate::chain::chain_mems;
use crate::collect::{project_raw_hits, MappingCandidate, MemCollectorConfig};
use crate::mem::Mem;

/// Exact base equality, case-insensitive, with `N` never matching (an `N` on
/// either side terminates the exact match, as in pufferfish's extension).
#[inline]
fn base_eq(a: u8, b: u8) -> bool {
    let ua = a.to_ascii_uppercase();
    ua == b.to_ascii_uppercase() && ua != b'N'
}

/// Extend a seed anchor to a maximal exact match against the **whole
/// reference** (a *reference MEM*): the match grows until a base mismatch or a
/// read/reference boundary, crossing unitig boundaries freely.
///
/// `query`/`ref_seq` are in the query frame: the caller passes the forward read
/// (forward group) or its reverse complement (reverse group), and `m`'s
/// coordinates are already expressed in that frame.
pub fn extend_mem(query: &[u8], ref_seq: &[u8], m: Mem) -> Mem {
    extend_mem_within(query, ref_seq, m, 0, ref_seq.len() as i32)
}

/// Extend a seed anchor to a maximal exact match whose reference span is
/// constrained to `[ref_lo, ref_hi)`.
///
/// With `ref_lo = 0`, `ref_hi = ref_seq.len()` this is an unconstrained
/// reference MEM ([`extend_mem`]). Passing the seed's **unitig** reference span
/// (`[base_pos, base_pos + contig_len)`) instead yields a true
/// **uni-MEM**: an exact match that cannot extend past the unitig boundary,
/// exactly as pufferfish's `expandHitEfficient` (which stops at `CONTIG_END`).
pub fn extend_mem_within(
    query: &[u8],
    ref_seq: &[u8],
    mut m: Mem,
    ref_lo: i32,
    ref_hi: i32,
) -> Mem {
    let qlen = query.len() as i32;
    let reflen = ref_seq.len() as i32;
    let hi = ref_hi.min(reflen);

    // Extend left along the shared diagonal, not past the unitig start.
    while m.read_start > 0
        && m.ref_start > ref_lo
        && base_eq(
            query[(m.read_start - 1) as usize],
            ref_seq[(m.ref_start - 1) as usize],
        )
    {
        m.read_start -= 1;
        m.ref_start -= 1;
        m.len += 1;
    }

    // Extend right along the shared diagonal, not past the unitig end.
    while m.read_end() < qlen
        && m.ref_end() < hi
        && base_eq(query[m.read_end() as usize], ref_seq[m.ref_end() as usize])
    {
        m.len += 1;
    }

    m
}

/// Project raw k-mer hits, extend each to a uni-MEM, dedup, and chain.
///
/// Identical to [`crate::collect::candidates_from_raw_hits`] except every anchor
/// is grown to a maximal exact match before chaining and colinear seeds that
/// extend to the same uni-MEM are collapsed. Pure given the reference bytes;
/// testable with hand-built hits + a [`RefProvider`].
pub fn candidates_from_raw_hits_unimems<R: RefProvider>(
    raw_hits: &[(i32, ProjectedHits<'_>)],
    encoding: &EntryEncoding,
    refs: &R,
    read: &[u8],
    k: i32,
    cfg: &MemCollectorConfig,
) -> Vec<MappingCandidate> {
    let read_len = read.len() as i32;
    let groups = project_raw_hits(raw_hits, encoding, read_len, k, cfg.max_hit_occ);
    let rc = revcomp(read);
    let mut chain_cfg = cfg.chain.clone();
    chain_cfg.seed_len = k;

    let mut candidates = Vec::new();
    for ((tid, is_fw), seeds) in groups {
        let ref_seq = refs.ref_seq(tid);
        let query: &[u8] = if is_fw { read } else { &rc };

        // Extend each seed to its uni-MEM; same-diagonal seeds extend to an
        // identical anchor and collapse via the (read_start, ref_start, len) key.
        let mut seen = HashSet::new();
        let mut unimems = Vec::with_capacity(seeds.len());
        for s in seeds {
            let e = extend_mem(query, ref_seq, s);
            if seen.insert((e.read_start, e.ref_start, e.len)) {
                unimems.push(e);
            }
        }

        for chain in chain_mems(&unimems, is_fw, &chain_cfg) {
            candidates.push(MappingCandidate {
                tid,
                is_fw,
                chain,
                read_len: read.len() as i32,
            });
        }
    }
    candidates
}

/// Collect reference-MEM mapping candidates for one read against a loaded index.
///
/// Drop-in alternative to [`crate::collect::collect_read_mems`] that seeds with
/// reference-extended MEMs instead of sparse fixed-`k` anchors. Needs the
/// reference sequences (`refs`) to extend against.
pub fn collect_read_unimems<'idx, R: RefProvider>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    refs: &R,
    read: &[u8],
    is_left: bool,
    cfg: &MemCollectorConfig,
) -> Vec<MappingCandidate> {
    let k = index.k();
    let read_len = read.len() as i32;
    if read_len < k as i32 {
        return Vec::new();
    }
    let encoding = index.encoding();

    dispatch_on_k!(k, K => {
        let mut query = PiscemStreamingQuery::<K>::new(index.dict());
        hs.get_raw_hits_sketch::<K>(read, &mut query, cfg.skipping, is_left);
        let raw = if is_left { hs.left_hits() } else { hs.right_hits() };
        candidates_from_raw_hits_unimems(raw, &encoding, refs, read, k as i32, cfg)
    })
}

/// Groups `(tid, is_fw, Mem)` anchors by `(tid, is_fw)` in O(N + G log G)
/// (N anchors, G distinct targets) without sorting the anchors themselves.
///
/// A direct `tid`-indexed slot table assigns each distinct target a counter
/// pair (reverse, forward) in one load per anchor; only the G distinct target
/// ids are then sorted, prefix offsets laid out in ascending `(tid, is_fw)`
/// order (`false < true`), and the anchors scattered into one flat buffer.
/// The resulting group order equals that of sorting by the `(tid, is_fw)` key.
/// Within a group anchors keep their input order — irrelevant downstream:
/// [`merge_same_diagonal`]'s output depends only on the group's set of anchors.
///
/// This replaced an LSD radix sort of the 20-byte tuples, whose 2–3 passes
/// (transcript ids up to ~400k on multi-isoform references such as GENCODE)
/// dominated mapping time. Unlike a linear-probe bucketing, the cost per anchor
/// is O(1) regardless of how many targets a repetitive read hits.
#[derive(Default)]
struct TidGrouper {
    /// `tid -> 1 + index into targets` for the current read; 0 = absent.
    /// Sized to the reference count on first use (4 bytes per target per
    /// thread); only touched entries are reset after each read.
    slot_of: Vec<u32>,
    /// `(tid, [reverse, forward])`: per-orientation counts, then (after
    /// [`group`](Self::group)) per-orientation group END offsets in `grouped`.
    /// Sorted by `tid` after `group`; the reverse group of `targets[i]` starts
    /// where `targets[i - 1]`'s forward group ends (0 for `i == 0`), and its
    /// forward group starts at its reverse end.
    targets: Vec<(u32, [u32; 2])>,
    /// Anchors laid out group-contiguously.
    grouped: Vec<Mem>,
}

impl TidGrouper {
    #[cold]
    fn grow(&mut self, t: usize, num_refs: usize) {
        // Allocate the full reference count once rather than growing (and
        // freeing) a large buffer repeatedly; geometric growth only as a
        // fallback for an implausible count or an id beyond it.
        const MAX_PREALLOC: usize = 1 << 26;
        let len = if t < num_refs && num_refs <= MAX_PREALLOC {
            num_refs
        } else {
            (t + 1).next_power_of_two()
        };
        self.slot_of.resize(len, 0);
    }

    /// Bucket `flat` by `(tid, is_fw)`; see the type docs for the layout.
    /// `num_refs` (the number of targets) only sizes the index.
    fn group(&mut self, flat: &[(u32, bool, Mem)], num_refs: usize) {
        self.targets.clear();
        self.grouped.clear();
        if flat.is_empty() {
            return;
        }
        // Counting pass: resolve each anchor's target slot and count it.
        for &(tid, fw, _) in flat {
            let t = tid as usize;
            if t >= self.slot_of.len() {
                self.grow(t, num_refs);
            }
            let mut s = self.slot_of[t];
            if s == 0 {
                self.targets.push((tid, [0, 0]));
                s = self.targets.len() as u32;
                self.slot_of[t] = s;
            }
            self.targets[(s - 1) as usize].1[fw as usize] += 1;
        }
        // Order the G distinct targets (ids are distinct, so unstable is exact)
        // and turn counts into start cursors.
        self.targets.sort_unstable_by_key(|&(tid, _)| tid);
        let mut off = 0u32;
        for (i, (tid, c)) in self.targets.iter_mut().enumerate() {
            self.slot_of[*tid as usize] = i as u32 + 1;
            let (rev, fwd) = (c[0], c[1]);
            c[0] = off;
            c[1] = off + rev;
            off += rev + fwd;
        }
        // Scatter pass; afterwards each cursor sits at its group's end.
        // Placeholder fill: every slot is overwritten below (counts sum to N).
        self.grouped.resize(flat.len(), flat[0].2);
        for &(tid, fw, m) in flat {
            let s = (self.slot_of[tid as usize] - 1) as usize;
            let cur = &mut self.targets[s].1[fw as usize];
            self.grouped[*cur as usize] = m;
            *cur += 1;
        }
        for &(tid, _) in &self.targets {
            self.slot_of[tid as usize] = 0;
        }
    }
}

/// Reused per-thread scratch for [`candidates_from_raw_hits_true_unimems`]: the
/// flat `(tid, is_fw, extended uni-MEM)` buffer, its `(tid, is_fw)` grouper,
/// and the merged-anchor output buffer (no per-read map or per-group `Vec`).
#[derive(Default)]
struct UnimemScratch {
    flat: Vec<(u32, bool, Mem)>,
    grouper: TidGrouper,
    merged: Vec<Mem>,
}

thread_local! {
    static UNIMEM_SCRATCH: std::cell::RefCell<UnimemScratch> =
        std::cell::RefCell::new(UnimemScratch::default());
}

/// Collapse same-diagonal overlapping/abutting uni-MEMs into maximal anchors.
///
/// Extending each raw k-mer hit independently yields several heavily-overlapping
/// uni-MEMs on one diagonal (one per surviving seed). The minimap2-style chaining
/// gain counts a successor's new bases by its *start* offset, which undercounts a
/// longer anchor's end extension — so the chainer scores the overlapping group
/// below the single longest anchor and splits it, deflating chain coverage. We
/// instead merge them here (as pufferfish's one-uni-MEM-per-diagonal seeding
/// does), so the chain sees a single maximal anchor. Anchors on the same diagonal
/// separated by a real gap (a SNP run) are left distinct for the chainer; only
/// overlapping/touching runs are merged. `ext` is consumed (sorted) and the
/// result written to `out`.
fn merge_same_diagonal(ext: &mut [Mem], out: &mut Vec<Mem>) {
    out.clear();
    if ext.is_empty() {
        return;
    }
    // diagonal = ref_start - read_start; collinear (gap-free) anchors share it.
    // The result depends only on the set of anchors, not their input order or
    // how sort ties break: per diagonal the loop below takes the union of the
    // intervals, and anchors with equal keys share their start, so whichever
    // comes first, the next merges into it.
    ext.sort_unstable_by_key(|m| (m.ref_start - m.read_start, m.read_start));
    let mut cur = ext[0];
    for &m in &ext[1..] {
        let same_diag = (m.ref_start - m.read_start) == (cur.ref_start - cur.read_start);
        if same_diag && m.read_start <= cur.read_end() {
            let new_end = cur.read_end().max(m.read_end());
            cur = Mem::new(cur.read_start, cur.ref_start, new_end - cur.read_start);
        } else {
            out.push(cur);
            cur = m;
        }
    }
    out.push(cur);
}

/// Project then extend each seed into a **true uni-MEM** (extension clamped to
/// the seed's unitig), dedup, and chain. Mirrors pufferfish's unitig-constrained
/// uni-MEM seeding: a junction-straddling match yields one uni-MEM per unitig,
/// which the chainer then stitches together. Contrast with
/// [`candidates_from_raw_hits_unimems`], whose reference MEMs bridge the
/// junction into a single anchor.
pub fn candidates_from_raw_hits_true_unimems<R: RefProvider>(
    raw_hits: &[(i32, ProjectedHits<'_>)],
    encoding: &EntryEncoding,
    refs: &R,
    read: &[u8],
    k: i32,
    cfg: &MemCollectorConfig,
) -> Vec<MappingCandidate> {
    let read_len = read.len() as i32;
    let rc = revcomp(read);
    let mut chain_cfg = cfg.chain.clone();
    chain_cfg.seed_len = k;

    UNIMEM_SCRATCH.with(|cell| {
        let UnimemScratch {
            flat,
            grouper,
            merged,
        } = &mut *cell.borrow_mut();
        flat.clear();
        // The previous (non-skipped) raw hit: (contig id, orientation on the
        // contig, per-orientation uni-MEM cache, span index of the first
        // occurrence of each orientation). See the duplicate check below.
        type Cache = [Option<(i32, i32, i32)>; 2];
        let mut prev: Option<(u32, bool, Cache, [usize; 2])> = None;
        // For each queried k-mer (one `ProjectedHits` = one unitig + all its
        // reference occurrences), extend the seed to its uni-MEM AT MOST ONCE PER
        // ORIENTATION, then project that uni-MEM to every occurrence.
        //
        // Why this is exact: a uni-MEM is a match against the *unitig*, whose
        // sequence is identical wherever it occurs, and extension is clamped to the
        // unitig — so the matched span is occurrence-independent. `decode_hit` gives
        // `is_fw = (contig_fw == contig_orientation)` with `contig_orientation` fixed
        // for this hit, so `is_fw` determines `contig_fw`; hence all occurrences with
        // the same `is_fw` share one contig-relative offset and we only differ the
        // reference start by the occurrence's `base_pos`. This replaces the previous
        // once-per-occurrence extension (the redundant cross-transcript re-walk).
        for (read_pos, phit) in raw_hits {
            if phit.num_hits() > cfg.max_hit_occ {
                continue;
            }
            let contig_len = phit.contig_len() as i32;
            let span = phit.ref_range();
            // Extend this hit's seed at occurrence `entry` to its uni-MEM, as
            // (read_start, contig_relative_start, len).
            let extend_at = |entry: u64| -> (bool, (i32, i32, i32)) {
                let tid = encoding.transcript_id(entry);
                let base_pos = encoding.pos(entry) as i32;
                let rp = phit.decode_hit(entry, encoding);
                let rs = if rp.is_fw {
                    *read_pos
                } else {
                    read_len - (*read_pos + k)
                };
                let query: &[u8] = if rp.is_fw { read } else { &rc };
                let seed = Mem::new(rs, rp.pos as i32, k);
                let e = extend_mem_within(
                    query,
                    refs.ref_seq(tid),
                    seed,
                    base_pos,
                    base_pos + contig_len,
                );
                (rp.is_fw, (e.read_start, e.ref_start - base_pos, e.len))
            };
            // Per-orientation cache of the extended uni-MEM (index 1 = forward),
            // and the span index of the occurrence it was computed at (always
            // the first occurrence with that orientation).
            let mut cache: Cache = [None, None];
            let mut first = [usize::MAX; 2];
            // Duplicate hit: consecutive k-mers on the same unitig (same contig,
            // same orientation on it) usually extend to the same uni-MEM. Then
            // every anchor this hit would push — same occurrences, same
            // `is_fw`, same `(read_start, crel, len)` — is already in `flat`, and
            // a repeated anchor is a no-op for `merge_same_diagonal` (an interval
            // union), so the hit is skipped. Each cached value is recomputed
            // exactly as the loop below would (same first occurrence), so the
            // comparison is exact; on a mismatch the values are reused.
            if let Some((pc, po, pcache, pfirst)) = &prev {
                if *pc == phit.contig_id() && *po == phit.hit_fw_on_contig() {
                    let mut same = true;
                    for slot in 0..2 {
                        if let Some(pv) = pcache[slot] {
                            let (_, v) = extend_at(span.get(pfirst[slot]));
                            cache[slot] = Some(v);
                            first[slot] = pfirst[slot];
                            same &= v == pv;
                        }
                    }
                    if same {
                        continue;
                    }
                }
            }
            for (idx, entry) in span.iter().enumerate() {
                let tid = encoding.transcript_id(entry);
                let base_pos = encoding.pos(entry) as i32;
                let is_fw = encoding.orientation(entry) == phit.hit_fw_on_contig();
                let slot = is_fw as usize;
                let (read_start, crel, len) = match cache[slot] {
                    Some(v) => v,
                    None => {
                        let (fw, v) = extend_at(entry);
                        debug_assert_eq!(fw, is_fw);
                        cache[slot] = Some(v);
                        first[slot] = idx;
                        v
                    }
                };
                flat.push((tid, is_fw, Mem::new(read_start, base_pos + crel, len)));
            }
            prev = Some((phit.contig_id(), phit.hit_fw_on_contig(), cache, first));
        }

        // Group by (tid, is_fw) in ascending key order, then merge + chain each.
        grouper.group(flat, refs.num_refs());
        let mut candidates = Vec::new();
        let mut start = 0usize;
        for &(tid, ends) in &grouper.targets {
            for (fw, end) in [(false, ends[0]), (true, ends[1])] {
                let end = end as usize;
                if end > start {
                    // Dedup duplicate uni-MEMs (same hit projected from several
                    // k-mers) and stitch same-diagonal uni-MEMs that abut across
                    // unitig junctions.
                    merge_same_diagonal(&mut grouper.grouped[start..end], merged);
                    for chain in chain_mems(merged, fw, &chain_cfg) {
                        candidates.push(MappingCandidate {
                            tid,
                            is_fw: fw,
                            chain,
                            read_len,
                        });
                    }
                }
                start = end;
            }
        }
        candidates
    })
}

/// Collect true uni-MEM mapping candidates (unitig-constrained extension) for
/// one read against a loaded index. See [`candidates_from_raw_hits_true_unimems`].
pub fn collect_read_true_unimems<'idx, R: RefProvider>(
    index: &'idx ReferenceIndex,
    hs: &mut HitSearcher<'idx>,
    refs: &R,
    read: &[u8],
    is_left: bool,
    cfg: &MemCollectorConfig,
) -> Vec<MappingCandidate> {
    let k = index.k();
    let read_len = read.len() as i32;
    if read_len < k as i32 {
        return Vec::new();
    }
    let encoding = index.encoding();

    dispatch_on_k!(k, K => {
        let mut query = PiscemStreamingQuery::<K>::new(index.dict());
        hs.get_raw_hits_sketch::<K>(read, &mut query, cfg.skipping, is_left);
        let raw = if is_left { hs.left_hits() } else { hs.right_hits() };
        candidates_from_raw_hits_true_unimems(raw, &encoding, refs, read, k as i32, cfg)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use piscem_rs::index::contig_table::{ContigTable, ContigTableBuilder};
    use piscem_rs::mapping::projected_hits::ProjectedHits;

    /// Reference grouping: stable sort by `(tid, is_fw)`, split into runs.
    fn sorted_groups(flat: &[(u32, bool, Mem)]) -> Vec<(u32, bool, Vec<Mem>)> {
        let mut v = flat.to_vec();
        v.sort_by_key(|&(t, f, _)| (t, f));
        let mut out: Vec<(u32, bool, Vec<Mem>)> = Vec::new();
        for (t, f, m) in v {
            match out.last_mut() {
                Some((lt, lf, ms)) if *lt == t && *lf == f => ms.push(m),
                _ => out.push((t, f, vec![m])),
            }
        }
        out
    }

    /// Groups as laid out by `TidGrouper` (after `group`).
    fn grouper_groups(g: &TidGrouper) -> Vec<(u32, bool, Vec<Mem>)> {
        let mut out = Vec::new();
        let mut start = 0usize;
        for &(tid, ends) in &g.targets {
            for (fw, end) in [(false, ends[0]), (true, ends[1])] {
                let end = end as usize;
                if end > start {
                    out.push((tid, fw, g.grouped[start..end].to_vec()));
                }
                start = end;
            }
        }
        assert_eq!(start, g.grouped.len());
        out
    }

    #[test]
    /// Anchors must end up grouped by (transcript, orientation) in ascending key
    /// order: the chainer is run per group, so a mis-grouping would mix
    /// unrelated anchors into one candidate, and the group order is the
    /// candidate order. Within a group the grouper keeps input order, so it must
    /// equal a stable sort exactly — across reads that reuse the scratch.
    fn tid_grouper_matches_stable_sort() {
        let mut g = TidGrouper::default();
        g.group(&[], 0);
        assert!(g.targets.is_empty() && g.grouped.is_empty());
        let mut x = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = |m: u64| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (x >> 33) % m
        };
        for round in 0..200 {
            let n = next(600) as usize;
            // Small and large id ranges (forces on-demand growth of the index).
            let tid_range = if round % 3 == 0 { 7 } else { 500_000 };
            let flat: Vec<(u32, bool, Mem)> = (0..n)
                .map(|_| {
                    (
                        next(tid_range) as u32,
                        next(2) == 1,
                        Mem::new(next(100) as i32, next(10_000) as i32, 1 + next(50) as i32),
                    )
                })
                .collect();
            // The size hint is only a hint: exercise growth past it.
            g.group(&flat, (tid_range / 2) as usize);
            assert_eq!(grouper_groups(&g), sorted_groups(&flat), "round {round}");
            assert!(g.slot_of.iter().all(|&s| s == 0), "index not reset");
        }
    }

    #[test]
    /// Overlapping anchors on one diagonal describe the same match and must merge
    /// into a single anchor.
    fn merge_collapses_same_diagonal_overlaps() {
        // The read-100000/ENST00000518938.1 geometry: three overlapping uni-MEMs
        // on diagonal 76 (one per surviving seed). They must collapse to a single
        // maximal anchor read[0,52) so the chain coverage matches the equivalent
        // fixed-k-mer chain (52), instead of the chainer splitting them.
        let mut ext = vec![
            Mem::new(0, 76, 34),
            Mem::new(4, 80, 34),
            Mem::new(8, 84, 44),
        ];
        let mut out = Vec::new();
        merge_same_diagonal(&mut ext, &mut out);
        assert_eq!(out.len(), 1, "same-diagonal overlaps must merge into one");
        assert_eq!(
            (out[0].read_start, out[0].ref_start, out[0].len),
            (0, 76, 52)
        );
    }

    #[test]
    /// But anchors on different diagonals, or with a genuine gap between them, must
    /// stay separate for the chainer to reason about.
    fn merge_keeps_distinct_diagonals_and_real_gaps() {
        // Same diagonal but a real read gap (a SNP run) stays distinct; a
        // different diagonal stays distinct — only overlapping/touching runs merge.
        let mut ext = vec![
            Mem::new(0, 100, 20),  // diag 100, read[0,20)
            Mem::new(10, 110, 20), // diag 100, overlaps prev -> merges to read[0,30)
            Mem::new(40, 140, 10), // diag 100, gap after read 30 -> distinct
            Mem::new(50, 200, 20), // diag 150 -> distinct
        ];
        let mut out = Vec::new();
        merge_same_diagonal(&mut ext, &mut out);
        assert_eq!(out.len(), 3);
        assert_eq!((out[0].read_start, out[0].len), (0, 30));
    }

    /// Minimal in-memory reference store for the extension tests.
    struct TestRefs {
        seqs: Vec<Vec<u8>>,
    }
    impl RefProvider for TestRefs {
        fn num_refs(&self) -> usize {
            self.seqs.len()
        }
        fn ref_seq(&self, tid: u32) -> &[u8] {
            &self.seqs[tid as usize]
        }
    }

    fn gen_seq(n: usize, seed: u64) -> Vec<u8> {
        const B: [u8; 4] = *b"ACGT";
        let mut x = seed;
        let mut s = Vec::with_capacity(n);
        for _ in 0..n {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            s.push(B[((x >> 33) & 3) as usize]);
        }
        s
    }

    #[test]
    /// Extension must reach the full extent of the exact match, which is the whole
    /// reason for extending.
    fn extend_grows_seed_to_full_exact_match() {
        // read == ref[50..130]; a length-31 seed at read 0 / ref 50 should grow
        // to cover the whole 80 bp exact match.
        let reference = gen_seq(300, 7);
        let read = reference[50..130].to_vec();
        let seed = Mem::new(0, 50, 31);
        let e = extend_mem(&read, &reference, seed);
        assert_eq!(e.read_start, 0);
        assert_eq!(e.ref_start, 50);
        assert_eq!(e.len, 80, "seed should extend to the full 80 bp match");
    }

    #[test]
    /// The uni-MEM variant must stop at the unitig edge, matching pufferfish, so
    /// one anchor per unitig is produced and the chainer stitches them.
    fn uni_mem_extension_stops_at_unitig_boundary() {
        // read == ref[50..130] (80 bp exact). A seed at read 0 / ref 50:
        // - reference MEM extends to the full 80 bp.
        // - uni-MEM clamped to a unitig spanning ref [50..100) stops at ref 100,
        //   i.e. 50 bp covered (read [0..50)).
        let reference = gen_seq(300, 7);
        let read = reference[50..130].to_vec();
        let seed = Mem::new(0, 50, 31);
        let ref_mem = extend_mem(&read, &reference, seed);
        assert_eq!(ref_mem.len, 80, "reference MEM crosses freely");
        let uni = extend_mem_within(&read, &reference, seed, 50, 100);
        assert_eq!(uni.ref_start, 50);
        assert_eq!(uni.ref_end(), 100, "uni-MEM must stop at the unitig end");
        assert_eq!(uni.len, 50);
    }

    #[test]
    /// Extension must stop at the first mismatch — an anchor is by definition an
    /// exact match.
    fn extend_stops_at_mismatch() {
        let mut reference = gen_seq(300, 9);
        let read = reference[50..130].to_vec();
        // Introduce a mismatch 10 bases to the right of the seed end.
        let break_at = 50 + 31 + 10;
        reference[break_at] = if reference[break_at] == b'A' {
            b'C'
        } else {
            b'A'
        };
        let e = extend_mem(&read, &reference, Mem::new(0, 50, 31));
        // right edge stops just before the mismatch; left edge reaches read 0.
        assert_eq!(e.read_start, 0);
        assert_eq!(e.ref_end(), break_at as i32);
    }

    /// One contig (id 0) occurring forward on ref 0 at position 50.
    fn fw_table_at(base: u32) -> ContigTable {
        let mut b = ContigTableBuilder::new(1, 10_000, 1);
        b.add_occurrence(0, 0, base, true);
        b.build()
    }

    #[test]
    /// Several seeds within one exact match must collapse to one anchor rather than
    /// being counted repeatedly as evidence.
    fn unimems_collapse_colinear_seeds_into_one_anchor() {
        // Reference where read (= ref[50..130]) matches exactly. Two forward
        // seeds at read offsets 0 and 10 both lie on the same diagonal and must
        // collapse to a single uni-MEM covering the whole read.
        let reference = gen_seq(300, 7);
        let read = reference[50..130].to_vec();
        let refs = TestRefs {
            seqs: vec![reference],
        };
        let table = fw_table_at(50);
        let enc = table.encoding();
        let k = 31;
        let contig_len = 200;
        let raw: Vec<(i32, ProjectedHits)> = [0i32, 10]
            .iter()
            .map(|&p| {
                (
                    p,
                    ProjectedHits::new(
                        0,
                        p as u32,
                        true,
                        contig_len,
                        0,
                        k as u32,
                        table.contig_entries(0),
                    ),
                )
            })
            .collect();

        let cfg = MemCollectorConfig::default();
        let cands = candidates_from_raw_hits_unimems(&raw, &enc, &refs, &read, k, &cfg);
        assert_eq!(cands.len(), 1);
        let c = &cands[0];
        assert!(c.is_fw);
        assert_eq!(
            c.chain.mems.len(),
            1,
            "colinear seeds should collapse to one uni-MEM"
        );
        assert_eq!(
            c.chain.covered_read_bases(),
            80,
            "uni-MEM covers the whole read"
        );
    }
}
