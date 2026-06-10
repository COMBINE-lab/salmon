//! Uni-MEM extraction: extend sparse k-mer anchors into maximal exact matches.
//!
//! The default [`collect`](crate::collect) path chains the *sparse, fixed-length*
//! k-mer anchors that piscem's skipping streaming query emits — one length-`k`
//! anchor per unitig transition. pufferfish/salmon instead seed with **uni-MEMs**:
//! each shared k-mer is extended along the unitig (here, against the reference)
//! into a maximal exact match before chaining, and colinear k-mers that belong to
//! the same exact stretch collapse into a single anchor.
//!
//! This module implements that uni-MEM seeding as an *additive* alternative to
//! the sparse-k-mer path (selected by [`MapConfig::use_unimems`]). It exists to
//! answer an experimental question: does the compacted-uni-MEM vs. sparse-k-mer
//! seed representation account for the read-placement differences vs. C++ salmon?
//! The extension is exact (no mismatches), in the same query frame the chainer
//! uses — the forward read for forward groups, its reverse complement for reverse
//! groups — so the produced anchors drop straight into [`chain_mems`].
//!
//! [`MapConfig::use_unimems`]: crate::mapper::MapConfig::use_unimems

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

/// Extend a seed anchor to a maximal exact match against `query`/`ref_seq`.
///
/// Both are in the query frame: the caller passes the forward read (forward
/// group) or its reverse complement (reverse group), and `m`'s coordinates are
/// already expressed in that frame. The match grows left and right while bases
/// agree and stay in bounds; the result is the uni-MEM containing the seed.
pub fn extend_mem(query: &[u8], ref_seq: &[u8], mut m: Mem) -> Mem {
    let qlen = query.len() as i32;
    let reflen = ref_seq.len() as i32;

    // Extend left along the shared diagonal.
    while m.read_start > 0
        && m.ref_start > 0
        && base_eq(
            query[(m.read_start - 1) as usize],
            ref_seq[(m.ref_start - 1) as usize],
        )
    {
        m.read_start -= 1;
        m.ref_start -= 1;
        m.len += 1;
    }

    // Extend right along the shared diagonal.
    while m.read_end() < qlen
        && m.ref_end() < reflen
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
            candidates.push(MappingCandidate { tid, is_fw, chain });
        }
    }
    candidates
}

/// Collect uni-MEM mapping candidates for one read against a loaded index.
///
/// Drop-in alternative to [`crate::collect::collect_read_mems`] that seeds with
/// extended uni-MEMs instead of sparse fixed-`k` anchors. Needs the reference
/// sequences (`refs`) to extend against.
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

#[cfg(test)]
mod tests {
    use super::*;
    use piscem_rs::index::contig_table::{ContigTable, ContigTableBuilder};
    use piscem_rs::mapping::projected_hits::ProjectedHits;

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
    fn extend_stops_at_mismatch() {
        let mut reference = gen_seq(300, 9);
        let read = reference[50..130].to_vec();
        // Introduce a mismatch 10 bases to the right of the seed end.
        let break_at = 50 + 31 + 10;
        reference[break_at] = if reference[break_at] == b'A' { b'C' } else { b'A' };
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
                    ProjectedHits::new(0, p as u32, true, contig_len, 0, k as u32, table.contig_entries(0)),
                )
            })
            .collect();

        let cfg = MemCollectorConfig::default();
        let cands = candidates_from_raw_hits_unimems(&raw, &enc, &refs, &read, k, &cfg);
        assert_eq!(cands.len(), 1);
        let c = &cands[0];
        assert!(c.is_fw);
        assert_eq!(c.chain.mems.len(), 1, "colinear seeds should collapse to one uni-MEM");
        assert_eq!(c.chain.covered_read_bases(), 80, "uni-MEM covers the whole read");
    }
}
