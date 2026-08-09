//! `salmon-map`: selective-alignment mapping for the salmon Rust port.
//!
//! # What "mapping" has to decide
//!
//! For each sequenced fragment, which transcripts could it have come from, and
//! how well does it fit each one? Answering that for hundreds of millions of
//! fragments against a hundred thousand transcripts is the bulk of a salmon run,
//! so the work is staged from cheap and approximate to expensive and exact:
//!
//! 1. **Seed.** Look up the fragment's k-mers in the index. Each hit says "this
//!    31-base window occurs here". Cheap, and most transcripts are eliminated
//!    outright because they share no k-mer at all.
//! 2. **Anchor.** Extend the k-mer hits into *MEMs*: maximal stretches that match
//!    exactly. A handful of long anchors is a much better summary of the
//!    evidence than a scatter of k-mer hits.
//! 3. **Chain.** Group anchors that could belong to one alignment — same
//!    transcript, consistent order, plausible gaps — into candidate chains. This
//!    is where a real placement is separated from coincidental shared k-mers.
//! 4. **Align.** Only for surviving candidates, compute an actual base-level
//!    alignment score. This is the expensive step, which is why so much is done
//!    to avoid reaching it.
//! 5. **Pair, score, filter.** Combine the two mates, discard placements that are
//!    too poor relative to the best, drop decoys, and emit weighted mappings.
//!
//! *Selective alignment* is the name for this whole approach: unlike
//! pseudoalignment it does verify candidates by alignment, but unlike a
//! traditional aligner it only ever aligns the few candidates chaining could not
//! rule out. [`sketch`] provides the alignment-free path for when even that is
//! more than a run needs.
//!
//! This crate implements the selective-alignment path: extracting MEM anchors
//! from piscem k-mer hits, chaining them into colinear mapping candidates,
//! recovering orphans, and validating candidates with alignment. It is built
//! incrementally:
//!
//! - [`mem`] — MEM anchor type and coverage helpers.
//! - [`chain`] — colinear MEM chaining DP (the mapping-candidate core).
//! - [`collect`] — MEM collection over piscem k-mer hits + per-target chaining.
//! - [`align`] — ksw2 validation of a chain against a reference window.
//! - [`pair`] — paired-end `joinReadsAndFilter`.
//! - [`score`] — per-transcript dedup, decoy filtering, and eq-class weighting.
//! - [`mapper`] — the assembled selective-alignment mapper (single & paired,
//!   with orphan recovery).
//! - [`sketch`] — the alignment-free pseudoalignment-only mode.
//! - `refprovider` (re-exported from `salmon-core`) — reference-sequence
//!   access for alignment.

pub mod align;
pub mod chain;
pub mod collect;
pub mod extend;
pub mod mapper;
pub mod mem;
pub mod pair;
pub mod score;
pub mod sketch;

pub use align::{align_chain, align_in_window, perfect_score, AlignConfig, Alignment};
pub use chain::{chain_mems, ChainConfig, MemChain};
pub use collect::{
    best_per_target, candidates_from_raw_hits, collect_read_mems, MappingCandidate,
    MemCollectorConfig,
};
pub use extend::{
    candidates_from_raw_hits_true_unimems, candidates_from_raw_hits_unimems,
    collect_read_true_unimems, collect_read_unimems, extend_mem, extend_mem_within,
};
pub use mapper::{
    debug_best_mapping, map_read_pair, map_read_pair_into, map_single_read, map_single_read_into,
    take_last_map_stats, DebugMapping, MapConfig, MapStats, SeedMode,
};
pub use mem::Mem;
pub use pair::{join_reads_and_filter, JointMapping, PairingConfig};
pub use salmon_core::RefProvider;
pub use score::{
    filter_sketch_decoys, finalize_mappings, DedupScratch, RawMapping, ScoreConfig, ScoredMapping,
};
pub use sketch::{
    map_read_pair_sketch, map_read_pair_sketch_into, map_single_read_sketch,
    map_single_read_sketch_into, SketchScratch,
};
