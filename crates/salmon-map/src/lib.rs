//! `salmon-map`: selective-alignment mapping for the salmon Rust port.
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
//! - [`refprovider`] — reference-sequence access for alignment.

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
    debug_best_mapping, map_read_pair, map_single_read, take_last_map_stats, DebugMapping,
    MapConfig, MapStats, SeedMode,
};
pub use mem::Mem;
pub use pair::{join_reads_and_filter, JointMapping, PairingConfig};
pub use salmon_core::RefProvider;
pub use score::{finalize_mappings, RawMapping, ScoreConfig, ScoredMapping};
pub use sketch::{map_read_pair_sketch, map_single_read_sketch};
