//! `salmon-model`: statistical models used during quantification.
//!
//! Currently provides the fragment-length distribution ([`fld`]) and automatic
//! library-type detection ([`libdetect`]). Bias models (sequence-specific, GC,
//! positional) and the alignment error model are added in later phases.

/// Fixed-point scale for deterministic bias-model mass accumulation. Bias
/// observed models sum per-fragment posterior masses (each in `[0,1]`) into
/// bins; an f64 `+=` is non-associative, so the accumulated model — and hence
/// bias correction — varies with the worker-thread fragment partition (thread
/// count). Accumulating `round_down(mass * BIAS_WEIGHT_SCALE)` as integers makes
/// the sum associative (order/thread-count independent). `2^20` keeps the
/// per-contribution resolution at ~1e-6 (ample for these coarse models) while a
/// bin total (~num_fragments × 2^20) stays far below `u64::MAX`.
pub const BIAS_WEIGHT_SCALE: f64 = (1u64 << 20) as f64;

/// Quantize a bias mass (`[0,∞)`, typically a `[0,1]` posterior) to the
/// fixed-point integer accumulator. Truncates (rounds toward zero) — cheaper
/// than `round()` and the ≤1-ULP downward bias is negligible and cancels under
/// the model's normalization.
#[inline]
pub fn bias_mass_to_fp(mass: f64) -> u64 {
    (mass * BIAS_WEIGHT_SCALE) as u64
}

pub mod bias;
pub mod dumps;
pub mod fld;
pub mod gcbias;
pub mod libdetect;
pub mod posbias;
pub mod seqbias;
pub mod spline;

pub use bias::{
    build_expected_pos, corrected_effective_length_full, positional_factor, BiasInputs,
};
pub use fld::FragLengthSource;
pub use fld::{
    ambig_frag_log_prob, smoothed_effective_length, DiscreteFld, FragmentLengthDistribution,
};
pub use gcbias::{
    build_expected_gc, gc_corrected_effective_length, gc_desc, gc_prefix, gc_ratio, GcFragModel,
    GcRank, GcStore, GcView, GC_SAMP_STRIDE,
};
pub use libdetect::{infer_format_from_counts, LibraryTypeDetector};
pub use posbias::{
    compute_length_quantiles, length_class_index, SimplePosBias, NUM_LENGTH_CLASSES, NUM_POS_BINS,
};
pub use seqbias::{build_expected, corrected_effective_length, LogBiasTable, SBModel};
