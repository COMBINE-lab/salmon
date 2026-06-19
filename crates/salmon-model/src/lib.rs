//! `salmon-model`: statistical models used during quantification.
//!
//! Currently provides the fragment-length distribution ([`fld`]) and automatic
//! library-type detection ([`libdetect`]). Bias models (sequence-specific, GC,
//! positional) and the alignment error model are added in later phases.

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
pub use fld::{ambig_frag_log_prob, smoothed_effective_length, FragmentLengthDistribution};
pub use gcbias::{
    build_expected_gc, gc_corrected_effective_length, gc_desc, gc_prefix, gc_ratio, GcFragModel,
    GcRank, GcStore, GcView, GC_SAMP_STRIDE,
};
pub use libdetect::LibraryTypeDetector;
pub use posbias::{
    compute_length_quantiles, length_class_index, SimplePosBias, NUM_LENGTH_CLASSES, NUM_POS_BINS,
};
pub use seqbias::{build_expected, corrected_effective_length, SBModel};
