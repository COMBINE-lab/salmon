//! `salmon-model`: statistical models used during quantification.
//!
//! Currently provides the fragment-length distribution ([`fld`]) and automatic
//! library-type detection ([`libdetect`]). Bias models (sequence-specific, GC,
//! positional) and the alignment error model are added in later phases.

pub mod fld;
pub mod gcbias;
pub mod libdetect;
pub mod seqbias;

pub use fld::{smoothed_effective_length, FragmentLengthDistribution};
pub use gcbias::{
    build_expected_gc, gc_corrected_effective_length, gc_desc, gc_prefix, gc_ratio, GcFragModel,
    GC_SAMP_STRIDE,
};
pub use libdetect::LibraryTypeDetector;
pub use seqbias::{build_expected, corrected_effective_length, SBModel};
