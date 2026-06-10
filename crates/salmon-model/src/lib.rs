//! `salmon-model`: statistical models used during quantification.
//!
//! Currently provides the fragment-length distribution ([`fld`]) and automatic
//! library-type detection ([`libdetect`]). Bias models (sequence-specific, GC,
//! positional) and the alignment error model are added in later phases.

pub mod fld;
pub mod libdetect;
pub mod seqbias;

pub use fld::FragmentLengthDistribution;
pub use libdetect::LibraryTypeDetector;
pub use seqbias::{build_expected, corrected_effective_length, SBModel};
