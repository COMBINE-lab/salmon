//! `salmon-core`: shared foundational types for the Rust port of salmon.
//!
//! This crate holds the vocabulary types used across the mapping, alignment,
//! modeling, and inference crates: library formats, mate status, the
//! [`Transcript`] record with its concurrent count/mass state, log-space math
//! helpers, and the shared error type. It deliberately has no heavy
//! dependencies so every other crate can depend on it cheaply.
//!
//! # How the workspace fits together
//!
//! A salmon run is a pipeline, and each stage is its own crate:
//!
//! ```text
//!   salmon-cli      command line, argument validation
//!        │
//!   salmon-index    build/load the reference index
//!        │
//!   salmon-map      place reads on transcripts (reads mode)
//!   salmon-align    read placements out of a BAM instead (alignment mode)
//!        │
//!   salmon-model    fragment length, sequence/GC/positional bias, error model
//!   salmon-eqclass  collapse fragments into equivalence classes
//!        │
//!   salmon-infer    EM / VBEM / Gibbs — turn classes into abundances
//!        │
//!   salmon-quant    drives all of the above and writes quant.sf
//!   salmon-rad      optional compact mapping-record file format
//! ```
//!
//! Everything above depends on this crate, and this crate depends on none of
//! them: it is the shared bottom of the stack, which is why it must stay light.

pub mod atomic;
pub mod diagnostics;
pub mod error;
pub mod genemap;
pub mod libtype;
pub mod mate;
pub mod math;
pub mod progress;
pub mod quantmerge;
pub mod refprovider;
pub mod timing;
pub mod transcript;

// Re-exports: the modules above are the real homes of these items, but naming
// them here lets other crates write `salmon_core::Transcript` instead of
// `salmon_core::transcript::Transcript`. It also lets a type move between
// modules later without breaking every caller.
pub use atomic::AtomicF64;
pub use diagnostics::{input_diagnostics, peak_rss_kb, Diagnostic};
pub use libtype::{compatible_paired, compatible_single, is_compatible, observed_paired_format};
pub use progress::ProgressCounters;
pub use refprovider::RefProvider;
pub use timing::PhaseTimer;

pub use error::{Result, SalmonError};
pub use libtype::{LibraryFormat, ReadOrientation, ReadStrandedness, ReadType};
pub use mate::MateStatus;
pub use transcript::Transcript;
