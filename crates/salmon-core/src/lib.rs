//! `salmon-core`: shared foundational types for the Rust port of salmon.
//!
//! This crate holds the vocabulary types used across the mapping, alignment,
//! modeling, and inference crates: library formats, mate status, the
//! [`Transcript`] record with its concurrent count/mass state, log-space math
//! helpers, and the shared error type. It deliberately has no heavy
//! dependencies so every other crate can depend on it cheaply.

pub mod atomic;
pub mod error;
pub mod genemap;
pub mod libtype;
pub mod mate;
pub mod math;
pub mod progress;
pub mod quantmerge;
pub mod refprovider;
pub mod transcript;

pub use atomic::AtomicF64;
pub use libtype::{compatible_paired, compatible_single, is_compatible, observed_paired_format};
pub use progress::ProgressCounters;
pub use refprovider::RefProvider;

pub use error::{Result, SalmonError};
pub use libtype::{LibraryFormat, ReadOrientation, ReadStrandedness, ReadType};
pub use mate::MateStatus;
pub use transcript::Transcript;
