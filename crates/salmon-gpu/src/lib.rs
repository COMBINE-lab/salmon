//! `salmon-gpu`: optional GPU backend for salmon's full-length selective
//! alignment.
//!
//! The hot path of mapping is the banded Smith-Waterman that validates each
//! candidate placement. In `--fullLengthAlignment` mode that is one banded
//! affine-gap DP per candidate, which is regular and independent across
//! candidates: ideal to batch onto a GPU. This crate provides
//!   * [`reference`] — a pure-Rust banded DP that is the bit-exact oracle for
//!     the GPU kernel (and a near line-by-line match of the WGSL shader), and
//!   * (under the `gpu` feature) the wgpu/WGSL backend itself, which runs the
//!     same DP on Metal or Vulkan and reproduces the reference scores.
//!
//! Scores are integer-valued, so a faithful backend is deterministic and
//! reproduces the CPU path bit-for-bit, preserving salmon's guarantee of
//! identical output regardless of where the work runs.

pub mod aligner;
pub mod reference;

#[cfg(feature = "gpu")]
pub mod gpu;

pub use aligner::RefAligner;
pub use reference::{banded_extz_score, banded_extz_score_dna5, dna5, BandedParams, NEG_INF};

#[cfg(feature = "gpu")]
pub use aligner::GpuAligner;
#[cfg(feature = "gpu")]
pub use gpu::{GpuContext, GpuTask};
