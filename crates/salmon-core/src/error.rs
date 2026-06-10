//! Shared error type for the salmon crates.

use thiserror::Error;

/// Errors surfaced by the core salmon crates.
#[derive(Error, Debug)]
pub enum SalmonError {
    #[error("unknown library format string: {0}")]
    UnknownLibraryFormat(String),

    #[error("invalid value for {field}: {value}")]
    InvalidValue { field: String, value: String },

    #[error("i/o error: {0}")]
    Io(#[from] std::io::Error),

    #[error("{0}")]
    Other(String),
}

/// Convenience result alias.
pub type Result<T> = std::result::Result<T, SalmonError>;
