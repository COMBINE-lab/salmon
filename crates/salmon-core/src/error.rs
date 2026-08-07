//! Shared error type for the salmon crates.
//!
//! Rust has no exceptions: a function that can fail returns a `Result`, which is
//! either `Ok(value)` or `Err(error)`, and the caller must deal with both. This
//! module defines the error side of that for the core crates.

use thiserror::Error;

/// Errors surfaced by the core salmon crates.
///
/// An `enum` here is a closed list of the things that can go wrong; callers can
/// match on the specific case (e.g. to retry only on I/O) rather than parsing a
/// message string.
///
/// The `#[error("…")]` lines are from the `thiserror` crate: they generate the
/// human-readable message for each case, so `{0}` is the first field of that
/// variant and `{field}` is a named one.
#[derive(Error, Debug)]
pub enum SalmonError {
    /// A `-l` library-type string that does not name a known format.
    #[error("unknown library format string: {0}")]
    UnknownLibraryFormat(String),

    /// A parameter outside its permitted range, named so the message can point
    /// at the offending option.
    #[error("invalid value for {field}: {value}")]
    InvalidValue { field: String, value: String },

    /// Anything the filesystem or OS reported. `#[from]` generates the
    /// conversion from `std::io::Error`, which is what lets call sites write
    /// `File::open(p)?` and have the error wrapped automatically.
    #[error("i/o error: {0}")]
    Io(#[from] std::io::Error),

    /// Escape hatch for one-off messages that do not deserve their own variant.
    #[error("{0}")]
    Other(String),
}

/// Convenience result alias, so signatures read `-> Result<T>` instead of
/// `-> std::result::Result<T, SalmonError>`.
pub type Result<T> = std::result::Result<T, SalmonError>;
