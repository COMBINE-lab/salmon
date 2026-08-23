//! `salmon-core`: shared foundational types for the Rust port of salmon.
//!
//! This crate holds the vocabulary types used across the mapping, alignment,
//! modeling, and inference crates: library formats, mate status, log-space math
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
pub mod compress;
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

// Re-exports: the modules above are the real homes of these items, but naming
// them here lets other crates write `salmon_core::MateStatus` instead of
// `salmon_core::mate::MateStatus`. It also lets a type move between
// modules later without breaking every caller.
pub use atomic::AtomicF64;
pub use diagnostics::{
    free_disk_bytes, input_diagnostics, peak_rss_kb, Diagnostic, MetaFieldSource, MissingMetaField,
};
pub use libtype::{
    compatible_paired, compatible_single, is_compatible, observed_paired_format,
    observed_single_format, summarize_lib_format_counts, LibFormatCountsArray, LibFormatCountsFile,
    LibFormatSummary, NUM_LIB_FORMATS,
};
pub use progress::ProgressCounters;
pub use refprovider::RefProvider;
pub use timing::PhaseTimer;

pub use error::{Result, SalmonError};
pub use libtype::{LibraryFormat, ReadOrientation, ReadStrandedness, ReadType};
pub use mate::MateStatus;

/// The output-row indices for a reference table laid out as
/// `[transcripts][decoys][shorts]`: quantified transcripts `[0, first_decoy)`,
/// the decoy block skipped, and sub-`k` "short" transcripts (reported with zero
/// counts) after it. `first_decoy_index == None` means no decoys — every
/// reference is a row.
///
/// Shared by the reads-mode writer and the RAD-requant writer so `quant.sf`
/// (and the files positionally aligned with it) contain the same row set in
/// every mode (#1140): decoys are references the index carries for mapping
/// specificity, never quantification targets.
pub fn quant_row_indices(
    total: usize,
    first_decoy_index: Option<usize>,
    num_decoys: usize,
) -> impl Iterator<Item = usize> {
    let fdi = first_decoy_index.unwrap_or(total).min(total);
    let decoy_end = match first_decoy_index {
        Some(_) if num_decoys > 0 => (fdi + num_decoys).min(total),
        _ => total,
    };
    (0..fdi).chain(decoy_end..total)
}

/// A set of reference sequences held as one contiguous buffer plus endpoint
/// offsets, handing out `&[u8]` views.
///
/// This is the shape the data already has on disk — salmon's index stores every
/// reference concatenated into a single blob alongside an offset table — and the
/// shape [`SalmonIndex::ref_seq`](../salmon_index/struct.SalmonIndex.html) has
/// always used internally. `RefSeqs` makes it a type so callers can pass the
/// whole set around without first exploding it into `Vec<Vec<u8>>`.
///
/// Against `Vec<Vec<u8>>` for `n` references this saves `n` heap allocations and
/// `n * 24` bytes of `Vec` headers, keeps the bases contiguous rather than
/// scattered, and — where the source is already concatenated — avoids copying
/// the bases at all. On a human transcriptome (~250k references) the header
/// overhead alone is ~6 MB before a single base is stored.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RefSeqs {
    /// Every reference's bases, back to back.
    seqs: Vec<u8>,
    /// `n + 1` endpoints delimiting each reference inside `seqs`; reference `i`
    /// occupies `offsets[i]..offsets[i + 1]`.
    offsets: Vec<u64>,
}

impl RefSeqs {
    /// Wrap an already-concatenated buffer and its offset table.
    ///
    /// This is the zero-copy path: both parts are taken by value and stored as
    /// they are.
    ///
    /// # Errors
    ///
    /// Returns an error if `offsets` is empty, is not non-decreasing, or runs
    /// past the end of `seqs` — any of which would otherwise surface later as a
    /// panicking slice or, worse, as a silently wrong reference.
    pub fn from_concatenated(seqs: Vec<u8>, offsets: Vec<u64>) -> Result<Self> {
        let bad = |m: String| SalmonError::Other(m);
        if offsets.is_empty() {
            return Err(bad(
                "reference offset table is empty; expected at least one endpoint".into(),
            ));
        }
        if offsets[0] != 0 {
            return Err(bad(format!(
                "reference offset table starts at {} rather than 0",
                offsets[0]
            )));
        }
        if !offsets.windows(2).all(|w| w[0] <= w[1]) {
            return Err(bad("reference offset table is not non-decreasing".into()));
        }
        let end = *offsets.last().expect("checked non-empty") as usize;
        if end > seqs.len() {
            return Err(bad(format!(
                "reference offset table ends at {end} but the sequence buffer is only {} bytes",
                seqs.len()
            )));
        }
        Ok(Self { seqs, offsets })
    }

    /// Build from individual sequences, concatenating them.
    ///
    /// For callers that genuinely produce one sequence at a time (assembling
    /// references from an annotation, or test fixtures). Prefer
    /// [`Self::from_concatenated`] when the bytes are already contiguous.
    pub fn from_sequences<I, S>(seqs: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<[u8]>,
    {
        let iter = seqs.into_iter();
        let (lo, _) = iter.size_hint();
        let mut buf = Vec::new();
        let mut offsets = Vec::with_capacity(lo + 1);
        offsets.push(0u64);
        for s in iter {
            buf.extend_from_slice(s.as_ref());
            offsets.push(buf.len() as u64);
        }
        Self { seqs: buf, offsets }
    }

    /// Number of references.
    pub fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    /// Whether there are no references at all.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The bases of reference `i`, or `None` if out of range.
    pub fn get(&self, i: usize) -> Option<&[u8]> {
        let s = *self.offsets.get(i)? as usize;
        let e = *self.offsets.get(i + 1)? as usize;
        self.seqs.get(s..e)
    }

    /// Iterate every reference's bases in order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = &[u8]> + '_ {
        (0..self.len()).map(move |i| {
            let s = self.offsets[i] as usize;
            let e = self.offsets[i + 1] as usize;
            &self.seqs[s..e]
        })
    }

    /// All references' bases, back to back.
    pub fn concatenated(&self) -> &[u8] {
        &self.seqs
    }

    /// The endpoint table; reference `i` spans `offsets[i]..offsets[i + 1]`.
    pub fn offsets(&self) -> &[u64] {
        &self.offsets
    }
}

impl std::ops::Index<usize> for RefSeqs {
    type Output = [u8];

    fn index(&self, i: usize) -> &[u8] {
        self.get(i)
            .unwrap_or_else(|| panic!("reference index {i} out of range (have {})", self.len()))
    }
}

impl<'a> IntoIterator for &'a RefSeqs {
    type Item = &'a [u8];
    type IntoIter = Box<dyn ExactSizeIterator<Item = &'a [u8]> + 'a>;

    fn into_iter(self) -> Self::IntoIter {
        Box::new(self.iter())
    }
}
