//! Access to reference sequences for alignment validation.
//!
//! The mapper needs the forward-strand bytes of a reference to align a read
//! against. piscem's index does not retain reference sequences, so they are
//! supplied through this trait (salmon-index persists and serves them). Keeping
//! it a trait decouples `salmon-map` from the concrete index/store.
//!
//! A "trait" in Rust is an interface: it lists methods a type must provide. Code
//! written against the trait works with any type that implements it, so the
//! mapper can be tested against a tiny in-memory stub and run against the real
//! on-disk index without changing a line.

/// Provides reference sequences (and decoy status) to the mapper.
pub trait RefProvider {
    /// Number of references.
    fn num_refs(&self) -> usize;

    /// Forward-strand sequence of reference `tid` (ASCII `ACGT`).
    ///
    /// Returns a borrowed slice rather than a copy: transcriptomes are hundreds
    /// of megabytes and this is called for every alignment.
    fn ref_seq(&self, tid: u32) -> &[u8];

    /// Whether reference `tid` is a decoy (genomic) sequence.
    ///
    /// Decoys are extra sequence (typically the genome) added to the index as a
    /// sink: a read that really came from an intron or an unannotated region
    /// will align better to the decoy than to any transcript, and is then
    /// discarded instead of being falsely credited to a transcript.
    ///
    /// This method has a body, so it is a *default*: an implementer that has no
    /// decoys can leave it out entirely.
    fn is_decoy(&self, _tid: u32) -> bool {
        false
    }
}
