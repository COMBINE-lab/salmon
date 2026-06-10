//! Access to reference sequences for alignment validation.
//!
//! The mapper needs the forward-strand bytes of a reference to align a read
//! against. piscem's index does not retain reference sequences, so they are
//! supplied through this trait (salmon-index persists and serves them). Keeping
//! it a trait decouples `salmon-map` from the concrete index/store.

/// Provides reference sequences (and decoy status) to the mapper.
pub trait RefProvider {
    /// Number of references.
    fn num_refs(&self) -> usize;

    /// Forward-strand sequence of reference `tid` (ASCII `ACGT`).
    fn ref_seq(&self, tid: u32) -> &[u8];

    /// Whether reference `tid` is a decoy (genomic) sequence. Reads that align
    /// best to a decoy are discarded. Defaults to no decoys.
    fn is_decoy(&self, _tid: u32) -> bool {
        false
    }
}
