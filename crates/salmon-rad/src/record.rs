//! The salmon bulk RAD record and its parsing context.
//!
//! This is the innermost layer of the format: the exact bytes written for one
//! fragment. The per-hit `(compressed_ori_ref, pos, frag_len)` triple is
//! byte-identical to piscem's `bulk_with_pos`; for the selective-alignment
//! profile salmon adds a per-hit `alignment_score`.
//!
//! On-disk record layout (little-endian), matching the order of the read- and
//! alignment-level [`TagSection`]s in
//! [`crate::schema::build_prelude`]:
//! ```text
//! na: u32
//! frag_map_type: u8           (read tag)
//! per alignment (na times):
//!   compressed_ori_ref: u32   (aln tag)
//!   pos: u32                  (aln tag)
//!   frag_len: u16             (aln tag)
//!   alignment_score: i32      (aln tag, SA profile only)
//! ```
//!
//! `na` ("number of alignments") comes first so a reader knows how many hit
//! blocks follow. "Little-endian" means the least significant byte of a
//! multi-byte integer is stored first; it is fixed by the format rather than
//! taken from the host, so files are portable across machines.

use std::io::{Read, Write};

use bio_types::strand::Strand;
use libradicl::rad_types::TagSection;
use libradicl::record::{MappedRecord, RecordContext};

use crate::{RadHit, RadProfile};

/// Parsing context for [`SalmonBulkRecord`]: which profile (and hence whether a
/// per-hit `alignment_score` is present).
///
/// The record layout is not self-describing per record — whether four extra
/// bytes follow each hit depends on the file's profile — so the decoder has to
/// carry that decision alongside it. That is what a "context" is here.
#[derive(Clone, Copy, Debug)]
pub struct SalmonBulkContext {
    pub profile: RadProfile,
}

impl RecordContext for SalmonBulkContext {
    /// Derive the context from the file's declared tag sections, i.e. work out
    /// how to decode records by reading the header.
    fn get_context_from_tag_section(
        ft: &TagSection,
        _rt: &TagSection,
        at: &TagSection,
    ) -> anyhow::Result<Self> {
        // salmon RAD is identified by its `rad_type` file tag (piscem uses
        // `known_rad_type`); require it so we don't mis-read a foreign RAD.
        if !ft.has_tag("rad_type") {
            anyhow::bail!("salmon bulk record context requires a \"rad_type\" file-level tag");
        }
        let profile = if at.has_tag("alignment_score") {
            RadProfile::SelectiveAlignment
        } else {
            RadProfile::Sketch
        };
        Ok(Self { profile })
    }
}

/// One fragment's worth of salmon mappings, as stored in a RAD chunk.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SalmonBulkRecord {
    /// Fragment-level mapping type (proper pair / orphan / single-end); see
    /// [`crate::frag_map_type`].
    pub frag_type: u8,
    /// Every place this fragment could have come from.
    pub hits: Vec<RadHit>,
    // cached transcript ids, to satisfy `MappedRecord::refs() -> &[u32]`.
    //
    // The trait hands out a borrowed slice, which means the ids must already
    // exist contiguously somewhere; they are interleaved inside `hits`, so a
    // flat copy is kept alongside.
    refs: Vec<u32>,
}

impl SalmonBulkRecord {
    /// Build a record from a fragment's `frag_type` and hits, deriving the
    /// cached `refs` slice.
    pub fn new(frag_type: u8, hits: Vec<RadHit>) -> Self {
        let refs = hits.iter().map(|h| h.tid).collect();
        Self {
            frag_type,
            hits,
            refs,
        }
    }
}

// `MappedRecord` is libradicl's interface for "a thing that can be read from and
// written to a RAD chunk". Implementing it lets libradicl's generic chunk
// machinery (compression, parallel chunk iteration) work with salmon's record.
impl MappedRecord for SalmonBulkRecord {
    type ParsingContext = SalmonBulkContext;
    type PeekResult = ();

    /// Peeking (inspecting a record without fully decoding it) is used by
    /// libradicl for barcode-style formats; bulk salmon needs nothing here.
    #[inline]
    fn peek_record(_buf: &[u8], _ctx: &Self::ParsingContext) {}

    /// Decode one record from a byte stream, following the layout documented at
    /// the top of this file.
    fn from_bytes_with_context<T: Read>(reader: &mut T, ctx: &Self::ParsingContext) -> Self {
        // Fixed-size scratch buffers, one per integer width, reused for every
        // field so the loop below performs no allocation.
        let mut u32b = [0u8; 4];
        let mut u16b = [0u8; 2];
        let mut u8b = [0u8; 1];

        // Header of the record: how many hits follow, then the fragment type.
        // `read_exact` fills the buffer completely or fails; `unwrap` is used
        // because the trait signature cannot return an error here, and a
        // truncated chunk is an unrecoverable format violation.
        reader.read_exact(&mut u32b).unwrap();
        let na = u32::from_le_bytes(u32b) as usize;
        reader.read_exact(&mut u8b).unwrap();
        let frag_type = u8b[0];

        // Pre-size the vector: we already know exactly how many hits there are.
        let mut hits = Vec::with_capacity(na);
        for _ in 0..na {
            reader.read_exact(&mut u32b).unwrap();
            // One u32 carries three fields; unpack it (see `RadHit`).
            let (tid, is_fw, mate_fw) = RadHit::decode_ori_ref(u32::from_le_bytes(u32b));
            reader.read_exact(&mut u32b).unwrap();
            let pos = u32::from_le_bytes(u32b);
            reader.read_exact(&mut u16b).unwrap();
            let frag_len = u16::from_le_bytes(u16b);
            // The score field exists only in the SA profile; in the sketch
            // profile there are no bytes to read and the score is left at 0.
            let score = if ctx.profile.has_scores() {
                reader.read_exact(&mut u32b).unwrap();
                i32::from_le_bytes(u32b)
            } else {
                0
            };
            hits.push(RadHit {
                tid,
                is_fw,
                mate_fw,
                pos,
                frag_len,
                score,
            });
        }
        Self::new(frag_type, hits)
    }

    /// Encode one record, the exact mirror of `from_bytes_with_context`. The two
    /// must stay in lockstep, and both must match the tag declaration order in
    /// [`crate::schema::build_prelude`].
    fn write<W: Write>(&self, writer: &mut W, ctx: &Self::ParsingContext) -> anyhow::Result<()> {
        writer.write_all(&(self.hits.len() as u32).to_le_bytes())?;
        writer.write_all(&self.frag_type.to_le_bytes())?;
        for h in &self.hits {
            writer.write_all(&h.compressed_ori_ref().to_le_bytes())?;
            writer.write_all(&h.pos.to_le_bytes())?;
            writer.write_all(&h.frag_len.to_le_bytes())?;
            if ctx.profile.has_scores() {
                writer.write_all(&h.score.to_le_bytes())?;
            }
        }
        Ok(())
    }

    /// A record with no hits is an unmapped fragment.
    fn is_empty(&self) -> bool {
        self.hits.is_empty()
    }

    /// Number of placements, used by readers to size their buffers.
    fn num_aln(&self) -> usize {
        self.hits.len()
    }

    /// The transcript ids this fragment touched (the cached flat copy).
    fn refs(&self) -> &[u32] {
        &self.refs
    }

    /// Whether any placement is on the given strand; libradicl uses this for
    /// strand-aware filtering. `Unknown` means "either strand will do", so it is
    /// satisfied by any hit at all.
    fn has_alignment_on_strand(&self, s: Strand) -> bool {
        match s {
            Strand::Unknown => !self.hits.is_empty(),
            Strand::Forward => self.hits.iter().any(|h| h.is_fw),
            Strand::Reverse => self.hits.iter().any(|h| !h.is_fw),
        }
    }
}
