//! The salmon bulk RAD record and its parsing context.
//!
//! The per-hit `(compressed_ori_ref, pos, frag_len)` triple is byte-identical to
//! piscem's `bulk_with_pos`; salmon adds a read-level `frag_name_hash` (after
//! `frag_map_type`) and, for the selective-alignment profile, a per-hit
//! `alignment_score`.
//!
//! On-disk record layout (little-endian), matching the order of the read- and
//! alignment-level [`TagSection`](libradicl::rad_types::TagSection)s in
//! [`crate::schema::build_prelude`]:
//! ```text
//! na: u32
//! frag_map_type: u8           (read tag)
//! frag_name_hash: u128        (read tag)
//! per alignment (na times):
//!   compressed_ori_ref: u32   (aln tag)
//!   pos: u32                  (aln tag)
//!   frag_len: u16             (aln tag)
//!   alignment_score: i32      (aln tag, SA profile only)
//! ```

use std::io::{Read, Write};

use bio_types::strand::Strand;
use libradicl::rad_types::TagSection;
use libradicl::record::{MappedRecord, RecordContext};

use crate::{RadHit, RadProfile};

/// Parsing context for [`SalmonBulkRecord`]: which profile (and hence whether a
/// per-hit `alignment_score` is present).
#[derive(Clone, Copy, Debug)]
pub struct SalmonBulkContext {
    pub profile: RadProfile,
}

impl RecordContext for SalmonBulkContext {
    fn get_context_from_tag_section(
        _ft: &TagSection,
        rt: &TagSection,
        at: &TagSection,
    ) -> anyhow::Result<Self> {
        // salmon RAD must carry the read-name hash; require it so we don't
        // silently mis-read a foreign RAD as salmon's.
        if !rt.has_tag("frag_name_hash") {
            anyhow::bail!(
                "salmon bulk record context requires a \"frag_name_hash\" read-level tag"
            );
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
    pub frag_type: u8,
    pub name_hash: u128,
    pub hits: Vec<RadHit>,
    // cached transcript ids, to satisfy `MappedRecord::refs() -> &[u32]`.
    refs: Vec<u32>,
}

impl SalmonBulkRecord {
    /// Build a record from a fragment's `frag_type`, read-name hash, and hits.
    pub fn new(frag_type: u8, name_hash: u128, hits: Vec<RadHit>) -> Self {
        let refs = hits.iter().map(|h| h.tid).collect();
        Self {
            frag_type,
            name_hash,
            hits,
            refs,
        }
    }
}

impl MappedRecord for SalmonBulkRecord {
    type ParsingContext = SalmonBulkContext;
    type PeekResult = ();

    #[inline]
    fn peek_record(_buf: &[u8], _ctx: &Self::ParsingContext) {}

    fn from_bytes_with_context<T: Read>(reader: &mut T, ctx: &Self::ParsingContext) -> Self {
        let mut u32b = [0u8; 4];
        let mut u16b = [0u8; 2];
        let mut u8b = [0u8; 1];
        let mut u128b = [0u8; 16];

        reader.read_exact(&mut u32b).unwrap();
        let na = u32::from_le_bytes(u32b) as usize;
        reader.read_exact(&mut u8b).unwrap();
        let frag_type = u8b[0];
        reader.read_exact(&mut u128b).unwrap();
        let name_hash = u128::from_le_bytes(u128b);

        let mut hits = Vec::with_capacity(na);
        for _ in 0..na {
            reader.read_exact(&mut u32b).unwrap();
            let (tid, is_fw, mate_fw) = RadHit::decode_ori_ref(u32::from_le_bytes(u32b));
            reader.read_exact(&mut u32b).unwrap();
            let pos = u32::from_le_bytes(u32b);
            reader.read_exact(&mut u16b).unwrap();
            let frag_len = u16::from_le_bytes(u16b);
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
        Self::new(frag_type, name_hash, hits)
    }

    fn write<W: Write>(&self, writer: &mut W, ctx: &Self::ParsingContext) -> anyhow::Result<()> {
        writer.write_all(&(self.hits.len() as u32).to_le_bytes())?;
        writer.write_all(&self.frag_type.to_le_bytes())?;
        writer.write_all(&self.name_hash.to_le_bytes())?;
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

    fn is_empty(&self) -> bool {
        self.hits.is_empty()
    }

    fn num_aln(&self) -> usize {
        self.hits.len()
    }

    fn refs(&self) -> &[u32] {
        &self.refs
    }

    fn has_alignment_on_strand(&self, s: Strand) -> bool {
        match s {
            Strand::Unknown => !self.hits.is_empty(),
            Strand::Forward => self.hits.iter().any(|h| h.is_fw),
            Strand::Reverse => self.hits.iter().any(|h| !h.is_fw),
        }
    }
}
