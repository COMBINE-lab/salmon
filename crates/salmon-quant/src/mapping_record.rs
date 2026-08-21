//! Format-neutral borrowed alignment records shared by SAM and BAM output.
//!
//! # Why this module exists
//!
//! salmon can write its per-fragment mappings as text (SAM) or as the equivalent
//! compressed binary (BAM). The two formats encode exactly the same information,
//! so if each writer derived flags, CIGARs, positions and tags for itself the two
//! would inevitably drift apart. Instead this module derives every field once, in
//! [`emit_fragment_records`], and hands each format a borrowed record to
//! serialize. The header is shared the same way, via [`header_text`].
//!
//! "Borrowed" means the record points at the caller's read name, sequence and
//! qualities rather than copying them, so emitting a mapping allocates nothing
//! beyond the reusable scratch in [`EmitScratch`].
//!
//! # Background: SAM records
//!
//! A SAM record describes one read's placement: which reference, at what
//! position, in which orientation, and how the read's bases line up with the
//! reference. The alignment shape is a CIGAR string.
//!
//! # Two grades of CIGAR
//!
//! Mapping scores candidate placements without keeping a base-level alignment
//! (see [`salmon_map::realize`] for why), so there are two ways to describe a
//! placement:
//!
//! * **Realized**, the placement is re-aligned with traceback, giving real `I`
//!   and `D` operations, an `NM` edit distance and an `MD` string. This is what
//!   gets written whenever the reference sequence is available, and it is what
//!   makes the output usable by tools that expect a genuine alignment.
//! * **Spoofed**, `<readLen>M` with soft clips where the read overhangs a
//!   transcript end. The fallback for when realization cannot run (an index
//!   without reference sequence) or cannot form an alignment.

use std::io;

use salmon_core::MateStatus;
use salmon_core::RefProvider as _;
use salmon_index::SalmonIndex;
use salmon_map::realize::{
    push_op, realize, CigarOp, CigarOpKind, RealizeScratch, RealizedAlignment,
};
use salmon_map::{AlignConfig, ScoredMapping};

// SAM FLAG bits (a bitfield in field 2 of every record). Each names one yes/no
// property of the alignment; a record's flags are these OR-ed together.
/// The read is part of a pair.
pub const PAIRED: u16 = 0x1;
/// Both mates placed consistently (right orientation, plausible distance).
pub const PROPER_PAIR: u16 = 0x2;
/// This read did not map anywhere.
pub const UNMAPPED: u16 = 0x4;
/// This read's mate did not map (an orphan).
pub const MATE_UNMAPPED: u16 = 0x8;
/// This read aligns to the reverse strand, so its stored sequence is the
/// reverse complement of what the sequencer produced.
pub const IS_RC: u16 = 0x10;
/// The mate aligns to the reverse strand.
pub const MATE_RC: u16 = 0x20;
/// This is read 1 of the pair.
pub const READ1: u16 = 0x40;
/// This is read 2 of the pair.
pub const READ2: u16 = 0x80;
/// An additional, non-primary placement of the same fragment. Exactly one
/// placement per fragment lacks this bit; downstream tools count primaries to
/// avoid counting a multi-mapping fragment several times.
pub const SECONDARY: u16 = 0x100;

/// One read's placement, with every field already in the form both writers need.
/// The `'a` lifetime ties the borrowed name, sequence, qualities and CIGAR to the
/// caller's buffers.
pub struct AlignmentRecord<'a> {
    pub name: &'a [u8],
    pub flags: u16,
    /// `None` for an unmapped record, which SAM writes as `*` and BAM as `-1`.
    pub reference_id: Option<usize>,
    pub position: i32,
    pub mapping_quality: u8,
    pub cigar: &'a [CigarOp],
    pub mate_reference_id: Option<usize>,
    pub mate_position: Option<i32>,
    pub template_length: i32,
    pub sequence: &'a [u8],
    /// Phred+33 base qualities in *read* order, as they came out of the FASTQ.
    /// `None` when the input carried none (a FASTA), which both formats write as
    /// "unavailable" rather than inventing values.
    pub qualities: Option<&'a [u8]>,
    /// Whether the stored sequence and qualities have to be reversed, because the
    /// read aligns to the reverse strand and SAM stores everything on the
    /// reference's forward strand.
    pub reverse_complement: bool,
    pub nh: usize,
    pub hi: usize,
    pub xt: u8,
    pub alignment_score: i32,
    /// The mapping's equivalence-class weight, written as `ZW`.
    pub weight: f32,
    /// `NM`, present only for a realized alignment.
    pub edit_distance: Option<u32>,
    /// `MD`, present only for a realized alignment.
    pub md: Option<&'a str>,
    /// `RG`, present when the run declared a read group.
    pub read_group: Option<&'a str>,
    /// `MC`, the mate's CIGAR. Lets a reader work out where the mate ends without
    /// having to find the mate's own record first, which is what every
    /// mate-aware tool otherwise has to sort the file to do.
    pub mate_cigar: Option<&'a [CigarOp]>,
    /// `MQ`, the mate's mapping quality, for the same reason.
    pub mate_mapping_quality: Option<u8>,
}

/// SAM version reported in `@HD VN:`. 1.6 is the version that defines the
/// `GO:query` grouping this output actually has.
const SAM_VERSION: &str = "1.6";

/// A read group: the `@RG` header line and the `ID` that every record's `RG` tag
/// refers back to.
///
/// Read groups are how a merged BAM keeps track of which sample, library and
/// sequencing run each record came from. salmon does not need one to quantify,
/// but anything that merges or genotypes salmon's output does.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadGroup {
    /// The `ID` field, repeated on every record as `RG:Z:<id>`.
    pub id: String,
    /// The full `@RG` line as it appears in the header, without the newline.
    pub header_line: String,
}

impl ReadGroup {
    /// Parse a read-group specification.
    ///
    /// The accepted spelling is the one bwa and STAR use: tab-separated
    /// `KEY:value` fields, with or without a leading `@RG`, and with a literal
    /// backslash-`t` accepted in place of a real tab so the whole thing survives a
    /// shell without quoting gymnastics. `ID` is mandatory, it is the key every
    /// record's `RG` tag points at, so a group without one would produce a header
    /// that no record could reference.
    pub fn parse(spec: &str) -> anyhow::Result<Self> {
        let normalized = spec.replace("\\t", "\t");
        let trimmed = normalized.trim();
        let body = trimmed
            .strip_prefix("@RG\t")
            .or_else(|| trimmed.strip_prefix("@RG "))
            .unwrap_or(trimmed)
            .trim();
        if body.is_empty() || body == "@RG" {
            anyhow::bail!("read group is empty; expected something like ID:sample1\\tSM:sample1");
        }
        let fields: Vec<&str> = body.split(['\t', ' ']).filter(|f| !f.is_empty()).collect();
        for field in &fields {
            let Some((tag, _)) = field.split_once(':') else {
                anyhow::bail!("read-group field {field:?} is not of the form KEY:value");
            };
            if tag.len() != 2 || !tag.bytes().all(|b| b.is_ascii_alphanumeric()) {
                anyhow::bail!("read-group field {field:?} does not start with a two-letter tag");
            }
        }
        let id = fields
            .iter()
            .find_map(|f| f.strip_prefix("ID:"))
            .ok_or_else(|| anyhow::anyhow!("read group has no ID: field"))?
            .to_string();
        if id.is_empty() {
            anyhow::bail!("read-group ID: is empty");
        }
        Ok(Self {
            id,
            header_line: format!("@RG\t{}", fields.join("\t")),
        })
    }
}

/// Everything about mapping output that is fixed for a whole run, rather than
/// derived per record.
#[derive(Default)]
pub struct RecordOptions<'a> {
    /// The read group to declare in the header and tag every record with.
    pub read_group: Option<&'a ReadGroup>,
    /// When set, placements are re-aligned to produce a real CIGAR plus `NM` and
    /// `MD`. `None` keeps the spoofed `<readLen>M` CIGAR and omits both tags.
    pub align_config: Option<&'a AlignConfig>,
    /// `@SQ UR:`, where the reference sequences came from, so a reader can find
    /// them again.
    pub reference_uri: Option<&'a str>,
}

impl RecordOptions<'_> {
    /// The reference name for `reference_id`, from whichever table this run's
    /// records index into.
    pub fn reference_name<'s>(&'s self, salmon: &'s SalmonIndex, id: usize) -> &'s str {
        salmon.ref_name(id)
    }

    /// The `(name, length)` reference table to declare in the header.
    pub fn reference_table(&self, salmon: &SalmonIndex) -> Vec<(String, u64)> {
        (0..salmon.num_refs())
            .map(|tid| (salmon.ref_name(tid).to_string(), salmon.ref_len(tid)))
            .collect()
    }
}

/// The `M5` checksum of one reference sequence.
///
/// `M5` is what ties a BAM to the exact reference bases it was aligned against:
/// two references with the same name and length but different sequence produce
/// different checksums, so a reader can tell them apart instead of silently
/// mixing them. It is also what CRAM uses to find the reference it needs.
///
/// The definition is fixed by the SAM spec: uppercase the bases, drop everything
/// that is not a printable character, then MD5 the result as lower-case hex.
fn sequence_checksum(sequence: &[u8], ambiguous: &[u32]) -> String {
    use md5::Digest as _;
    let mut hasher = md5::Md5::new();
    let mut block = [0u8; 8192];
    let mut filled = 0;
    for (offset, &base) in sequence.iter().enumerate() {
        if base.is_ascii_whitespace() || !base.is_ascii_graphic() {
            continue;
        }
        // The checksum has to identify the reference the user supplied, which is
        // the one `NM` and `MD` now describe. Where indexing replaced a base, hash
        // the `N` that was there rather than the substitute, or the header would
        // name a sequence nothing else in the file refers to.
        block[filled] = if ambiguous.binary_search(&(offset as u32)).is_ok() {
            b'N'
        } else {
            base.to_ascii_uppercase()
        };
        filled += 1;
        if filled == block.len() {
            hasher.update(block);
            filled = 0;
        }
    }
    hasher.update(&block[..filled]);
    hasher
        .finalize()
        .iter()
        .fold(String::with_capacity(32), |mut hex, byte| {
            use std::fmt::Write as _;
            let _ = write!(hex, "{byte:02x}");
            hex
        })
}

/// Reusable per-worker buffers for record emission.
///
/// Realizing an alignment needs somewhere to put the CIGAR, the `MD` string and
/// the aligner's own workspace. Keeping them here means a worker allocates them
/// once rather than once per record.
#[derive(Default)]
pub struct EmitScratch {
    realize_scratch: RealizeScratch,
    /// Realized alignments for the two mates of the fragment being emitted.
    realized: [RealizedAlignment; 2],
    /// Spoofed CIGARs, used when realization is off or fails.
    spoofed: [Vec<CigarOp>; 2],
}

/// The header text shared by both mapping-output formats: `@HD`, one `@SQ` per
/// reference (including decoys, matching salmon's full ref table), the optional
/// `@RG`, and `@PG`.
///
/// This is the single source of truth for the header, mirroring
/// [`emit_fragment_records`] for record bodies. SAM writes the text verbatim;
/// BAM stores it as the `text` block and then repeats the reference table in
/// binary form (see `bam::encode_header`).
///
/// `SO:unsorted GO:query` is the honest description of what salmon produces: no
/// coordinate order at all, but all records for one fragment together. Declaring
/// the grouping lets a reader pair mates in a single streaming pass instead of
/// assuming it has to sort first.
pub fn header_text(salmon: &SalmonIndex, command: &str, options: &RecordOptions<'_>) -> String {
    use std::fmt::Write as _;
    // ~48 B/@SQ line for typical transcript names, plus the @HD, @RG and @PG lines.
    let mut text = String::with_capacity(salmon.num_refs() * 48 + command.len() + 128);
    let _ = writeln!(text, "@HD\tVN:{SAM_VERSION}\tSO:unsorted\tGO:query");
    let references = options.reference_table(salmon);
    // Hashing a whole transcriptome is a few hundred megabytes of MD5, which is
    // enough to be noticeable at the start of a short run, and it is one
    // independent hash per reference. Run them in parallel so the header costs
    // wall-clock proportional to the largest transcript rather than to the sum.
    let checksums: Vec<String> = {
        use rayon::prelude::*;
        (0..references.len())
            .into_par_iter()
            .map(|index| {
                let sequence = salmon.ref_seq(index as u32);
                if sequence.is_empty() {
                    String::new()
                } else {
                    sequence_checksum(sequence, salmon.ambiguous_offsets(index))
                }
            })
            .collect()
    };
    for (index, (name, length)) in references.into_iter().enumerate() {
        let _ = write!(text, "@SQ\tSN:{name}\tLN:{length}");
        if let Some(checksum) = checksums.get(index).filter(|c| !c.is_empty()) {
            let _ = write!(text, "\tM5:{checksum}");
        }
        if let Some(uri) = options.reference_uri {
            let _ = write!(text, "\tUR:{uri}");
        }
        text.push('\n');
    }
    if let Some(read_group) = options.read_group {
        let _ = writeln!(text, "{}", read_group.header_line);
    }
    let _ = writeln!(
        text,
        "@PG\tID:salmon\tPN:salmon\tVN:{}\tCL:{}",
        crate::output::SALMON_VERSION,
        command
    );
    // State the one ordering guarantee in the file itself. `GO:query` says the
    // records are grouped by name, but not that nothing else is ordered, and a
    // reader that assumes more will be wrong in ways that are hard to trace.
    let _ = writeln!(
        text,
        "@CO\tRecords for one fragment are contiguous. No other order is imposed: \
         fragments appear in whichever order the mapping threads finished them."
    );
    text
}

/// The read name to report: the FASTQ header up to the first space or tab, with
/// a trailing `/1` or `/2` mate suffix removed.
///
/// Both mates of a pair must carry the *same* name for a SAM/BAM reader to pair
/// them up, but many FASTQ files distinguish them with that suffix, so it has to
/// go.
pub fn read_name(id: &[u8]) -> &[u8] {
    let end = id
        .iter()
        .position(|&b| b == b' ' || b == b'\t')
        .unwrap_or(id.len());
    let name = &id[..end];
    if name.len() > 2 && name[name.len() - 2] == b'/' && matches!(name[name.len() - 1], b'1' | b'2')
    {
        &name[..name.len() - 2]
    } else {
        name
    }
}

/// Complement one base (unknown bases become `N`), used when emitting a
/// reverse-strand read's sequence.
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

/// The qualities to write for this record, or `None` when there are none to
/// write.
///
/// Base qualities are not optional output: whenever the input carries them they
/// go into the record, which is why there is no flag to turn them on. There are
/// only two ways to end up without them, and they are not the same thing:
///
/// * The input never had any (a FASTA). Nothing is wrong and nothing is said.
/// * A quality string disagrees in length with its sequence. The FASTQ reader
///   already rejects that outright, so this is a guard against some other input
///   path rather than a case users meet: writing such a record would produce one
///   no reader accepts. It warns rather than dropping them quietly, because a
///   guard that fires silently is one nobody finds out about.
pub fn usable_qualities<'a>(record: &AlignmentRecord<'a>) -> Option<&'a [u8]> {
    let qualities = record.qualities?;
    if qualities.len() == record.sequence.len() {
        return (!qualities.is_empty()).then_some(qualities);
    }
    static WARNED: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
    if !WARNED.swap(true, std::sync::atomic::Ordering::Relaxed) {
        tracing::warn!(
            read = %String::from_utf8_lossy(record.name),
            qualities = qualities.len(),
            bases = record.sequence.len(),
            "read has a quality string of a different length than its sequence; \
             writing its record without qualities (reported once)"
        );
    }
    None
}

/// The `MAPQ` to report for a fragment with `nh` placements.
///
/// SAM's `MAPQ` is a phred-scaled probability that the placement is wrong, and
/// salmon does not compute one: its whole design is to keep every plausible
/// placement and let the EM apportion them. What it does know is *how many*
/// placements a fragment has, which is the same thing STAR reports, so this uses
/// STAR's mapping: a uniquely placed fragment gets 255 and the value falls off
/// as the fragment becomes more ambiguous. Matching STAR matters because the
/// RNA-seq tools downstream of a salmon BAM are the ones already tuned to it: to
/// all of them `-q 255` means "uniquely placed".
///
/// The previous hardcoded `1` made every record look ambiguous, so any such
/// filter discarded the whole file. It is the one placeholder in this output
/// that costs nothing to replace: the count is already in hand when the records
/// are written, so no alignment work is implied (COMBINE-lab/salmon#1140,
/// #1141).
pub fn mapping_quality(nh: usize) -> u8 {
    match nh {
        // No placements is not a case a mapped record can reach (the writer
        // skips a fragment with none), so this arm exists for the record kinds
        // that have no placement by construction: an unmapped read, once such
        // records are written. 255 there would claim certainty about a
        // placement that does not exist; 0 is what samtools convention expects
        // on an unmapped record.
        0 => 0,
        1 => 255,
        2 => 3,
        3..=4 => 1,
        _ => 0,
    }
}

/// The CIGAR and clamped position for a read placed at `pos` on a transcript of
/// length `txp_len`, without a base-level alignment.
///
/// A read can hang off either end of a transcript, the mapper places it by seed
/// position, and the fragment may genuinely extend past the annotated boundary.
/// SAM cannot express a negative position or an alignment past the reference end,
/// so the overhanging bases are soft-clipped and the position clamped into range.
fn overhang_cigar(ops: &mut Vec<CigarOp>, pos: i32, read_len: i32, txp_len: i32) -> i32 {
    ops.clear();
    let read_len_u = read_len.max(0) as u32;
    if pos + read_len < 0 {
        push_op(ops, CigarOpKind::SoftClip, read_len_u);
        0
    } else if pos < 0 {
        let matched = (read_len + pos).max(0) as u32;
        push_op(ops, CigarOpKind::SoftClip, read_len_u - matched);
        push_op(ops, CigarOpKind::Match, matched);
        0
    } else if pos > txp_len {
        push_op(ops, CigarOpKind::SoftClip, read_len_u);
        pos
    } else if pos + read_len > txp_len {
        let matched = (txp_len - pos).max(0) as u32;
        push_op(ops, CigarOpKind::Match, matched);
        push_op(ops, CigarOpKind::SoftClip, read_len_u - matched);
        pos
    } else {
        push_op(ops, CigarOpKind::Match, read_len_u);
        pos
    }
}

/// One mate's placement after CIGAR derivation: everything that differs between
/// the realized and spoofed paths, resolved.
struct Placed {
    position: i32,
    edit_distance: Option<u32>,
    /// The realized alignment's own score, when there is one. `AS` reports this
    /// in preference to the mapper's, so that the score, the CIGAR and `NM`
    /// describe the same alignment.
    score: Option<i32>,
    /// Which of [`EmitScratch`]'s two mate slots holds this placement's CIGAR.
    slot: usize,
    /// Whether that CIGAR is in `realized` (`true`) or `spoofed` (`false`).
    realized: bool,
}

/// Derive one mate's CIGAR and start position, realizing the alignment when the
/// run asked for it and the reference sequence is there to align against.
///
/// Falls back to the spoofed CIGAR silently: a missing reference sequence is a
/// property of the index, not an error, and a record with a spoofed CIGAR is
/// still a valid record.
#[allow(clippy::too_many_arguments)]
fn place(
    scratch: &mut EmitScratch,
    slot: usize,
    options: &RecordOptions<'_>,
    salmon: &SalmonIndex,
    tid: usize,
    sequence: &[u8],
    forward: bool,
    pos: i32,
    txp_len: i32,
) -> Placed {
    if let Some(cfg) = options.align_config {
        let ref_seq = salmon.ref_seq(tid as u32);
        if !ref_seq.is_empty() && !sequence.is_empty() {
            let EmitScratch {
                realize_scratch,
                realized,
                ..
            } = scratch;
            if realize(
                sequence,
                forward,
                ref_seq,
                salmon.ambiguous_offsets(tid),
                pos,
                cfg,
                realize_scratch,
                &mut realized[slot],
            ) {
                return Placed {
                    position: realized[slot].ref_start,
                    edit_distance: Some(realized[slot].edit_distance),
                    score: Some(realized[slot].score),
                    slot,
                    realized: true,
                };
            }
        }
    }
    let position = overhang_cigar(
        &mut scratch.spoofed[slot],
        pos,
        sequence.len() as i32,
        txp_len,
    );
    Placed {
        position,
        edit_distance: None,
        score: None,
        slot,
        realized: false,
    }
}

/// The CIGAR and `MD` a [`Placed`] refers to, borrowed out of the scratch.
fn borrow_cigar<'s>(scratch: &'s EmitScratch, placed: &Placed) -> (&'s [CigarOp], Option<&'s str>) {
    if placed.realized {
        let realized = &scratch.realized[placed.slot];
        (&realized.ops, Some(realized.md.as_str()))
    } else {
        (&scratch.spoofed[placed.slot], None)
    }
}

/// Emit borrowed records immediately; no record, name, sequence or tag allocation
/// is retained between callback invocations.
///
/// This is the single source of truth for record bodies: one call per fragment
/// produces every record that fragment implies (two per placement for a proper
/// pair, one otherwise), and the caller's closure serializes them as SAM text or
/// BAM binary.
#[allow(clippy::too_many_arguments)]
pub fn emit_fragment_records(
    salmon: &SalmonIndex,
    r1_id: &[u8],
    r1_seq: &[u8],
    r1_qual: Option<&[u8]>,
    r2: Option<(&[u8], &[u8], Option<&[u8]>)>,
    maps: &[ScoredMapping],
    options: &RecordOptions<'_>,
    scratch: &mut EmitScratch,
    emit_record: impl FnMut(&AlignmentRecord<'_>) -> io::Result<()>,
) -> io::Result<()> {
    let name1 = read_name(r1_id);
    let (name2, r2_seq, r2_qual) = r2.map_or((name1, &[][..], None), |(id, seq, qual)| {
        (read_name(id), seq, qual)
    });
    let nh = maps.len();
    let mapq = mapping_quality(nh);
    let read_group = options.read_group.map(|rg| rg.id.as_str());
    let mut emit = emit_record;

    for (index, mapping) in maps.iter().enumerate() {
        // The first placement is primary; the rest are marked secondary so a
        // downstream counter does not count this fragment once per placement.
        let secondary = if index == 0 { 0 } else { SECONDARY };
        // HI is 1-based, unlike the loop index.
        let hi = index + 1;
        let tid = mapping.tid as usize;
        let txp_len = salmon.ref_len(tid) as i32;
        // salmon's XT tag. In practice it is always `T`: a placement on a decoy
        // means the read probably came from no transcript at all, so scoring
        // drops it long before a record is written, and a decoy-aware run was
        // confirmed to emit no `D` at all. The branch stays because the tag is
        // defined over the reference table, not over what currently survives
        // filtering, and a future policy that kept decoy placements would need
        // it to be right.
        let xt = if salmon.is_decoy(mapping.tid) {
            b'D'
        } else {
            b'T'
        };
        let weight = mapping.weight as f32;
        match mapping.status {
            MateStatus::PairedEndPaired => {
                let p1 = place(
                    scratch,
                    0,
                    options,
                    salmon,
                    tid,
                    r1_seq,
                    mapping.is_fw,
                    mapping.r1_pos,
                    txp_len,
                );
                let p2 = place(
                    scratch,
                    1,
                    options,
                    salmon,
                    tid,
                    r2_seq,
                    mapping.r2_fw,
                    mapping.r2_pos,
                    txp_len,
                );
                // TLEN is the observed template (fragment) length, truncated at
                // the transcript end for the same reason positions are clamped.
                let mut fragment_length = mapping.fragment_len;
                let min_pos = mapping.r1_pos.min(mapping.r2_pos);
                if min_pos + fragment_length > txp_len {
                    fragment_length = txp_len - min_pos;
                }
                // SAM signs TLEN: positive on the leftmost mate, negative on the
                // rightmost, so the two records sum to zero.
                let tlen1 = if mapping.r1_pos <= mapping.r2_pos {
                    fragment_length
                } else {
                    -fragment_length
                };
                let mut f1 = PAIRED | PROPER_PAIR | READ1 | secondary;
                let mut f2 = PAIRED | PROPER_PAIR | READ2 | secondary;
                // Each mate's own strand sets its IS_RC bit and the *other*
                // record's MATE_RC bit, so the pair stays mutually consistent.
                if !mapping.is_fw {
                    f1 |= IS_RC;
                    f2 |= MATE_RC;
                }
                if !mapping.r2_fw {
                    f2 |= IS_RC;
                    f1 |= MATE_RC;
                }
                let (cigar1, md1) = borrow_cigar(scratch, &p1);
                let (cigar2, md2) = borrow_cigar(scratch, &p2);
                emit(&AlignmentRecord {
                    name: name1,
                    flags: f1,
                    reference_id: Some(tid),
                    position: p1.position,
                    mapping_quality: mapq,
                    cigar: cigar1,
                    mate_reference_id: Some(tid),
                    mate_position: Some(p2.position),
                    template_length: tlen1,
                    sequence: r1_seq,
                    qualities: r1_qual,
                    reverse_complement: !mapping.is_fw,
                    nh,
                    hi,
                    xt,
                    alignment_score: p1.score.unwrap_or(mapping.r1_score),
                    weight,
                    edit_distance: p1.edit_distance,
                    md: md1,
                    read_group,
                    mate_cigar: Some(cigar2),
                    mate_mapping_quality: Some(mapq),
                })?;
                emit(&AlignmentRecord {
                    name: name2,
                    flags: f2,
                    reference_id: Some(tid),
                    position: p2.position,
                    mapping_quality: mapq,
                    cigar: cigar2,
                    mate_reference_id: Some(tid),
                    mate_position: Some(p1.position),
                    template_length: -tlen1,
                    sequence: r2_seq,
                    qualities: r2_qual,
                    reverse_complement: !mapping.r2_fw,
                    nh,
                    hi,
                    xt,
                    alignment_score: p2.score.unwrap_or(mapping.score - mapping.r1_score),
                    weight,
                    edit_distance: p2.edit_distance,
                    md: md2,
                    read_group,
                    mate_cigar: Some(cigar1),
                    mate_mapping_quality: Some(mapq),
                })?;
            }
            MateStatus::SingleEnd => {
                let p = place(
                    scratch,
                    0,
                    options,
                    salmon,
                    tid,
                    r1_seq,
                    mapping.is_fw,
                    mapping.r1_pos,
                    txp_len,
                );
                let (cigar, md) = borrow_cigar(scratch, &p);
                emit(&AlignmentRecord {
                    name: name1,
                    flags: secondary | if mapping.is_fw { 0 } else { IS_RC },
                    reference_id: Some(tid),
                    position: p.position,
                    mapping_quality: mapq,
                    cigar,
                    mate_reference_id: None,
                    mate_position: None,
                    template_length: 0,
                    sequence: r1_seq,
                    qualities: r1_qual,
                    reverse_complement: !mapping.is_fw,
                    nh,
                    hi,
                    xt,
                    alignment_score: p.score.unwrap_or(mapping.r1_score),
                    weight,
                    edit_distance: p.edit_distance,
                    md,
                    read_group,
                    mate_cigar: None,
                    mate_mapping_quality: None,
                })?;
            }
            // Orphan: one mate placed, the other not. One record is emitted, for
            // whichever mate mapped, flagged so a reader knows the pair is
            // incomplete rather than that this read is single-end.
            MateStatus::PairedEndLeft | MateStatus::PairedEndRight => {
                let left = mapping.status == MateStatus::PairedEndLeft;
                let (name, sequence, qualities, position, forward, score, read_flag) = if left {
                    (
                        name1,
                        r1_seq,
                        r1_qual,
                        mapping.r1_pos,
                        mapping.is_fw,
                        mapping.r1_score,
                        READ1,
                    )
                } else {
                    (
                        name2,
                        r2_seq,
                        r2_qual,
                        mapping.r2_pos,
                        mapping.r2_fw,
                        mapping.score - mapping.r1_score,
                        READ2,
                    )
                };
                let p = place(
                    scratch, 0, options, salmon, tid, sequence, forward, position, txp_len,
                );
                let (cigar, md) = borrow_cigar(scratch, &p);
                emit(&AlignmentRecord {
                    name,
                    flags: PAIRED
                        | MATE_UNMAPPED
                        | read_flag
                        | secondary
                        | if forward { 0 } else { IS_RC },
                    reference_id: Some(tid),
                    position: p.position,
                    mapping_quality: mapq,
                    cigar,
                    // The mate did not map and no record is written for it, so
                    // naming a mate position would be a dangling reference.
                    mate_reference_id: None,
                    mate_position: None,
                    template_length: 0,
                    sequence,
                    qualities,
                    reverse_complement: !forward,
                    nh,
                    hi,
                    xt,
                    alignment_score: p.score.unwrap_or(score),
                    weight,
                    edit_distance: p.edit_distance,
                    md,
                    read_group,
                    mate_cigar: None,
                    mate_mapping_quality: None,
                })?;
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The read group accepts the spelling bwa and STAR users already type, with
    /// or without the `@RG` prefix and with escaped tabs.
    #[test]
    fn read_group_parses_the_conventional_spellings() {
        let expected = ReadGroup {
            id: "s1".into(),
            header_line: "@RG\tID:s1\tSM:sample".into(),
        };
        for spec in [
            "ID:s1\tSM:sample",
            "ID:s1\\tSM:sample",
            "@RG\\tID:s1\\tSM:sample",
            "@RG\tID:s1\tSM:sample",
        ] {
            assert_eq!(ReadGroup::parse(spec).unwrap(), expected, "{spec:?}");
        }
    }

    /// A read group without an `ID` is rejected: every record's `RG` tag has to
    /// point at one, so a header line without it could never be referenced.
    #[test]
    fn read_group_without_an_id_is_rejected() {
        for spec in ["SM:sample", "", "@RG", "notafield"] {
            assert!(
                ReadGroup::parse(spec).is_err(),
                "{spec:?} should be rejected"
            );
        }
    }

    /// MAPQ has to distinguish a uniquely placed fragment from an ambiguous one;
    /// the old constant 1 said "ambiguous" for every record ever written.
    #[test]
    fn mapping_quality_follows_the_star_convention() {
        assert_eq!(mapping_quality(1), 255);
        assert_eq!(mapping_quality(2), 3);
        assert_eq!(mapping_quality(3), 1);
        assert_eq!(mapping_quality(4), 1);
        assert_eq!(mapping_quality(5), 0);
        assert_eq!(mapping_quality(100), 0);
        // An unmapped record has no placement to be confident about.
        assert_eq!(mapping_quality(0), 0);
    }

    /// The spoofed CIGAR still has to describe the whole read, clip included, at
    /// both transcript ends and when the read misses the transcript entirely.
    #[test]
    fn spoofed_cigar_always_covers_the_read() {
        let mut ops = Vec::new();
        for (pos, read_len, txp_len) in [
            (0, 100, 1000),
            (-20, 100, 1000),
            (-200, 100, 1000),
            (950, 100, 1000),
            (1200, 100, 1000),
        ] {
            let reported = overhang_cigar(&mut ops, pos, read_len, txp_len);
            let covered: u32 = ops
                .iter()
                .filter(|op| op.kind.consumes_read())
                .map(|op| op.len)
                .sum();
            assert_eq!(covered as i32, read_len, "pos={pos}");
            assert!(reported >= 0, "pos={pos} reported a negative position");
        }
    }
}
