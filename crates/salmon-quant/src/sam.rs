//! `--writeMappings` SAM output.
//!
//! SAM is the plain-text alignment format: one tab-separated line per placed
//! read, after a header describing the references. Every field here is derived by
//! [`crate::mapping_record`], which BAM output shares, so this file is purely the
//! text encoding.

use std::fmt::Write as _;
use std::io::Write;
use std::sync::Mutex;

use salmon_index::SalmonIndex;
use salmon_map::ScoredMapping;

use crate::mapping_record::{self, AlignmentRecord, EmitScratch, RecordOptions};

/// Shared SAM sink. Worker threads format whole blocks of records into their own
/// string buffers and hand them over; the mutex is taken once per block rather
/// than once per record, so it is never contended in practice.
pub struct SamWriter {
    inner: Mutex<Box<dyn Write + Send>>,
}

impl SamWriter {
    pub fn create(
        path: &std::path::Path,
        salmon: &SalmonIndex,
        cmd: &str,
        options: &RecordOptions<'_>,
    ) -> anyhow::Result<Self> {
        use anyhow::Context as _;
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                std::fs::create_dir_all(parent)
                    .with_context(|| format!("creating SAM output dir {}", parent.display()))?;
            }
        }
        let file = std::fs::File::create(path)
            .with_context(|| format!("creating SAM output {}", path.display()))?;
        let mut writer: Box<dyn Write + Send> =
            Box::new(std::io::BufWriter::with_capacity(1 << 20, file));
        // Same header text BAM embeds, so the two formats cannot drift.
        writer.write_all(mapping_record::header_text(salmon, cmd, options).as_bytes())?;
        Ok(Self {
            inner: Mutex::new(writer),
        })
    }

    /// Append one worker's accumulated block of records.
    ///
    /// Blocks land in whatever order threads finish, so records for one fragment
    /// are contiguous but no order is imposed across fragments, as documented
    /// for `--writeMappings`.
    pub fn write_block(&self, block: &str) -> std::io::Result<()> {
        if block.is_empty() {
            return Ok(());
        }
        self.inner.lock().unwrap().write_all(block.as_bytes())
    }

    pub fn flush(&self) -> std::io::Result<()> {
        self.inner.lock().unwrap().flush()
    }
}

/// Write the read's bases, reverse-complementing them when the alignment is on
/// the reverse strand: SAM always stores the sequence as it appears on the
/// reference's forward strand.
fn write_sequence(buf: &mut String, record: &AlignmentRecord<'_>) {
    if record.sequence.is_empty() {
        buf.push('*');
        return;
    }
    if record.reverse_complement {
        for &base in record.sequence.iter().rev() {
            buf.push(mapping_record::complement(base) as char);
        }
    } else {
        for &base in record.sequence {
            buf.push(base as char);
        }
    }
}

/// Write `QUAL`: the Phred+33 characters as the FASTQ carried them, reversed to
/// match a reverse-strand `SEQ` so that base *i* of one still describes base *i*
/// of the other. `*` when the input had no qualities (a FASTA), which is SAM's
/// spelling for "not available".
///
/// A quality string of the wrong length is dropped rather than written: readers
/// reject a record whose `QUAL` and `SEQ` disagree, so a truncated FASTQ record
/// must not be able to produce an unreadable file.
fn write_qualities(buf: &mut String, record: &AlignmentRecord<'_>) {
    match mapping_record::usable_qualities(record) {
        Some(qualities) if record.reverse_complement => {
            for &q in qualities.iter().rev() {
                buf.push(q as char);
            }
        }
        Some(qualities) => {
            for &q in qualities {
                buf.push(q as char);
            }
        }
        None => buf.push('*'),
    }
}

/// Serialize one record as a SAM line: the eleven mandatory tab-separated fields
/// followed by salmon's optional tags.
fn write_record(buf: &mut String, salmon: &SalmonIndex, record: &AlignmentRecord<'_>) {
    let name = String::from_utf8_lossy(record.name);
    // QNAME, FLAG, RNAME, POS, MAPQ. SAM positions are 1-based, ours are 0-based,
    // and an unmapped record has neither a reference nor a position.
    let _ = write!(buf, "{name}\t{}\t", record.flags);
    match record.reference_id {
        Some(tid) => {
            let _ = write!(
                buf,
                "{}\t{}\t{}",
                salmon.ref_name(tid),
                record.position + 1,
                record.mapping_quality
            );
        }
        None => {
            let _ = write!(buf, "*\t0\t{}", record.mapping_quality);
        }
    }
    buf.push('\t');
    if record.cigar.is_empty() {
        buf.push('*');
    } else {
        for op in record.cigar {
            let _ = write!(buf, "{}{}", op.len, op.kind.letter());
        }
    }
    // RNEXT, PNEXT, TLEN. `=` is SAM's shorthand for "same reference as this
    // record", which is always the case for a proper pair here.
    if let (Some(mate_id), Some(mate_position)) = (record.mate_reference_id, record.mate_position) {
        let mate_reference = if Some(mate_id) == record.reference_id {
            "="
        } else {
            salmon.ref_name(mate_id)
        };
        let _ = write!(
            buf,
            "\t{mate_reference}\t{}\t{}",
            mate_position + 1,
            record.template_length
        );
    } else {
        let _ = write!(buf, "\t*\t0\t{}", record.template_length);
    }
    buf.push('\t');
    write_sequence(buf, record);
    buf.push('\t');
    write_qualities(buf, record);
    write_tags(buf, record);
    buf.push('\n');
}

/// Append the optional tags.
///
/// An unmapped record carries none of the placement tags: `NH`/`HI` describe a
/// set of placements it does not have, and `AS`/`XT` describe an alignment that
/// was never made. Only the read group, which is a property of the read rather
/// than of any placement, survives.
fn write_tags(buf: &mut String, record: &AlignmentRecord<'_>) {
    {
        // NH = number of placements for this fragment, HI = which one this is,
        // XT = transcript/decoy, AS = score, ZW = the mapping's EM weight.
        let _ = write!(
            buf,
            "\tNH:i:{}\tHI:i:{}\tXT:A:{}\tAS:i:{}\tZW:f:{}",
            record.nh, record.hi, record.xt as char, record.alignment_score, record.weight
        );
        if let Some(edit_distance) = record.edit_distance {
            let _ = write!(buf, "\tNM:i:{edit_distance}");
        }
        if let Some(md) = record.md {
            let _ = write!(buf, "\tMD:Z:{md}");
        }
        if let Some(mate_mapping_quality) = record.mate_mapping_quality {
            let _ = write!(buf, "\tMQ:i:{mate_mapping_quality}");
        }
        if let Some(mate_cigar) = record.mate_cigar {
            if !mate_cigar.is_empty() {
                buf.push_str("\tMC:Z:");
                for op in mate_cigar {
                    let _ = write!(buf, "{}{}", op.len, op.kind.letter());
                }
            }
        }
    }
    if let Some(read_group) = record.read_group {
        let _ = write!(buf, "\tRG:Z:{read_group}");
    }
}

/// Append every SAM record implied by one fragment's placements.
#[allow(clippy::too_many_arguments)]
pub fn write_fragment(
    buf: &mut String,
    salmon: &SalmonIndex,
    r1_id: &[u8],
    r1_seq: &[u8],
    r1_qual: Option<&[u8]>,
    r2: Option<(&[u8], &[u8], Option<&[u8]>)>,
    maps: &[ScoredMapping],
    options: &RecordOptions<'_>,
    scratch: &mut EmitScratch,
) {
    mapping_record::emit_fragment_records(
        salmon,
        r1_id,
        r1_seq,
        r1_qual,
        r2,
        maps,
        options,
        scratch,
        |record| {
            write_record(buf, salmon, record);
            Ok(())
        },
    )
    .expect("writing to String is infallible");
}
