//! `--writeMappings` SAM output.

use std::fmt::Write as _;
use std::io::Write;
use std::sync::Mutex;

use salmon_index::SalmonIndex;
use salmon_map::ScoredMapping;

use crate::mapping_record::{self, AlignmentRecord, CigarKind};

const SAM_VERSION: &str = "1.0";

pub struct SamWriter {
    inner: Mutex<Box<dyn Write + Send>>,
}

impl SamWriter {
    pub fn create(path: &std::path::Path, salmon: &SalmonIndex, cmd: &str) -> anyhow::Result<Self> {
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
        writeln!(writer, "@HD\tVN:{SAM_VERSION}\tSO:unknown")?;
        for tid in 0..salmon.num_refs() {
            writeln!(
                writer,
                "@SQ\tSN:{}\tLN:{}",
                salmon.ref_name(tid),
                salmon.ref_len(tid)
            )?;
        }
        writeln!(
            writer,
            "@PG\tID:salmon\tPN:salmon\tVN:{}\tCL:{}",
            crate::output::SALMON_VERSION,
            cmd
        )?;
        Ok(Self {
            inner: Mutex::new(writer),
        })
    }

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

fn write_sequence(buf: &mut String, record: &AlignmentRecord<'_>) {
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

fn write_record(buf: &mut String, salmon: &SalmonIndex, record: &AlignmentRecord<'_>) {
    let name = String::from_utf8_lossy(record.name);
    let reference = salmon.ref_name(record.reference_id);
    let _ = write!(
        buf,
        "{name}\t{}\t{reference}\t{}\t{}",
        record.flags,
        record.position + 1,
        record.mapping_quality
    );
    buf.push('\t');
    for op in record.cigar.as_slice() {
        let kind = match op.kind {
            CigarKind::Match => 'M',
            CigarKind::SoftClip => 'S',
        };
        let _ = write!(buf, "{}{kind}", op.len);
    }
    if let (Some(mate_id), Some(mate_position)) = (record.mate_reference_id, record.mate_position) {
        let mate_reference = if mate_id == record.reference_id {
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
    let _ = writeln!(
        buf,
        "\t*\tNH:i:{}\tHI:i:{}\tXT:A:{}\tAS:i:{}",
        record.nh, record.hi, record.xt as char, record.alignment_score
    );
}

pub fn write_fragment(
    buf: &mut String,
    salmon: &SalmonIndex,
    r1_id: &[u8],
    r1_seq: &[u8],
    r2: Option<(&[u8], &[u8])>,
    maps: &[ScoredMapping],
) {
    mapping_record::emit_fragment_records(salmon, r1_id, r1_seq, r2, maps, |record| {
        write_record(buf, salmon, record);
        Ok(())
    })
    .expect("writing to String is infallible");
}
