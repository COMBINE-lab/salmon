//! Replicate a BAM N times into a larger BAM for align-mode timing benchmarks.
//! Each copy's records stay contiguous, so per-fragment name grouping (by
//! adjacency) is preserved — copies are simply counted as additional fragments.
//!
//! Usage: replicate_bam <in.bam> <out.bam> <n>

use noodles_bam as bam;
use noodles_sam::alignment::io::Write as _;

fn main() -> anyhow::Result<()> {
    let mut args = std::env::args().skip(1);
    let input = args.next().expect("in.bam");
    let output = args.next().expect("out.bam");
    let n: usize = args.next().expect("n").parse()?;

    let mut reader = bam::io::Reader::new(std::fs::File::open(&input)?);
    let header = reader.read_header()?;
    let recs: Vec<_> = reader
        .record_bufs(&header)
        .collect::<std::io::Result<Vec<_>>>()?;
    eprintln!("read {} records; writing {} copies", recs.len(), n);

    let mut writer = bam::io::Writer::new(std::io::BufWriter::new(std::fs::File::create(&output)?));
    writer.write_header(&header)?;
    for _ in 0..n {
        for r in &recs {
            writer.write_alignment_record(&header, r)?;
        }
    }
    writer.finish(&header)?;
    eprintln!("wrote {} records total", recs.len() * n);
    Ok(())
}
