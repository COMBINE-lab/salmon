//! Dump the deduplicated transcriptome held by a salmon index as FASTA, so a
//! simulator (polyester) draws reads from exactly the transcript set that will
//! be quantified — no off-by-one between the simulated truth and the index.
//!
//! Usage: dump_txome <index_dir> <out.fa>

use salmon_index::SalmonIndex;
use std::io::{BufWriter, Write};

fn main() {
    let mut args = std::env::args().skip(1);
    let index_dir = args.next().expect("index dir");
    let out_path = args.next().expect("out fasta");
    let idx = SalmonIndex::load(index_dir).expect("load index");
    let n = idx.num_refs();
    let f = std::fs::File::create(&out_path).expect("create out");
    let mut w = BufWriter::new(f);
    for tid in 0..n {
        let name = idx.ref_name(tid);
        let seq = idx.ref_seq(tid as u32);
        writeln!(w, ">{name}").unwrap();
        // 70-col wrap (FASTA convention; keeps downstream parsers happy)
        for chunk in seq.chunks(70) {
            w.write_all(chunk).unwrap();
            w.write_all(b"\n").unwrap();
        }
    }
    w.flush().unwrap();
    eprintln!("wrote {n} transcripts to {out_path}");
}
