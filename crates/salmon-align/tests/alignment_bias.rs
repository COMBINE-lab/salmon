//! Smoke coverage for bias correction through the direct alignment-mode entry
//! point. RAD has more exhaustive invariance coverage; this test ensures the
//! shared correction remains wired into BAM/SAM quantification.

use salmon_align::{quantify_alignments, AlignQuantOptions};
use std::io::Write;
use std::path::{Path, PathBuf};

fn write_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let sam = dir.join("bias.sam");
    let fasta = dir.join("transcripts.fa");
    let names = ["txA", "txB", "txC"];
    let mut file = std::fs::File::create(&sam).unwrap();
    writeln!(file, "@HD\tVN:1.6\tSO:queryname").unwrap();
    for name in names {
        writeln!(file, "@SQ\tSN:{name}\tLN:3000").unwrap();
    }

    let mut fragment = 0;
    for (tid, count) in [(0usize, 200), (1, 150), (2, 250)] {
        for i in 0..count {
            let left = 1 + (i * 7) % 2500;
            let right = left + 200;
            writeln!(
                file,
                "r{fragment}\t99\t{}\t{left}\t60\t100M\t=\t{right}\t300\t*\t*",
                names[tid]
            )
            .unwrap();
            writeln!(
                file,
                "r{fragment}\t147\t{}\t{right}\t60\t100M\t=\t{left}\t-300\t*\t*",
                names[tid]
            )
            .unwrap();
            fragment += 1;
        }
    }

    let ref_seqs: Vec<Vec<u8>> = (0..names.len())
        .map(|tid| (0..3000).map(|pos| b"ACGT"[(pos + tid) % 4]).collect())
        .collect();
    let mut file = std::fs::File::create(&fasta).unwrap();
    for (name, seq) in names.iter().zip(&ref_seqs) {
        writeln!(file, ">{name}\n{}", String::from_utf8_lossy(seq)).unwrap();
    }
    (sam, fasta)
}

#[test]
/// All three bias corrections at once through the alignment entry point: the
/// corrections are shared with the reads path, so this only has to prove they are
/// still wired up and produce a sane result, not re-verify the models.
fn all_bias_alignment_mode_smoke() {
    let tmp = tempfile::tempdir().unwrap();
    let (sam, fasta) = write_fixture(tmp.path());
    let out = tmp.path().join("quant");
    let mut opts = AlignQuantOptions::new(sam, out);
    opts.lib_type = "IU".to_string();
    opts.no_error_model = true;
    opts.transcripts = Some(fasta);
    opts.seq_bias = true;
    opts.gc_bias = true;
    opts.pos_bias = true;

    let result = rayon::ThreadPoolBuilder::new()
        .num_threads(4)
        .build()
        .unwrap()
        .install(|| quantify_alignments(&opts).expect("all-bias alignment quantification"));

    assert_eq!(result.num_mapped, 600);
    assert!(result.eff_lengths.iter().all(|value| value.is_finite()));
    assert!(!result.bias_dump.obs_gc.is_empty());
    assert!(!result.bias_dump.exp_gc.is_empty());
    assert!(!result.bias_dump.obs5_seq.is_empty());
    assert!(!result.bias_dump.exp5_seq.is_empty());
    assert!(!result.bias_dump.obs5_pos.is_empty());
    assert!(!result.bias_dump.exp5_pos.is_empty());
}
