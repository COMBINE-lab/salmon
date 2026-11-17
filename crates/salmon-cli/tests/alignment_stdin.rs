//! `salmon quant -a /dev/stdin` fed from a pipe (`samtools view -b ... | salmon
//! quant -a /dev/stdin`), SAM and BAM, on every alignment driver: the default
//! deterministic path, `--online`, and genome projection. A pipe has no
//! extension and cannot be reopened, so this needs the input opened once and its
//! format read from its bytes.
#![cfg(unix)]

use noodles_sam::alignment::io::Write as _;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

const SALMON: &str = env!("CARGO_BIN_EXE_salmon");

/// A name-grouped transcriptome SAM of inward pairs, unique and ambiguous.
fn write_txome_sam(dir: &Path) -> PathBuf {
    let mut text = String::from("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:txA\tLN:3000\n");
    text.push_str("@SQ\tSN:txB\tLN:3000\n");
    for i in 0..60 {
        let (l, r) = (1 + i * 20, 201 + i * 20);
        let txs: &[&str] = match i % 3 {
            0 => &["txA"],
            1 => &["txB"],
            _ => &["txA", "txB"],
        };
        for (k, tx) in txs.iter().enumerate() {
            let sf = if k > 0 { 256 } else { 0 };
            text.push_str(&format!(
                "p{i}\t{}\t{tx}\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:-{}\n",
                99 | sf,
                k * 3
            ));
            text.push_str(&format!(
                "p{i}\t{}\t{tx}\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:0\n",
                147 | sf
            ));
        }
    }
    let path = dir.join("aln.sam");
    std::fs::write(&path, text).unwrap();
    path
}

/// A name-grouped single-end genome SAM and the GTF of its two transcripts.
fn write_genome_sam(dir: &Path) -> (PathBuf, PathBuf) {
    let mut gtf = String::new();
    for (tx, start, end, strand) in [("txA", 1001, 3000, '+'), ("txB", 5001, 7000, '-')] {
        let attrs = format!("gene_id \"g{tx}\"; transcript_id \"{tx}\";");
        for kind in ["transcript", "exon"] {
            gtf.push_str(&format!(
                "chr1\tt\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n"
            ));
        }
    }
    let mut sam = String::from("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:10000\n");
    for i in 0..150 {
        let (flag, p) = if i % 3 == 0 { (16, 5001) } else { (0, 1001) };
        sam.push_str(&format!(
            "r{i}\t{flag}\tchr1\t{}\t60\t100M\t*\t0\t0\t*\t*\tAS:i:0\n",
            p + i
        ));
    }
    let (sam_path, gtf_path) = (dir.join("g.sam"), dir.join("a.gtf"));
    std::fs::write(&sam_path, sam).unwrap();
    std::fs::write(&gtf_path, gtf).unwrap();
    (sam_path, gtf_path)
}

fn sam_to_bam(sam: &Path, out: &Path) -> PathBuf {
    let mut reader =
        noodles_sam::io::Reader::new(std::io::BufReader::new(std::fs::File::open(sam).unwrap()));
    let header = reader.read_header().unwrap();
    let mut writer = noodles_bam::io::Writer::new(std::fs::File::create(out).unwrap());
    writer.write_header(&header).unwrap();
    for rec in reader.records() {
        writer
            .write_alignment_record(&header, &rec.unwrap())
            .unwrap();
    }
    writer.finish(&header).unwrap();
    out.to_path_buf()
}

/// Run `salmon quant -a <input|/dev/stdin> <args> -o <out>`; with `piped`, the
/// input's bytes go to salmon's stdin through a pipe. Bounded, so a regression
/// that waits on the pipe fails instead of hanging the suite. Returns quant.sf.
fn quant(input: &Path, piped: bool, args: &[&str], out: &Path) -> Vec<u8> {
    let mut cmd = Command::new(SALMON);
    cmd.args(["quant", "-a"])
        .arg(if piped {
            Path::new("/dev/stdin")
        } else {
            input
        })
        .args(args)
        .arg("-o")
        .arg(out)
        .stdout(Stdio::null())
        .stderr(Stdio::piped());
    if piped {
        cmd.stdin(Stdio::piped());
    }
    let mut child = cmd.spawn().unwrap();
    if piped {
        let mut stdin = child.stdin.take().unwrap();
        let bytes = std::fs::read(input).unwrap();
        std::thread::spawn(move || {
            let _ = stdin.write_all(&bytes);
        });
    }
    let start = Instant::now();
    while child.try_wait().unwrap().is_none() {
        if start.elapsed() > Duration::from_secs(120) {
            child.kill().unwrap();
            panic!("salmon hung on -a {}", input.display());
        }
        std::thread::sleep(Duration::from_millis(20));
    }
    let o = child.wait_with_output().unwrap();
    assert!(
        o.status.success(),
        "salmon quant -a {} {args:?} (piped={piped}) failed:\n{}",
        input.display(),
        String::from_utf8_lossy(&o.stderr)
    );
    std::fs::read(out.join("quant.sf")).unwrap()
}

#[test]
fn dev_stdin_from_a_pipe_quantifies_like_the_file() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sam = write_txome_sam(d);
    let bam = sam_to_bam(&sam, &d.join("aln.bam"));
    let (gsam, gtf) = write_genome_sam(d);
    let gbam = sam_to_bam(&gsam, &d.join("g.bam"));
    let gtf = gtf.to_str().unwrap();
    // `-l A` throughout: its pairedness peek reads ahead into the records.
    let runs: [(&str, &Path, Vec<&str>); 6] = [
        ("sam", &sam, vec!["-l", "A", "-p", "2"]),
        ("bam", &bam, vec!["-l", "A", "-p", "2"]),
        ("sam-online", &sam, vec!["-l", "A", "-p", "1", "--online"]),
        ("bam-online", &bam, vec!["-l", "A", "-p", "1", "--online"]),
        (
            "gsam",
            &gsam,
            vec!["-l", "A", "-p", "2", "--annotation", gtf],
        ),
        (
            "gbam",
            &gbam,
            vec!["-l", "A", "-p", "2", "--annotation", gtf],
        ),
    ];
    for (tag, input, args) in &runs {
        let want = quant(input, false, args, &d.join(format!("{tag}-file")));
        let got = quant(input, true, args, &d.join(format!("{tag}-pipe")));
        assert_eq!(
            String::from_utf8_lossy(&want),
            String::from_utf8_lossy(&got),
            "{tag}: -a /dev/stdin from a pipe differs from the file"
        );
    }
}
