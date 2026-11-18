//! Alignment input from a pipe. Every alignment path makes one pass over the
//! records and opens its input once, so a named pipe, `<(samtools view -b ...)`
//! or `/dev/stdin` fed from a pipe is read exactly as the file it carries would
//! be. SAM and BAM are told apart by content, since a pipe's name says nothing.
#![cfg(unix)]

use noodles_sam::alignment::io::Write as _;
use salmon_align::{
    project_genome_bam_to_rad, quantify_alignments, quantify_rad, write_alignment_rad,
    AlignQuantOptions, GenomeProjectOptions,
};
use salmon_rad::ChunkCodec;
use std::path::{Path, PathBuf};
use std::time::Duration;

/// A name-grouped two-transcript SAM of inward pairs: unique to txA, unique to
/// txB, and ambiguous (txA primary + txB secondary).
fn write_txome_sam(dir: &Path) -> PathBuf {
    let mut text = String::from("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:txA\tLN:5000\n");
    text.push_str("@SQ\tSN:txB\tLN:5000\n");
    let mut pair = |name: &str, tx: &str, p: usize, secondary: bool| {
        let (l, r) = (p + 1, p + 201);
        let sf = if secondary { 256 } else { 0 };
        text.push_str(&format!(
            "{name}\t{}\t{tx}\t{l}\t60\t100M\t=\t{r}\t300\t*\t*\tAS:i:{}\n",
            99 | sf,
            -(p as i64 % 7)
        ));
        text.push_str(&format!(
            "{name}\t{}\t{tx}\t{r}\t60\t100M\t=\t{l}\t-300\t*\t*\tAS:i:0\n",
            147 | sf
        ));
    };
    for i in 0..40 {
        pair(&format!("a{i}"), "txA", 1 + i * 100, false);
    }
    for i in 0..20 {
        pair(&format!("b{i}"), "txB", 1 + i * 100, false);
    }
    for i in 0..30 {
        pair(&format!("amb{i}"), "txA", 2001 + i * 50, false);
        pair(&format!("amb{i}"), "txB", 2001 + i * 50, true);
    }
    let path = dir.join("aln.sam");
    std::fs::write(&path, text).unwrap();
    path
}

/// A name-grouped single-end genome SAM over two annotated transcripts on
/// opposite strands, plus the GTF describing them.
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
    for i in 0..200 {
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

/// The same records as a BAM, at `out`.
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

/// A named pipe at `fifo` that a detached thread fills with `src`'s bytes.
fn fifo_of(src: &Path, fifo: &Path) -> PathBuf {
    assert!(std::process::Command::new("mkfifo")
        .arg(fifo)
        .status()
        .unwrap()
        .success());
    let (src, dst) = (src.to_path_buf(), fifo.to_path_buf());
    std::thread::spawn(move || {
        // A refused input closes the read end early; the writer just stops.
        let _ = std::io::copy(
            &mut std::fs::File::open(src).unwrap(),
            &mut std::fs::File::create(dst).unwrap(),
        );
    });
    fifo.to_path_buf()
}

/// Run `f` on another thread, failing the test rather than hanging the suite if
/// a pipe is waited on that nobody will ever write to.
fn bounded<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> T {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let _ = tx.send(f());
    });
    rx.recv_timeout(Duration::from_secs(60))
        .expect("salmon hung on an alignment pipe")
}

fn opts(input: &Path, out: &Path, lib: &str) -> AlignQuantOptions {
    std::fs::create_dir_all(out).unwrap();
    let mut o = AlignQuantOptions::new(input.to_path_buf(), out.to_path_buf());
    o.lib_type = lib.to_string();
    o.no_error_model = true;
    o
}

/// The default (deterministic) path: write the RAD, quantify it. Returns the
/// RAD and the quant.sf bytes.
fn deterministic(input: &Path, out: &Path, lib: &str) -> (Vec<u8>, Vec<u8>) {
    let (input, out, lib) = (input.to_path_buf(), out.to_path_buf(), lib.to_string());
    bounded(move || {
        let o = opts(&input, &out, &lib);
        let rad = out.join("aln.rad");
        write_alignment_rad(&o, &rad, ChunkCodec::None).unwrap();
        quantify_rad(&o, &rad).unwrap();
        (
            std::fs::read(&rad).unwrap(),
            std::fs::read(out.join("quant.sf")).unwrap(),
        )
    })
}

/// The `--online` path. Returns quant.sf.
fn online(input: &Path, out: &Path, lib: &str) -> Vec<u8> {
    let (input, out, lib) = (input.to_path_buf(), out.to_path_buf(), lib.to_string());
    bounded(move || {
        quantify_alignments(&opts(&input, &out, &lib)).unwrap();
        std::fs::read(out.join("quant.sf")).unwrap()
    })
}

/// Genome projection to RAD, then quantification. Returns RAD and quant.sf.
fn projected(input: &Path, gtf: &Path, out: &Path, lib: &str) -> (Vec<u8>, Vec<u8>) {
    let (input, gtf, out, lib) = (
        input.to_path_buf(),
        gtf.to_path_buf(),
        out.to_path_buf(),
        lib.to_string(),
    );
    bounded(move || {
        let gp = GenomeProjectOptions {
            bam: input.clone(),
            annotation: gtf,
            genome_fasta: None,
            lib_type: lib.clone(),
            junc_miss_discount: None,
            bias: false,
            fld_mean: 250.0,
            fld_sd: 25.0,
            fld_max: 1000,
            bias_seed_em_iters: 0,
            em: Default::default(),
            rad_codec: ChunkCodec::None,
        };
        let rad = out.join("g.rad");
        std::fs::create_dir_all(&out).unwrap();
        project_genome_bam_to_rad(&gp, &rad).unwrap();
        quantify_rad(&opts(&input, &out, &lib), &rad).unwrap();
        (
            std::fs::read(&rad).unwrap(),
            std::fs::read(out.join("quant.sf")).unwrap(),
        )
    })
}

/// Both drivers over the transcriptome input, SAM and BAM, explicit and `-l A`
/// (whose pairedness peek reads ahead into the records): a FIFO gives the same
/// RAD and the same abundances as the file it carries.
#[test]
fn a_fifo_alignment_input_reads_like_its_file() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sam = write_txome_sam(d);
    let bam = sam_to_bam(&sam, &d.join("aln.bam"));
    for (fmt, file) in [("sam", &sam), ("bam", &bam)] {
        for lib in ["IU", "A"] {
            let tag = format!("{fmt}-{lib}");
            let pipe = |n: &str| fifo_of(file, &d.join(format!("{tag}-{n}.fifo")));
            let want = deterministic(file, &d.join(format!("{tag}-file")), lib);
            let got = deterministic(&pipe("det"), &d.join(format!("{tag}-fifo")), lib);
            assert!(
                want.0 == got.0,
                "{tag}: RAD from a FIFO differs from the file's"
            );
            assert_eq!(want.1, got.1, "{tag}: quant.sf from a FIFO differs");
            let want = online(file, &d.join(format!("{tag}-ofile")), lib);
            let got = online(&pipe("online"), &d.join(format!("{tag}-ofifo")), lib);
            assert_eq!(want, got, "{tag}: --online quant.sf from a FIFO differs");
        }
    }
}

/// The genome-projection driver, SAM and BAM, under `-l A`.
#[test]
fn a_fifo_genome_alignment_input_projects_like_its_file() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let (sam, gtf) = write_genome_sam(d);
    let bam = sam_to_bam(&sam, &d.join("g.bam"));
    for (fmt, file) in [("sam", &sam), ("bam", &bam)] {
        let want = projected(file, &gtf, &d.join(format!("{fmt}-file")), "A");
        let fifo = fifo_of(file, &d.join(format!("{fmt}.fifo")));
        let got = projected(&fifo, &gtf, &d.join(format!("{fmt}-fifo")), "A");
        assert!(want.0 == got.0, "{fmt}: projected RAD from a FIFO differs");
        assert_eq!(want.1, got.1, "{fmt}: quant.sf from a FIFO differs");
    }
}

/// SAM vs BAM comes from the bytes, not the name: a BAM called `.sam` and a SAM
/// called `.bam` (or with no extension, like `/dev/fd/63`) both read correctly.
#[test]
fn alignment_format_is_read_from_content_not_name() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sam = write_txome_sam(d);
    let bam = sam_to_bam(&sam, &d.join("aln.bam"));
    let want = deterministic(&sam, &d.join("ref"), "IU").1;
    let sam_as_bam = d.join("really_sam.bam");
    std::fs::copy(&sam, &sam_as_bam).unwrap();
    let bam_as_sam = d.join("really_bam.sam");
    std::fs::copy(&bam, &bam_as_sam).unwrap();
    let bare = d.join("noext");
    std::fs::copy(&bam, &bare).unwrap();
    for p in [&sam_as_bam, &bam_as_sam, &bare] {
        let name = p.file_name().unwrap().to_string_lossy().into_owned();
        assert_eq!(
            deterministic(p, &d.join(format!("q-{name}")), "IU").1,
            want,
            "{name}"
        );
    }
}

/// `--errorModel` on the deterministic path is the one pass that reads the
/// records twice, so it refuses a pipe at once instead of hanging on it.
#[test]
fn the_two_pass_error_model_refuses_a_fifo() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let sam = write_txome_sam(d);
    let fa = d.join("t.fa");
    let seq = "ACGT".repeat(1250);
    std::fs::write(&fa, format!(">txA\n{seq}\n>txB\n{seq}\n")).unwrap();
    let fifo = fifo_of(&sam, &d.join("em.fifo"));
    let out = d.join("q");
    let err = bounded(move || {
        let mut o = opts(&fifo, &out, "IU");
        o.no_error_model = false;
        o.deterministic_error_model = true;
        o.transcripts = Some(fa);
        write_alignment_rad(&o, &out.join("aln.rad"), ChunkCodec::None)
            .map(|_| ())
            .map_err(|e| format!("{e:#}"))
    })
    .expect_err("a FIFO under --errorModel must be refused");
    assert!(err.contains("regular file"), "{err}");
}

/// Unreadable input still fails with a message that names it.
#[test]
fn missing_or_empty_alignment_input_is_reported() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let out = d.join("q");
    let run = |p: &Path| {
        write_alignment_rad(&opts(p, &out, "IU"), &out.join("a.rad"), ChunkCodec::None)
            .map(|_| ())
            .map_err(|e| format!("{e:#}"))
            .expect_err("must fail")
    };
    let missing = d.join("nope.bam");
    let err = run(&missing);
    assert!(err.contains("opening alignment input"), "{err}");
    let empty = d.join("empty.bam");
    std::fs::write(&empty, "").unwrap();
    let err = run(&empty);
    assert!(err.contains("is empty"), "{err}");
}
