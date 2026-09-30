//! `salmon index --unspliced` / `salmon quant` on an unspliced index, end to
//! end through the binary: option wiring and validation.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const SALMON: &str = env!("CARGO_BIN_EXE_salmon");

fn sequence(mut seed: u64, len: usize) -> String {
    (0..len)
        .map(|_| {
            seed = seed
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            "ACGT".as_bytes()[(seed >> 33) as usize & 3] as char
        })
        .collect()
}

struct Fx {
    txome: PathBuf,
    genome: PathBuf,
    gtf: PathBuf,
}

/// One two-exon gene on a 4 kb chromosome, its transcript, and the GTF.
fn fixture(dir: &Path) -> Fx {
    let chr = sequence(7, 4000);
    let tx = format!("{}{}", &chr[500..1000], &chr[2000..2500]);
    let txome = dir.join("txome.fa");
    std::fs::write(&txome, format!(">T1\n{tx}\n")).unwrap();
    let genome = dir.join("genome.fa");
    std::fs::write(&genome, format!(">chr1\n{chr}\n")).unwrap();
    let gtf = dir.join("ann.gtf");
    std::fs::write(
        &gtf,
        "chr1\tt\texon\t501\t1000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
         chr1\tt\texon\t2001\t2500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n",
    )
    .unwrap();
    Fx { txome, genome, gtf }
}

fn salmon(args: &[&std::ffi::OsStr]) -> Output {
    Command::new(SALMON).args(args).output().unwrap()
}

fn os(s: &str) -> &std::ffi::OsStr {
    s.as_ref()
}

#[test]
fn read_length_sets_the_flank_and_the_table_is_written() {
    let tmp = tempfile::tempdir().unwrap();
    let fx = fixture(tmp.path());
    let idx = tmp.path().join("idx");
    let out = salmon(&[
        os("index"),
        os("-t"),
        fx.txome.as_os_str(),
        os("-i"),
        idx.as_os_str(),
        os("-p"),
        os("1"),
        os("--unspliced"),
        os("intron"),
        os("--genome"),
        fx.genome.as_os_str(),
        os("--gtf"),
        fx.gtf.as_os_str(),
        os("--readLength"),
        os("51"),
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let info: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(idx.join("info.json")).unwrap()).unwrap();
    assert_eq!(info["unspliced"]["mode"], "intron");
    assert_eq!(info["unspliced"]["flank"], 50);
    assert_eq!(
        std::fs::read_to_string(idx.join("t2g_3col.tsv")).unwrap(),
        "T1\tG1\tS\nG1-I\tG1\tU\n"
    );
}

#[test]
fn unspliced_options_are_validated() {
    let tmp = tempfile::tempdir().unwrap();
    let fx = fixture(tmp.path());
    let idx = tmp.path().join("idx");
    let base = [
        os("index"),
        os("-t"),
        fx.txome.as_os_str(),
        os("-i"),
        idx.as_os_str(),
        os("-p"),
        os("1"),
    ];
    let run = |extra: &[&std::ffi::OsStr]| {
        let mut a = base.to_vec();
        a.extend_from_slice(extra);
        let o = salmon(&a);
        assert!(!o.status.success(), "{extra:?} should fail");
        String::from_utf8_lossy(&o.stderr).into_owned()
    };
    // --genome and --gtf go together
    run(&[
        os("--unspliced"),
        os("intron"),
        os("--genome"),
        fx.genome.as_os_str(),
    ]);
    // --flank and --readLength are two spellings of one setting
    run(&[os("--flank"), os("10"), os("--readLength"), os("51")]);
    // a mode clap does not know
    run(&[os("--unspliced"), os("splici")]);
    // intron mode needs a flank
    let e = run(&[
        os("--unspliced"),
        os("intron"),
        os("--genome"),
        fx.genome.as_os_str(),
        os("--gtf"),
        fx.gtf.as_os_str(),
    ]);
    assert!(e.contains("flank"), "{e}");
    // annotation inputs without a mode
    let e = run(&[
        os("--genome"),
        fx.genome.as_os_str(),
        os("--gtf"),
        fx.gtf.as_os_str(),
    ]);
    assert!(e.contains("without --unspliced"), "{e}");
}

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect()
}

/// Inward 2x75 pairs of 250-base fragments: `n_spliced` from the mature
/// transcript, `n_pre` from the intron (1000..2000) of the fixture gene.
fn reads(dir: &Path, n_spliced: usize, n_pre: usize) -> (PathBuf, PathBuf) {
    let chr = sequence(7, 4000);
    let tx = format!("{}{}", &chr[500..1000], &chr[2000..2500]);
    let intron = &chr[1000..2000];
    let (mut r1, mut r2) = (String::new(), String::new());
    let q = "I".repeat(75);
    let mut id = 0;
    for (src, n) in [(tx.as_str(), n_spliced), (intron, n_pre)] {
        for i in 0..n {
            let start = (i * 37) % (src.len() - 250);
            let frag = &src[start..start + 250];
            r1.push_str(&format!("@f{id}\n{}\n+\n{q}\n", &frag[..75]));
            r2.push_str(&format!("@f{id}\n{}\n+\n{q}\n", revcomp(&frag[175..])));
            id += 1;
        }
    }
    let (p1, p2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    std::fs::write(&p1, r1).unwrap();
    std::fs::write(&p2, r2).unwrap();
    (p1, p2)
}

fn build_unspliced_index(dir: &Path, fx: &Fx) -> PathBuf {
    let idx = dir.join("idx");
    let out = salmon(&[
        os("index"),
        os("-t"),
        fx.txome.as_os_str(),
        os("-i"),
        idx.as_os_str(),
        os("-p"),
        os("1"),
        os("--unspliced"),
        os("intron"),
        os("--genome"),
        fx.genome.as_os_str(),
        os("--gtf"),
        fx.gtf.as_os_str(),
        os("--readLength"),
        os("75"),
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    idx
}

fn num_reads(quant_sf: &Path) -> f64 {
    std::fs::read_to_string(quant_sf)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| l.rsplit('\t').next().unwrap().parse::<f64>().unwrap())
        .sum()
}

/// Both inference paths write `quant.genes.usa.tsv` and the `unspliced` block
/// of `meta_info.json`, and the three columns add up to `quant.sf`.
#[test]
fn quant_writes_the_spliced_unspliced_summary_on_both_paths() {
    let tmp = tempfile::tempdir().unwrap();
    let fx = fixture(tmp.path());
    let idx = build_unspliced_index(tmp.path(), &fx);
    let (r1, r2) = reads(tmp.path(), 300, 200);
    for (label, extra) in [("deterministic", None), ("online", Some("--online"))] {
        let out = tmp.path().join(label);
        let mut args = vec![
            os("quant"),
            os("-i"),
            idx.as_os_str(),
            os("-l"),
            os("A"),
            os("-1"),
            r1.as_os_str(),
            os("-2"),
            r2.as_os_str(),
            os("-p"),
            os("2"),
            os("-o"),
            out.as_os_str(),
        ];
        if let Some(e) = extra {
            args.push(os(e));
        }
        let o = salmon(&args);
        assert!(
            o.status.success(),
            "{label}: {}",
            String::from_utf8_lossy(&o.stderr)
        );

        let usa = std::fs::read_to_string(out.join("quant.genes.usa.tsv")).unwrap();
        let mut lines = usa.lines();
        assert_eq!(lines.next(), Some("Name\tSpliced\tUnspliced\tAmbiguous"));
        let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
        assert_eq!(row[0], "G1");
        let (s, u, a): (f64, f64, f64) = (
            row[1].parse().unwrap(),
            row[2].parse().unwrap(),
            row[3].parse().unwrap(),
        );
        let total = num_reads(&out.join("quant.sf"));
        assert!(
            (s + u + a - total).abs() < 0.01 * total,
            "{label}: {s}+{u}+{a} vs {total}"
        );
        // Every simulated intronic fragment lies inside the intron target, and
        // fragments that stay within the flank are compatible with both forms.
        assert!(u > 150.0 && u < 205.0, "{label}: unspliced {u}");
        assert!(s > 150.0, "{label}: spliced {s}");

        let meta: serde_json::Value = serde_json::from_str(
            &std::fs::read_to_string(out.join("aux_info/meta_info.json")).unwrap(),
        )
        .unwrap();
        let m = &meta["unspliced"];
        assert_eq!(m["mode"], "intron", "{label}");
        assert_eq!(m["num_unspliced_targets"], 1);
        let f = m["unspliced_fraction"].as_f64().unwrap();
        assert!((f - u / (s + u + a)).abs() < 1e-3, "{label}: {f}");
    }
}

/// With the genome as decoys, intergenic fragments go to the decoy and form
/// the fourth category of the run-level split; intronic ones tie between the
/// intron target and the genome and stay unspliced.
#[test]
fn decoy_fragments_are_the_fourth_category() {
    let tmp = tempfile::tempdir().unwrap();
    let fx = fixture(tmp.path());
    let chr = sequence(7, 4000);
    let gentrome = tmp.path().join("gentrome.fa");
    std::fs::write(
        &gentrome,
        format!(
            "{}>chr1\n{chr}\n",
            std::fs::read_to_string(&fx.txome).unwrap()
        ),
    )
    .unwrap();
    let decoys = tmp.path().join("decoys.txt");
    std::fs::write(&decoys, "chr1\n").unwrap();
    let idx = tmp.path().join("idx");
    let o = salmon(&[
        os("index"),
        os("-t"),
        gentrome.as_os_str(),
        os("-d"),
        decoys.as_os_str(),
        os("-i"),
        idx.as_os_str(),
        os("-p"),
        os("1"),
        os("--unspliced"),
        os("intron"),
        os("--genome"),
        fx.genome.as_os_str(),
        os("--gtf"),
        fx.gtf.as_os_str(),
        os("--readLength"),
        os("75"),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));

    // 300 spliced + 200 intronic fragments, then 100 intergenic ones
    let (r1, r2) = reads(tmp.path(), 300, 200);
    let q = "I".repeat(75);
    let (mut a, mut b) = (
        std::fs::read_to_string(&r1).unwrap(),
        std::fs::read_to_string(&r2).unwrap(),
    );
    let inter = &chr[3000..4000];
    for i in 0..100 {
        let start = (i * 7) % (inter.len() - 250);
        let frag = &inter[start..start + 250];
        a.push_str(&format!("@g{i}\n{}\n+\n{q}\n", &frag[..75]));
        b.push_str(&format!("@g{i}\n{}\n+\n{q}\n", revcomp(&frag[175..])));
    }
    std::fs::write(&r1, a).unwrap();
    std::fs::write(&r2, b).unwrap();

    let out = tmp.path().join("q");
    let o = salmon(&[
        os("quant"),
        os("-i"),
        idx.as_os_str(),
        os("-l"),
        os("A"),
        os("-1"),
        r1.as_os_str(),
        os("-2"),
        r2.as_os_str(),
        os("-p"),
        os("2"),
        os("-o"),
        out.as_os_str(),
    ]);
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let meta: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(out.join("aux_info/meta_info.json")).unwrap(),
    )
    .unwrap();
    let m = &meta["unspliced"];
    let decoy = m["num_decoy_fragments"].as_f64().unwrap();
    assert_eq!(decoy, meta["num_decoy_fragments"].as_f64().unwrap());
    assert!((90.0..=100.0).contains(&decoy), "decoy fragments {decoy}");
    let w = &m["fractions_with_decoys"];
    let f = |k: &str| w[k].as_f64().unwrap();
    let sum = f("spliced") + f("unspliced") + f("ambiguous") + f("decoy");
    assert!((sum - 1.0).abs() < 1e-9, "{w}");
    // the intronic fragments were not lost to the decoy
    let u = m["num_unspliced_fragments"].as_f64().unwrap();
    assert!(u > 150.0, "unspliced {u}");
}
