//! Total-RNA simulation: a known mix of spliced and pre-mRNA fragments,
//! quantified against a decoy-aware gentrome index with and without
//! `--unspliced intron`.
//!
//! The situation it reproduces: every gene has a `retained_intron` isoform,
//! and half of the fragments come from unspliced pre-mRNA. Against the plain
//! gentrome, a pre-mRNA fragment inside the retained intron scores the same on
//! that isoform as on the genome decoy, the tie goes to the transcript, and the
//! retained-intron isoform soaks up intronic signal. With intron targets in the
//! index those fragments have a better home, and the EM shares them by
//! abundance. The test asserts that the retained-intron excess shrinks and the
//! spliced isoform estimates get closer to the truth; run it with
//! `--nocapture` to print the metrics.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

const SALMON: &str = env!("CARGO_BIN_EXE_salmon");
const READ_LEN: usize = 75;
const NUM_GENES: usize = 12;
const GENE_SPACING: usize = 10_000;
const GENE_OFFSET: usize = 5_000;
/// Exons relative to the gene start: E1, E2, E3, E4.
const EXONS: [(usize, usize); 4] = [(0, 300), (2300, 2500), (4500, 4700), (6700, 7100)];

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

fn sequence(rng: &mut Rng, len: usize) -> String {
    (0..len)
        .map(|_| "ACGT".as_bytes()[rng.below(4)] as char)
        .collect()
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

struct Isoform {
    name: String,
    gene: usize,
    /// exons in genome coordinates, ascending
    exons: Vec<(usize, usize)>,
    retained_intron: bool,
    fragments: usize,
}

struct Sim {
    gentrome: PathBuf,
    decoys: PathBuf,
    genome: PathBuf,
    gtf: PathBuf,
    r1: PathBuf,
    r2: PathBuf,
    /// true spliced fragments per isoform
    truth: HashMap<String, usize>,
    retained: Vec<String>,
    pre_mrna_fragments: usize,
}

fn simulate(dir: &Path) -> Sim {
    let mut rng = Rng(0x5eed);
    let chr = sequence(&mut rng, GENE_OFFSET * 2 + NUM_GENES * GENE_SPACING);

    let mut isoforms = Vec::new();
    for g in 0..NUM_GENES {
        let o = GENE_OFFSET + g * GENE_SPACING;
        let ex = |i: usize| (o + EXONS[i].0, o + EXONS[i].1);
        let models = [
            ("A", vec![ex(0), ex(1), ex(2), ex(3)], false, 200 + 50 * g),
            ("B", vec![ex(0), ex(2), ex(3)], false, 100 + (g % 3) * 80),
            (
                "RI",
                vec![(ex(0).0, ex(1).1), ex(2), ex(3)],
                true,
                30 + (g % 4) * 20,
            ),
        ];
        for (suffix, exons, retained_intron, fragments) in models {
            isoforms.push(Isoform {
                name: format!("T{g}{suffix}"),
                gene: g,
                exons,
                retained_intron,
                fragments,
            });
        }
    }
    let minus = |g: usize| g % 2 == 1;
    let mature = |iso: &Isoform| {
        let s: String = iso.exons.iter().map(|&(a, b)| &chr[a..b]).collect();
        if minus(iso.gene) {
            revcomp(&s)
        } else {
            s
        }
    };

    // references: transcripts, then the chromosome as the decoy
    let mut fa = String::new();
    for iso in &isoforms {
        fa.push_str(&format!(">{}\n{}\n", iso.name, mature(iso)));
    }
    fa.push_str(&format!(">chr1\n{chr}\n"));
    let gentrome = dir.join("gentrome.fa");
    std::fs::write(&gentrome, fa).unwrap();
    let decoys = dir.join("decoys.txt");
    std::fs::write(&decoys, "chr1\n").unwrap();
    let genome = dir.join("genome.fa");
    std::fs::write(&genome, format!(">chr1\n{chr}\n")).unwrap();
    let mut gtf = String::new();
    for iso in &isoforms {
        let strand = if minus(iso.gene) { '-' } else { '+' };
        let ttype = if iso.retained_intron {
            "retained_intron"
        } else {
            "protein_coding"
        };
        for &(a, b) in &iso.exons {
            gtf.push_str(&format!(
                "chr1\tsim\texon\t{}\t{b}\t.\t{strand}\t.\tgene_id \"G{}\"; transcript_id \"{}\"; transcript_type \"{ttype}\";\n",
                a + 1,
                iso.gene,
                iso.name
            ));
        }
    }
    let gtf_path = dir.join("sim.gtf");
    std::fs::write(&gtf_path, gtf).unwrap();

    // fragments: 200..300 bases, unstranded, 2x75 inward pairs
    let (mut r1, mut r2) = (String::new(), String::new());
    let qual = "I".repeat(READ_LEN);
    let mut id = 0usize;
    let mut emit = |src: &str, rng: &mut Rng, id: &mut usize| {
        let flen = 200 + rng.below(101);
        let start = rng.below(src.len() - flen + 1);
        let mut frag = src[start..start + flen].to_string();
        if rng.below(2) == 1 {
            frag = revcomp(&frag);
        }
        r1.push_str(&format!("@f{id}\n{}\n+\n{qual}\n", &frag[..READ_LEN]));
        r2.push_str(&format!(
            "@f{id}\n{}\n+\n{qual}\n",
            revcomp(&frag[flen - READ_LEN..])
        ));
        *id += 1;
    };
    let mut truth = HashMap::new();
    let mut per_gene = [0usize; NUM_GENES];
    for iso in &isoforms {
        let m = mature(iso);
        for _ in 0..iso.fragments {
            emit(&m, &mut rng, &mut id);
        }
        truth.insert(iso.name.clone(), iso.fragments);
        per_gene[iso.gene] += iso.fragments;
    }
    // as many pre-mRNA fragments as mature ones, spread over the gene body
    let mut pre_mrna_fragments = 0;
    for (g, &n) in per_gene.iter().enumerate() {
        let o = GENE_OFFSET + g * GENE_SPACING;
        let body = &chr[o..o + EXONS[3].1];
        for _ in 0..n {
            emit(body, &mut rng, &mut id);
        }
        pre_mrna_fragments += n;
    }
    let (p1, p2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    std::fs::write(&p1, r1).unwrap();
    std::fs::write(&p2, r2).unwrap();

    Sim {
        gentrome,
        decoys,
        genome,
        gtf: gtf_path,
        r1: p1,
        r2: p2,
        truth,
        retained: isoforms
            .iter()
            .filter(|i| i.retained_intron)
            .map(|i| i.name.clone())
            .collect(),
        pre_mrna_fragments,
    }
}

fn run(args: &[&std::ffi::OsStr]) {
    let o = Command::new(SALMON).args(args).output().unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
}

fn os(s: &str) -> &std::ffi::OsStr {
    s.as_ref()
}

/// Build a gentrome index (plus unspliced targets of `mode` in `layout`) and
/// quantify the simulated reads; returns NumReads by target name.
fn index_and_quant(
    sim: &Sim,
    dir: &Path,
    unspliced: Option<&str>,
    layout: &str,
) -> HashMap<String, f64> {
    let idx = dir.join("idx");
    let mut args = vec![
        os("index"),
        os("-t"),
        sim.gentrome.as_os_str(),
        os("-d"),
        sim.decoys.as_os_str(),
        os("-i"),
        idx.as_os_str(),
        os("-p"),
        os("2"),
    ];
    if let Some(mode) = unspliced {
        args.extend([
            os("--unspliced"),
            os(mode),
            os("--genome"),
            sim.genome.as_os_str(),
            os("--gtf"),
            sim.gtf.as_os_str(),
            os("--readLength"),
            os("75"),
            os("--unsplicedLayout"),
            os(layout),
        ]);
    }
    run(&args);
    if unspliced.is_some() {
        let info = std::fs::read_to_string(idx.join("info.json")).unwrap();
        assert!(
            info.contains(&format!("\"layout\": \"{layout}\"")),
            "{info}"
        );
    }
    let out = dir.join("quant");
    run(&[
        os("quant"),
        os("-i"),
        idx.as_os_str(),
        os("-l"),
        os("A"),
        os("-1"),
        sim.r1.as_os_str(),
        os("-2"),
        sim.r2.as_os_str(),
        os("-p"),
        os("4"),
        os("-o"),
        out.as_os_str(),
    ]);
    std::fs::read_to_string(out.join("quant.sf"))
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[4].parse().unwrap())
        })
        .collect()
}

struct Metrics {
    retained_excess: f64,
    mean_rel_error: f64,
    pearson: f64,
    unspliced_targets: f64,
}

fn metrics(sim: &Sim, est: &HashMap<String, f64>) -> Metrics {
    let names: Vec<&String> = sim.truth.keys().collect();
    let t: Vec<f64> = names.iter().map(|n| sim.truth[*n] as f64).collect();
    let e: Vec<f64> = names.iter().map(|n| est[*n]).collect();
    let mean_rel_error = t
        .iter()
        .zip(&e)
        .map(|(t, e)| (e - t).abs() / t)
        .sum::<f64>()
        / t.len() as f64;
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    let (mt, me) = (mean(&t), mean(&e));
    let cov: f64 = t.iter().zip(&e).map(|(a, b)| (a - mt) * (b - me)).sum();
    let var = |v: &[f64], m: f64| v.iter().map(|x| (x - m).powi(2)).sum::<f64>();
    let pearson = cov / (var(&t, mt) * var(&e, me)).sqrt();
    let retained_excess = sim
        .retained
        .iter()
        .map(|n| est[n] - sim.truth[n] as f64)
        .sum();
    let unspliced_targets = est
        .iter()
        .filter(|(n, _)| n.contains("-I"))
        .map(|(_, v)| v)
        .sum();
    Metrics {
        retained_excess,
        mean_rel_error,
        pearson,
        unspliced_targets,
    }
}

#[test]
fn intron_targets_stop_retained_intron_isoforms_absorbing_pre_mrna() {
    let tmp = tempfile::tempdir().unwrap();
    let sim = simulate(tmp.path());
    let dirs =
        ["gentrome", "intron", "premrna", "intron_seq", "premrna_seq"].map(|d| tmp.path().join(d));
    for d in &dirs {
        std::fs::create_dir_all(d).unwrap();
    }
    // With the genome as decoys the default layout is projection.
    let a = metrics(&sim, &index_and_quant(&sim, &dirs[0], None, "auto"));
    let b_est = index_and_quant(&sim, &dirs[1], Some("intron"), "projection");
    let c_est = index_and_quant(&sim, &dirs[2], Some("premrna"), "projection");
    let (b, c) = (metrics(&sim, &b_est), metrics(&sim, &c_est));
    // The sequence layout indexes the same targets' k-mers instead: the
    // projection is meant to find exactly the alignments it finds.
    for (est, dir, mode) in [(&b_est, &dirs[3], "intron"), (&c_est, &dirs[4], "premrna")] {
        let seq = index_and_quant(&sim, dir, Some(mode), "sequence");
        let total: f64 = seq.values().sum();
        let moved: f64 = seq.iter().map(|(n, v)| (v - est[n]).abs()).sum();
        eprintln!(
            "{mode}: projection vs sequence layout, {:.4}% of fragment mass moved",
            100.0 * moved / total
        );
        assert!(
            moved < 0.005 * total,
            "{mode}: {moved} of {total} moved between layouts"
        );
    }
    let ri_truth: usize = sim.retained.iter().map(|n| sim.truth[n]).sum();
    eprintln!(
        "simulated: {} spliced + {} pre-mRNA fragments; retained-intron isoforms carry {} true fragments",
        sim.truth.values().sum::<usize>(),
        sim.pre_mrna_fragments,
        ri_truth
    );
    eprintln!("metric                          gentrome     +intron    +premrna");
    eprintln!(
        "retained-intron excess (frags)  {:>8.1}  {:>10.1}  {:>10.1}",
        a.retained_excess, b.retained_excess, c.retained_excess
    );
    eprintln!(
        "mean |rel. error|, isoforms     {:>8.3}  {:>10.3}  {:>10.3}",
        a.mean_rel_error, b.mean_rel_error, c.mean_rel_error
    );
    eprintln!(
        "Pearson r, isoforms vs truth    {:>8.4}  {:>10.4}  {:>10.4}",
        a.pearson, b.pearson, c.pearson
    );
    eprintln!(
        "fragments on unspliced targets  {:>8.1}  {:>10.1}  {:>10.1}",
        a.unspliced_targets, b.unspliced_targets, c.unspliced_targets
    );

    // The plain gentrome hands the retained-intron isoforms a large excess of
    // intronic fragments; the intron targets remove most of it. What remains
    // is mostly pre-mRNA fragments that lie wholly in exons: they are
    // indistinguishable from mature ones in any index.
    assert!(
        a.retained_excess > 0.5 * ri_truth as f64,
        "the simulation should reproduce the problem: excess {}",
        a.retained_excess
    );
    assert!(
        b.retained_excess < 0.4 * a.retained_excess,
        "excess {} vs {}",
        b.retained_excess,
        a.retained_excess
    );
    assert!(b.mean_rel_error < a.mean_rel_error);
    assert!(b.pearson > a.pearson);
    assert!(b.unspliced_targets > 0.5 * sim.pre_mrna_fragments as f64);
    // Gene bodies absorb the intronic signal too.
    assert!(c.retained_excess < 0.4 * a.retained_excess);
    assert!(c.pearson > a.pearson);
}
