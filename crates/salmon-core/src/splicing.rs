//! Splicing status of the targets of an index built with
//! `salmon index --unspliced`, and the gene-level spliced / unspliced /
//! ambiguous summary written next to `quant.sf`.
//!
//! The index records every quantified target's gene and status in
//! `t2g_3col.tsv` (`target<TAB>gene<TAB>S|U`, pyroe's layout). This module
//! parses that table ([`SpliceAnnotation`]), aligns it to a run's reference
//! numbering ([`SpliceTable`]), and holds and writes the per-gene summary
//! ([`SplicingSummary`], `quant.genes.usa.tsv`, [`UnsplicedMeta`]). The
//! summary itself is computed from the equivalence classes by
//! `salmon_infer::splicing_summary`.

use std::collections::HashMap;
use std::io::{self, Write};
use std::path::Path;

use serde::Serialize;

/// The parsed `t2g_3col.tsv` of an unspliced index.
#[derive(Debug, Clone)]
pub struct SpliceAnnotation {
    /// `intron` or `premrna` (from the index's `info.json`)
    pub mode: String,
    /// target -> (gene, is unspliced)
    rows: HashMap<String, (String, bool)>,
}

fn invalid(msg: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

impl SpliceAnnotation {
    /// Parse `t2g_3col.tsv` text; `origin` names the file in errors.
    pub fn parse(mode: &str, text: &str, origin: &str) -> io::Result<Self> {
        let mut rows = HashMap::new();
        for (i, line) in text.lines().enumerate() {
            let mut f = line.split('\t');
            match (f.next(), f.next(), f.next()) {
                (Some(t), Some(g), Some(st @ ("S" | "U"))) => {
                    rows.insert(t.to_string(), (g.to_string(), st == "U"));
                }
                _ => {
                    return Err(invalid(format!(
                        "{origin}:{}: malformed line {line:?} (expected target, gene, S|U)",
                        i + 1
                    )))
                }
            }
        }
        Ok(Self {
            mode: mode.to_string(),
            rows,
        })
    }

    /// Align the annotation to a run's references: `names[i]` is reference
    /// `i`, and `quantified` lists the reference ids that are quantified (not
    /// decoys). Every quantified reference must be in the annotation.
    pub fn table<S: AsRef<str>>(
        &self,
        names: &[S],
        quantified: impl Iterator<Item = usize>,
    ) -> io::Result<SpliceTable> {
        let n = names.len();
        let mut genes: Vec<String> = Vec::new();
        let mut gene_ix: HashMap<&str, u32> = HashMap::new();
        let mut gene_of = vec![SpliceTable::NO_GENE; n];
        let mut unspliced = vec![false; n];
        for t in quantified {
            let name = names[t].as_ref();
            let (g, is_u) = self.rows.get(name).ok_or_else(|| {
                invalid(format!(
                    "target {name:?} is missing from the index's t2g_3col.tsv"
                ))
            })?;
            let gi = *gene_ix.entry(g.as_str()).or_insert_with(|| {
                genes.push(g.clone());
                (genes.len() - 1) as u32
            });
            gene_of[t] = gi;
            unspliced[t] = *is_u;
        }
        Ok(SpliceTable {
            mode: self.mode.clone(),
            genes,
            gene_of,
            unspliced,
        })
    }
}

/// Gene and splicing status of every reference of a run.
#[derive(Debug, Clone)]
pub struct SpliceTable {
    /// `intron` or `premrna`
    pub mode: String,
    /// gene names, in reference order of first appearance
    pub genes: Vec<String>,
    /// per reference id: index into `genes`, or [`SpliceTable::NO_GENE`] for a
    /// reference that is not quantified (a decoy)
    pub gene_of: Vec<u32>,
    /// per reference id: whether it is an unspliced target
    pub unspliced: Vec<bool>,
}

impl SpliceTable {
    /// `gene_of` value of references with no gene (decoys).
    pub const NO_GENE: u32 = u32::MAX;

    /// Number of (spliced, unspliced) quantified targets.
    pub fn target_counts(&self) -> (usize, usize) {
        let mut n = (0, 0);
        for (g, &u) in self.gene_of.iter().zip(&self.unspliced) {
            if *g != Self::NO_GENE {
                if u {
                    n.1 += 1;
                } else {
                    n.0 += 1;
                }
            }
        }
        n
    }
}

/// Per-gene spliced / unspliced / ambiguous fragment counts.
#[derive(Debug, Clone)]
pub struct SplicingSummary {
    /// index build mode (`intron` or `premrna`)
    pub mode: String,
    pub genes: Vec<String>,
    pub spliced: Vec<f64>,
    pub unspliced: Vec<f64>,
    pub ambiguous: Vec<f64>,
    pub num_spliced_targets: usize,
    pub num_unspliced_targets: usize,
    /// fragments dropped because their best alignment was to a decoy; they
    /// belong to no gene, but are the fourth category of the run-level split
    /// (the genome-derived share an unspliced index did not claim)
    pub decoy_fragments: f64,
}

impl SplicingSummary {
    /// Total (spliced, unspliced, ambiguous) fragment mass over all genes.
    pub fn totals(&self) -> (f64, f64, f64) {
        (
            self.spliced.iter().sum(),
            self.unspliced.iter().sum(),
            self.ambiguous.iter().sum(),
        )
    }
}

/// File name of the per-gene summary in the output directory.
pub const USA_FILE: &str = "quant.genes.usa.tsv";

/// Write `quant.genes.usa.tsv`: per gene, the fragments filed as spliced,
/// unspliced and ambiguous, with `quant.sf`'s NumReads precision.
pub fn write_usa_tsv(path: &Path, sp: &SplicingSummary, sig_digits: usize) -> io::Result<()> {
    let mut w = io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(w, "Name\tSpliced\tUnspliced\tAmbiguous")?;
    for (i, g) in sp.genes.iter().enumerate() {
        writeln!(
            w,
            "{g}\t{:.*}\t{:.*}\t{:.*}",
            sig_digits, sp.spliced[i], sig_digits, sp.unspliced[i], sig_digits, sp.ambiguous[i]
        )?;
    }
    w.flush()
}

/// Run-level split of assigned *and* decoy fragments, as fractions of their
/// sum.
#[derive(Debug, Clone, Serialize)]
pub struct FractionsWithDecoys {
    pub spliced: f64,
    pub unspliced: f64,
    pub ambiguous: f64,
    pub decoy: f64,
}

/// The `unspliced` block of `meta_info.json`.
#[derive(Debug, Clone, Serialize)]
pub struct UnsplicedMeta {
    pub mode: String,
    pub num_spliced_targets: usize,
    pub num_unspliced_targets: usize,
    pub num_genes: usize,
    /// assigned fragments filed as spliced / unspliced / ambiguous
    pub num_spliced_fragments: f64,
    pub num_unspliced_fragments: f64,
    pub num_ambiguous_fragments: f64,
    /// fragments dropped as decoy-dominated (no gene)
    pub num_decoy_fragments: f64,
    /// spliced / unspliced / ambiguous as fractions of the assigned fragments
    pub spliced_fraction: f64,
    pub unspliced_fraction: f64,
    pub ambiguous_fraction: f64,
    /// the four categories as fractions of assigned + decoy fragments, so
    /// libraries compare whatever share their decoys took
    pub fractions_with_decoys: FractionsWithDecoys,
}

impl From<&SplicingSummary> for UnsplicedMeta {
    fn from(sp: &SplicingSummary) -> Self {
        let (s, u, a) = sp.totals();
        let total = s + u + a;
        let frac = |x: f64| if total > 0.0 { x / total } else { 0.0 };
        let d = sp.decoy_fragments;
        let all = total + d;
        let frac_all = |x: f64| if all > 0.0 { x / all } else { 0.0 };
        Self {
            mode: sp.mode.clone(),
            num_spliced_targets: sp.num_spliced_targets,
            num_unspliced_targets: sp.num_unspliced_targets,
            num_genes: sp.genes.len(),
            num_spliced_fragments: s,
            num_unspliced_fragments: u,
            num_ambiguous_fragments: a,
            num_decoy_fragments: d,
            spliced_fraction: frac(s),
            unspliced_fraction: frac(u),
            ambiguous_fraction: frac(a),
            fractions_with_decoys: FractionsWithDecoys {
                spliced: frac_all(s),
                unspliced: frac_all(u),
                ambiguous: frac_all(a),
                decoy: frac_all(d),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_follows_the_run_numbering_and_skips_decoys() {
        let a = SpliceAnnotation::parse("intron", "T2\tG2\tS\nT1\tG1\tS\nG1-I\tG1\tU\n", "t2g")
            .unwrap();
        let names = ["T1", "T2", "G1-I", "chr1"];
        let t = a.table(&names, 0..3).unwrap();
        assert_eq!(t.genes, ["G1", "G2"]);
        assert_eq!(t.gene_of, [0, 1, 0, SpliceTable::NO_GENE]);
        assert_eq!(t.unspliced, [false, false, true, false]);
        assert_eq!(t.target_counts(), (2, 1));
        // a quantified target the table does not know
        assert!(a.table(&names, 0..4).is_err());
    }

    #[test]
    fn malformed_rows_are_rejected() {
        assert!(SpliceAnnotation::parse("intron", "T1\tG1\tX\n", "t2g").is_err());
        assert!(SpliceAnnotation::parse("intron", "T1\tG1\n", "t2g").is_err());
    }

    #[test]
    fn meta_fractions_sum_to_one() {
        let sp = SplicingSummary {
            mode: "intron".into(),
            genes: vec!["A".into(), "B".into()],
            spliced: vec![3.0, 1.0],
            unspliced: vec![4.0, 0.0],
            ambiguous: vec![2.0, 0.0],
            num_spliced_targets: 2,
            num_unspliced_targets: 1,
            decoy_fragments: 10.0,
        };
        let m = UnsplicedMeta::from(&sp);
        let w = &m.fractions_with_decoys;
        assert!((w.spliced + w.unspliced + w.ambiguous + w.decoy - 1.0).abs() < 1e-12);
        assert!((w.decoy - 0.5).abs() < 1e-12);
        assert_eq!(m.num_unspliced_fragments, 4.0);
        assert!(
            (m.spliced_fraction + m.unspliced_fraction + m.ambiguous_fraction - 1.0).abs() < 1e-12
        );
        assert!((m.unspliced_fraction - 0.4).abs() < 1e-12);
    }
}
