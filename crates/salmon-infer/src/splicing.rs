//! Gene-level spliced / unspliced / ambiguous summary for an index built with
//! `salmon index --unspliced` (`quant.genes.usa.tsv` and the `unspliced`
//! block of `meta_info.json`; the types live in `salmon_core::splicing`).
//!
//! # Definition
//!
//! This is alevin-fry's USA split carried over to bulk EM counts. Each
//! equivalence class's fragments are shared between genes by the final EM
//! responsibilities (`alpha_t * w_t / sum`, the same E-step that produced
//! `quant.sf`'s NumReads). A gene's share of a class is then filed by what the
//! class says about that gene:
//!
//! * **spliced**: the class holds only spliced targets of the gene;
//! * **unspliced**: only unspliced targets;
//! * **ambiguous**: both, i.e. the fragments are compatible with the mature and
//!   the unspliced molecule alike (exonic fragments of an intron-bearing gene
//!   under `premrna`, retained-intron or flank fragments under `intron`).
//!
//! So `Spliced + Unspliced + Ambiguous` of a gene is the sum of its targets'
//! NumReads, up to the EM's convergence tolerance.

use salmon_core::splicing::{SpliceTable, SplicingSummary};

use crate::PackedEqClasses;

/// Split every class's fragments between genes by the final abundances and
/// file each gene's share as spliced, unspliced or ambiguous (see the module
/// docs). `decoy_fragments` (the run's decoy-dominated fragments) is carried
/// into the summary as the run-level fourth category.
pub fn splicing_summary(
    p: &PackedEqClasses,
    alphas: &[f64],
    table: &SpliceTable,
    decoy_fragments: u64,
) -> SplicingSummary {
    let ng = table.genes.len();
    let (mut s, mut u, mut a) = (vec![0.0; ng], vec![0.0; ng], vec![0.0; ng]);
    // (gene, responsibility mass, has spliced, has unspliced) for one class;
    // classes are small, so a linear scan beats a map.
    let mut per_gene: Vec<(u32, f64, bool, bool)> = Vec::new();
    for ci in 0..p.num_classes() {
        let (lo, hi) = (p.starts[ci] as usize, p.starts[ci + 1] as usize);
        let labels = &p.labels[lo..hi];
        let w = &p.combined[lo..hi];
        let denom: f64 = labels
            .iter()
            .zip(w)
            .map(|(&t, &w)| alphas[t as usize] * w)
            .sum();
        if denom <= 0.0 || !denom.is_finite() {
            // every member truncated: this mass is unassigned in quant.sf too
            continue;
        }
        let scale = p.counts[ci] as f64 / denom;
        per_gene.clear();
        for (&t, &w) in labels.iter().zip(w) {
            let g = table.gene_of[t as usize];
            if g == SpliceTable::NO_GENE {
                continue;
            }
            let r = alphas[t as usize] * w * scale;
            let is_u = table.unspliced[t as usize];
            match per_gene.iter_mut().find(|e| e.0 == g) {
                Some(e) => {
                    e.1 += r;
                    e.2 |= !is_u;
                    e.3 |= is_u;
                }
                None => per_gene.push((g, r, !is_u, is_u)),
            }
        }
        for &(g, r, has_s, has_u) in &per_gene {
            let g = g as usize;
            match (has_s, has_u) {
                (true, true) => a[g] += r,
                (false, true) => u[g] += r,
                _ => s[g] += r,
            }
        }
    }
    let (num_spliced_targets, num_unspliced_targets) = table.target_counts();
    SplicingSummary {
        mode: table.mode.clone(),
        genes: table.genes.clone(),
        spliced: s,
        unspliced: u,
        ambiguous: a,
        num_spliced_targets,
        num_unspliced_targets,
        decoy_fragments: decoy_fragments as f64,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use salmon_eqclass::{EquivalenceClassBuilder, TranscriptGroup};

    fn packed(classes: &[(Vec<u32>, u64)], num_txps: usize) -> PackedEqClasses {
        let b = EquivalenceClassBuilder::new();
        for (txps, count) in classes {
            b.add_group(
                TranscriptGroup::new(txps.clone()),
                vec![1.0; txps.len()],
                *count,
            );
        }
        let mut eq = b.finish();
        eq.update_eff_lengths(&vec![1.0; num_txps]);
        PackedEqClasses::from_collapsed(&eq, num_txps)
    }

    #[test]
    fn classes_are_filed_by_what_they_say_about_each_gene() {
        // gene A: t0 (S), t1 (S), t2 (U); gene B: t3 (S), t4 (U); t5 decoy
        let table = SpliceTable {
            mode: "intron".into(),
            genes: vec!["A".into(), "B".into()],
            gene_of: vec![0, 0, 0, 1, 1, SpliceTable::NO_GENE],
            unspliced: vec![false, false, true, false, true, false],
        };
        let p = packed(
            &[
                (vec![0], 10),    // A spliced
                (vec![0, 1], 6),  // A spliced (two isoforms)
                (vec![2], 4),     // A unspliced
                (vec![1, 2], 8),  // A ambiguous (retained intron)
                (vec![2, 4], 12), // A unspliced + B unspliced, split 1:3
                (vec![3], 5),     // B spliced
            ],
            6,
        );
        let alphas = [1.0, 1.0, 1.0, 1.0, 3.0, 0.0];
        let sum = splicing_summary(&p, &alphas, &table, 0);
        let close = |x: f64, y: f64| (x - y).abs() < 1e-9;
        assert!(close(sum.spliced[0], 16.0), "{:?}", sum.spliced);
        assert!(close(sum.unspliced[0], 4.0 + 3.0), "{:?}", sum.unspliced);
        assert!(close(sum.ambiguous[0], 8.0), "{:?}", sum.ambiguous);
        assert!(close(sum.spliced[1], 5.0));
        assert!(close(sum.unspliced[1], 9.0));
        assert!(close(sum.ambiguous[1], 0.0));
        let (s, u, a) = sum.totals();
        assert!(close(s + u + a, 45.0));
        assert_eq!((sum.num_spliced_targets, sum.num_unspliced_targets), (3, 2));
    }

    #[test]
    fn fully_truncated_classes_are_left_unassigned() {
        let table = SpliceTable {
            mode: "premrna".into(),
            genes: vec!["A".into()],
            gene_of: vec![0, 0],
            unspliced: vec![false, true],
        };
        let p = packed(&[(vec![0, 1], 7), (vec![1], 3)], 2);
        let sum = splicing_summary(&p, &[5.0, 0.0], &table, 0);
        // the first class holds both forms of A: ambiguous, whatever the split
        assert_eq!(sum.totals(), (0.0, 0.0, 7.0));
    }
}
