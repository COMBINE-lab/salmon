//! Transcript→gene mapping and gene-level aggregation (`-g/--geneMap` →
//! `quant.genes.sf`), matching salmon's `aggregateEstimatesToGeneLevel`.
//!
//! Gene NumReads and TPM are the sums over the gene's transcripts; gene Length
//! and EffectiveLength are the TPM-weighted means of the transcript (effective)
//! lengths, falling back to the unweighted mean when the gene's total TPM is 0.

use std::collections::{BTreeMap, HashMap};
use std::io::{self, BufRead, Write};
use std::path::Path;

/// Parse a transcript→gene map. A `.gtf`/`.gff`/`.gff3` file is parsed for the
/// `transcript_id` and `gene_id` attributes (GTF `key "value"` and GFF3
/// `key=value` syntaxes); anything else is read as a TSV whose first two
/// whitespace-separated columns are `transcript` and `gene`.
pub fn read_transcript_gene_map(path: &Path) -> io::Result<HashMap<String, String>> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let is_gtf = matches!(ext.as_str(), "gtf" | "gff" | "gff3");
    let reader = io::BufReader::new(std::fs::File::open(path)?);
    let mut map = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if is_gtf {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 9 {
                continue;
            }
            if let (Some(t), Some(g)) = (
                extract_attr(cols[8], "transcript_id"),
                extract_attr(cols[8], "gene_id"),
            ) {
                map.entry(t).or_insert(g);
            }
        } else {
            let mut it = trimmed.split_whitespace();
            if let (Some(t), Some(g)) = (it.next(), it.next()) {
                map.insert(t.to_string(), g.to_string());
            }
        }
    }

    if map.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "no transcript->gene mappings parsed from {}",
                path.display()
            ),
        ));
    }
    Ok(map)
}

/// Extract attribute `key` from a GTF/GFF column-9 string. Each `;`-separated
/// entry is `key "value"` (GTF) or `key=value` (GFF3); the value's surrounding
/// quotes are stripped.
fn extract_attr(attrs: &str, key: &str) -> Option<String> {
    for entry in attrs.split(';') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        let (k, v) = if let Some(eq) = entry.find('=') {
            (entry[..eq].trim(), entry[eq + 1..].trim())
        } else if let Some(sp) = entry.find(char::is_whitespace) {
            (entry[..sp].trim(), entry[sp + 1..].trim())
        } else {
            continue;
        };
        if k == key {
            let v = v.trim_matches('"');
            if !v.is_empty() {
                return Some(v.to_string());
            }
        }
    }
    None
}

/// A match rate at or below this is reported as a problem rather than a note.
///
/// The failure this guards against is total: a versioned/unversioned identifier
/// mismatch matches *nothing*, and pairing an annotation with the FASTA it was
/// built from matches essentially everything, so the two regimes are far apart
/// and the exact cut is not delicate. Half is chosen because below it the
/// majority of what was quantified is missing from `quant.genes.sf`, which is
/// worth a look in any workflow, while deliberately partial annotations — a
/// single-chromosome GTF, a spike-in-only map — stay well clear of being
/// mistaken for correct. It gates a warning only; nothing fails.
pub const LOW_MATCH_RATE: f64 = 0.5;

/// What gene-level aggregation actually did, for the caller to report.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct GeneQuantSummary {
    /// Quantified transcripts that found an entry in the gene map.
    pub matched: usize,
    /// Quantified transcripts with no entry; omitted from `quant.genes.sf`.
    pub unmapped: usize,
    /// Gene rows actually written — not the number of genes in the gene map.
    pub genes_written: usize,
}

impl GeneQuantSummary {
    /// Quantified transcripts considered.
    pub fn total(&self) -> usize {
        self.matched + self.unmapped
    }

    /// Fraction of quantified transcripts that found a gene. An empty transcript
    /// set is vacuously fully matched, so it does not trip the low-rate warning.
    pub fn match_rate(&self) -> f64 {
        if self.total() == 0 {
            1.0
        } else {
            self.matched as f64 / self.total() as f64
        }
    }

    /// Whether the caller should warn about how little the gene map matched.
    pub fn match_rate_is_low(&self) -> bool {
        self.match_rate() <= LOW_MATCH_RATE
    }
}

/// Strip a trailing `.<digits>` version suffix, as tximport's `ignoreTxVersion`
/// does. Identifiers whose last dot-separated part is not purely numeric are
/// returned unchanged, which leaves GENCODE's `…_PAR_Y` suffixes alone.
pub fn strip_tx_version(id: &str) -> &str {
    match id.rsplit_once('.') {
        Some((base, ver)) if !ver.is_empty() && ver.bytes().all(|b| b.is_ascii_digit()) => base,
        _ => id,
    }
}

/// How many of `names` would match if a version suffix were ignored on both
/// sides. Purely diagnostic: it explains a poor match rate and never changes
/// what is written.
pub fn matches_ignoring_tx_version(names: &[String], gene_map: &HashMap<String, String>) -> usize {
    let stripped: HashMap<&str, &str> = gene_map
        .iter()
        .map(|(t, g)| (strip_tx_version(t), g.as_str()))
        .collect();
    names
        .iter()
        .filter(|n| stripped.contains_key(strip_tx_version(n)))
        .count()
}

/// Aggregate transcript-level estimates to gene level and write `quant.genes.sf`.
///
/// Transcripts absent from `gene_map` are omitted from the output. Note this
/// differs from C++ salmon, whose `aggregateEstimatesToGeneLevel` emitted an
/// unmatched transcript as its own single-transcript gene and warned once per
/// transcript; see #1075 for the discussion of which behaviour to keep.
///
/// The returned [`GeneQuantSummary`] is what the caller should report — in
/// particular `genes_written`, which is the number of rows in the file and not
/// the number of genes in the gene map.
#[allow(clippy::too_many_arguments)]
pub fn write_gene_quant(
    out_path: &Path,
    names: &[String],
    lengths: &[u32],
    eff_lengths: &[f64],
    tpm: &[f64],
    counts: &[f64],
    gene_map: &HashMap<String, String>,
) -> io::Result<GeneQuantSummary> {
    // Group transcript indices by gene (gene name order is sorted for determinism).
    let mut genes: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    let mut summary = GeneQuantSummary::default();
    for (i, name) in names.iter().enumerate() {
        match gene_map.get(name) {
            Some(g) => {
                genes.entry(g.as_str()).or_default().push(i);
                summary.matched += 1;
            }
            None => summary.unmapped += 1,
        }
    }

    // smallest positive double, matching salmon's `minTPM = denorm_min`.
    let min_tpm = f64::from_bits(1);
    let mut w = io::BufWriter::new(std::fs::File::create(out_path)?);
    writeln!(w, "Name\tLength\tEffectiveLength\tTPM\tNumReads")?;
    for (gene, idxs) in &genes {
        let total_tpm: f64 = idxs.iter().map(|&i| tpm[i]).sum();
        let total_reads: f64 = idxs.iter().map(|&i| counts[i]).sum();
        let (mut g_len, mut g_eff) = (0.0f64, 0.0f64);
        if total_tpm > min_tpm {
            for &i in idxs {
                let frac = tpm[i] / total_tpm;
                g_len += lengths[i] as f64 * frac;
                g_eff += eff_lengths[i] * frac;
            }
        } else {
            let frac = 1.0 / idxs.len() as f64;
            for &i in idxs {
                g_len += lengths[i] as f64 * frac;
                g_eff += eff_lengths[i] * frac;
            }
        }
        writeln!(
            w,
            "{gene}\t{g_len:.3}\t{g_eff:.3}\t{total_tpm:.6}\t{total_reads:.3}"
        )?;
    }
    w.flush()?;
    summary.genes_written = genes.len();
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(case: &str) -> std::path::PathBuf {
        let dir =
            std::env::temp_dir().join(format!("salmon_genemap_{}_{case}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// The reference human bulk RNA-seq setup: an index built from the Ensembl
    /// cDNA FASTA, whose headers carry the version inline, against an Ensembl
    /// GTF, which splits it into a separate `transcript_version` attribute.
    /// The join is a literal string compare, so nothing matches at all.
    #[test]
    fn versioned_quant_names_match_no_unversioned_gene_map_entry() {
        let gtf = "\
1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000290825\"; gene_version \"2\"; transcript_id \"ENST00000456328\"; transcript_version \"2\";
1\thavana\ttranscript\t12010\t13670\t.\t+\t.\tgene_id \"ENSG00000223972\"; gene_version \"6\"; transcript_id \"ENST00000450305\"; transcript_version \"3\";
";
        let gtf_path = scratch("nomatch").join("annotation.gtf");
        std::fs::write(&gtf_path, gtf).unwrap();
        let map = read_transcript_gene_map(&gtf_path).unwrap();
        assert_eq!(map.len(), 2);
        assert!(
            map.contains_key("ENST00000456328"),
            "GTF keys are unversioned"
        );

        // ...but quant.sf names them as the cDNA FASTA does.
        let names: Vec<String> = vec!["ENST00000456328.2".into(), "ENST00000450305.3".into()];
        let out = scratch("nomatch").join("quant.genes.sf");
        let summary = write_gene_quant(
            &out,
            &names,
            &[1000, 2000],
            &[800.0, 1800.0],
            &[10.0, 20.0],
            &[100.0, 200.0],
            &map,
        )
        .unwrap();

        assert_eq!(summary.unmapped, 2, "every transcript should fail to match");
        assert_eq!(summary.matched, 0);
        assert_eq!(summary.genes_written, 0);
        assert_eq!(summary.match_rate(), 0.0);
        assert!(summary.match_rate_is_low());
        let body = std::fs::read_to_string(&out).unwrap();
        assert_eq!(
            body.lines().count(),
            1,
            "quant.genes.sf should be header-only: {body:?}"
        );

        // The diagnosis the warning is built on: the identifiers differ only by
        // the version suffix, so ignoring it would match everything.
        assert_eq!(matches_ignoring_tx_version(&names, &map), 2);
    }

    #[test]
    fn summary_counts_genes_written_not_gene_map_size() {
        // The gene map knows three genes; only one of them is quantified.
        let mut gm = HashMap::new();
        gm.insert("t1".to_string(), "G".to_string());
        gm.insert("t2".to_string(), "H".to_string());
        gm.insert("t3".to_string(), "I".to_string());
        let names = vec!["t1".to_string(), "zzz".to_string()];
        let out = scratch("written").join("quant.genes.sf");
        let summary = write_gene_quant(
            &out,
            &names,
            &[10, 20],
            &[8.0, 18.0],
            &[1.0, 2.0],
            &[3.0, 4.0],
            &gm,
        )
        .unwrap();
        assert_eq!(summary.matched, 1);
        assert_eq!(summary.unmapped, 1);
        assert_eq!(summary.genes_written, 1, "one row written, not 3 map genes");
        assert_eq!(summary.total(), 2);
        assert_eq!(summary.match_rate(), 0.5);
        assert!(summary.match_rate_is_low(), "50% is at the threshold");
    }

    #[test]
    fn full_match_does_not_trip_the_low_rate_warning() {
        let mut gm = HashMap::new();
        gm.insert("t1".to_string(), "G".to_string());
        gm.insert("t2".to_string(), "G".to_string());
        let names = vec!["t1".to_string(), "t2".to_string()];
        let out = scratch("full").join("quant.genes.sf");
        let summary = write_gene_quant(
            &out,
            &names,
            &[10, 20],
            &[8.0, 18.0],
            &[1.0, 2.0],
            &[3.0, 4.0],
            &gm,
        )
        .unwrap();
        assert_eq!(summary.match_rate(), 1.0);
        assert!(!summary.match_rate_is_low());
        assert_eq!(summary.genes_written, 1);
        // Nothing to suggest when everything already matches.
        assert_eq!(matches_ignoring_tx_version(&names, &gm), 2);
    }

    #[test]
    fn version_stripping_leaves_non_numeric_suffixes_alone() {
        assert_eq!(strip_tx_version("ENST00000456328.2"), "ENST00000456328");
        assert_eq!(strip_tx_version("ENST00000456328.12"), "ENST00000456328");
        assert_eq!(strip_tx_version("ENST00000456328"), "ENST00000456328");
        // GENCODE PAR_Y entries must survive intact.
        assert_eq!(
            strip_tx_version("ENST00000456328.2_PAR_Y"),
            "ENST00000456328.2_PAR_Y"
        );
        // a name that is dotted but not versioned
        assert_eq!(strip_tx_version("gene.symbol"), "gene.symbol");
        assert_eq!(strip_tx_version("t."), "t.");
    }

    /// An empty transcript set must not be reported as a mismatch.
    #[test]
    fn empty_transcript_set_is_not_low() {
        let s = GeneQuantSummary::default();
        assert_eq!(s.match_rate(), 1.0);
        assert!(!s.match_rate_is_low());
    }

    #[test]
    fn gtf_and_gff_attr_extraction() {
        let gtf = r#"gene_id "ENSG1"; transcript_id "ENST1"; gene_name "A";"#;
        assert_eq!(extract_attr(gtf, "transcript_id").as_deref(), Some("ENST1"));
        assert_eq!(extract_attr(gtf, "gene_id").as_deref(), Some("ENSG1"));
        let gff = "ID=ENST1;gene_id=ENSG1;Parent=g1";
        assert_eq!(extract_attr(gff, "gene_id").as_deref(), Some("ENSG1"));
        // a key that is a substring of another must not false-match
        assert_eq!(extract_attr("havana_gene_id \"X\";", "gene_id"), None);
    }

    #[test]
    fn aggregation_sums_and_weights() {
        // two transcripts of gene G (lengths 100/200, TPM 30/10), one of gene H.
        let names = vec!["t1".into(), "t2".into(), "t3".into()];
        let lengths = vec![100u32, 200, 300];
        let eff = vec![80.0, 180.0, 280.0];
        let tpm = vec![30.0, 10.0, 5.0];
        let counts = vec![300.0, 100.0, 50.0];
        let mut gm = HashMap::new();
        gm.insert("t1".into(), "G".into());
        gm.insert("t2".into(), "G".into());
        gm.insert("t3".into(), "H".into());
        let dir = std::env::temp_dir().join(format!("gq_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("quant.genes.sf");
        let summary = write_gene_quant(&p, &names, &lengths, &eff, &tpm, &counts, &gm).unwrap();
        assert_eq!(summary.unmapped, 0);
        assert_eq!(summary.matched, 3);
        assert_eq!(summary.genes_written, 2);
        let body = std::fs::read_to_string(&p).unwrap();
        let g_line = body.lines().find(|l| l.starts_with("G\t")).unwrap();
        let f: Vec<&str> = g_line.split('\t').collect();
        // TPM-weighted length = 100*0.75 + 200*0.25 = 125; reads = 400; tpm = 40
        assert_eq!(f[1], "125.000");
        assert_eq!(f[3], "40.000000");
        assert_eq!(f[4], "400.000");
    }
}
