//! Transcript→gene mapping and gene-level aggregation (`-g/--geneMap` →
//! `quant.genes.sf`), matching salmon's `aggregateEstimatesToGeneLevel`.
//!
//! # Why aggregate at all
//!
//! salmon quantifies *transcripts* (isoforms), but many analyses want one number
//! per *gene*. A gene typically has several isoforms, so gene-level output is a
//! roll-up of its transcripts.
//!
//! Counts and TPM are additive, so they simply sum. Lengths are not: a gene has
//! no single length, and averaging its isoforms equally would be wrong when one
//! isoform dominates expression. So gene Length and EffectiveLength are the
//! TPM-weighted means of the transcript (effective) lengths — the length of the
//! "average molecule actually present" — falling back to the unweighted mean when
//! the gene's total TPM is 0 and there is nothing to weight by.
//!
//! ("TPM" is transcripts per million: a within-sample relative abundance that
//! already accounts for length, so all TPMs in a sample sum to one million.)

use std::collections::{BTreeMap, HashMap};
use std::io::{self, BufRead, Write};
use std::path::Path;

/// Parse a transcript→gene map. A `.gtf`/`.gff`/`.gff3` file is parsed for the
/// `transcript_id` and `gene_id` attributes (GTF `key "value"` and GFF3
/// `key=value` syntaxes); anything else is read as a TSV whose first two
/// whitespace-separated columns are `transcript` and `gene`.
///
/// GTF/GFF are the standard annotation formats: tab-separated, one feature per
/// line, with a free-form attribute string in the ninth column. The simple
/// two-column TSV is the escape hatch for annotations salmon cannot parse.
pub fn read_transcript_gene_map(path: &Path) -> io::Result<HashMap<String, String>> {
    // Format is chosen by file extension; there is no reliable content sniff for
    // "GTF vs arbitrary TSV" since both are tab-separated text.
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
        // Blank lines and `#` comment/pragma lines appear in real annotations.
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if is_gtf {
            let cols: Vec<&str> = line.split('\t').collect();
            // A well-formed feature line has 9 columns; skip anything shorter
            // rather than failing the whole file.
            if cols.len() < 9 {
                continue;
            }
            if let (Some(t), Some(g)) = (
                extract_attr(cols[8], "transcript_id"),
                extract_attr(cols[8], "gene_id"),
            ) {
                // `or_insert` keeps the first mapping seen: an annotation
                // repeats a transcript on every one of its exon lines, and they
                // all name the same gene.
                map.entry(t).or_insert(g);
            }
        } else {
            let mut it = trimmed.split_whitespace();
            if let (Some(t), Some(g)) = (it.next(), it.next()) {
                map.insert(t.to_string(), g.to_string());
            }
        }
    }

    // An empty map almost always means the wrong file or the wrong format was
    // passed, and would otherwise surface much later as an empty gene table.
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
///
/// Both syntaxes are handled by one parser because real files mix conventions,
/// and because the caller only knows the file's extension, not its dialect.
fn extract_attr(attrs: &str, key: &str) -> Option<String> {
    for entry in attrs.split(';') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        // Split into key and value at the first `=` (GFF3) or the first space
        // (GTF). Anything with neither separator is not an attribute.
        let (k, v) = if let Some(eq) = entry.find('=') {
            (entry[..eq].trim(), entry[eq + 1..].trim())
        } else if let Some(sp) = entry.find(char::is_whitespace) {
            (entry[..sp].trim(), entry[sp + 1..].trim())
        } else {
            continue;
        };
        // Exact key comparison, never a substring test: `havana_gene_id` must
        // not satisfy a request for `gene_id` (see the test below).
        if k == key {
            let v = v.trim_matches('"');
            if !v.is_empty() {
                return Some(v.to_string());
            }
        }
    }
    None
}

/// Aggregate transcript-level estimates to gene level and write `quant.genes.sf`.
/// Transcripts absent from `gene_map` are skipped (salmon's behavior); the count
/// of skipped transcripts is returned for the caller to report.
///
/// All the input slices are parallel arrays: index `i` describes the same
/// transcript in every one of them.
#[allow(clippy::too_many_arguments)]
pub fn write_gene_quant(
    out_path: &Path,
    names: &[String],
    lengths: &[u32],
    eff_lengths: &[f64],
    tpm: &[f64],
    counts: &[f64],
    gene_map: &HashMap<String, String>,
) -> io::Result<usize> {
    // Group transcript indices by gene (gene name order is sorted for determinism).
    //
    // A `BTreeMap` iterates in sorted key order, unlike a `HashMap`; that is the
    // whole reason it is used here, so the output file is byte-identical between
    // runs.
    let mut genes: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    let mut unmapped = 0usize;
    for (i, name) in names.iter().enumerate() {
        match gene_map.get(name) {
            Some(g) => genes.entry(g.as_str()).or_default().push(i),
            None => unmapped += 1,
        }
    }

    // smallest positive double, matching salmon's `minTPM = denorm_min`.
    // Used as the "is the total effectively zero" threshold below.
    let min_tpm = f64::from_bits(1);
    let mut w = io::BufWriter::new(std::fs::File::create(out_path)?);
    writeln!(w, "Name\tLength\tEffectiveLength\tTPM\tNumReads")?;
    for (gene, idxs) in &genes {
        // Additive quantities: straight sums over the gene's transcripts.
        let total_tpm: f64 = idxs.iter().map(|&i| tpm[i]).sum();
        let total_reads: f64 = idxs.iter().map(|&i| counts[i]).sum();
        let (mut g_len, mut g_eff) = (0.0f64, 0.0f64);
        if total_tpm > min_tpm {
            // Weight each isoform's length by its share of the gene's expression.
            for &i in idxs {
                let frac = tpm[i] / total_tpm;
                g_len += lengths[i] as f64 * frac;
                g_eff += eff_lengths[i] * frac;
            }
        } else {
            // Nothing expressed: weighting by zero would give 0/0, so fall back
            // to a plain average over the isoforms.
            let frac = 1.0 / idxs.len() as f64;
            for &i in idxs {
                g_len += lengths[i] as f64 * frac;
                g_eff += eff_lengths[i] * frac;
            }
        }
        // Fixed decimal places matching salmon's output precision exactly, so
        // downstream tools that diff the two see no difference.
        writeln!(
            w,
            "{gene}\t{g_len:.3}\t{g_eff:.3}\t{total_tpm:.6}\t{total_reads:.3}"
        )?;
    }
    w.flush()?;
    Ok(unmapped)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Both annotation dialects must parse, and — the subtle one — a key that
    /// merely *ends with* the requested key must not match.
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

    /// Counts and TPM sum; length is TPM-weighted, not a plain average. The
    /// numbers are chosen so the weighted answer (125) differs clearly from the
    /// unweighted one (150).
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
        let unmapped = write_gene_quant(&p, &names, &lengths, &eff, &tpm, &counts, &gm).unwrap();
        assert_eq!(unmapped, 0);
        let body = std::fs::read_to_string(&p).unwrap();
        let g_line = body.lines().find(|l| l.starts_with("G\t")).unwrap();
        let f: Vec<&str> = g_line.split('\t').collect();
        // TPM-weighted length = 100*0.75 + 200*0.25 = 125; reads = 400; tpm = 40
        assert_eq!(f[1], "125.000");
        assert_eq!(f[3], "40.000000");
        assert_eq!(f[4], "400.000");
    }
}
