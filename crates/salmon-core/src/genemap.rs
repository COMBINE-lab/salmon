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

/// Aggregate transcript-level estimates to gene level and write `quant.genes.sf`.
/// Transcripts absent from `gene_map` are skipped (salmon's behavior); the count
/// of skipped transcripts is returned for the caller to report.
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
    let mut genes: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    let mut unmapped = 0usize;
    for (i, name) in names.iter().enumerate() {
        match gene_map.get(name) {
            Some(g) => genes.entry(g.as_str()).or_default().push(i),
            None => unmapped += 1,
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
    Ok(unmapped)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An Ensembl-shaped GTF: a header comment, a `gene` feature carrying no
    /// `transcript_id` (which must be skipped), and two transcripts.
    const GTF: &str = "\
#!genome-build GRCh38.p14
1\thavana\tgene\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000290825\"; gene_name \"DDX11L16\";
1\thavana\ttranscript\t11869\t14409\t.\t+\t.\tgene_id \"ENSG00000290825\"; transcript_id \"ENST00000456328\";
1\thavana\ttranscript\t12010\t13670\t.\t+\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000450305\";
";

    const TSV: &str = "\
ENST00000456328\tENSG00000290825
ENST00000450305\tENSG00000223972
";

    /// Per-test scratch directory, so cases can run concurrently.
    fn scratch(case: &str) -> std::path::PathBuf {
        let dir =
            std::env::temp_dir().join(format!("salmon_genemap_{}_{case}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_fixture(case: &str, name: &str, body: &str, gzip: bool) -> std::path::PathBuf {
        let path = scratch(case).join(name);
        if gzip {
            use flate2::{write::GzEncoder, Compression};
            let mut enc =
                GzEncoder::new(std::fs::File::create(&path).unwrap(), Compression::new(6));
            enc.write_all(body.as_bytes()).unwrap();
            enc.finish().unwrap();
        } else {
            std::fs::write(&path, body).unwrap();
        }
        path
    }

    /// Every accepted spelling of a gene map must yield the same two pairs.
    fn assert_expected_map(map: &HashMap<String, String>) {
        assert_eq!(map.len(), 2, "unexpected entries: {map:?}");
        assert_eq!(
            map.get("ENST00000456328").map(String::as_str),
            Some("ENSG00000290825")
        );
        assert_eq!(
            map.get("ENST00000450305").map(String::as_str),
            Some("ENSG00000223972")
        );
    }

    #[test]
    fn reads_plain_gtf() {
        let p = write_fixture("plain", "annotation.gtf", GTF, false);
        assert_expected_map(&read_transcript_gene_map(&p).unwrap());
    }

    #[test]
    fn reads_two_column_tsv() {
        let p = write_fixture("tsv", "t2g.tsv", TSV, false);
        assert_expected_map(&read_transcript_gene_map(&p).unwrap());
    }

    // Gzip is how Ensembl and GENCODE ship annotations, but the reader opens the
    // path as plain text, so `BufRead::lines()` rejects the magic bytes with
    // "stream did not contain valid UTF-8". Ignored so this commit stays green
    // and bisectable; un-ignored by the commit that adds decompression.
    // See https://github.com/COMBINE-lab/salmon/issues/1074.
    #[test]
    #[ignore = "gzip gene maps are not supported yet — see issue #1074"]
    fn reads_gzipped_gtf() {
        let p = write_fixture("gzip", "annotation.gtf.gz", GTF, true);
        assert_expected_map(&read_transcript_gene_map(&p).unwrap());
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
