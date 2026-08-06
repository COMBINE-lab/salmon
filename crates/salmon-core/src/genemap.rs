//! Transcript→gene mapping and gene-level aggregation (`-g/--geneMap` →
//! `quant.genes.sf`), matching salmon's `aggregateEstimatesToGeneLevel`.
//!
//! Gene NumReads and TPM are the sums over the gene's transcripts; gene Length
//! and EffectiveLength are the TPM-weighted means of the transcript (effective)
//! lengths, falling back to the unweighted mean when the gene's total TPM is 0.

use flate2::read::MultiGzDecoder;
use std::collections::{BTreeMap, HashMap};
use std::io::{self, BufRead, Write};
use std::path::Path;

/// Parse a transcript→gene map. A `.gtf`/`.gff`/`.gff3` file is parsed for the
/// `transcript_id` and `gene_id` attributes (GTF `key "value"` and GFF3
/// `key=value` syntaxes); anything else is read as a TSV whose first two
/// whitespace-separated columns are `transcript` and `gene`.
///
/// The file may be gzip-compressed — as Ensembl and GENCODE ship it. Compression
/// is detected from the content, not the name, so a compressed file that kept a
/// bare `.gtf` name is read correctly too.
pub fn read_transcript_gene_map(path: &Path) -> io::Result<HashMap<String, String>> {
    let is_gtf = has_gtf_extension(path);
    let reader = open_maybe_gzip(path)?;
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

/// Whether the path names a GTF/GFF, ignoring a trailing compression suffix:
/// the last extension of `annotation.gtf.gz` is `gz`, but the content is GTF.
/// Getting this wrong is not a cosmetic error — a GTF handed to the TSV branch
/// parses into `{sequence name: source}` pairs, which is a non-empty map that
/// matches no transcript.
fn has_gtf_extension(path: &Path) -> bool {
    let ext = |p: &Path| {
        p.extension()
            .and_then(|e| e.to_str())
            .map(str::to_ascii_lowercase)
    };
    let outer = ext(path);
    let effective = match outer.as_deref() {
        Some("gz") => ext(Path::new(path.file_stem().unwrap_or_default())),
        _ => outer,
    };
    matches!(effective.as_deref(), Some("gtf" | "gff" | "gff3"))
}

/// Open `path`, transparently decompressing it when it is gzipped.
///
/// Detection is by the `1f 8b` magic bytes rather than the file name, so a
/// compressed annotation that was renamed without a `.gz` suffix still works.
/// `MultiGzDecoder` handles multi-member files — concatenating gzip members is
/// legal and some published annotations are built that way; a single-member
/// decoder would silently stop at the first member's end.
fn open_maybe_gzip(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut reader = io::BufReader::new(std::fs::File::open(path)?);
    // `fill_buf` peeks without consuming, so the magic bytes stay in the stream
    // for whichever reader we hand back.
    let gzipped = matches!(reader.fill_buf()?, [0x1f, 0x8b, ..]);
    if gzipped {
        Ok(Box::new(io::BufReader::new(MultiGzDecoder::new(reader))))
    } else {
        Ok(Box::new(reader))
    }
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

    /// Gzip is how Ensembl and GENCODE ship annotations. This also pins the
    /// GTF-vs-TSV decision: the last extension here is `gz`, and routing the
    /// decompressed GTF to the TSV branch would yield `{"1": "havana"}` instead.
    #[test]
    fn reads_gzipped_gtf() {
        let p = write_fixture("gzip", "annotation.gtf.gz", GTF, true);
        assert_expected_map(&read_transcript_gene_map(&p).unwrap());
    }

    #[test]
    fn reads_gzipped_tsv() {
        let p = write_fixture("gziptsv", "t2g.tsv.gz", TSV, true);
        assert_expected_map(&read_transcript_gene_map(&p).unwrap());
    }

    /// Compression is sniffed from the magic bytes, so a gzipped file that kept
    /// a bare `.gtf` name is still read as a compressed GTF.
    #[test]
    fn reads_gzip_content_under_a_plain_name() {
        let p = write_fixture("misnamed", "annotation.gtf", GTF, true);
        assert_expected_map(&read_transcript_gene_map(&p).unwrap());
    }

    /// A multi-member gzip stream (concatenated members) must be read in full,
    /// not truncated at the first member's end.
    #[test]
    fn reads_multi_member_gzip() {
        use flate2::{write::GzEncoder, Compression};
        let path = scratch("multimember").join("annotation.gtf.gz");
        let mut f = std::fs::File::create(&path).unwrap();
        // split the fixture so each half becomes its own gzip member
        let (head, tail) = GTF.split_at(GTF.find("1\thavana\ttranscript\t12010").unwrap());
        for part in [head, tail] {
            let mut enc = GzEncoder::new(Vec::new(), Compression::new(6));
            enc.write_all(part.as_bytes()).unwrap();
            f.write_all(&enc.finish().unwrap()).unwrap();
        }
        f.flush().unwrap();
        assert_expected_map(&read_transcript_gene_map(&path).unwrap());
    }

    #[test]
    fn gtf_extension_detection_ignores_compression_suffix() {
        for name in ["a.gtf", "a.gtf.gz", "a.GTF.GZ", "a.gff3", "a.gff.gz"] {
            assert!(has_gtf_extension(Path::new(name)), "{name}");
        }
        for name in ["t2g.tsv", "t2g.tsv.gz", "t2g", "t2g.txt.gz"] {
            assert!(!has_gtf_extension(Path::new(name)), "{name}");
        }
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
