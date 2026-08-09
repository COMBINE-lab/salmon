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
/// The file may be compressed — gzip is how Ensembl and GENCODE ship it, and
/// `niffler` also covers BGZF, bzip2, xz and zstd. Compression is detected from
/// the content, not the name, so a compressed file that kept a bare `.gtf` name
/// is read correctly too.
pub fn read_transcript_gene_map(path: &Path) -> io::Result<HashMap<String, String>> {
    let is_gtf = has_gtf_extension(path);
    let reader = open_maybe_compressed(path)?;
    let mut map = HashMap::new();

    for line in reader.lines() {
        let line = line.map_err(|e| annotate_read_error(path, e))?;
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

/// Name the file in a read failure, and add a hint when the bytes are not
/// readable as text. `BufRead::lines()` reports only "stream did not contain
/// valid UTF-8", which says neither which file nor what to do about it; past
/// the format sniff, that error means the input is neither text nor any
/// compression format we recognise. The original message is kept so a
/// corrupt-gzip failure still reads as one.
fn annotate_read_error(path: &Path, err: io::Error) -> io::Error {
    let hint = if err.kind() == io::ErrorKind::InvalidData {
        " (expected plain text, or an annotation compressed with gzip, BGZF, \
         bzip2, xz or zstd)"
    } else {
        ""
    };
    io::Error::new(
        err.kind(),
        format!("reading gene map {}: {err}{hint}", path.display()),
    )
}

/// Compression suffixes stripped before deciding GTF-vs-TSV, one per format
/// `niffler` decodes. Content detection is the reader's job; this is only about
/// the *name*, and a name that ends in a codec suffix says nothing about the
/// format underneath it.
const COMPRESSION_SUFFIXES: [&str; 6] = ["gz", "bgz", "bgzf", "bz2", "xz", "zst"];

/// Whether the path names a GTF/GFF, ignoring a trailing compression suffix:
/// the last extension of `annotation.gtf.gz` is `gz`, but the content is GTF.
/// Getting this wrong is not a cosmetic error — a GTF handed to the TSV branch
/// parses into `{sequence name: source}` pairs, which is a non-empty map that
/// matches no transcript, so it fails silently rather than loudly.
fn has_gtf_extension(path: &Path) -> bool {
    let ext = |p: &Path| {
        p.extension()
            .and_then(|e| e.to_str())
            .map(str::to_ascii_lowercase)
    };
    let outer = ext(path);
    let effective = match outer.as_deref() {
        Some(e) if COMPRESSION_SUFFIXES.contains(&e) => {
            ext(Path::new(path.file_stem().unwrap_or_default()))
        }
        _ => outer,
    };
    matches!(effective.as_deref(), Some("gtf" | "gff" | "gff3"))
}

/// Open `path`, transparently decompressing it when it is compressed.
///
/// Delegates to [`crate::compress::open_maybe_compressed`] so the gene map is
/// decoded by exactly the same rules as every other user-supplied input, and
/// adds this module's path-annotated error context.
fn open_maybe_compressed(path: &Path) -> io::Result<Box<dyn BufRead>> {
    match crate::compress::open_maybe_compressed(path) {
        Ok(r) => Ok(r),
        Err(e) => Err(annotate_read_error(path, e)),
    }
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
    /// Quantified transcripts with no entry. They are still written, each as its
    /// own single-transcript gene, so they contribute rows to `quant.genes.sf`
    /// and are listed by name in the unmatched-transcript file.
    pub unmatched: usize,
    /// Gene rows actually written — not the number of genes in the gene map.
    /// Counts the `unmatched` fallback rows as well as the joined genes.
    pub genes_written: usize,
}

impl GeneQuantSummary {
    /// Quantified transcripts considered.
    pub fn total(&self) -> usize {
        self.matched + self.unmatched
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
/// A transcript with no entry in `gene_map` is emitted as its own
/// single-transcript gene, keyed by the transcript's own name — C++ salmon's
/// `aggregateEstimatesToGeneLevel` behaviour. No abundance is dropped on the way
/// to gene level, so downstream tools (`tximport` in particular, which has its
/// own transcript-to-gene policy) see everything that was quantified and decide
/// for themselves what to do with the leftovers.
///
/// What C++ got wrong was the reporting, not the aggregation: it warned once per
/// unmatched transcript, so a failed join buried the signal under half a million
/// lines. Here the names go to `unmatched_out` instead, and the caller reports
/// the counts in [`GeneQuantSummary`] once.
///
/// Note that the fallback keys on the transcript name, so an unmatched
/// transcript whose name collides with a real gene name in the map merges into
/// that gene's row. C++ had the same property; it needs a gene named exactly
/// like an unmapped transcript, which does not occur in practice for the
/// Ensembl/GENCODE/RefSeq identifier schemes.
///
/// `unmatched_out`, when given, receives the unmatched transcript names as JSON
/// (see [`write_unmatched_list`]) in `names` order, and is removed when nothing
/// is unmatched so a stale list from an earlier run cannot be mistaken for a
/// current one. The returned
/// [`GeneQuantSummary`] is what the caller should report — in particular
/// `genes_written`, which is the number of rows in the file and not the number
/// of genes in the gene map.
///
/// With `ignore_tx_version`, identifiers are compared with a trailing `.<digits>`
/// suffix removed on both sides (`--ignoreTxVersion`). Off by default: silently
/// normalising identifiers is how a wrong join goes unnoticed in the first place.
#[allow(clippy::too_many_arguments)]
pub fn write_gene_quant(
    out_path: &Path,
    names: &[String],
    lengths: &[u32],
    eff_lengths: &[f64],
    tpm: &[f64],
    counts: &[f64],
    gene_map: &HashMap<String, String>,
    ignore_tx_version: bool,
    unmatched_out: Option<&Path>,
) -> io::Result<GeneQuantSummary> {
    // Re-key the map once when versions are ignored, rather than per lookup.
    let stripped: Option<HashMap<&str, &str>> = ignore_tx_version.then(|| {
        gene_map
            .iter()
            .map(|(t, g)| (strip_tx_version(t), g.as_str()))
            .collect()
    });

    // Group transcript indices by gene (gene name order is sorted for determinism).
    //
    // A `BTreeMap` iterates in sorted key order, unlike a `HashMap`; that is the
    // whole reason it is used here, so the output file is byte-identical between
    // runs.
    let mut genes: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    let mut summary = GeneQuantSummary::default();
    let mut unmatched: Vec<&str> = Vec::new();
    for (i, name) in names.iter().enumerate() {
        let gene = match &stripped {
            Some(m) => m.get(strip_tx_version(name)).copied(),
            None => gene_map.get(name).map(String::as_str),
        };
        match gene {
            Some(g) => {
                genes.entry(g).or_default().push(i);
                summary.matched += 1;
            }
            // No entry: the transcript becomes its own gene, so its abundance
            // still reaches the output. Recorded by name for the caller's list.
            None => {
                genes.entry(name.as_str()).or_default().push(i);
                summary.unmatched += 1;
                unmatched.push(name.as_str());
            }
        }
    }

    if let Some(list_path) = unmatched_out {
        write_unmatched_list(list_path, &unmatched)?;
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
    summary.genes_written = genes.len();
    Ok(summary)
}

/// JSON key holding the unmatched identifiers. Named for what the list *is*
/// rather than for the file it sits in, so a consumer that has the parsed object
/// but not the path still knows what it is looking at.
const UNMATCHED_JSON_KEY: &str = "unmatched_transcripts";

/// Write the names of transcripts the gene map did not cover, in `quant.sf`
/// order, as a one-key JSON object:
///
/// ```json
/// {
///   "unmatched_transcripts": [
///     "ENST00000456328.2",
///     "ENST00000450305.3"
///   ]
/// }
/// ```
///
/// JSON rather than one-per-line text so it parses without a bespoke reader,
/// matching the other non-primary artifacts under `aux_info/`. Wrapped in an
/// object rather than written as a bare array so the file can gain fields later
/// without breaking consumers written against it today.
///
/// An empty list removes the file rather than writing an empty one: the run that
/// produced it is the only thing that can be said to be current, and a stale
/// list left over from a previous run in the same output directory would claim
/// failures this run did not have.
fn write_unmatched_list(path: &Path, unmatched: &[&str]) -> io::Result<()> {
    if unmatched.is_empty() {
        return match std::fs::remove_file(path) {
            Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(()),
            other => other,
        };
    }
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    // `to_string_pretty` puts one identifier per line, so the file stays as
    // greppable as the text version it replaces.
    let doc = serde_json::json!({ UNMATCHED_JSON_KEY: unmatched });
    let body = serde_json::to_string_pretty(&doc).map_err(io::Error::other)?;
    let mut w = io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(w, "{body}")?;
    w.flush()
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
        let list = scratch("nomatch").join("genemap_unmatched_txps.json");
        let summary = write_gene_quant(
            &out,
            &names,
            &[1000, 2000],
            &[800.0, 1800.0],
            &[10.0, 20.0],
            &[100.0, 200.0],
            &map,
            false,
            Some(&list),
        )
        .unwrap();

        assert_eq!(
            summary.unmatched, 2,
            "every transcript should fail to match"
        );
        assert_eq!(summary.matched, 0);
        assert_eq!(summary.match_rate(), 0.0);
        assert!(summary.match_rate_is_low());

        // Nothing is dropped: each unmatched transcript stands in as its own
        // gene, so the abundance still reaches gene level for tximport to
        // re-aggregate under its own policy.
        assert_eq!(summary.genes_written, 2);
        let body = std::fs::read_to_string(&out).unwrap();
        assert_eq!(body.lines().count(), 3, "header + 2 rows: {body:?}");
        let row = body
            .lines()
            .find(|l| l.starts_with("ENST00000456328.2\t"))
            .unwrap_or_else(|| panic!("rows are keyed by transcript name: {body:?}"));
        assert_eq!(row.split('\t').nth(4), Some("100.000"), "reads preserved");
        // Mass is conserved through the join even when nothing joins.
        let reads: f64 = body
            .lines()
            .skip(1)
            .map(|l| l.split('\t').nth(4).unwrap().parse::<f64>().unwrap())
            .sum();
        assert_eq!(reads, 300.0);

        // And the failure is recorded where it outlives the terminal. Asserted
        // through a JSON parse, not on the raw text, so the test pins the
        // document a consumer sees rather than the formatting around it.
        let listed: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&list).unwrap()).unwrap();
        assert_eq!(
            listed[UNMATCHED_JSON_KEY],
            serde_json::json!(["ENST00000456328.2", "ENST00000450305.3"]),
            "listed under the documented key, in quant.sf order"
        );
        assert_eq!(
            listed.as_object().unwrap().len(),
            1,
            "one key, so a consumer can match on it exactly"
        );

        // The diagnosis the warning is built on: the identifiers differ only by
        // the version suffix, so ignoring it would match everything.
        assert_eq!(matches_ignoring_tx_version(&names, &map), 2);
    }

    /// An unmatched transcript named exactly like a gene in the map merges into
    /// that gene's row rather than adding one. C++ behaved the same way; this
    /// pins it so the collision is a known property and not a surprise.
    #[test]
    fn a_fallback_row_colliding_with_a_real_gene_name_merges_into_it() {
        let mut gm = HashMap::new();
        gm.insert("t1".to_string(), "G".to_string());
        let names = vec!["t1".to_string(), "G".to_string()];
        let out = scratch("collide").join("quant.genes.sf");
        let s = write_gene_quant(
            &out,
            &names,
            &[10, 20],
            &[8.0, 18.0],
            &[1.0, 2.0],
            &[3.0, 4.0],
            &gm,
            false,
            None,
        )
        .unwrap();
        assert_eq!(s.matched, 1);
        assert_eq!(s.unmatched, 1);
        assert_eq!(s.genes_written, 1, "both land on the row named G");
        let body = std::fs::read_to_string(&out).unwrap();
        let row = body.lines().find(|l| l.starts_with("G\t")).unwrap();
        assert_eq!(row.split('\t').nth(4), Some("7.000"), "3 + 4 reads");
    }

    /// A clean run must not leave the previous run's failure list behind: an
    /// empty list is the absence of the file, not an empty file.
    #[test]
    fn a_stale_unmatched_list_is_removed_when_everything_matches() {
        let mut gm = HashMap::new();
        gm.insert("t1".to_string(), "G".to_string());
        let out = scratch("stale").join("quant.genes.sf");
        let list = scratch("stale").join("genemap_unmatched_txps.json");
        std::fs::write(&list, "leftover_from_a_previous_run\n").unwrap();

        let s = write_gene_quant(
            &out,
            &["t1".to_string()],
            &[10],
            &[8.0],
            &[1.0],
            &[3.0],
            &gm,
            false,
            Some(&list),
        )
        .unwrap();
        assert_eq!(s.unmatched, 0);
        assert!(!list.exists(), "stale list must not survive a clean run");
    }

    /// The same input as the mismatch case, with `--ignoreTxVersion` on: the
    /// join succeeds and the gene rows appear.
    #[test]
    fn ignore_tx_version_matches_versioned_names_to_an_unversioned_map() {
        let mut gm = HashMap::new();
        gm.insert("ENST00000456328".to_string(), "ENSG00000290825".to_string());
        gm.insert("ENST00000450305".to_string(), "ENSG00000223972".to_string());
        let names = vec![
            "ENST00000456328.2".to_string(),
            "ENST00000450305.3".to_string(),
        ];
        let out = scratch("ignorever").join("quant.genes.sf");

        let off = write_gene_quant(
            &out,
            &names,
            &[1000, 2000],
            &[800.0, 1800.0],
            &[10.0, 20.0],
            &[100.0, 200.0],
            &gm,
            false,
            None,
        )
        .unwrap();
        assert_eq!(off.matched, 0, "default must not strip versions");
        // Two rows either way; the flag changes what they are named after, and
        // that is the whole difference the warnings have to convey.
        assert_eq!(off.unmatched, 2);
        assert_eq!(off.genes_written, 2);
        let off_body = std::fs::read_to_string(&out).unwrap();
        assert!(
            off_body.contains("ENST00000456328.2\t"),
            "rows keyed by transcript when the join fails: {off_body:?}"
        );

        let on = write_gene_quant(
            &out,
            &names,
            &[1000, 2000],
            &[800.0, 1800.0],
            &[10.0, 20.0],
            &[100.0, 200.0],
            &gm,
            true,
            None,
        )
        .unwrap();
        assert_eq!(on.matched, 2);
        assert_eq!(on.unmatched, 0);
        assert_eq!(on.genes_written, 2);
        assert!(!on.match_rate_is_low());

        let body = std::fs::read_to_string(&out).unwrap();
        assert_eq!(body.lines().count(), 3, "header + 2 genes: {body:?}");
        // Counts must survive the join unchanged.
        let row = body
            .lines()
            .find(|l| l.starts_with("ENSG00000290825\t"))
            .unwrap();
        assert_eq!(row.split('\t').nth(4), Some("100.000"));
    }

    /// Stripping applies to both sides, so a versioned gene map joins to
    /// unversioned quant names too.
    #[test]
    fn ignore_tx_version_strips_the_gene_map_side_as_well() {
        let mut gm = HashMap::new();
        gm.insert("ENST00000456328.2".to_string(), "ENSG1".to_string());
        let names = vec!["ENST00000456328".to_string()];
        let out = scratch("ignorever2").join("quant.genes.sf");
        let s =
            write_gene_quant(&out, &names, &[10], &[8.0], &[1.0], &[3.0], &gm, true, None).unwrap();
        assert_eq!(s.matched, 1);
        assert_eq!(s.genes_written, 1);
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
            false,
            None,
        )
        .unwrap();
        assert_eq!(summary.matched, 1);
        assert_eq!(summary.unmatched, 1);
        // Two rows (gene G, plus `zzz` standing in for itself) — not the 3 genes
        // the map knows about. Reporting the map's size is what let a
        // header-only file be announced as a success in #1075.
        assert_eq!(summary.genes_written, 2, "rows written, not 3 map genes");
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
            false,
            None,
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
        let summary =
            write_gene_quant(&p, &names, &lengths, &eff, &tpm, &counts, &gm, false, None).unwrap();
        assert_eq!(summary.unmatched, 0);
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
