//! `salmon quantmerge`: merge one column across several samples' `quant.sf`
//! (or `quant.genes.sf`) into a single target × sample matrix TSV.
//!
//! # What this is for
//!
//! A salmon run produces one `quant.sf` per sample. Almost every downstream
//! analysis (differential expression, clustering, plotting) instead wants a
//! single matrix: one row per transcript or gene, one column per sample, filled
//! with a single quantity such as TPM. This module does exactly that join.
//!
//! Values are copied through as text rather than parsed and reformatted, so the
//! merged matrix contains the same digits the individual files did.

use std::collections::{HashMap, HashSet};
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

/// Columns in a `quant.sf` row; see [`MergeColumn::index`].
const QUANT_COLUMNS: usize = 5;

/// Which `quant.sf` column to merge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MergeColumn {
    /// Reference length.
    Len,
    /// Effective length.
    Elen,
    /// Transcripts per million.
    Tpm,
    /// Estimated fragment count.
    NumReads,
}

impl MergeColumn {
    /// Parse a column name (case-insensitive), accepting salmon's aliases.
    ///
    /// The generous alias list exists because users reasonably type any of
    /// `elen`, `effLen`, `EffectiveLength` for the same column.
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_ascii_uppercase().as_str() {
            "LEN" | "LENGTH" => Some(Self::Len),
            "ELEN" | "ELENGTH" | "EFFLEN" | "EFFLENGTH" | "EFFECTIVELEN" | "EFFECTIVELENGTH" => {
                Some(Self::Elen)
            }
            "TPM" => Some(Self::Tpm),
            "NUMREADS" | "NUM_READS" | "NREADS" => Some(Self::NumReads),
            _ => None,
        }
    }

    /// Column index in a `quant.sf` row (Name=0, Length=1, EffectiveLength=2,
    /// TPM=3, NumReads=4). The layout is fixed by salmon's output format.
    fn index(self) -> usize {
        match self {
            Self::Len => 1,
            Self::Elen => 2,
            Self::Tpm => 3,
            Self::NumReads => 4,
        }
    }
}

/// Merge the chosen `column` of each sample's quant file into `out_path`. The
/// output has a `Name` header plus one column per sample (`names`), then one row
/// per target (in first-seen order across the samples), with `missing` filling
/// targets absent from a sample. The per-sample value strings are taken verbatim
/// from each quant file. Returns the number of missing entries written.
///
/// **Why targets can be missing.** Samples quantified against different index
/// versions, or gene files from annotations that disagree, will not have exactly
/// the same target set. Rather than failing or silently dropping rows, the union
/// of targets is emitted and the gaps are filled with the caller's `missing`
/// placeholder (typically `NA`), which is both visible and easy to filter.
///
/// Damage is different from absence, so these are errors, reported before the
/// output is created: an empty or repeated sample name, a row with too few
/// columns, and a target appearing twice in one file.
pub fn merge_quants(
    quant_files: &[PathBuf],
    names: &[String],
    column: MergeColumn,
    missing: &str,
    out_path: &Path,
) -> io::Result<usize> {
    // Caller-side invariant: one display name per input file.
    assert_eq!(quant_files.len(), names.len());
    let invalid = |msg: String| io::Error::new(io::ErrorKind::InvalidData, msg);
    let mut seen_names = HashSet::new();
    for n in names {
        if n.is_empty() {
            return Err(invalid(
                "a sample has an empty name; pass --names to name the samples".to_string(),
            ));
        }
        if !seen_names.insert(n) {
            return Err(invalid(format!(
                "two samples are both named \"{n}\"; pass --names to tell them apart"
            )));
        }
    }
    let nsamp = quant_files.len();
    let col = column.index();

    // `order` preserves first-seen target order for the output; `values` holds
    // one slot per sample for each target. Two structures because a `HashMap`
    // has no meaningful iteration order.
    let mut order: Vec<String> = Vec::new();
    let mut values: HashMap<String, Vec<Option<String>>> = HashMap::new();

    for (s, qf) in quant_files.iter().enumerate() {
        let reader = io::BufReader::new(
            // Rewrap the error with the path, so a failure names the file that
            // caused it rather than just "no such file".
            std::fs::File::open(qf)
                .map_err(|e| io::Error::new(e.kind(), format!("opening {}: {e}", qf.display())))?,
        );
        for (lineno, line) in reader.lines().enumerate() {
            let line = line?;
            if lineno == 0 {
                continue; // header
            }
            if line.trim().is_empty() {
                continue;
            }
            let toks: Vec<&str> = line.split('\t').collect();
            // `NA` means "not in this sample", so a damaged row must not become one.
            if toks.len() < QUANT_COLUMNS {
                return Err(invalid(format!(
                    "{} line {}: expected {QUANT_COLUMNS} columns, found {}",
                    qf.display(),
                    lineno + 1,
                    toks.len()
                )));
            }
            let name = toks[0];
            // First time this target is seen anywhere: record its position in
            // the output order and create its all-missing row.
            let entry = values.entry(name.to_string()).or_insert_with(|| {
                order.push(name.to_string());
                vec![None; nsamp]
            });
            if entry[s].is_some() {
                return Err(invalid(format!(
                    "{} line {}: target {name} appears twice",
                    qf.display(),
                    lineno + 1
                )));
            }
            entry[s] = Some(toks[col].to_string());
        }
    }

    // Create the output directory if the caller asked for a nested path. The
    // empty check covers a bare filename, whose parent is `""`.
    if let Some(parent) = out_path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    let mut w = io::BufWriter::new(std::fs::File::create(out_path)?);
    write!(w, "Name")?;
    for n in names {
        write!(w, "\t{n}")?;
    }
    writeln!(w)?;

    let mut n_missing = 0usize;
    for name in &order {
        write!(w, "{name}")?;
        for v in &values[name] {
            match v {
                Some(s) => write!(w, "\t{s}")?,
                None => {
                    n_missing += 1;
                    write!(w, "\t{missing}")?;
                }
            }
        }
        writeln!(w)?;
    }
    w.flush()?;
    Ok(n_missing)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Every documented alias must resolve, and an unknown name must be
    /// rejected rather than defaulting to some column.
    #[test]
    fn column_aliases() {
        assert_eq!(MergeColumn::parse("tpm"), Some(MergeColumn::Tpm));
        assert_eq!(MergeColumn::parse("NumReads"), Some(MergeColumn::NumReads));
        assert_eq!(
            MergeColumn::parse("effectiveLength"),
            Some(MergeColumn::Elen)
        );
        assert_eq!(MergeColumn::parse("len"), Some(MergeColumn::Len));
        assert_eq!(MergeColumn::parse("bogus"), None);
    }

    /// Samples with different target sets must union, keep first-seen row order,
    /// and fill the gaps both ways (t1 missing from s2, t3 missing from s1).
    #[test]
    fn merges_with_missing() {
        let dir = std::env::temp_dir().join(format!("qm_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let s1 = dir.join("s1.sf");
        let s2 = dir.join("s2.sf");
        std::fs::write(
            &s1,
            "Name\tLength\tEffectiveLength\tTPM\tNumReads\nt1\t100\t80.0\t10.5\t5\nt2\t200\t180.0\t20.0\t9\n",
        )
        .unwrap();
        // s2 is missing t1, has a new t3
        std::fs::write(
            &s2,
            "Name\tLength\tEffectiveLength\tTPM\tNumReads\nt2\t200\t180.0\t30.0\t12\nt3\t50\t30.0\t1.0\t1\n",
        )
        .unwrap();
        let out = dir.join("merged.tsv");
        let n = merge_quants(
            &[s1, s2],
            &["A".into(), "B".into()],
            MergeColumn::Tpm,
            "NA",
            &out,
        )
        .unwrap();
        let body = std::fs::read_to_string(&out).unwrap();
        let lines: Vec<&str> = body.lines().collect();
        assert_eq!(lines[0], "Name\tA\tB");
        // first-seen order: t1, t2, t3
        assert_eq!(lines[1], "t1\t10.5\tNA");
        assert_eq!(lines[2], "t2\t20.0\t30.0");
        assert_eq!(lines[3], "t3\tNA\t1.0");
        assert_eq!(n, 2);
    }

    const HEADER: &str = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n";

    /// Merge one file per `bodies` entry under `names`, in a directory of its own.
    fn merge_in(tag: &str, bodies: &[&str], names: &[&str]) -> (io::Result<usize>, PathBuf) {
        let dir = std::env::temp_dir().join(format!("qm_{}_{tag}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let files: Vec<PathBuf> = bodies
            .iter()
            .enumerate()
            .map(|(i, b)| {
                let f = dir.join(format!("{i}.sf"));
                std::fs::write(&f, format!("{HEADER}{b}")).unwrap();
                f
            })
            .collect();
        let names: Vec<String> = names.iter().map(|n| n.to_string()).collect();
        let out = dir.join("merged.tsv");
        (
            merge_quants(&files, &names, MergeColumn::Tpm, "NA", &out),
            out,
        )
    }

    /// Two columns with one name cannot be told apart downstream, and an empty
    /// one (`--quants .`) has no name at all.
    #[test]
    fn duplicate_or_empty_sample_names_are_refused() {
        let row = "t1\t100\t80.0\t10.5\t5\n";
        let (r, out) = merge_in("dup", &[row, row], &["quant", "quant"]);
        let e = r.unwrap_err().to_string();
        assert!(e.contains("quant"), "{e}");
        assert!(!out.exists());
        let (r, out) = merge_in("empty", &[row], &[""]);
        assert!(r.is_err());
        assert!(!out.exists());
    }

    /// A damaged row used to become `NA`, which reads as "target not in sample".
    #[test]
    fn a_short_row_is_refused_and_named() {
        let (r, out) = merge_in("short", &["t1\t100\t80.0\n"], &["A"]);
        let e = r.unwrap_err().to_string();
        assert!(e.contains("0.sf") && e.contains("line 2"), "{e}");
        assert!(!out.exists());
        // Short only in a column not being merged is still a damaged file.
        let (r, _) = merge_in("short2", &["t1\t100\t80.0\t10.5\n"], &["A"]);
        assert!(r.is_err());
    }

    /// The later row used to overwrite the earlier one silently.
    #[test]
    fn a_target_repeated_in_one_sample_is_refused() {
        let body = "t1\t100\t80.0\t10.5\t5\nt1\t100\t80.0\t99.0\t7\n";
        let (r, out) = merge_in("rep", &[body], &["A"]);
        let e = r.unwrap_err().to_string();
        assert!(e.contains("t1") && e.contains("line 3"), "{e}");
        assert!(!out.exists());
    }
}
