//! `salmon quantmerge`: merge one column across several samples' `quant.sf`
//! (or `quant.genes.sf`) into a single target × sample matrix TSV.

use std::collections::HashMap;
use std::io::{self, BufRead, Write};
use std::path::{Path, PathBuf};

/// Which `quant.sf` column to merge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MergeColumn {
    Len,
    Elen,
    Tpm,
    NumReads,
}

impl MergeColumn {
    /// Parse a column name (case-insensitive), accepting salmon's aliases.
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
    /// TPM=3, NumReads=4).
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
pub fn merge_quants(
    quant_files: &[PathBuf],
    names: &[String],
    column: MergeColumn,
    missing: &str,
    out_path: &Path,
) -> io::Result<usize> {
    assert_eq!(quant_files.len(), names.len());
    let nsamp = quant_files.len();
    let col = column.index();

    let mut order: Vec<String> = Vec::new();
    let mut values: HashMap<String, Vec<Option<String>>> = HashMap::new();

    for (s, qf) in quant_files.iter().enumerate() {
        let reader = io::BufReader::new(std::fs::File::open(qf).map_err(|e| {
            io::Error::new(e.kind(), format!("opening {}: {e}", qf.display()))
        })?);
        for (lineno, line) in reader.lines().enumerate() {
            let line = line?;
            if lineno == 0 {
                continue; // header
            }
            if line.trim().is_empty() {
                continue;
            }
            let toks: Vec<&str> = line.split('\t').collect();
            if toks.len() <= col {
                continue;
            }
            let name = toks[0];
            let entry = values.entry(name.to_string()).or_insert_with(|| {
                order.push(name.to_string());
                vec![None; nsamp]
            });
            entry[s] = Some(toks[col].to_string());
        }
    }

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

    #[test]
    fn column_aliases() {
        assert_eq!(MergeColumn::parse("tpm"), Some(MergeColumn::Tpm));
        assert_eq!(MergeColumn::parse("NumReads"), Some(MergeColumn::NumReads));
        assert_eq!(MergeColumn::parse("effectiveLength"), Some(MergeColumn::Elen));
        assert_eq!(MergeColumn::parse("len"), Some(MergeColumn::Len));
        assert_eq!(MergeColumn::parse("bogus"), None);
    }

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
}
