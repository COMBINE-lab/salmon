//! Unspliced (intronic / pre-mRNA) targets for total-RNA quantification.
//!
//! # Why
//!
//! Total (ribo-depleted) RNA-seq libraries carry a large share of nascent,
//! not-yet-spliced RNA: typically half or more of the fragments fall in introns.
//! With a transcripts-only index those fragments are either lost or, worse,
//! forced onto whichever annotated isoform contains the intronic sequence,
//! which is usually a `retained_intron` isoform. A genome decoy does not fix
//! this: a fragment only goes to the decoy when its genome score is *strictly*
//! greater than its best transcript score, so an intronic fragment that lies
//! inside a retained-intron isoform ties and stays on the isoform.
//!
//! The remedy used by the single-cell tools (alevin-fry's *splici* /
//! *spliceu* references, He et al. 2022) is to index the unspliced sequence of
//! every gene next to its spliced transcripts, so pre-mRNA fragments have a
//! target of their own and the EM apportions the shared (retained-intron)
//! evidence between the two by abundance. This module derives those unspliced
//! targets from a genome FASTA and a GTF, inside `salmon index`, following the
//! conventions of `pyroe make-splici` so the resulting tables are familiar.
//!
//! # Coordinates
//!
//! GTF coordinates are 1-based and inclusive; everything here is 0-based and
//! half-open (`[start, end)`), converted once at parse time.

use std::collections::HashMap;
use std::io::BufRead;
use std::path::Path;

use anyhow::{bail, Context, Result};
use salmon_core::genemap::extract_attr;

/// A genomic interval, 0-based half-open `[start, end)`.
pub type Interval = (u64, u64);

/// Genomic strand of a gene.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Strand {
    Plus,
    Minus,
}

/// One gene's exon structure on one sequence and strand.
///
/// A `gene_id` that the annotation places on two sequences (GENCODE's PAR
/// genes before the `_PAR_Y` suffix, for instance) yields two models: they
/// share the id but not the coordinates.
#[derive(Debug, Clone)]
pub struct GeneModel {
    pub gene_id: String,
    pub seqname: String,
    pub strand: Strand,
    /// Exons of each transcript, sorted by start.
    pub transcripts: Vec<Vec<Interval>>,
}

/// The parts of a GTF the unspliced targets need.
#[derive(Debug, Default)]
pub struct Annotation {
    /// Gene models in order of first appearance in the GTF.
    pub genes: Vec<GeneModel>,
    /// `transcript_id` -> `gene_id`. When the GTF carries the version in its
    /// own attribute (Ensembl: `transcript_id "ENST…"; transcript_version "2"`)
    /// the versioned id `ENST….2` is recorded too, since that is how Ensembl
    /// names the records of its cDNA FASTA.
    pub tx_to_gene: HashMap<String, String>,
    /// Exon records skipped because their strand was neither `+` nor `-`.
    pub unstranded_exons: usize,
}

impl Annotation {
    /// The gene of transcript `name`, trying the name as given and then with a
    /// trailing `.<digits>` version removed.
    pub fn gene_of(&self, name: &str) -> Option<&str> {
        self.tx_to_gene
            .get(name)
            .or_else(|| {
                self.tx_to_gene
                    .get(salmon_core::genemap::strip_tx_version(name))
            })
            .map(String::as_str)
    }
}

/// Read the exon records of a GTF (optionally compressed).
pub fn read_gtf(path: &Path) -> Result<Annotation> {
    let reader = salmon_core::compress::open_maybe_compressed(path)
        .with_context(|| format!("opening GTF {}", path.display()))?;
    parse_gtf(reader).with_context(|| format!("parsing GTF {}", path.display()))
}

/// Parse the `exon` records of a GTF. Only exon lines are read: a gene's
/// introns and body are both derived from its transcripts' exons, which is
/// what `pyroe make-splici` does too, and it spares us from trusting `gene` /
/// `transcript` lines that not every annotation provides.
pub fn parse_gtf<R: BufRead>(reader: R) -> Result<Annotation> {
    let mut ann = Annotation::default();
    // (gene_id, seqname, strand) -> index into `ann.genes`
    let mut gene_index: HashMap<(String, String, Strand), usize> = HashMap::new();
    // per gene model: transcript_id -> index into its `transcripts`
    let mut tx_index: Vec<HashMap<String, usize>> = Vec::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.splitn(9, '\t').collect();
        if cols.len() < 9 || cols[2] != "exon" {
            continue;
        }
        let attrs = cols[8];
        let (Some(gene_id), Some(tx_id)) = (
            extract_attr(attrs, "gene_id"),
            extract_attr(attrs, "transcript_id"),
        ) else {
            continue;
        };
        let parse = |s: &str| -> Result<u64> {
            s.parse::<u64>()
                .with_context(|| format!("line {}: invalid coordinate {s:?}", lineno + 1))
        };
        let (start, end) = (parse(cols[3])?, parse(cols[4])?);
        if start == 0 || end < start {
            bail!(
                "line {}: invalid exon coordinates {start}..{end} (GTF is 1-based, start <= end)",
                lineno + 1
            );
        }
        let strand = match cols[6] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            _ => {
                ann.unstranded_exons += 1;
                continue;
            }
        };

        if let Some(v) = extract_attr(attrs, "transcript_version") {
            if salmon_core::genemap::strip_tx_version(&tx_id) == tx_id {
                ann.tx_to_gene
                    .entry(format!("{tx_id}.{v}"))
                    .or_insert_with(|| gene_id.clone());
            }
        }
        ann.tx_to_gene
            .entry(tx_id.clone())
            .or_insert_with(|| gene_id.clone());

        let key = (gene_id, cols[0].to_string(), strand);
        let gi = match gene_index.get(&key) {
            Some(&gi) => gi,
            None => {
                ann.genes.push(GeneModel {
                    gene_id: key.0.clone(),
                    seqname: key.1.clone(),
                    strand,
                    transcripts: Vec::new(),
                });
                tx_index.push(HashMap::new());
                gene_index.insert(key, ann.genes.len() - 1);
                ann.genes.len() - 1
            }
        };
        let gene = &mut ann.genes[gi];
        let ti = *tx_index[gi].entry(tx_id).or_insert_with(|| {
            gene.transcripts.push(Vec::new());
            gene.transcripts.len() - 1
        });
        // 1-based inclusive -> 0-based half-open
        gene.transcripts[ti].push((start - 1, end));
    }

    for gene in &mut ann.genes {
        for exons in &mut gene.transcripts {
            exons.sort_unstable();
        }
    }
    if ann.genes.is_empty() {
        bail!(
            "no stranded exon records carrying both gene_id and transcript_id were found \
             (a GTF is required; GFF3 is not supported for unspliced targets)"
        );
    }
    Ok(ann)
}

/// Sort and merge overlapping or book-ended intervals.
pub fn merge_intervals(mut v: Vec<Interval>) -> Vec<Interval> {
    v.sort_unstable();
    let mut out: Vec<Interval> = Vec::with_capacity(v.len());
    for (s, e) in v {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

/// A transcript's introns: the gaps between its (sorted) exons. Overlapping or
/// abutting exons leave no gap.
fn transcript_introns(exons: &[Interval], out: &mut Vec<Interval>) {
    let mut reach: Option<u64> = None;
    for &(s, e) in exons {
        if let Some(r) = reach {
            if s > r {
                out.push((r, s));
            }
        }
        reach = Some(reach.map_or(e, |r| r.max(e)));
    }
}

/// The union of the introns of every transcript of `gene`, merged.
///
/// This is the *splici* definition: an interval is intronic if it is an intron
/// of at least one transcript, even when another transcript of the same gene
/// (typically a retained-intron isoform) covers it with an exon. Those shared
/// stretches are exactly the ambiguity the EM has to resolve, and they are
/// resolved using the evidence each target has elsewhere.
pub fn merged_introns(gene: &GeneModel) -> Vec<Interval> {
    let mut v = Vec::new();
    for exons in &gene.transcripts {
        transcript_introns(exons, &mut v);
    }
    merge_intervals(v)
}

/// The gene body (first exon start to last exon end) when the gene has at
/// least one intron; `None` for a gene whose transcripts are all single-exon
/// and contiguous, whose "pre-mRNA" would be a copy of its mature transcript.
pub fn gene_body(gene: &GeneModel) -> Option<Interval> {
    if merged_introns(gene).is_empty() {
        return None;
    }
    let start = gene
        .transcripts
        .iter()
        .filter_map(|t| t.first())
        .map(|e| e.0)
        .min()?;
    let end = gene.transcripts.iter().flatten().map(|e| e.1).max()?;
    Some((start, end))
}

/// Which unspliced targets `salmon index --unspliced` adds next to the
/// spliced transcripts.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UnsplicedMode {
    /// No unspliced targets (the default; the index is unchanged).
    #[default]
    None,
    /// Merged introns of each gene, extended by a flank on both sides
    /// (`pyroe make-splici`).
    Intron,
    /// The full gene body, first exon start to last exon end
    /// (`pyroe make-spliceu`).
    Premrna,
}

impl UnsplicedMode {
    /// Name used on the command line and in `info.json`.
    pub fn as_str(self) -> &'static str {
        match self {
            UnsplicedMode::None => "none",
            UnsplicedMode::Intron => "intron",
            UnsplicedMode::Premrna => "premrna",
        }
    }
}

impl std::str::FromStr for UnsplicedMode {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "none" => Ok(UnsplicedMode::None),
            "intron" => Ok(UnsplicedMode::Intron),
            "premrna" => Ok(UnsplicedMode::Premrna),
            other => Err(format!(
                "unknown unspliced mode {other:?} (expected none, intron or premrna)"
            )),
        }
    }
}

/// How unspliced targets are held in the index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UnsplicedLayout {
    /// `projection` when every target lies on a sequence given as a decoy
    /// (the genome is in the gentrome), `sequence` otherwise.
    #[default]
    Auto,
    /// Targets carry no k-mers; alignments to the genome decoy inside a
    /// target's interval are projected onto it at quant time. The k-mer index
    /// and mapping cost stay those of the gentrome.
    Projection,
    /// Target sequences are indexed like transcripts (works without decoys).
    Sequence,
}

impl UnsplicedLayout {
    /// Name used on the command line and in `info.json`.
    pub fn as_str(self) -> &'static str {
        match self {
            UnsplicedLayout::Auto => "auto",
            UnsplicedLayout::Projection => "projection",
            UnsplicedLayout::Sequence => "sequence",
        }
    }
}

impl std::str::FromStr for UnsplicedLayout {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "auto" => Ok(UnsplicedLayout::Auto),
            "projection" => Ok(UnsplicedLayout::Projection),
            "sequence" => Ok(UnsplicedLayout::Sequence),
            other => Err(format!(
                "unknown unspliced layout {other:?} (expected auto, projection or sequence)"
            )),
        }
    }
}

/// Suffix that marks an unspliced target, followed by an ordinal from the
/// second interval of a gene on (`G-I`, `G-I1`, `G-I2`, ...). This is the
/// naming `pyroe make-splici` / `make-spliceu` use, so `t2g_3col.tsv` and
/// downstream tooling written for alevin-fry read these indices as they are.
pub const UNSPLICED_SUFFIX: &str = "-I";

/// One unspliced target: a genomic interval of one gene, read on the gene's
/// strand.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnsplicedTarget {
    pub name: String,
    pub gene_id: String,
    pub seqname: String,
    pub strand: Strand,
    /// 0-based half-open. `end` may exceed the sequence length when a flank
    /// runs off the end of a chromosome; extraction clips it.
    pub start: u64,
    pub end: u64,
}

/// Derive the unspliced targets of every gene, named deterministically.
///
/// * `Intron`: the gene's merged introns, each extended by `flank` bases on
///   both sides (clipped at 0 here and at the sequence end on extraction), then
///   merged again so flanks that meet across a short exon become one target.
///   With `flank = L - 1` for reads of length `L`, any read with at least one
///   intronic base lies wholly inside a target, so its unspliced placement
///   scores as well as its genome placement and it is not lost to a decoy.
/// * `Premrna`: the gene body; `flank` is ignored.
///
/// Genes without an intron get no target in either mode (their unspliced form
/// would duplicate a mature transcript). Names are `<gene_id>-I`, then
/// `<gene_id>-I1`, `-I2`, ... over a gene's intervals in annotation order, then
/// by start, across all the sequences a `gene_id` appears on.
pub fn unspliced_targets(
    ann: &Annotation,
    mode: UnsplicedMode,
    flank: u64,
) -> Vec<UnsplicedTarget> {
    let mut out = Vec::new();
    let mut ordinal: HashMap<&str, usize> = HashMap::new();
    for gene in &ann.genes {
        let intervals = match mode {
            UnsplicedMode::None => return out,
            UnsplicedMode::Intron => merge_intervals(
                merged_introns(gene)
                    .into_iter()
                    .map(|(s, e)| (s.saturating_sub(flank), e.saturating_add(flank)))
                    .collect(),
            ),
            UnsplicedMode::Premrna => gene_body(gene).into_iter().collect(),
        };
        for (start, end) in intervals {
            let n = ordinal.entry(gene.gene_id.as_str()).or_insert(0);
            let name = if *n == 0 {
                format!("{}{UNSPLICED_SUFFIX}", gene.gene_id)
            } else {
                format!("{}{UNSPLICED_SUFFIX}{n}", gene.gene_id)
            };
            *n += 1;
            out.push(UnsplicedTarget {
                name,
                gene_id: gene.gene_id.clone(),
                seqname: gene.seqname.clone(),
                strand: gene.strand,
                start,
                end,
            });
        }
    }
    out
}

/// What [`write_unspliced_fasta`] did.
#[derive(Debug, Default)]
pub struct ExtractStats {
    /// targets written to the FASTA
    pub written: usize,
    /// targets whose flank ran past the end of their sequence and was clipped
    pub clipped: usize,
    /// targets dropped because less than half of their bases are ACGT (for
    /// example a gene on a hard-masked PAR copy of chrY): non-ACGT bases would
    /// otherwise be replaced by random ones at indexing time
    pub dropped_low_acgt: Vec<String>,
    /// targets dropped because their sequence is absent from the genome FASTA
    pub dropped_missing_seq: usize,
    /// annotation sequence names absent from the genome FASTA
    pub missing_seqnames: Vec<String>,
    /// every written target with its final (clipped) genomic interval, in
    /// output order
    pub targets: Vec<WrittenTarget>,
    /// digest of each genome sequence that carries targets (see
    /// [`sequence_digest`]), to check that the decoys are the same sequences
    pub seq_digests: HashMap<String, u64>,
}

/// A target as written: its genomic interval after clipping at the sequence
/// end, 0-based half-open.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WrittenTarget {
    pub name: String,
    pub seqname: String,
    pub start: u64,
    pub end: u64,
    pub minus: bool,
}

/// Case-insensitive digest of a sequence (xxh3 over its uppercased bytes),
/// used to check that a genome FASTA and the index's decoys hold the same
/// sequences before alignments to the latter are projected with coordinates
/// from the former.
pub fn sequence_digest(seq: &[u8]) -> u64 {
    let mut h = xxhash_rust::xxh3::Xxh3::new();
    let mut buf = [0u8; 8192];
    for chunk in seq.chunks(buf.len()) {
        let b = &mut buf[..chunk.len()];
        for (o, &c) in b.iter_mut().zip(chunk) {
            *o = c.to_ascii_uppercase();
        }
        h.update(b);
    }
    h.digest()
}

/// Complement of an (uppercase) base; anything but ACGT becomes `N`.
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// Stream `genome` one sequence at a time and write the sequence of every
/// target to `out` (uppercased, reverse-complemented for minus-strand genes).
///
/// Only one chromosome is held in memory at a time. Records are written in
/// genome FASTA order, then by start, end and name, so the output (and hence
/// the index) is a deterministic function of the inputs. Sequences are matched
/// on the first whitespace-delimited token of the FASTA header, the GTF
/// `seqname` convention of both GENCODE (`chr1`) and Ensembl (`1`).
pub fn write_unspliced_fasta(
    genome: &Path,
    targets: &[UnsplicedTarget],
    out: &Path,
) -> Result<ExtractStats> {
    use std::io::Write as _;
    let mut by_seq: HashMap<&str, Vec<&UnsplicedTarget>> = HashMap::new();
    for t in targets {
        by_seq.entry(t.seqname.as_str()).or_default().push(t);
    }
    for v in by_seq.values_mut() {
        v.sort_by(|a, b| (a.start, a.end, &a.name).cmp(&(b.start, b.end, &b.name)));
    }
    let mut w = std::io::BufWriter::new(
        std::fs::File::create(out).with_context(|| format!("creating {}", out.display()))?,
    );
    let mut stats = ExtractStats::default();
    let mut seen_seqs: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut reader = needletail::parse_fastx_file(genome)
        .with_context(|| format!("opening genome {}", genome.display()))?;
    let mut buf: Vec<u8> = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec.context("reading genome FASTA record")?;
        let id = rec.id();
        let end = id
            .iter()
            .position(|b| b.is_ascii_whitespace())
            .unwrap_or(id.len());
        let seqname = String::from_utf8_lossy(&id[..end]).into_owned();
        if !seen_seqs.insert(seqname.clone()) {
            bail!(
                "genome sequence {seqname:?} appears twice in {}",
                genome.display()
            );
        }
        let Some(ts) = by_seq.get(seqname.as_str()) else {
            continue;
        };
        let seq = rec.seq();
        let len = seq.len() as u64;
        stats
            .seq_digests
            .insert(seqname.clone(), sequence_digest(&seq));
        for t in ts {
            if t.start >= len {
                stats.dropped_missing_seq += 1;
                continue;
            }
            let stop = t.end.min(len);
            if stop < t.end {
                stats.clipped += 1;
            }
            let slice = &seq[t.start as usize..stop as usize];
            buf.clear();
            match t.strand {
                Strand::Plus => buf.extend(slice.iter().map(u8::to_ascii_uppercase)),
                Strand::Minus => buf.extend(
                    slice
                        .iter()
                        .rev()
                        .map(|b| complement(b.to_ascii_uppercase())),
                ),
            }
            let acgt = buf
                .iter()
                .filter(|&&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
                .count();
            if acgt * 2 < buf.len() {
                stats.dropped_low_acgt.push(t.name.clone());
                continue;
            }
            w.write_all(b">")?;
            w.write_all(t.name.as_bytes())?;
            w.write_all(b"\n")?;
            w.write_all(&buf)?;
            w.write_all(b"\n")?;
            stats.written += 1;
            stats.targets.push(WrittenTarget {
                name: t.name.clone(),
                seqname: seqname.clone(),
                start: t.start,
                end: stop,
                minus: t.strand == Strand::Minus,
            });
        }
    }
    w.flush()?;
    let mut missing: Vec<&str> = by_seq
        .iter()
        .filter(|(s, _)| !seen_seqs.contains(**s))
        .map(|(s, ts)| {
            stats.dropped_missing_seq += ts.len();
            *s
        })
        .collect();
    missing.sort_unstable();
    stats.missing_seqnames = missing.into_iter().map(str::to_string).collect();
    if stats.written == 0 && !targets.is_empty() {
        bail!(
            "none of the {} unspliced targets could be extracted from {}: do the GTF \
             sequence names match the genome FASTA headers?",
            targets.len(),
            genome.display()
        );
    }
    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Synthetic GTF covering the awkward cases: two transcripts with a
    /// retained intron, a single-exon gene, genes overlapping on opposite
    /// strands, a mitochondrial gene, and a GENCODE `_PAR_Y` gene.
    pub(crate) const GTF: &str = "\
##description: synthetic
chr1\tt\tgene\t101\t1000\t.\t+\t.\tgene_id \"G1\";
chr1\tt\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1a\";
chr1\tt\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1a\";
chr1\tt\texon\t601\t1000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1a\";
chr1\tt\texon\t101\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1ri\"; transcript_type \"retained_intron\";
chr1\tt\texon\t601\t700\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1ri\";
chr1\tt\texon\t1201\t1500\t.\t+\t.\tgene_id \"G2\"; transcript_id \"T2\";
chr1\tt\texon\t901\t950\t.\t-\t.\tgene_id \"G3\"; transcript_id \"T3\";
chr1\tt\texon\t1101\t1300\t.\t-\t.\tgene_id \"G3\"; transcript_id \"T3\";
chrM\tt\texon\t1\t500\t.\t+\t.\tgene_id \"MT1\"; transcript_id \"TM1\";
chrY\tt\texon\t11\t60\t.\t+\t.\tgene_id \"G4_PAR_Y\"; transcript_id \"T4_PAR_Y\";
chrY\tt\texon\t121\t200\t.\t+\t.\tgene_id \"G4_PAR_Y\"; transcript_id \"T4_PAR_Y\";
chr1\tt\texon\t5\t9\t.\t.\t.\tgene_id \"GX\"; transcript_id \"TX\";
";

    fn gene<'a>(a: &'a Annotation, id: &str) -> &'a GeneModel {
        a.genes.iter().find(|g| g.gene_id == id).unwrap()
    }

    #[test]
    fn parses_exons_into_gene_models() {
        let a = parse_gtf(GTF.as_bytes()).unwrap();
        let ids: Vec<&str> = a.genes.iter().map(|g| g.gene_id.as_str()).collect();
        assert_eq!(ids, ["G1", "G2", "G3", "MT1", "G4_PAR_Y"]);
        let g1 = gene(&a, "G1");
        assert_eq!(g1.strand, Strand::Plus);
        assert_eq!(g1.transcripts.len(), 2);
        // 1-based inclusive 101..200 -> 0-based half-open 100..200
        assert_eq!(g1.transcripts[0], vec![(100, 200), (300, 400), (600, 1000)]);
        assert_eq!(gene(&a, "G3").strand, Strand::Minus);
        assert_eq!(a.unstranded_exons, 1);
        assert_eq!(a.gene_of("T1ri"), Some("G1"));
        assert_eq!(a.gene_of("T4_PAR_Y"), Some("G4_PAR_Y"));
        assert_eq!(a.gene_of("nope"), None);
    }

    #[test]
    fn transcript_version_attribute_is_honoured() {
        let gtf = "1\te\texon\t1\t10\t.\t+\t.\tgene_id \"ENSG1\"; gene_version \"3\"; \
                   transcript_id \"ENST1\"; transcript_version \"2\";\n";
        let a = parse_gtf(gtf.as_bytes()).unwrap();
        assert_eq!(a.gene_of("ENST1.2"), Some("ENSG1"));
        assert_eq!(a.gene_of("ENST1"), Some("ENSG1"));
        // a version the GTF does not know still matches with the version dropped
        assert_eq!(a.gene_of("ENST1.9"), Some("ENSG1"));
    }

    #[test]
    fn introns_are_the_union_over_transcripts() {
        let a = parse_gtf(GTF.as_bytes()).unwrap();
        // T1a: introns 200..300 and 400..600; T1ri: 400..600 (its first exon
        // retains 200..300). The union keeps the retained intron.
        assert_eq!(merged_introns(gene(&a, "G1")), vec![(200, 300), (400, 600)]);
        assert_eq!(gene_body(gene(&a, "G1")), Some((100, 1000)));
        assert_eq!(merged_introns(gene(&a, "G3")), vec![(950, 1100)]);
        assert_eq!(merged_introns(gene(&a, "G4_PAR_Y")), vec![(60, 120)]);
    }

    #[test]
    fn single_exon_genes_have_no_unspliced_form() {
        let a = parse_gtf(GTF.as_bytes()).unwrap();
        for id in ["G2", "MT1"] {
            assert!(merged_introns(gene(&a, id)).is_empty(), "{id}");
            assert_eq!(gene_body(gene(&a, id)), None, "{id}");
        }
    }

    #[test]
    fn merge_handles_overlap_and_abutment() {
        assert_eq!(
            merge_intervals(vec![(10, 20), (0, 5), (5, 8), (15, 30), (40, 41)]),
            vec![(0, 8), (10, 30), (40, 41)]
        );
        assert!(merge_intervals(Vec::new()).is_empty());
    }

    #[test]
    fn overlapping_exons_leave_no_intron() {
        let g = GeneModel {
            gene_id: "g".into(),
            seqname: "c".into(),
            strand: Strand::Plus,
            transcripts: vec![vec![(0, 50), (40, 100), (100, 120), (150, 200)]],
        };
        assert_eq!(merged_introns(&g), vec![(120, 150)]);
    }

    #[test]
    fn rejects_annotation_without_exons() {
        let err =
            parse_gtf("chr1\tt\tgene\t1\t10\t.\t+\t.\tgene_id \"G\";\n".as_bytes()).unwrap_err();
        assert!(
            err.to_string().contains("no stranded exon records"),
            "{err}"
        );
    }

    #[test]
    fn rejects_bad_coordinates() {
        let gtf = "chr1\tt\texon\t0\t10\t.\t+\t.\tgene_id \"G\"; transcript_id \"T\";\n";
        assert!(parse_gtf(gtf.as_bytes()).is_err());
    }

    fn targets(mode: UnsplicedMode, flank: u64) -> Vec<(String, String, u64, u64)> {
        let a = parse_gtf(GTF.as_bytes()).unwrap();
        unspliced_targets(&a, mode, flank)
            .into_iter()
            .map(|t| (t.name, t.seqname, t.start, t.end))
            .collect()
    }

    fn row(n: &str, c: &str, s: u64, e: u64) -> (String, String, u64, u64) {
        (n.to_string(), c.to_string(), s, e)
    }

    #[test]
    fn intron_targets_without_flank() {
        assert_eq!(
            targets(UnsplicedMode::Intron, 0),
            vec![
                row("G1-I", "chr1", 200, 300),
                row("G1-I1", "chr1", 400, 600),
                row("G3-I", "chr1", 950, 1100),
                row("G4_PAR_Y-I", "chrY", 60, 120),
            ]
        );
    }

    #[test]
    fn flanks_extend_clip_at_zero_and_merge_across_short_exons() {
        // flank 50: G1's introns become 150..350 and 350..650, which meet
        // across the 100-base exon 300..400 and merge; G4's flank would run
        // before position 10 only at flank > 60, so use 70 to see the clip.
        assert_eq!(
            targets(UnsplicedMode::Intron, 50),
            vec![
                row("G1-I", "chr1", 150, 650),
                row("G3-I", "chr1", 900, 1150),
                row("G4_PAR_Y-I", "chrY", 10, 170),
            ]
        );
        let t = targets(UnsplicedMode::Intron, 70);
        assert_eq!(t[2], row("G4_PAR_Y-I", "chrY", 0, 190));
    }

    #[test]
    fn premrna_targets_are_gene_bodies_of_intron_bearing_genes() {
        assert_eq!(
            targets(UnsplicedMode::Premrna, 1000),
            vec![
                row("G1-I", "chr1", 100, 1000),
                row("G3-I", "chr1", 900, 1300),
                row("G4_PAR_Y-I", "chrY", 10, 200),
            ]
        );
        assert!(targets(UnsplicedMode::None, 10).is_empty());
    }

    #[test]
    fn a_gene_on_two_sequences_gets_distinct_names() {
        let gtf = "\
chrX\tt\texon\t1\t10\t.\t+\t.\tgene_id \"P\"; transcript_id \"T\";
chrX\tt\texon\t21\t30\t.\t+\t.\tgene_id \"P\"; transcript_id \"T\";
chrY\tt\texon\t1\t10\t.\t+\t.\tgene_id \"P\"; transcript_id \"T\";
chrY\tt\texon\t21\t30\t.\t+\t.\tgene_id \"P\"; transcript_id \"T\";
";
        let a = parse_gtf(gtf.as_bytes()).unwrap();
        let names: Vec<String> = unspliced_targets(&a, UnsplicedMode::Intron, 0)
            .into_iter()
            .map(|t| format!("{}@{}", t.name, t.seqname))
            .collect();
        assert_eq!(names, ["P-I@chrX", "P-I1@chrY"]);
    }

    #[test]
    fn mode_round_trips_through_its_name() {
        for m in [
            UnsplicedMode::None,
            UnsplicedMode::Intron,
            UnsplicedMode::Premrna,
        ] {
            assert_eq!(m.as_str().parse::<UnsplicedMode>().unwrap(), m);
        }
        assert!("splici".parse::<UnsplicedMode>().is_err());
    }

    #[test]
    fn extracts_stranded_uppercase_sequences_in_genome_order() {
        let dir = tempfile::tempdir().unwrap();
        let genome = dir.path().join("g.fa");
        // chrY first in the FASTA: output follows genome order, not GTF order.
        // chrY's second half is hard-masked (N), like a PAR copy.
        let chr1: String = (0..1600).map(|i| b"ACGTTGCA"[i % 8] as char).collect();
        let chr_y = format!("{}{}", "acgt".repeat(10), "N".repeat(160));
        std::fs::write(
            &genome,
            format!(
                ">chrY  Y chromosome\n{chr_y}\n>chr1 1\n{}\n{}\n>chrM\nACGT\n",
                &chr1[..800],
                &chr1[800..]
            ),
        )
        .unwrap();
        let a = parse_gtf(GTF.as_bytes()).unwrap();
        let mut ts = unspliced_targets(&a, UnsplicedMode::Intron, 0);
        // a target on a sequence the genome lacks, and one running off chr1's end
        ts.push(UnsplicedTarget {
            name: "Z-I".into(),
            gene_id: "Z".into(),
            seqname: "chrZ".into(),
            strand: Strand::Plus,
            start: 0,
            end: 10,
        });
        ts.push(UnsplicedTarget {
            name: "E-I".into(),
            gene_id: "E".into(),
            seqname: "chr1".into(),
            strand: Strand::Plus,
            start: 1590,
            end: 1700,
        });
        let out = dir.path().join("u.fa");
        let st = write_unspliced_fasta(&genome, &ts, &out).unwrap();
        // G4_PAR_Y-I (60..120 on chrY) is all N: dropped.
        assert_eq!(st.dropped_low_acgt, ["G4_PAR_Y-I"]);
        assert_eq!(st.missing_seqnames, ["chrZ"]);
        assert_eq!(st.dropped_missing_seq, 1);
        assert_eq!(st.clipped, 1);
        assert_eq!(st.written, 4);
        let text = std::fs::read_to_string(&out).unwrap();
        let recs: Vec<(&str, &str)> = text
            .lines()
            .collect::<Vec<_>>()
            .chunks(2)
            .map(|c| (c[0], c[1]))
            .collect();
        let names: Vec<&str> = recs.iter().map(|r| r.0).collect();
        assert_eq!(names, [">G1-I", ">G1-I1", ">G3-I", ">E-I"]);
        assert_eq!(recs[0].1, &chr1[200..300]);
        assert_eq!(recs[3].1, &chr1[1590..]);
        // G3 is on the minus strand: reverse complement of 950..1100
        let rc: String = chr1[950..1100]
            .bytes()
            .rev()
            .map(|b| complement(b) as char)
            .collect();
        assert_eq!(recs[2].1, rc);
    }

    #[test]
    fn extraction_fails_when_no_sequence_name_matches() {
        let dir = tempfile::tempdir().unwrap();
        let genome = dir.path().join("g.fa");
        std::fs::write(&genome, ">1\nACGTACGTACGT\n").unwrap();
        let a = parse_gtf(GTF.as_bytes()).unwrap();
        let ts = unspliced_targets(&a, UnsplicedMode::Intron, 0);
        let err = write_unspliced_fasta(&genome, &ts, &dir.path().join("u.fa")).unwrap_err();
        assert!(err.to_string().contains("sequence names"), "{err}");
    }
}
