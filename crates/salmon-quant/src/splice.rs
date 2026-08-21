//! Projection of transcriptome alignments into spliced genome alignments.
//!
//! # Why this exists
//!
//! salmon aligns reads to transcripts, and a transcript has no introns: a read
//! sitting across an exon junction is contiguous in transcript coordinates, so
//! the CIGAR that describes it never contains an `N`. That is the correct
//! description of the alignment salmon actually made, but it is not the
//! description a genome browser, a junction counter, or anything else built
//! around STAR's output can read.
//!
//! This module converts one into the other. Given an annotation, each transcript
//! is a list of exon blocks on a chromosome; a transcript-coordinate alignment is
//! cut at those exon boundaries and the gaps between them reappear as `N`
//! operations. The result is a spliced genome alignment of exactly the kind a
//! genome aligner would have produced, derived from the transcript evidence
//! salmon already has.
//!
//! # Reverse-strand transcripts
//!
//! Half of all transcripts run antisense to the genome, and for those the
//! transcript coordinate axis points the opposite way from the genomic one. A
//! read at transcript position 0 is then at the *rightmost* genomic base of the
//! transcript, and the alignment has to be turned around: the CIGAR is reversed,
//! the leftmost genomic position becomes the one the last transcript base landed
//! on, and the record's strand flips (its stored sequence and qualities are
//! reverse-complemented once more). [`Projection::flipped`] reports that flip so
//! the record writer can apply it to `SEQ`, `QUAL` and the `FLAG` bits together.
//!
//! # `MD` after projection
//!
//! `MD` describes the reference bases along the alignment, so it has to follow
//! the alignment through projection rather than be carried across it. Two things
//! change: a reverse-strand transcript's string is reverse-complemented with
//! everything else, and a deletion the projection cuts in two becomes two
//! deletions, which `MD` spells with an empty match between them. Both are done
//! by re-deriving the string from the projected operations, which gets the
//! grouping right by construction.
//!
//! # What is deliberately not done here
//!
//! Two transcripts of the same gene often imply the *same* genomic alignment for
//! a read. This module projects each placement independently and leaves the
//! duplicates in place, tagged as secondary exactly as they were in transcript
//! space, rather than silently collapsing them: which of them to keep is a
//! quantification decision, and salmon's answer to it is `quant.sf`, not the
//! BAM.

use std::collections::HashMap;
use std::path::Path;

use anyhow::Context as _;
use salmon_map::realize::{push_op, CigarOp, CigarOpKind};

/// One exon, as a genomic interval plus where it starts along the transcript.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ExonBlock {
    /// 0-based genomic start of the exon.
    genome_start: u32,
    /// Transcript coordinate of this exon's first base.
    transcript_start: u32,
    len: u32,
}

/// A transcript's exon structure, in transcript order.
///
/// "Transcript order" means exon 0 contains transcript coordinate 0, which for a
/// reverse-strand transcript is the exon with the *highest* genomic coordinates.
#[derive(Debug, Clone)]
struct TranscriptModel {
    chromosome: usize,
    forward: bool,
    exons: Vec<ExonBlock>,
}

impl TranscriptModel {
    fn length(&self) -> u32 {
        self.exons.last().map_or(0, |e| e.transcript_start + e.len)
    }

    /// The exon containing transcript coordinate `pos`, or the last one when
    /// `pos` runs past the transcript end.
    fn exon_at(&self, pos: u32) -> usize {
        match self
            .exons
            .binary_search_by(|e| e.transcript_start.cmp(&pos))
        {
            Ok(index) => index,
            // `binary_search_by` reports where `pos` would be inserted, so the
            // exon containing it is the one before that.
            Err(index) => index.saturating_sub(1),
        }
    }

    /// The genomic coordinate of transcript position `pos`.
    ///
    /// On the forward strand this counts up from the exon's genomic start; on the
    /// reverse strand it counts down from the exon's genomic end, because the two
    /// axes run opposite ways.
    fn genome_position(&self, pos: u32) -> u32 {
        let exon = &self.exons[self.exon_at(pos)];
        let offset = pos.saturating_sub(exon.transcript_start);
        if self.forward {
            exon.genome_start + offset
        } else {
            (exon.genome_start + exon.len).saturating_sub(offset + 1)
        }
    }
}

/// The result of projecting one placement onto the genome.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Projection {
    /// Index into [`GenomeProjector::chromosomes`].
    pub chromosome: usize,
    /// Leftmost genomic base the alignment consumes, 0-based.
    pub position: i32,
    /// Whether the alignment had to be turned around, because the transcript is
    /// on the reverse strand. The record's strand bits, `SEQ` and `QUAL` all have
    /// to be flipped to match.
    pub flipped: bool,
}

/// Transcript exon structures for every reference in the salmon index, plus the
/// chromosome table the projected records refer to.
pub struct GenomeProjector {
    /// Placements dropped because their reference has no exon structure. Counted
    /// so the run can say how much of the mapping output the annotation cost,
    /// instead of quietly producing a file that disagrees with `quant.sf`.
    dropped: std::sync::atomic::AtomicU64,
    /// `(name, length)` per chromosome, in the order they appear in `@SQ`.
    chromosomes: Vec<(String, u64)>,
    /// One entry per salmon reference id; `None` for a reference the annotation
    /// does not describe (a decoy, or a transcript absent from the GTF).
    models: Vec<Option<TranscriptModel>>,
}

impl GenomeProjector {
    /// Build the projector for `reference_names` (the salmon index's reference
    /// table, in index order) from an annotation and a genome FASTA.
    ///
    /// The genome is needed only for the chromosome *lengths* that `@SQ` has to
    /// declare: an annotation says where exons are but never how long the
    /// chromosome they sit on is, and a header that understates a chromosome
    /// length produces a BAM that readers reject. A FASTA index (`.fai`) is used
    /// when one is present, so the common case costs one small file read rather
    /// than a pass over the genome.
    pub fn build(
        reference_names: &[&str],
        annotation: &Path,
        genome: &Path,
    ) -> anyhow::Result<Self> {
        let chromosome_lengths = chromosome_lengths(genome)?;
        let transcripts = bramble_rs::annotation::load_transcripts(annotation)
            .with_context(|| format!("reading annotation {}", annotation.display()))?;

        // Only the chromosomes the annotation actually places transcripts on go
        // into the header, in the genome's own order, so the reference table stays
        // as small as the output needs while keeping a stable, familiar ordering.
        let mut used: HashMap<&str, usize> = HashMap::new();
        let mut chromosomes: Vec<(String, u64)> = Vec::new();
        for (name, length) in &chromosome_lengths {
            if transcripts.iter().any(|t| t.seqname == *name) {
                used.insert(name.as_str(), chromosomes.len());
                chromosomes.push((name.clone(), *length));
            }
        }

        let by_id: HashMap<&str, &bramble_rs::annotation::Transcript> =
            transcripts.iter().map(|t| (t.id.as_str(), t)).collect();
        let mut missing_chromosome = 0usize;
        let models = reference_names
            .iter()
            .map(|name| {
                let transcript = by_id.get(name)?;
                let Some(&chromosome) = used.get(transcript.seqname.as_str()) else {
                    missing_chromosome += 1;
                    return None;
                };
                Some(model_for(transcript, chromosome))
            })
            .collect::<Vec<_>>();

        let described = models.iter().filter(|m| m.is_some()).count();
        anyhow::ensure!(
            described > 0,
            "no transcript in the index appears in {}: check that the annotation matches the \
             reference the index was built from",
            annotation.display()
        );
        if missing_chromosome > 0 {
            tracing::warn!(
                transcripts = missing_chromosome,
                "annotation places transcripts on sequences absent from the genome FASTA; those \
                 placements will not be projected"
            );
        }
        tracing::info!(
            transcripts = described,
            total = reference_names.len(),
            chromosomes = chromosomes.len(),
            "genome projection ready"
        );
        Ok(Self {
            dropped: std::sync::atomic::AtomicU64::new(0),
            chromosomes,
            models,
        })
    }

    /// How many placements were dropped for want of an exon structure.
    pub fn dropped(&self) -> u64 {
        self.dropped.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// The `(name, length)` table projected records index into.
    pub fn chromosomes(&self) -> &[(String, u64)] {
        &self.chromosomes
    }

    /// Whether the transcript at `tid` runs along the genome's forward strand,
    /// or `None` if the annotation does not describe it. This is what `XS`
    /// reports for a spliced alignment.
    pub fn strand(&self, tid: usize) -> Option<bool> {
        Some(self.models.get(tid)?.as_ref()?.forward)
    }

    /// Whether reference `tid` can be projected at all.
    pub fn describes(&self, tid: usize) -> bool {
        self.models.get(tid).is_some_and(|m| m.is_some())
    }

    /// Project a transcript-coordinate alignment onto the genome, writing the
    /// genomic CIGAR into `out`.
    ///
    /// Returns `None` when the reference has no exon structure (a decoy, or a
    /// transcript the annotation does not mention), in which case the caller has
    /// nothing genomic to write for this placement.
    #[allow(clippy::too_many_arguments)]
    pub fn project(
        &self,
        tid: usize,
        transcript_position: i32,
        ops: &[CigarOp],
        md: Option<(&str, &mut Vec<MdBase>, &mut String)>,
        out: &mut Vec<CigarOp>,
    ) -> Option<Projection> {
        out.clear();
        let Some(Some(model)) = self.models.get(tid) else {
            self.dropped
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            return None;
        };
        if model.exons.is_empty() || transcript_position.max(0) as u32 >= model.length() {
            self.dropped
                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            return None;
        }
        let start = transcript_position.max(0) as u32;

        // Walk the transcript CIGAR, cutting every reference-consuming operation
        // at exon boundaries and filling the gap between consecutive exons with an
        // N. Operations that do not consume the reference pass straight through.
        let mut position = start;
        let mut exon_index = model.exon_at(start);
        for op in ops {
            if !op.kind.consumes_reference() {
                push_op(out, op.kind, op.len);
                continue;
            }
            let mut remaining = op.len;
            while remaining > 0 {
                // A junction can fall between two operations just as easily as
                // inside one: an operation ending exactly on an exon boundary
                // leaves the next one starting in the following exon. Checking
                // only mid-operation let those alignments cross an intron with no
                // `N` at all, so the bases after the boundary silently claimed
                // the wrong reference.
                let current = model.exon_at(position);
                while exon_index < current {
                    let exon = &model.exons[exon_index];
                    let next = model.exons.get(exon_index + 1)?;
                    push_op(out, CigarOpKind::Skip, intron_length(model, exon, next));
                    exon_index += 1;
                }
                let exon = &model.exons[exon_index];
                let exon_end = exon.transcript_start + exon.len;
                if position >= exon_end {
                    // Past the final exon: the rest of the operation has no
                    // genomic home, so the alignment cannot be projected.
                    return None;
                }
                let take = remaining.min(exon_end - position);
                push_op(out, op.kind, take);
                position += take;
                remaining -= take;
            }
        }
        if out.is_empty() {
            return None;
        }

        // The genomic extremes of the alignment. On the forward strand the first
        // transcript base is leftmost; on the reverse strand the last one is.
        let last = position.saturating_sub(1).max(start);
        let position = if model.forward {
            model.genome_position(start)
        } else {
            model.genome_position(last)
        } as i32;
        // MD is re-derived here, while the operations are still in transcript
        // order and line up with the transcript's reference bases one for one.
        if let Some((transcript_md, scratch, md_out)) = md {
            reproject_md(transcript_md, out, scratch, md_out);
            if !model.forward {
                let forward_md = std::mem::take(md_out);
                reverse_complement_md(&forward_md, md_out);
            }
        }
        if !model.forward {
            out.reverse();
        }
        Some(Projection {
            chromosome: model.chromosome,
            position,
            flipped: !model.forward,
        })
    }
}

/// One reference base, as `MD` describes it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MdBase {
    /// The read agrees with the reference here.
    Match,
    /// The read differs; the reference base is carried.
    Mismatch(u8),
    /// The reference base is absent from the read.
    Deleted(u8),
}

/// Expand an `MD` string into one entry per reference base.
///
/// `MD` is written as runs, which is compact but positional: nothing in it says
/// which CIGAR operation a given base belongs to. Expanding it makes the bases
/// addressable, so they can be re-emitted against a different operation list.
fn expand_md(md: &str, out: &mut Vec<MdBase>) {
    out.clear();
    let mut digits = 0usize;
    let mut chars = md.chars().peekable();
    while let Some(c) = chars.next() {
        if let Some(digit) = c.to_digit(10) {
            digits = digits * 10 + digit as usize;
            continue;
        }
        out.extend(std::iter::repeat_n(MdBase::Match, digits));
        digits = 0;
        if c == '^' {
            while chars.peek().is_some_and(|n| n.is_ascii_alphabetic()) {
                out.push(MdBase::Deleted(chars.next().unwrap() as u8));
            }
        } else {
            out.push(MdBase::Mismatch(c as u8));
        }
    }
    out.extend(std::iter::repeat_n(MdBase::Match, digits));
}

/// Re-emit `MD` against the projected operations.
///
/// The reference bases are unchanged by projection, but their grouping is not: a
/// deletion the projection cuts in two is two deletions, and `MD` spells those
/// `^G0^G`, not `^GG`. Deriving the string from the transcript's operations and
/// carrying it across produced the merged form, which no longer matches the CIGAR
/// beside it. Walking the projected operations instead gets the grouping right by
/// construction, and the empty match between the halves falls out of the normal
/// rule that a run length separates every pair of tokens.
///
/// `ops` must be in transcript order, before any reversal for strand.
fn reproject_md(transcript_md: &str, ops: &[CigarOp], scratch: &mut Vec<MdBase>, out: &mut String) {
    use std::fmt::Write as _;
    expand_md(transcript_md, scratch);
    out.clear();
    let mut next = 0usize;
    let mut run = 0u32;
    for op in ops {
        match op.kind {
            CigarOpKind::Match => {
                for _ in 0..op.len {
                    match scratch.get(next) {
                        Some(MdBase::Mismatch(base)) => {
                            let _ = write!(out, "{run}");
                            out.push(*base as char);
                            run = 0;
                        }
                        _ => run += 1,
                    }
                    next += 1;
                }
            }
            CigarOpKind::Deletion => {
                let _ = write!(out, "{run}");
                run = 0;
                out.push('^');
                for _ in 0..op.len {
                    if let Some(MdBase::Deleted(base)) = scratch.get(next) {
                        out.push(*base as char);
                    }
                    next += 1;
                }
            }
            // An intron consumes genome but no transcript, so it takes no
            // reference bases from the stream. It still ends the current token,
            // which is what a following deletion needs to start a fresh one.
            CigarOpKind::Skip => {}
            CigarOpKind::Insertion | CigarOpKind::SoftClip => {}
        }
    }
    let _ = write!(out, "{run}");
}

/// Rewrite an `MD` string for an alignment that has been turned around.
///
/// `MD` describes the reference bases along the alignment, so when a reverse
/// strand transcript's alignment is flipped into genome order the string has to
/// be flipped with it: the tokens reverse, and every base in them is
/// complemented. Leaving it alone produces an `MD` that disagrees with the very
/// `SEQ` and CIGAR it sits beside, which is a mistake no reader can catch and
/// `samtools calmd` catches immediately.
///
/// A deletion's bases reverse within the deletion as well, since they too are
/// read along the reference.
pub fn reverse_complement_md(md: &str, out: &mut String) {
    out.clear();
    // Split into tokens: match run lengths, single mismatched bases, and
    // `^`-prefixed deleted stretches.
    let mut tokens: Vec<String> = Vec::new();
    let mut digits = String::new();
    let mut chars = md.chars().peekable();
    while let Some(c) = chars.next() {
        if c.is_ascii_digit() {
            digits.push(c);
            continue;
        }
        if !digits.is_empty() {
            tokens.push(std::mem::take(&mut digits));
        }
        if c == '^' {
            let mut deleted = String::from("^");
            while chars.peek().is_some_and(|n| n.is_ascii_alphabetic()) {
                deleted.push(complement(chars.next().unwrap() as u8) as char);
            }
            // The deleted bases were pushed in transcript order; reverse them.
            let reversed: String = deleted[1..].chars().rev().collect();
            tokens.push(format!("^{reversed}"));
        } else {
            tokens.push((complement(c as u8) as char).to_string());
        }
    }
    if !digits.is_empty() {
        tokens.push(digits);
    }
    for token in tokens.iter().rev() {
        out.push_str(token);
    }
}

/// Complement one base, for [`reverse_complement_md`].
fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        other => other,
    }
}

/// The genomic gap between two consecutive exons, whichever way the transcript
/// runs.
fn intron_length(model: &TranscriptModel, exon: &ExonBlock, next: &ExonBlock) -> u32 {
    if model.forward {
        next.genome_start
            .saturating_sub(exon.genome_start + exon.len)
    } else {
        exon.genome_start
            .saturating_sub(next.genome_start + next.len)
    }
}

/// Turn a bramble transcript into exon blocks in transcript order, which for a
/// reverse-strand transcript means genomic order reversed.
fn model_for(
    transcript: &bramble_rs::annotation::Transcript,
    chromosome: usize,
) -> TranscriptModel {
    let forward = transcript.strand != '-';
    // bramble stores a GTF exon as `[start_1, end_1 + 1)`: the length is right,
    // but the interval is still anchored at the 1-based start, so every
    // coordinate is one too high. Subtract the one here, once, rather than
    // letting a systematic off-by-one into every projected position. Exon
    // lengths are unaffected, which is why the introns between them come out
    // correct either way and only the reported positions gave it away.
    let mut exons: Vec<(u32, u32)> = transcript
        .exons
        .iter()
        .map(|e| (e.start.saturating_sub(1), e.end.saturating_sub(e.start)))
        .filter(|(_, len)| *len > 0)
        .collect();
    exons.sort_unstable();
    if !forward {
        exons.reverse();
    }
    let mut transcript_start = 0;
    let exons = exons
        .into_iter()
        .map(|(genome_start, len)| {
            let block = ExonBlock {
                genome_start,
                transcript_start,
                len,
            };
            transcript_start += len;
            block
        })
        .collect();
    TranscriptModel {
        chromosome,
        forward,
        exons,
    }
}

/// Chromosome names and lengths, from a FASTA index when there is one and from
/// the FASTA itself otherwise.
fn chromosome_lengths(genome: &Path) -> anyhow::Result<Vec<(String, u64)>> {
    let index_path = genome.with_extension(format!(
        "{}.fai",
        genome
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or_default()
    ));
    if index_path.exists() {
        let text = std::fs::read_to_string(&index_path)
            .with_context(|| format!("reading FASTA index {}", index_path.display()))?;
        let lengths = parse_fai(&text);
        if !lengths.is_empty() {
            return Ok(lengths);
        }
    }
    scan_fasta_lengths(genome)
}

/// Parse the first two columns of a `.fai`: name and length.
fn parse_fai(text: &str) -> Vec<(String, u64)> {
    text.lines()
        .filter_map(|line| {
            let mut fields = line.split('\t');
            let name = fields.next()?;
            let length = fields.next()?.parse().ok()?;
            (!name.is_empty()).then(|| (name.to_string(), length))
        })
        .collect()
}

/// Measure every sequence in a FASTA by streaming it. Only reached when no `.fai`
/// is available; a genome is large, so this counts bytes rather than holding
/// sequences.
fn scan_fasta_lengths(genome: &Path) -> anyhow::Result<Vec<(String, u64)>> {
    use std::io::BufRead as _;
    let file = std::fs::File::open(genome)
        .with_context(|| format!("reading genome FASTA {}", genome.display()))?;
    let reader = std::io::BufReader::with_capacity(1 << 20, file);
    let mut lengths: Vec<(String, u64)> = Vec::new();
    for line in reader.lines() {
        let line = line.with_context(|| format!("reading {}", genome.display()))?;
        if let Some(header) = line.strip_prefix('>') {
            let name = header.split_whitespace().next().unwrap_or_default();
            lengths.push((name.to_string(), 0));
        } else if let Some(last) = lengths.last_mut() {
            last.1 += line.trim_end().len() as u64;
        }
    }
    anyhow::ensure!(
        !lengths.is_empty(),
        "genome FASTA {} contains no sequences",
        genome.display()
    );
    Ok(lengths)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A three-exon transcript on the forward strand: exons at 100..150,
    /// 300..350 and 500..600, so transcript coordinates 0..50, 50..100, 100..200.
    fn forward_model() -> TranscriptModel {
        TranscriptModel {
            chromosome: 0,
            forward: true,
            exons: vec![
                ExonBlock {
                    genome_start: 100,
                    transcript_start: 0,
                    len: 50,
                },
                ExonBlock {
                    genome_start: 300,
                    transcript_start: 50,
                    len: 50,
                },
                ExonBlock {
                    genome_start: 500,
                    transcript_start: 100,
                    len: 100,
                },
            ],
        }
    }

    /// The same exons, on a transcript transcribed the other way, so transcript
    /// coordinate 0 is the last genomic base of the 500..600 exon.
    fn reverse_model() -> TranscriptModel {
        TranscriptModel {
            chromosome: 0,
            forward: false,
            exons: vec![
                ExonBlock {
                    genome_start: 500,
                    transcript_start: 0,
                    len: 100,
                },
                ExonBlock {
                    genome_start: 300,
                    transcript_start: 100,
                    len: 50,
                },
                ExonBlock {
                    genome_start: 100,
                    transcript_start: 150,
                    len: 50,
                },
            ],
        }
    }

    fn projector(model: TranscriptModel) -> GenomeProjector {
        GenomeProjector {
            dropped: std::sync::atomic::AtomicU64::new(0),
            chromosomes: vec![("chr1".into(), 1_000_000)],
            models: vec![Some(model)],
        }
    }

    fn cigar(ops: &[CigarOp]) -> String {
        ops.iter()
            .map(|op| format!("{}{}", op.len, op.kind.letter()))
            .collect()
    }

    fn matched(len: u32) -> Vec<CigarOp> {
        vec![CigarOp::new(CigarOpKind::Match, len)]
    }

    /// A read wholly inside one exon needs no intron: it projects to a plain
    /// match at the corresponding genomic offset.
    #[test]
    fn a_read_inside_one_exon_gets_no_skip() {
        let projector = projector(forward_model());
        let mut ops = Vec::new();
        let projection = projector
            .project(0, 10, &matched(20), None, &mut ops)
            .unwrap();
        assert_eq!(cigar(&ops), "20M");
        assert_eq!(projection.position, 110);
        assert!(!projection.flipped);
    }

    /// A read across a junction is where the whole module earns its keep: the
    /// intron between the two exons has to appear as an N of exactly its genomic
    /// length.
    #[test]
    fn a_junction_read_gains_an_n_operation() {
        let projector = projector(forward_model());
        let mut ops = Vec::new();
        // Transcript 30..70 spans the first junction 20 bases in.
        let projection = projector
            .project(0, 30, &matched(40), None, &mut ops)
            .unwrap();
        // 20 bases to the end of exon 1, an intron of 300 - 150 = 150, then 20 more.
        assert_eq!(cigar(&ops), "20M150N20M");
        assert_eq!(projection.position, 130);
    }

    /// A read crossing two junctions produces two introns, and the short middle
    /// exon appears in full between them.
    #[test]
    fn a_read_over_two_junctions_gains_two_skips() {
        let projector = projector(forward_model());
        let mut ops = Vec::new();
        let projection = projector
            .project(0, 40, &matched(70), None, &mut ops)
            .unwrap();
        assert_eq!(cigar(&ops), "10M150N50M150N10M");
        assert_eq!(projection.position, 140);
    }

    /// A junction falling exactly between two operations is still a junction. A
    /// deletion ending on an exon boundary used to leave the following match
    /// crossing the intron with no `N`, so it claimed reference bases from the
    /// wrong side of it.
    #[test]
    fn a_junction_between_two_operations_still_gets_a_skip() {
        let projector = projector(forward_model());
        let mut ops = Vec::new();
        // Exon 1 is transcript 0..50. Align 47 bases, delete 3 (ending exactly on
        // the boundary), then match 10 more, which belong in the next exon.
        let transcript_cigar = vec![
            CigarOp::new(CigarOpKind::Match, 47),
            CigarOp::new(CigarOpKind::Deletion, 3),
            CigarOp::new(CigarOpKind::Match, 10),
        ];
        let projection = projector
            .project(0, 0, &transcript_cigar, None, &mut ops)
            .unwrap();
        assert_eq!(cigar(&ops), "47M3D150N10M");
        assert_eq!(projection.position, 100);
    }

    /// Soft clips and insertions do not consume the reference, so they must pass
    /// through the projection untouched rather than being cut at an exon edge.
    #[test]
    fn clips_and_insertions_survive_projection() {
        let projector = projector(forward_model());
        let mut ops = Vec::new();
        let transcript_cigar = vec![
            CigarOp::new(CigarOpKind::SoftClip, 5),
            CigarOp::new(CigarOpKind::Match, 15),
            CigarOp::new(CigarOpKind::Insertion, 3),
            CigarOp::new(CigarOpKind::Match, 15),
        ];
        let projection = projector
            .project(0, 40, &transcript_cigar, None, &mut ops)
            .unwrap();
        assert_eq!(cigar(&ops), "5S10M150N5M3I15M");
        assert_eq!(projection.position, 140);
    }

    /// A reverse-strand transcript runs against the genomic axis, so the CIGAR
    /// comes out reversed, the reported position is the alignment's leftmost
    /// genomic base, and the record has to be flipped to match.
    #[test]
    fn reverse_strand_transcripts_are_turned_around() {
        let projector = projector(reverse_model());
        let mut ops = Vec::new();
        // Transcript 90..130 crosses the junction between the 500..600 exon and
        // the 300..350 one.
        let projection = projector
            .project(0, 90, &matched(40), None, &mut ops)
            .unwrap();
        assert!(projection.flipped);
        // 10 bases at the bottom of the 500..600 exon, the 150-base intron, then
        // 30 in the 300..350 exon, read left to right in genome order.
        assert_eq!(cigar(&ops), "30M150N10M");
        // The alignment's last transcript base (129) sits at 350 - 30 = 320.
        assert_eq!(projection.position, 320);
    }

    /// Reverse-strand transcript coordinates must map to the mirror-image genomic
    /// positions, including at the exon boundaries where an off-by-one is easiest
    /// to make.
    #[test]
    fn reverse_strand_positions_mirror_the_transcript() {
        let model = reverse_model();
        assert_eq!(model.genome_position(0), 599);
        assert_eq!(model.genome_position(99), 500);
        assert_eq!(model.genome_position(100), 349);
        assert_eq!(model.genome_position(149), 300);
        assert_eq!(model.genome_position(150), 149);
        assert_eq!(model.genome_position(199), 100);
    }

    /// Whatever the projection does to the reference side, it must not change how
    /// many read bases the CIGAR accounts for.
    #[test]
    fn projection_preserves_the_read_length() {
        for model in [forward_model(), reverse_model()] {
            let projector = projector(model);
            let mut ops = Vec::new();
            for start in (0..190).step_by(7) {
                for len in [10u32, 45, 80] {
                    if start + len > 200 {
                        continue;
                    }
                    let transcript_cigar = vec![
                        CigarOp::new(CigarOpKind::SoftClip, 4),
                        CigarOp::new(CigarOpKind::Match, len),
                    ];
                    let Some(_) =
                        projector.project(0, start as i32, &transcript_cigar, None, &mut ops)
                    else {
                        continue;
                    };
                    let read_bases: u32 = ops
                        .iter()
                        .filter(|op| op.kind.consumes_read())
                        .map(|op| op.len)
                        .sum();
                    assert_eq!(read_bases, len + 4, "start={start} len={len}");
                }
            }
        }
    }

    /// A reference the annotation says nothing about (a decoy) has no projection,
    /// and must be reported as such rather than guessed at.
    #[test]
    fn an_undescribed_reference_does_not_project() {
        let projector = GenomeProjector {
            dropped: std::sync::atomic::AtomicU64::new(0),
            chromosomes: vec![("chr1".into(), 1000)],
            models: vec![None],
        };
        let mut ops = Vec::new();
        assert!(!projector.describes(0));
        assert_eq!(projector.project(0, 0, &matched(10), None, &mut ops), None);
    }

    /// bramble hands back GTF exons anchored at the 1-based start, so the
    /// conversion to genomic coordinates has to subtract the one. Getting this
    /// wrong shifts every projected record by a base while leaving intron
    /// lengths correct, which is exactly the kind of error that survives a
    /// casual look at the output.
    #[test]
    fn gtf_exons_become_zero_based_genomic_blocks() {
        let transcript = bramble_rs::annotation::Transcript {
            id: "t".into(),
            seqname: "chr1".into(),
            strand: '+',
            exons: vec![
                // GTF `101 300` and `801 1000`, as bramble reports them.
                bramble_rs::annotation::Exon {
                    start: 101,
                    end: 301,
                },
                bramble_rs::annotation::Exon {
                    start: 801,
                    end: 1001,
                },
            ],
        };
        let model = model_for(&transcript, 0);
        assert_eq!(model.exons[0].genome_start, 100);
        assert_eq!(model.exons[0].len, 200);
        assert_eq!(model.exons[1].genome_start, 800);
        assert_eq!(model.exons[1].len, 200);
        // Transcript coordinate 0 is the exon's first genomic base, and the base
        // after the first exon is the second exon's first base.
        assert_eq!(model.genome_position(0), 100);
        assert_eq!(model.genome_position(199), 299);
        assert_eq!(model.genome_position(200), 800);
    }

    /// The reverse-strand `MD` bug that a mismatch-free fixture cannot show: with
    /// no mismatches the string is just a run length, which is unchanged by
    /// reversing it. Only a mismatch reveals that the whole string has to be
    /// reverse-complemented alongside the CIGAR and the sequence.
    #[test]
    fn md_is_reverse_complemented_with_the_alignment() {
        let mut out = String::new();
        for (forward, reversed) in [
            ("71T28", "28A71"),
            ("50", "50"),
            ("0C1T0", "0A1G0"),
            ("10^ACGT5", "5^ACGT10"),
            ("3A2^GG4C1", "1G4^CC2T3"),
        ] {
            reverse_complement_md(forward, &mut out);
            assert_eq!(out, reversed, "reversing {forward}");
            // Reversing twice returns the original, which is what makes it a
            // faithful change of orientation rather than a rewrite.
            let once = out.clone();
            reverse_complement_md(&once, &mut out);
            assert_eq!(out, forward, "round-tripping {forward}");
        }
    }

    /// A deletion the projection has to cut in two is two deletions, and `MD`
    /// says so with an empty match between them. Carrying the transcript's
    /// string across produced the merged spelling, which disagrees with the
    /// CIGAR beside it.
    #[test]
    fn a_deletion_split_by_an_intron_is_two_md_deletions() {
        let projector = projector(forward_model());
        let mut ops = Vec::new();
        let mut scratch = Vec::new();
        let mut md = String::new();
        // Exon 1 is transcript 0..50: match 48, delete 4 across the boundary,
        // then match 10.
        let transcript_cigar = vec![
            CigarOp::new(CigarOpKind::Match, 48),
            CigarOp::new(CigarOpKind::Deletion, 4),
            CigarOp::new(CigarOpKind::Match, 10),
        ];
        projector
            .project(
                0,
                0,
                &transcript_cigar,
                Some(("48^ACGT10", &mut scratch, &mut md)),
                &mut ops,
            )
            .unwrap();
        assert_eq!(cigar(&ops), "48M2D150N2D10M");
        assert_eq!(md, "48^AC0^GT10");
    }

    /// Re-emitting must leave an unsplit alignment's `MD` exactly as it was.
    #[test]
    fn reprojecting_an_intact_alignment_changes_nothing() {
        let projector = projector(forward_model());
        let mut ops = Vec::new();
        let mut scratch = Vec::new();
        let mut md = String::new();
        let transcript_cigar = vec![
            CigarOp::new(CigarOpKind::Match, 20),
            CigarOp::new(CigarOpKind::Deletion, 2),
            CigarOp::new(CigarOpKind::Match, 15),
        ];
        projector
            .project(
                0,
                5,
                &transcript_cigar,
                Some(("12A7^GG15", &mut scratch, &mut md)),
                &mut ops,
            )
            .unwrap();
        assert_eq!(cigar(&ops), "20M2D15M");
        assert_eq!(md, "12A7^GG15");
    }

    /// A `.fai` is the cheap path to chromosome lengths, and only its first two
    /// columns matter.
    #[test]
    fn fai_parsing_takes_name_and_length() {
        let text = "chr1\t248956422\t112\t60\t61\nchr2\t242193529\t253404903\t60\t61\n";
        assert_eq!(
            parse_fai(text),
            vec![
                ("chr1".to_string(), 248_956_422),
                ("chr2".to_string(), 242_193_529)
            ]
        );
    }
}
