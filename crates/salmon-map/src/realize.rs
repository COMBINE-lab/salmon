//! Base-level alignment realization for SAM/BAM output.
//!
//! # Why this is separate from [`crate::align`]
//!
//! Mapping only ever needs a *score*: is this candidate placement good enough to
//! keep? [`crate::align`] therefore runs ksw2 in `KSW_EZ_SCORE_ONLY` mode, which
//! skips the traceback matrix entirely, the single largest cost in a banded DP.
//! That is the right trade for the hot path, where the aligner runs once per
//! candidate across every read in the library.
//!
//! Writing a SAM or BAM record needs more than a score. A record has to say *how*
//! the read lines up: where each insertion and deletion falls (the CIGAR), how
//! many bases differ from the reference (`NM`), and which reference bases they
//! were (`MD`). None of that survives a score-only DP.
//!
//! So this module re-derives the alignment, with traceback, for the placements
//! that actually get written. That is a much smaller set than the candidates the
//! mapper scored: only the placements that survived filtering, and only when
//! `--writeMappings`/`--writeBam` is on. A run without mapping output never calls
//! into this module at all, so the mapping hot path is untouched.
//!
//! # What "realizing" a placement means
//!
//! The mapper hands over an approximate start position on a transcript. This
//! module aligns the read against a window starting there and reports:
//!
//! * the CIGAR, with real `I`/`D` operations where the read has indels and `S`
//!   where it overhangs a transcript end,
//! * `ref_start`, the leftmost reference base the alignment actually consumes,
//! * `NM`, the edit distance (mismatches + inserted + deleted bases),
//! * `MD`, the string that lets a reader reconstruct the reference from the read.
//!
//! # Gap placement
//!
//! ksw2 is run *without* `KSW_EZ_RIGHT` here, unlike the scoring path. That flag
//! right-aligns gaps; SAM convention (and every downstream variant caller) wants
//! them left-aligned, so the CIGAR this module produces is the conventional one
//! even though the score it implies is identical either way.

use crate::align::AlignConfig;

/// A CIGAR operation code, using the BAM numeric encoding so that both writers
/// and the genome projector can share one vocabulary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum CigarOpKind {
    /// `M`: consumes read and reference (match *or* mismatch, SAM's `M` does not
    /// distinguish them; `MD`/`NM` carry that information).
    Match = 0,
    /// `I`: bases present in the read, absent from the reference.
    Insertion = 1,
    /// `D`: bases present in the reference, absent from the read.
    Deletion = 2,
    /// `N`: reference skipped without the read crossing it. Only produced by the
    /// genome projector, for introns; a transcriptome alignment never has one.
    Skip = 3,
    /// `S`: bases present in the read but not aligned, and still stored in `SEQ`.
    SoftClip = 4,
}

impl CigarOpKind {
    /// The SAM letter for this operation.
    pub fn letter(self) -> char {
        match self {
            Self::Match => 'M',
            Self::Insertion => 'I',
            Self::Deletion => 'D',
            Self::Skip => 'N',
            Self::SoftClip => 'S',
        }
    }

    /// Whether this operation advances the position in the read.
    pub fn consumes_read(self) -> bool {
        matches!(self, Self::Match | Self::Insertion | Self::SoftClip)
    }

    /// Whether this operation advances the position on the reference.
    pub fn consumes_reference(self) -> bool {
        matches!(self, Self::Match | Self::Deletion | Self::Skip)
    }
}

/// One CIGAR operation: a kind and how many bases it covers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarOp {
    pub kind: CigarOpKind,
    pub len: u32,
}

impl CigarOp {
    pub fn new(kind: CigarOpKind, len: u32) -> Self {
        Self { kind, len }
    }
}

/// Append `op`, merging it into the previous operation when they are the same
/// kind. Zero-length operations are dropped: they are legal to *produce* while
/// stitching an alignment together but illegal to write into a record.
pub fn push_op(ops: &mut Vec<CigarOp>, kind: CigarOpKind, len: u32) {
    if len == 0 {
        return;
    }
    match ops.last_mut() {
        Some(last) if last.kind == kind => last.len += len,
        _ => ops.push(CigarOp::new(kind, len)),
    }
}

/// A realized alignment: everything a SAM/BAM record needs beyond the fields the
/// mapper already knows.
///
/// Buffers are owned so one instance can be reused across every record a worker
/// writes; [`RealizedAlignment::clear`] resets it without freeing.
#[derive(Debug, Default, Clone)]
pub struct RealizedAlignment {
    /// CIGAR operations, in reference order (5'→3' on the forward strand).
    pub ops: Vec<CigarOp>,
    /// Leftmost reference position the alignment consumes, 0-based.
    pub ref_start: i32,
    /// Edit distance: mismatches plus inserted plus deleted bases (SAM `NM`).
    pub edit_distance: u32,
    /// SAM `MD`: matched run lengths, mismatched reference bases, and `^`-prefixed
    /// deleted reference stretches.
    pub md: String,
    /// Alignment score of the realized alignment under the configured affine
    /// scoring. Reported separately from the mapper's `AS` so the two never get
    /// confused.
    pub score: i32,
}

impl RealizedAlignment {
    pub fn clear(&mut self) {
        self.ops.clear();
        self.ref_start = 0;
        self.edit_distance = 0;
        self.md.clear();
        self.score = 0;
    }
}

/// Scratch buffers for [`realize`], so realizing a record allocates nothing after
/// the first call on a given worker.
#[derive(Default)]
pub struct RealizeScratch {
    aligner: ksw2rs::Aligner,
    query_codes: Vec<u8>,
    target_codes: Vec<u8>,
    /// The read in reference-forward orientation.
    oriented: Vec<u8>,
    /// Candidate starts from the probe search, reused across records.
    candidates: Vec<i32>,
}

/// DNA5 code for a base: `A,C,G,T` → `0..3`, anything else → `4`.
fn dna5(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

/// 5x5 DNA5 scoring matrix, matching [`crate::align`]'s.
fn dna5_matrix(cfg: &AlignConfig) -> [i8; 25] {
    let mut mat = [-cfg.mismatch_pen; 25];
    for i in 0..5 {
        mat[i * 5 + i] = cfg.match_score;
    }
    mat[24] = 0;
    mat
}

/// Two bases count as a match for `MD`/`NM` only if they are the same known base.
/// An `N` on either side is a mismatch, which is what samtools' own `calmd` does.
fn bases_match(read: u8, reference: u8) -> bool {
    let (r, t) = (read.to_ascii_uppercase(), reference.to_ascii_uppercase());
    r == t && matches!(r, b'A' | b'C' | b'G' | b'T')
}

/// Realize the placement of `read` at approximately `approx_pos` on `ref_seq`.
///
/// `is_forward` says whether the read aligns to the reference's forward strand;
/// when it does not, the reverse complement is what gets aligned, which is also
/// what SAM stores in `SEQ`. `approx_pos` is the mapper's implied start and may
/// fall outside the reference: a negative position means the read overhangs the
/// 5' end, and a position plus read length beyond the reference means it
/// overhangs the 3' end. Both become soft clips.
///
/// Returns `false`, leaving `out` cleared, when no alignment can be formed at all
/// (an empty read or reference, or a window that lands entirely off the
/// reference). Callers fall back to the spoofed CIGAR in that case.
#[allow(clippy::too_many_arguments)]
pub fn realize(
    read: &[u8],
    is_forward: bool,
    ref_seq: &[u8],
    ambiguous: &[u32],
    approx_pos: i32,
    cfg: &AlignConfig,
    scratch: &mut RealizeScratch,
    out: &mut RealizedAlignment,
) -> bool {
    out.clear();
    if read.is_empty() || ref_seq.is_empty() {
        return false;
    }

    // SEQ is stored on the reference's forward strand, so that is the orientation
    // every coordinate below is expressed in.
    scratch.oriented.clear();
    if is_forward {
        scratch.oriented.extend_from_slice(read);
    } else {
        scratch
            .oriented
            .extend(read.iter().rev().map(|&b| complement(b)));
    }
    let query = std::mem::take(&mut scratch.oriented);
    let realized = realize_anchored(&query, ref_seq, ambiguous, approx_pos, cfg, scratch, out);
    scratch.oriented = query;
    realized
}

/// Mismatches when the read is laid down ungapped at `pos`, or `None` when it
/// would not fit inside the reference there.
fn ungapped_mismatches(query: &[u8], ref_seq: &[u8], pos: i32) -> Option<i32> {
    if pos < 0 || pos as usize + query.len() > ref_seq.len() {
        return None;
    }
    let reference = &ref_seq[pos as usize..pos as usize + query.len()];
    Some(
        query
            .iter()
            .zip(reference)
            .filter(|(&read, &base)| !bases_match(read, base))
            .count() as i32,
    )
}

/// Whether an ungapped alignment with `mismatches` is provably optimal: no
/// alignment containing a gap can score better. See [`realize_oriented`] for the
/// derivation.
fn ungapped_is_optimal(mismatches: i32, cfg: &AlignConfig) -> bool {
    mismatches * (cfg.match_score + cfg.mismatch_pen) as i32
        <= cfg.gap_open_pen as i32 + cfg.gap_extend_pen as i32
}

/// Length of the exact probe used to locate a misplaced read. Long enough to be
/// specific over the few hundred bases searched, short enough that a 1% error
/// rate usually leaves at least one of the two probes clean.
const PROBE: usize = 20;

/// Candidate starts for `query` near `approx_pos`, found by locating an exact
/// probe from the read inside the surrounding window.
///
/// Comparing the whole read at every offset answers the question but costs the
/// read length at each one, and the sweep has to reach a read length in either
/// direction. Searching for one exact substring costs the window instead of the
/// window times the read, and every occurrence of it names a start directly.
///
/// Two probes are taken, a third and two thirds of the way in, because a single
/// probe landing on a sequencing error finds nothing. Ends are avoided: they are
/// where reads are worst and where a genuine indel is most likely to sit.
fn probe_candidates(query: &[u8], ref_seq: &[u8], approx_pos: i32, reach: i32, out: &mut Vec<i32>) {
    out.clear();
    if query.len() < PROBE * 2 {
        return;
    }
    let low = (approx_pos - reach).max(0) as usize;
    let high = ((approx_pos + reach) as usize + query.len()).min(ref_seq.len());
    if high <= low + PROBE {
        return;
    }
    let window = &ref_seq[low..high];
    for offset in [query.len() / 3, query.len() * 2 / 3] {
        let Some(probe) = query.get(offset..offset + PROBE) else {
            continue;
        };
        for (index, candidate) in window.windows(PROBE).enumerate() {
            if candidate == probe {
                let start = low as i32 + index as i32 - offset as i32;
                if start >= 0 && !out.contains(&start) {
                    out.push(start);
                }
            }
        }
    }
}

/// Find where the read actually starts, when it plainly does not start where the
/// mapper said, and report how many mismatches it has there.
///
/// The mapper's position comes from a seed chain, and a mismatch near the read's
/// 5' end pushes that chain rightwards. A read can also be anchored so close to a
/// transcript end that it hangs off it, which turns into a soft-clipped record
/// that quietly throws bases away when the read often belongs further in and
/// matches in full.
///
/// The scan is gated on the read looking badly placed to begin with. A handful of
/// sequencing errors is not evidence of a wrong anchor, and relocating a read that
/// is already where it belongs would cost more than it saves.
fn repair_anchor(
    query: &[u8],
    ref_seq: &[u8],
    approx_pos: i32,
    cfg: &AlignConfig,
    candidates: &mut Vec<i32>,
) -> Option<(i32, i32)> {
    // A read that does not even fit inside the reference where it was placed is
    // as badly anchored as one that mismatches everywhere, and is the case most
    // worth looking at. Treating "does not fit" as "matches nothing" is what lets
    // the search run there at all.
    let here = ungapped_mismatches(query, ref_seq, approx_pos).unwrap_or(query.len() as i32);
    // An eighth of the read differing is far past what sequencing error explains,
    // and squarely in "this is not where the read goes" territory.
    if here <= (query.len() as i32) / 8 {
        return None;
    }
    // How far the anchor can be wrong is not a question about indel size, which is
    // what `indel_margin` bounds. A seed chain on a repeat or a paralogous stretch
    // can imply a start most of a read length away.
    let reach = query.len() as i32 + cfg.indel_margin as i32;
    probe_candidates(query, ref_seq, approx_pos, reach, candidates);
    let mut best = (here, approx_pos);
    let consider = |best: &mut (i32, i32), candidate: i32| {
        if let Some(mismatches) = ungapped_mismatches(query, ref_seq, candidate) {
            if mismatches < best.0 {
                *best = (mismatches, candidate);
            }
        }
    };
    for &candidate in candidates.iter() {
        consider(&mut best, candidate);
    }
    // Both probes can land on a sequencing error, and a read shorter than two
    // probes has none to offer. That leaves a read the cheap search cannot place,
    // and the alternative to placing it is a record with no base-level alignment
    // at all, so it is worth the exhaustive sweep. This runs on a fraction of a
    // percent of records, which is what makes it affordable.
    if best.1 == approx_pos {
        for delta in 1..=reach {
            consider(&mut best, approx_pos - delta);
            consider(&mut best, approx_pos + delta);
            if ungapped_is_optimal(best.0, cfg) {
                break;
            }
        }
    }
    (best.1 != approx_pos).then_some((best.1, best.0))
}

/// Realize at `approx_pos`, repairing the anchor first when the read plainly does
/// not sit there.
///
/// ksw2's extension pins the read's first base to the window's first base, so the
/// alignment can grow rightwards but never shift left. That is fine when the
/// mapper's position is the read's true start, and wrong when it is not: the DP
/// then has no way to express "the read starts earlier" except by opening a
/// leading insertion, which is both a worse alignment and a false one.
///
/// [`repair_anchor`] settles that before the DP runs, by comparison rather than
/// alignment. What remains afterwards is a safety net for the cases it declines
/// to act on: a leading insertion of `n` bases is itself the statement that the
/// read began `n` bases earlier, so shift by that and realize again.
#[allow(clippy::too_many_arguments)]
fn realize_anchored(
    query: &[u8],
    ref_seq: &[u8],
    ambiguous: &[u32],
    approx_pos: i32,
    cfg: &AlignConfig,
    scratch: &mut RealizeScratch,
    out: &mut RealizedAlignment,
) -> bool {
    let mut anchor = approx_pos;
    if !realize_oriented(query, ref_seq, ambiguous, approx_pos, cfg, scratch, out) {
        return false;
    }
    // A repaired anchor is a hypothesis, not a verdict. Reads with a genuine
    // deletion also look badly placed to an ungapped scan, since every base past
    // the gap is shifted, so the scan can talk itself into relocating a read that
    // was where it belonged. Realizing at both anchors and keeping the better
    // score makes the repair unable to do harm: the worst it can do is cost one
    // extra alignment on a read that already looked wrong.
    let mut candidates = std::mem::take(&mut scratch.candidates);
    let repaired = repair_anchor(query, ref_seq, approx_pos, cfg, &mut candidates);
    scratch.candidates = candidates;
    if let Some((repaired, mismatches)) = repaired {
        let mut alternative = RealizedAlignment::default();
        if realize_oriented(
            query,
            ref_seq,
            ambiguous,
            repaired,
            cfg,
            scratch,
            &mut alternative,
        ) && (alternative.score > out.score || ungapped_is_optimal(mismatches, cfg))
        {
            *out = alternative;
            anchor = repaired;
        }
    }
    // Refuse to publish an alignment the mapper itself would have rejected.
    //
    // When the anchor is wrong by more than the search could reach, the banded DP
    // still returns *an* alignment: a sawtooth of one-base insertions and
    // deletions that walks the read diagonally across the reference. It is a
    // valid CIGAR and a completely false claim, and it contradicts the `AS` on
    // the same record, which came from the mapper and says the read aligned well.
    // `min_accepted_score` is the mapper's own bar for believing a placement, so
    // an alignment below it is one this module has no business asserting. Falling
    // back to the spoofed CIGAR says "no base-level alignment here", which is at
    // least true.
    //
    // The bar is set by the bases the alignment actually claims, not by the whole
    // read: a read hanging off a transcript end is soft-clipped, and its clipped
    // bases are not evidence of anything either way.
    let claimed: u32 = out
        .ops
        .iter()
        .filter(|op| matches!(op.kind, CigarOpKind::Match | CigarOpKind::Insertion))
        .map(|op| op.len)
        .sum();
    if out.score < crate::align::min_accepted_score(claimed as usize, cfg) {
        out.clear();
        return false;
    }
    let leading_insertion = match out.ops.first() {
        Some(op) if op.kind == CigarOpKind::Insertion => op.len as i32,
        _ => 0,
    };
    if leading_insertion > 0 && anchor - leading_insertion >= 0 {
        let mut retried = RealizedAlignment::default();
        std::mem::swap(out, &mut retried);
        if realize_oriented(
            query,
            ref_seq,
            ambiguous,
            anchor - leading_insertion,
            cfg,
            scratch,
            out,
        ) && out.score >= retried.score
        {
            return true;
        }
        // The retry did not help, so keep the original rather than a worse one.
        std::mem::swap(out, &mut retried);
    }
    // An insertion at either edge is a legal CIGAR but a meaningless claim: those
    // bases are not aligned to anything, which is what a soft clip says.
    clip_edge_insertions(out);
    true
}

/// Turn insertions at the very start or end of the alignment into soft clips.
///
/// Both consume read bases and no reference, so the CIGAR stays consistent
/// either way, but only the clip is honest: an insertion says "these bases are
/// present in the read and absent from the reference at this point", which at
/// the edge of an alignment is a claim about nothing. It also inflates `NM`,
/// so the edit distance gives those bases back.
fn clip_edge_insertions(out: &mut RealizedAlignment) {
    let mut clipped = 0;
    let mut clip = |ops: &mut [CigarOp], index: usize| {
        if let Some(op) = ops.get_mut(index) {
            if op.kind == CigarOpKind::Insertion {
                op.kind = CigarOpKind::SoftClip;
                clipped += op.len;
            }
        }
    };
    // The edge is the first operation, or the second when a clip already leads.
    let leading = usize::from(out.ops.first().map(|op| op.kind) == Some(CigarOpKind::SoftClip));
    clip(&mut out.ops, leading);
    let last = out.ops.len().saturating_sub(1);
    let trailing = if out.ops.last().map(|op| op.kind) == Some(CigarOpKind::SoftClip) {
        last.saturating_sub(1)
    } else {
        last
    };
    if trailing > leading {
        clip(&mut out.ops, trailing);
    }
    out.edit_distance = out.edit_distance.saturating_sub(clipped);
    // Converting next to an existing clip leaves two clips in a row.
    let mut merged: Vec<CigarOp> = Vec::with_capacity(out.ops.len());
    for op in &out.ops {
        push_op(&mut merged, op.kind, op.len);
    }
    out.ops = merged;
}

/// Record the ungapped alignment: one `M` run over the whole read, with `NM` and
/// `MD` derived from a single pass over the bases.
fn finish_ungapped(
    query: &[u8],
    ref_seq: &[u8],
    window_start: i32,
    query_len: i32,
    score: i32,
    out: &mut RealizedAlignment,
) {
    out.ref_start = window_start;
    out.score = score;
    push_op(&mut out.ops, CigarOpKind::Match, query_len as u32);
    fill_md_and_edit_distance(query, ref_seq, &[], out);
}

/// The body of [`realize`], with the query already in reference orientation.
#[allow(clippy::too_many_arguments)]
fn realize_oriented(
    query: &[u8],
    ref_seq: &[u8],
    ambiguous: &[u32],
    approx_pos: i32,
    cfg: &AlignConfig,
    scratch: &mut RealizeScratch,
    out: &mut RealizedAlignment,
) -> bool {
    let ref_len = ref_seq.len() as i32;
    let query_len = query.len() as i32;

    // A read hanging off the 5' end contributes a leading soft clip; alignment
    // starts at reference 0 with whatever bases are left.
    let (leading_clip, window_start) = if approx_pos < 0 {
        ((-approx_pos).min(query_len), 0)
    } else {
        (0, approx_pos.min(ref_len))
    };
    let aligned_query = &query[leading_clip as usize..];
    if aligned_query.is_empty() || window_start >= ref_len {
        // Nothing of the read reaches the reference: it is all clip, which is not
        // a record worth spoofing a CIGAR for.
        return false;
    }

    // Most records do not need a dynamic program at all. A read placed where the
    // mapper put it usually lines up base for base, differing only by
    // substitutions, and for that alignment the CIGAR is a single `M` run whose
    // `NM` and `MD` fall out of one pass over the bases.
    //
    // The question is when that ungapped alignment is not merely plausible but
    // *optimal*, since anything less would make the output depend on a shortcut.
    // Affine scoring answers it: any alignment containing a gap pays at least
    // `gap_open + gap_extend`, and the very best it can do with the rest is match
    // everywhere, so its score cannot exceed `perfect - (gap_open +
    // gap_extend)`. The ungapped alignment scores `perfect - mismatches *
    // (match + mismatch)`. Whenever the second is at least the first, no gapped
    // alignment can win, and the DP would only spend time rediscovering that.
    //
    // At salmon's defaults (match 2, mismatch 4, open 6, extend 2) the rule
    // admits up to one mismatch, which covers most reads at realistic error
    // rates. `ungapped_score` is kept for the second shortcut below.
    // A reference whose bases the build had to replace cannot answer "does this
    // read match here" from the stored sequence alone, so the shortcuts that rest
    // on that answer are skipped and the full walk decides. This costs nothing in
    // practice: it is four references in the human transcriptome.
    let ungapped_mismatch_count = (leading_clip == 0 && ambiguous.is_empty())
        .then(|| ungapped_mismatches(query, ref_seq, window_start))
        .flatten();
    let ungapped = ungapped_mismatch_count.map(|mismatches| {
        query_len * cfg.match_score as i32
            - mismatches * (cfg.match_score + cfg.mismatch_pen) as i32
    });
    if let (Some(mismatches), Some(ungapped_score)) = (ungapped_mismatch_count, ungapped) {
        if ungapped_is_optimal(mismatches, cfg) {
            finish_ungapped(query, ref_seq, window_start, query_len, ungapped_score, out);
            return true;
        }
    }

    // Give the aligner room past the read length so a net deletion in the read
    // still has reference to consume.
    let window_end =
        (window_start + aligned_query.len() as i32 + cfg.indel_margin as i32).min(ref_len);
    let window = &ref_seq[window_start as usize..window_end as usize];

    let matrix = dna5_matrix(cfg);
    encode_into(aligned_query, &mut scratch.query_codes);
    encode_into(window, &mut scratch.target_codes);

    // Pass 1, score only: find how much of the window the best full-query
    // alignment consumes (`mqe_t`). ksw2's extension anchors at (0, 0), which is
    // exactly the "read starts here" constraint we want.
    let (query_end_target, query_end_reachable, target_end_query, query_end_score) = {
        let input = ksw2rs::Extz2Input {
            query: &scratch.query_codes,
            target: &scratch.target_codes,
            m: 5,
            mat: &matrix,
            q: cfg.gap_open_pen,
            e: cfg.gap_extend_pen,
            // salmon's own `dpBandwidth`. It used to have to cover
            // `indel_margin` as well, because the anchor correction read its
            // offset off a leading insertion and could not see one wider than
            // the band. `repair_anchor` now settles the anchor by comparison
            // before the DP runs, so the band only has to cover genuine indels,
            // and a banded DP is linear in its width.
            w: cfg.bandwidth,
            zdrop: -1,
            end_bonus: 0,
            flag: ksw2rs::KSW_EZ_SCORE_ONLY,
        };
        let ez = scratch.aligner.align(&input);
        (ez.mqe_t, ez.mqe > ksw2rs::KSW_NEG_INF, ez.mte_q, ez.mqe)
    };

    // The DP has now priced the optimal alignment. When that price is the one the
    // ungapped alignment already pays, the traceback would only redraw a single
    // `M` run, so skip it: filling and walking the traceback matrix is the
    // expensive half of a banded alignment, and this is the case where it buys
    // nothing. Reads with several mismatches but no indel land here.
    if let Some(ungapped_score) = ungapped {
        // The optimal full-query alignment ends exactly `query_len` bases into
        // the window (so it has no net indel) and scores what the ungapped one
        // scores: the ungapped alignment is therefore an optimal one, and the
        // conventional choice among any ties.
        if query_end_reachable
            && query_end_target == query_len - 1
            && query_end_score == ungapped_score
        {
            finish_ungapped(query, ref_seq, window_start, query_len, ungapped_score, out);
            return true;
        }
    }

    // How much window the traceback pass should align against, and how much of
    // the read's tail is soft-clipped.
    //
    // The read reaching its end inside the window is the normal case. It stops
    // being the right answer when the read is simply longer than the reference
    // left to align against: the DP can still "reach the query end" by dumping
    // every leftover base into an insertion at the last reference position, which
    // scores worse than clipping and describes an alignment nobody believes. A
    // read overhanging a transcript's 3' end is exactly that case, so when the
    // query cannot fit in the window at all, clip the tail instead.
    let overhangs_reference_end = aligned_query.len() > window.len();
    let (target_len, trailing_clip) =
        if query_end_reachable && query_end_target >= 0 && !overhangs_reference_end {
            ((query_end_target + 1).min(window.len() as i32), 0)
        } else if target_end_query >= 0 {
            (
                window.len() as i32,
                aligned_query.len() as i32 - (target_end_query + 1),
            )
        } else {
            return false;
        };
    if target_len <= 0 || trailing_clip < 0 {
        return false;
    }
    let aligned_query = &aligned_query[..aligned_query.len() - trailing_clip as usize];
    if aligned_query.is_empty() {
        return false;
    }

    // Pass 2, with traceback: a global alignment of the read (minus clips) to the
    // window prefix picked above. No `KSW_EZ_RIGHT`, so gaps come out
    // left-aligned as SAM expects.
    encode_into(aligned_query, &mut scratch.query_codes);
    let target = &window[..target_len as usize];
    encode_into(target, &mut scratch.target_codes);
    let (packed, score) = {
        let input = ksw2rs::Extz2Input {
            query: &scratch.query_codes,
            target: &scratch.target_codes,
            m: 5,
            mat: &matrix,
            q: cfg.gap_open_pen,
            e: cfg.gap_extend_pen,
            w: cfg
                .bandwidth
                .max((aligned_query.len() as i32 - target_len).abs() + 4),
            zdrop: -1,
            end_bonus: 0,
            flag: 0,
        };
        let ez = scratch.aligner.align(&input);
        (ez.cigar.clone(), ez.score)
    };
    if packed.is_empty() {
        return false;
    }

    out.ref_start = window_start;
    out.score = score;
    push_op(&mut out.ops, CigarOpKind::SoftClip, leading_clip as u32);
    for entry in &packed {
        let len = entry >> 4;
        let kind = match entry & 0xf {
            0 => CigarOpKind::Match,
            1 => CigarOpKind::Insertion,
            2 => CigarOpKind::Deletion,
            other => {
                debug_assert!(false, "ksw2 emitted an unexpected CIGAR op {other}");
                return false;
            }
        };
        push_op(&mut out.ops, kind, len);
    }
    push_op(&mut out.ops, CigarOpKind::SoftClip, trailing_clip as u32);

    // A leading or trailing deletion is a legal DP result but a meaningless
    // record: it claims reference bases outside the read's span. Trim them and
    // move the start position accordingly.
    trim_edge_deletions(out);
    fill_md_and_edit_distance(query, ref_seq, ambiguous, out);
    true
}

/// Drop `D`/`N` operations at either end of the CIGAR, shifting `ref_start` past
/// a leading one.
fn trim_edge_deletions(out: &mut RealizedAlignment) {
    while let Some(first) = out.ops.first().copied() {
        // A leading clip is fine; look past it for a deletion.
        let index = if first.kind == CigarOpKind::SoftClip {
            1
        } else {
            0
        };
        match out.ops.get(index).copied() {
            Some(op) if matches!(op.kind, CigarOpKind::Deletion | CigarOpKind::Skip) => {
                out.ref_start += op.len as i32;
                out.ops.remove(index);
            }
            _ => break,
        }
    }
    while let Some(&last) = out.ops.last() {
        let index = if last.kind == CigarOpKind::SoftClip {
            out.ops.len().saturating_sub(2)
        } else {
            out.ops.len() - 1
        };
        match out.ops.get(index).copied() {
            Some(op) if matches!(op.kind, CigarOpKind::Deletion | CigarOpKind::Skip) => {
                out.ops.remove(index);
            }
            _ => break,
        }
    }
}

/// Walk the finished CIGAR against read and reference to fill in `MD` and `NM`.
///
/// `MD` alternates: a count of reference bases that matched, then either the
/// reference base at a mismatch or `^` followed by the deleted reference bases. A
/// count is emitted between every pair of those, including zero-length ones, which
/// is what makes the string unambiguously parseable.
pub(crate) fn fill_md_and_edit_distance(
    query: &[u8],
    ref_seq: &[u8],
    ambiguous: &[u32],
    out: &mut RealizedAlignment,
) {
    use std::fmt::Write as _;

    let mut read_pos = 0usize;
    let mut ref_pos = out.ref_start as usize;
    let mut matched_run = 0u32;
    let mut edit_distance = 0u32;

    for op in &out.ops {
        match op.kind {
            CigarOpKind::Match => {
                for _ in 0..op.len {
                    let read_base = query.get(read_pos).copied().unwrap_or(b'N');
                    // A base the build replaced is reported as the `N` it was, not
                    // as the substitute the index carries: the record describes the
                    // reference the user supplied.
                    let ref_base = if ambiguous.binary_search(&(ref_pos as u32)).is_ok() {
                        b'N'
                    } else {
                        ref_seq.get(ref_pos).copied().unwrap_or(b'N')
                    };
                    if bases_match(read_base, ref_base) {
                        matched_run += 1;
                    } else {
                        let _ = write!(out.md, "{matched_run}");
                        out.md.push(ref_base.to_ascii_uppercase() as char);
                        matched_run = 0;
                        edit_distance += 1;
                    }
                    read_pos += 1;
                    ref_pos += 1;
                }
            }
            CigarOpKind::Insertion => {
                edit_distance += op.len;
                read_pos += op.len as usize;
            }
            CigarOpKind::Deletion => {
                let _ = write!(out.md, "{matched_run}");
                matched_run = 0;
                out.md.push('^');
                for _ in 0..op.len {
                    let ref_base = if ambiguous.binary_search(&(ref_pos as u32)).is_ok() {
                        b'N'
                    } else {
                        ref_seq.get(ref_pos).copied().unwrap_or(b'N')
                    };
                    out.md.push(ref_base.to_ascii_uppercase() as char);
                    ref_pos += 1;
                }
                edit_distance += op.len;
            }
            // An intron is not an edit: the read simply does not cross it, and MD
            // is defined over the aligned blocks only.
            CigarOpKind::Skip => ref_pos += op.len as usize,
            CigarOpKind::SoftClip => read_pos += op.len as usize,
        }
    }
    let _ = write!(out.md, "{matched_run}");
    out.edit_distance = edit_distance;
}

fn encode_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(seq.len());
    out.extend(seq.iter().map(|&b| dna5(b)));
}

/// How many reference bases a CIGAR spans.
pub fn reference_span(ops: &[CigarOp]) -> i32 {
    ops.iter()
        .filter(|op| op.kind.consumes_reference())
        .map(|op| op.len as i32)
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> AlignConfig {
        AlignConfig::default()
    }

    /// A deterministic sequence with no repeats long enough to give the anchor
    /// repair a second equally good home for a read. A periodic reference would
    /// make these tests assert on an arbitrary choice between ties.
    fn nonrepeating(len: usize) -> Vec<u8> {
        let mut state: u64 = 0x2545F4914F6CDD1D;
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[((state >> 33) % 4) as usize]
            })
            .collect()
    }

    fn realize_at(read: &[u8], reference: &[u8], pos: i32) -> RealizedAlignment {
        let mut scratch = RealizeScratch::default();
        let mut out = RealizedAlignment::default();
        assert!(realize(
            read,
            true,
            reference,
            &[],
            pos,
            &config(),
            &mut scratch,
            &mut out
        ));
        out
    }

    fn cigar_string(out: &RealizedAlignment) -> String {
        out.ops
            .iter()
            .map(|op| format!("{}{}", op.len, op.kind.letter()))
            .collect()
    }

    /// A read that matches the reference exactly is one `M` run, no edits, and an
    /// `MD` that is just its length.
    #[test]
    fn exact_match_has_no_edits() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let out = realize_at(&reference[4..24], reference, 4);
        assert_eq!(cigar_string(&out), "20M");
        assert_eq!(out.ref_start, 4);
        assert_eq!(out.edit_distance, 0);
        assert_eq!(out.md, "20");
    }

    /// A single substitution stays one `M` run, SAM's `M` covers mismatches, and
    /// shows up in `NM` and as a reference base in `MD`.
    #[test]
    fn substitution_is_reported_in_md_not_cigar() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut read = reference[4..24].to_vec();
        read[5] = b'A'; // reference has C here
        let out = realize_at(&read, reference, 4);
        assert_eq!(cigar_string(&out), "20M");
        assert_eq!(out.edit_distance, 1);
        assert_eq!(out.md, "5C14");
    }

    /// A deleted reference stretch becomes a `D` operation, is counted in `NM`,
    /// and is spelled out after a `^` in `MD`.
    #[test]
    fn deletion_produces_a_d_operation_and_caret_in_md() {
        let reference = b"ACGTTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA";
        let mut read = Vec::new();
        read.extend_from_slice(&reference[2..12]);
        read.extend_from_slice(&reference[15..27]);
        let out = realize_at(&read, reference, 2);
        assert_eq!(out.ref_start, 2);
        assert_eq!(cigar_string(&out), "10M3D12M");
        assert_eq!(out.edit_distance, 3);
        assert!(
            out.md.contains('^'),
            "MD should mark the deletion: {}",
            out.md
        );
    }

    /// Bases inserted in the read become an `I` operation and count toward `NM`,
    /// but leave no trace in `MD`, `MD` describes the reference.
    #[test]
    fn insertion_produces_an_i_operation() {
        let reference = b"ACGTTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA";
        let mut read = Vec::new();
        read.extend_from_slice(&reference[2..14]);
        read.extend_from_slice(b"GTA");
        read.extend_from_slice(&reference[14..26]);
        let out = realize_at(&read, reference, 2);
        assert_eq!(cigar_string(&out), "12M3I12M");
        assert_eq!(out.edit_distance, 3);
        assert_eq!(out.md, "24");
    }

    /// A read hanging off the 5' end of the transcript is clipped, not placed at a
    /// negative position.
    #[test]
    fn overhang_at_the_start_becomes_a_leading_clip() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut read = b"TTTT".to_vec();
        read.extend_from_slice(&reference[0..16]);
        let out = realize_at(&read, reference, -4);
        assert_eq!(out.ref_start, 0);
        assert_eq!(cigar_string(&out), "4S16M");
        assert_eq!(out.edit_distance, 0);
    }

    /// A read running past the 3' end is clipped there instead of claiming
    /// reference bases that do not exist.
    #[test]
    fn overhang_at_the_end_becomes_a_trailing_clip() {
        let reference = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut read = reference[24..].to_vec();
        read.extend_from_slice(b"TTTTTT");
        let out = realize_at(&read, reference, 24);
        assert_eq!(out.ref_start, 24);
        assert_eq!(cigar_string(&out), "8M6S");
    }

    /// A reverse-strand read is realized against the forward reference, matching
    /// what SAM stores in `SEQ`.
    #[test]
    fn reverse_strand_reads_are_realized_on_the_forward_strand() {
        let reference = b"ACGTACGTTTGGCCAAGGTTACGTACGTACGT";
        let forward = &reference[6..26];
        let read: Vec<u8> = forward.iter().rev().map(|&b| complement(b)).collect();
        let mut scratch = RealizeScratch::default();
        let mut out = RealizedAlignment::default();
        assert!(realize(
            &read,
            false,
            reference,
            &[],
            6,
            &config(),
            &mut scratch,
            &mut out
        ));
        assert_eq!(cigar_string(&out), "20M");
        assert_eq!(out.edit_distance, 0);
        assert_eq!(out.md, "20");
    }

    /// The mapper's position comes from a seed chain, and a mismatch near the
    /// read's 5' end pushes that chain rightwards. Anchoring the read there
    /// leaves the DP no way to say "it started earlier" except with a leading
    /// insertion, which is a worse alignment and a false one. The anchor has to
    /// be corrected from the alignment itself.
    #[test]
    fn a_late_anchor_is_pulled_back_instead_of_becoming_an_insertion() {
        let reference = b"ACGTTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAAACGTACGTTTGGCCAA";
        let mut read = reference[4..44].to_vec();
        read[3] = complement(read[3]);
        // Hand it a position 8 bases past the read's true start.
        let out = realize_at(&read, reference, 12);
        assert_eq!(cigar_string(&out), "40M");
        assert_eq!(out.ref_start, 4);
        assert_eq!(out.edit_distance, 1);
    }

    /// The anchor repair has to see offsets the DP's band cannot, which is the
    /// whole reason it exists: it slides the read by comparison, so a start that
    /// is wrong by more than the bandwidth is still found.
    #[test]
    fn a_badly_placed_read_is_slid_back_further_than_the_band() {
        let cfg = config();
        assert!(
            cfg.indel_margin as i32 > cfg.bandwidth,
            "this test is only meaningful when the offset can exceed the band"
        );
        let reference = nonrepeating(400);
        let read = &reference[60..160];
        // Anchored 25 bases late, well past the bandwidth of 15.
        assert_eq!(
            repair_anchor(read, &reference, 85, &cfg, &mut Vec::new()).map(|(p, _)| p),
            Some(60)
        );
        let out = realize_at(read, &reference, 85);
        assert_eq!(cigar_string(&out), "100M");
        assert_eq!(out.ref_start, 60);
        assert_eq!(out.edit_distance, 0);
    }

    /// A read that is where it belongs must not be slid. Sequencing errors are
    /// not evidence of a wrong anchor, and scanning for one would cost more than
    /// the alignment it replaces.
    #[test]
    fn a_correctly_placed_read_is_left_alone() {
        let cfg = config();
        let reference = nonrepeating(400);
        let mut read = reference[60..160].to_vec();
        for i in [10usize, 41, 77] {
            read[i] = complement(read[i]);
        }
        assert_eq!(
            repair_anchor(&read, &reference, 60, &cfg, &mut Vec::new()),
            None
        );
        let out = realize_at(&read, &reference, 60);
        assert_eq!(cigar_string(&out), "100M");
        assert_eq!(out.ref_start, 60);
        assert_eq!(out.edit_distance, 3);
    }

    /// A seed chain on a repeat can imply a start most of a read away from the
    /// truth. The old search reached only `indel_margin`, which is a bound on
    /// indel size and has nothing to do with how wrong an anchor can be, so the
    /// banded DP was left to describe a read that was nowhere near where it was
    /// placed. It obliged, with a sawtooth of one-base indels walking the read
    /// diagonally across the reference: a valid CIGAR and a false claim.
    #[test]
    fn an_anchor_wrong_by_most_of_a_read_is_still_found() {
        let reference = nonrepeating(600);
        let read = &reference[200..300];
        // 69 bases late, far past indel_margin.
        let out = realize_at(read, &reference, 269);
        assert_eq!(cigar_string(&out), "100M");
        assert_eq!(out.ref_start, 200);
        assert_eq!(out.edit_distance, 0);
    }

    /// When no anchor produces an alignment the mapper would have accepted, the
    /// honest answer is no base-level alignment at all. Publishing the DP's
    /// best-effort sawtooth would contradict the `AS` on the same record.
    #[test]
    fn an_unplaceable_read_is_refused_rather_than_described() {
        let reference = nonrepeating(600);
        let read = nonrepeating(100);
        let mut scratch = RealizeScratch::default();
        let mut out = RealizedAlignment::default();
        assert!(
            !realize(
                &read,
                true,
                &reference,
                &[],
                250,
                &config(),
                &mut scratch,
                &mut out
            ),
            "an unrelated read should not be given a CIGAR"
        );
        assert!(out.ops.is_empty());
    }

    /// An insertion at the very edge of an alignment claims nothing and inflates
    /// NM; a soft clip says the same thing honestly.
    #[test]
    fn edge_insertions_become_soft_clips() {
        let mut out = RealizedAlignment {
            ops: vec![
                CigarOp::new(CigarOpKind::Insertion, 3),
                CigarOp::new(CigarOpKind::Match, 20),
                CigarOp::new(CigarOpKind::Insertion, 2),
            ],
            edit_distance: 5,
            ..RealizedAlignment::default()
        };
        clip_edge_insertions(&mut out);
        assert_eq!(cigar_string(&out), "3S20M2S");
        assert_eq!(out.edit_distance, 0);
    }

    /// A base the index had to replace is reported as the `N` it was in the
    /// input, not as the substitute the index stores. Otherwise the record
    /// describes a reference the user never supplied, and disagrees with the
    /// `M5` on its own `@SQ` line.
    #[test]
    fn a_replaced_base_is_reported_as_n() {
        let reference = nonrepeating(200);
        let read = &reference[50..150];
        // Offset 80 of the reference held a non-ACGT base at index time.
        let mut scratch = RealizeScratch::default();
        let mut out = RealizedAlignment::default();
        assert!(realize(
            read,
            true,
            &reference,
            &[80],
            50,
            &config(),
            &mut scratch,
            &mut out
        ));
        assert_eq!(cigar_string(&out), "100M");
        // The read agrees with the stored base, but the reference base is unknown,
        // so it counts as an edit and MD names it.
        assert_eq!(out.edit_distance, 1);
        assert_eq!(out.md, "30N69");
    }

    /// `MD` and `NM` must be consistent with the CIGAR: every `M` base is either
    /// counted in a run or named as a mismatch, so the run lengths plus the named
    /// bases add back up to the aligned length.
    #[test]
    fn md_covers_every_aligned_reference_base() {
        let reference = b"ACGTTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA";
        let mut read = reference[3..30].to_vec();
        read[4] = complement(read[4]);
        read[19] = complement(read[19]);
        let out = realize_at(&read, reference, 3);
        let mut covered = 0usize;
        let mut digits = String::new();
        let mut chars = out.md.chars().peekable();
        while let Some(c) = chars.next() {
            if c.is_ascii_digit() {
                digits.push(c);
                continue;
            }
            covered += digits.parse::<usize>().unwrap_or(0);
            digits.clear();
            if c == '^' {
                while chars.peek().is_some_and(|n| n.is_ascii_alphabetic()) {
                    chars.next();
                    covered += 1;
                }
            } else {
                covered += 1;
            }
        }
        covered += digits.parse::<usize>().unwrap_or(0);
        assert_eq!(covered, reference_span(&out.ops) as usize);
        assert_eq!(out.edit_distance, 2);
    }

    /// Realizing must never produce a CIGAR whose read length disagrees with the
    /// read: that is the single check every SAM validator makes first.
    #[test]
    fn cigar_read_length_always_matches_the_read() {
        let reference = b"ACGTTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAAACGTACGT";
        for start in 0..20i32 {
            for len in [12usize, 25, 31] {
                let end = (start as usize + len).min(reference.len());
                if end <= start as usize {
                    continue;
                }
                let read = &reference[start as usize..end];
                let mut scratch = RealizeScratch::default();
                let mut out = RealizedAlignment::default();
                if !realize(
                    read,
                    true,
                    reference,
                    &[],
                    start,
                    &config(),
                    &mut scratch,
                    &mut out,
                ) {
                    continue;
                }
                let consumed: u32 = out
                    .ops
                    .iter()
                    .filter(|op| op.kind.consumes_read())
                    .map(|op| op.len)
                    .sum();
                assert_eq!(consumed as usize, read.len(), "start={start} len={len}");
            }
        }
    }
}
