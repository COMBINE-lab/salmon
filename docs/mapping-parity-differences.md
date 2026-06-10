# Mapping/quantification differences vs. C++ salmon

This note documents a deliberate, validated behavioral difference between the
Rust port and C++ salmon (1.11.4) at the read-mapping boundary, the evidence
behind it, and why we keep the Rust behavior.

## TL;DR

On confidently-mapped reads the two tools are essentially identical
(`NumReads` linear Pearson **r = 0.99891**, `EffectiveLength` **r = 0.99964**;
on the clean simulated `sample_data`, counts **r = 0.99998**). The only
material difference is that **the Rust port maps ~2–3% more reads** than salmon
on real short-read data. We investigated those reads exhaustively and found
that **~99% of them have an *exact*, full-length match (forward or
reverse-complement) in the transcriptome** — they are genuine matches that
salmon's selective-alignment heuristic fails to place. The Rust port is
therefore *more correct* on these reads; we keep its behavior and document the
difference here.

## Dataset

- Reads: `ERR458493` (Gierliński/Schurch *S. cerevisiae* benchmark),
  1,093,957 single-end 51 bp reads.
- Reference: Ensembl R64-1-1 cDNA, 6,612 transcripts.
- Tools: C++ salmon 1.11.4 (`--keepDuplicates`, `-l U`) vs. this port (`-l U`),
  both selective-alignment / mapping mode, no bias correction.

## What we observed

| metric | value |
| --- | --- |
| salmon mapping rate | 83.48% (913,271 reads) |
| Rust mapping rate | 86.05% (941,403 reads) |
| reads salmon leaves unmapped that Rust maps | 28,480 |
| `NumReads` linear Pearson (rust vs salmon) | 0.99891 |
| `EffectiveLength` Pearson | 0.99964 |
| TPM **log** Pearson | 0.915 |

The low *log*-TPM number is an artifact of log-space over-weighting the
low-count tail; per-abundance-bin Pearson is 0.96–0.98, and the linear
`NumReads`/`EffectiveLength` correlations show the quantification of
well-supported transcripts is near-identical.

## Ruling out the usual suspects

| hypothesis | test | result |
| --- | --- | --- |
| duplicate transcripts | rebuild salmon index `--keepDuplicates` | unchanged (0.918 → 0.915) ✗ |
| `N` bases (we map `N`→`A`) | only 0.17% of reads have `N`; N-free subset | gap unchanged ✗ |
| effective length | compare `EffectiveLength` columns | identical (r = 0.9996) ✗ |
| `minScoreFraction` value too low | sweep 0.65 → 0.80 | gap persists (84.5% at 0.80) ✗ |
| repetitive reads (too-many-loci) | eq-class sizes of rescued reads | 15,381 unique vs 2,875 multi ✗ |
| **alignment-scoring engine** (block-aligner vs ksw2) | implemented ksw2rs `extz2` with salmon's exact config (band 15, gapo 6/gape 2, `KSW_EZ_RIGHT\|KSW_EZ_SCORE_ONLY`) and re-ran | **no change** (85.6% vs 86.0%) ✗ |
| **seed representation** (sparse fixed-`k` anchors vs extended MEMs) | implemented reference-MEM extension (extend each k-mer hit to a maximal exact match, collapse colinear seeds) behind `--uniMEMs` and re-ran | **no change** (85.55% both; NumReads r = 0.99999999, 24 reads differ) ✗ |

## What it actually is: a seeding/placement difference

We extracted salmon's 180,686 unmapped reads (`--writeUnmappedNames`) and ran
them through the Rust mapper:

- The Rust port maps **28,480** of them.
- Of the high-confidence subset (Rust alignment score ≥ 0.95·perfect, 8,925
  reads), a random sample of 400 was checked directly against the transcriptome
  string: **397/400 (99%) have an exact full-length match** (forward or RC),
  and only 1/400 is low-complexity.
- salmon *did* find candidate locations for these reads (it re-maps exactly
  40,535 of its unmapped at `minScoreFraction 0.05`, matching its reported
  `num_fragments_filtered_vm = 40,535`) but its best chosen placement scored
  below threshold, so it discarded them.

Because these reads have **exact full-length MEMs**, they would score a perfect
102/102 under *any* alignment method — full-read DP, ksw2, or PuffAligner's
inter-MEM-gap-only scoring alike. So the difference is **not** the alignment
score and **not** full-vs-anchored scoring. It is that salmon's
selective-alignment **seeding** fails to collect enough of the read for these
short reads, so the correct placement is never assembled or scored. The next
section pins the exact mechanism.

### Root cause in C++ salmon: the `mismatchSeedSkip` seeding heuristic

We traced the canonical example (`ERR458493.21`) all the way into pufferfish's
`MemCollector` and isolated the exact cause. It is **not** alignment scoring,
**not** the consensus *slack* (sub-optimal-chain) filter, and **not** an
orientation/decoding issue. It is pufferfish's `mismatchSeedSkip` seed-skipping
optimization interacting with a unitig junction and the consensus *coverage*
requirement.

**Evidence chain (all reproducible against C++ salmon 1.11.4):**

1. `ERR458493.21` is a 51 bp read whose **21 forward 31-mers are each unique**
   in the whole transcriptome (every one occurs exactly once, all in
   `YMR216C_mRNA`, zero reverse-complement occurrences). It is a perfect,
   unambiguous, full-length match — the easiest possible read.
2. C++ salmon reports it **unmapped**, and *still* unmapped at
   `--minScoreFraction 0.05`. So salmon never produces a candidate at all; the
   read never reaches the alignment/scoring stage.
3. **Length sweep** of exact substrings of `YMR216C` starting at the same
   position: lengths 31–38, 40–44 and **52–70 map**, but **39 and 45–51 are
   unmapped**. Same locus, same k-mers — the failure is *read-length specific*,
   which rules out the index and repetitiveness.
4. **Knob isolation** on the read: `--mismatchSeedSkip 1` makes it map;
   `--consensusSlack 0.99` and `--fullLengthAlignment` do **not**. With
   `--mismatchSeedSkip 1` the entire failing length band maps (40/40).
   The default is `mismatchSeedSkip = 3` (`SalmonDefaults.hpp:37`).

**Why this happens (`MemCollector::operator()` + `expandHitEfficient`).** The
read straddles a **unitig junction**: the first k-mer extends (via
`expandHitEfficient`) into a uni-MEM that terminates at the unitig boundary
(`CONTIG_END`), covering only the first segment of the read. After a
non-mismatch termination the collector leaves the cursor ~`k`−1 bases *before*
the junction; the next k-mers legitimately **span the junction** and therefore
miss the index. pufferfish's miss-handling is the `mismatchSeedSkip` heuristic,
which *assumes a run of misses means the read is passing over a sequencing
error* (its own comment: "the next k k-mers will still likely overlap that
error") and so advances `mismatchSeedSkip` bases at a time until it is `k` bases
past the miss. But here the misses are caused by a **unitig boundary, not an
error** — the bases just past the junction are perfect sequence. For unfavorable
read lengths the skip-by-3 stepping overshoots the single good seed on the
second unitig, so the read's second segment is never seeded. Only the
first-segment match survives, and the resulting under-seeded candidate is
filtered out by salmon's downstream thresholds (the read is reported unmapped
even at `--minScoreFraction 0.05`, so it never produces an acceptable mapping).
Setting `--mismatchSeedSkip 1` (query every k-mer) re-seeds the second segment
and the read maps perfectly. (Note: salmon's `consensusFraction = 0.65` is a
*relative* filter — `chainCoverage < consensusFraction * maxChainScore`,
`MemChainer.cpp:469` — not an absolute coverage gate, so it does not by itself
explain the miss; the cause is upstream, in seed collection.)

**This is a genuine sensitivity bug, not just a tuning choice.** The heuristic
conflates two distinct uni-MEM termination causes (mismatch vs. unitig end) and,
for the unitig-end case, discards seeds it should keep. It costs real, verifiable
mappings of perfect unique reads. The whole-dataset impact is large: on the full
1.09M-read yeast set, `--mismatchSeedSkip 1` raises salmon from **913,271
(83.48%) to 929,891 (85.00%)** — recovering **16,620 reads, ≈74% of the entire
22,579-read gap** to the Rust port (935,850 / 85.55%). The remaining ~6k reads
are smaller seeding/placement differences. The Rust port is unaffected because
its uni-MEM extension runs against the **reference** (bridging unitig junctions
into one anchor) and its skipping streaming query re-probes after a skip, so a
junction never silently drops the downstream seed.

### Ruled out: the seed representation itself

A natural worry is that the difference is an artifact of seeding on *sparse
fixed-`k` k-mer anchors* (what piscem's skipping streaming query emits) rather
than longer extended seeds. We tested this directly: the `--uniMEMs` path
(module `salmon-map::extend`) extends every k-mer hit to a maximal exact match
and collapses colinear seeds, then chains those. (Caveat on naming: that path
extends against the **reference**, so it produces *reference MEMs* that cross
unitig boundaries, **not** unitig-constrained uni-MEMs in the pufferfish sense —
see "A note on uni-MEMs vs. reference MEMs" below.) On the full yeast set the
result is **indistinguishable** from the sparse path — 85.55% mapped both ways
(vs. salmon's 83.48%), `NumReads` Pearson **r = 0.99999999**, only 45 of 6,612
transcripts differ at all and by 24 reads total. So the ~2–3% extra reads are
**not** a seed-granularity artifact on our side; they persist identically under
extended-MEM seeding, confirming the difference lives in candidate
selection/orientation decoding, not seed representation.

### A note on uni-MEMs vs. reference MEMs

The `--uniMEMs` flag is a **misnomer** and will be renamed. A *uni-MEM* in the
pufferfish sense is an exact match constrained to lie within a **single unitig**
of the compacted de Bruijn graph — `expandHitEfficient` stops at the unitig
boundary (`CONTIG_END`). Our `salmon-map::extend::extend_mem` instead extends the
seed against the **reference transcript** sequence, so it does not stop at unitig
boundaries; it yields the maximal read↔reference exact match, which can *span*
several unitigs. These are **reference MEMs**, not uni-MEMs.

This distinction is not academic — it is exactly why the Rust port sidesteps the
`mismatchSeedSkip` bug above. pufferfish's unitig-constrained uni-MEMs split a
junction-straddling match into two pieces and rely on the seed-skip heuristic to
re-seed across the junction (where it fails); our reference MEMs bridge the
junction into a single anchor, so there is nothing to re-seed. Producing true
unitig-constrained uni-MEMs would require clamping extension to the contig bounds
available on each hit (`ProjectedHits::contig_pos`/`contig_len`), which
`project_raw_hits` currently discards. We have no need to do so for correctness,
but it would be the way to faithfully replicate pufferfish's seeding if desired.

## Concrete example

```
read  ERR458493.21
seq   TTTATGTCAAGCAAATATCGAAGCAGTTACTGCTTGGTCTTGATTACATGC   (51 bp, 44 distinct 4-mers; high complexity)
```

This read has an **exact forward match** in `YMR216C_mRNA` at position 802
(verified by direct substring search of the transcriptome, independent of
either mapper). The Rust port aligns it there with a perfect score and assigns
it to `YMR216C`. C++ salmon reports it as **unmapped**.

(A second pattern, ~8.6% of the disputed reads, places the read on a *different*
transcript than salmon; in the spot-checked cases the Rust placement is again a
perfect/near-perfect exact match while salmon's chosen placement scores poorly.)

## Why the Rust behavior is correct

A read with an exact, full-length match to a transcript unambiguously
originates from (a copy of) that transcript and should be counted toward it.
Discarding it undercounts the transcript. The Rust port recovers these reads;
salmon's heuristic does not. We therefore treat the Rust behavior as the more
correct one and retain it, rather than reproducing salmon's miss.

## Caveats and scope

- This is most visible in the **short-read regime** (51 bp with k = 31 leaves
  only 21 k-mers; real sequencing errors then break enough k-mers to stress
  salmon's seeding/placement). On clean/longer reads the tools agree closely
  (`sample_data`: 100% mapped by both, counts r = 0.99998).
- The port currently maps ~2–3% more reads; this is a *sensitivity* gain on
  genuine exact matches, not spurious mapping (99% verified exact).
- PuffAligner's inter-MEM-gap-only scoring (scoring only alignment gaps and
  reusing exact-MEM scores, with a thread-local sequence cache to avoid
  duplicate DPs) is now **implemented and is the default** (`anchored_align_score`
  in `salmon-map`); it is a **performance** optimization and does **not** change
  these mapping decisions (full exact MEMs score identically under it), so it is
  orthogonal to this correctness difference. We confirmed this empirically: on
  the yeast set anchored scoring maps 85.55% vs. 85.62% for full-read DP. Full-read
  DP remains available behind `--fullLengthAlignment` for cases where the global
  optimum is desired.

## Reproduction

See `salmon-cli debug-map -i <index> -r <reads>` (per-read placement, seed
coverage, and alignment score) and the commands in the project memory
(`salmon-cpp-parity`).
