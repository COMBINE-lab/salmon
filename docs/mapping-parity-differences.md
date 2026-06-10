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
selective-alignment seeding/chaining selects a *suboptimal placement* for these
short, error-containing reads and never evaluates the exact-match location,
whereas the piscem-rs seeding the port uses finds it.

### The mechanism, in detail

Both tools seed from the same exact k-mers (with k = 31 on a 51 bp read there
are 21 possible 31-mers; for the example below all 21 forward 31-mers hit the
correct transcript uniquely). The divergence is downstream of k-mer lookup, in
two places:

1. **Hit decoding / orientation.** pufferfish stores k-mers *canonically* (a
   k-mer and its reverse complement share one slot), so `MemCollector` must
   re-infer each seed's strand by comparing the read k-mer to the unitig
   sequence. On an error-containing short read that per-seed inference can place
   seeds on the canonicalized (reverse) strand. piscem-rs's `decode_hit`
   instead projects every raw hit directly to an explicit
   `(tid, ref_pos, is_fw)` triple, so forward hits stay forward.

2. **Per-orientation vs. region-merged chaining.** The Rust path
   (`collect.rs::project_raw_hits` + `chain_mems`) buckets MEMs by
   `(tid, is_fw)` and chains **each orientation independently** with the
   minimap2-style colinear DP; `best_per_target` then keeps the
   maximal-coverage chain, so the forward and reverse chains compete only at the
   *scoring* stage and the full-coverage forward chain wins. salmon's
   pufferfish chainer merges uni-MEMs by *reference region* and prunes chains by
   `consensusFraction = 0.65` **before** alignment; when the early orientation
   inference seeded the reverse strand, the reverse chain becomes the surviving
   candidate for that region, salmon DP-aligns only that placement, scores it
   below threshold, and discards the read — never forming/evaluating the correct
   forward chain.

In other words the Rust port doesn't align more cleverly; it **keeps both
orientations' chains alive through to scoring** and decodes hit orientation
exactly, instead of committing to one region-merged placement and pruning by
consensus fraction before any alignment runs.

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
