# salmon 2.10.0

A performance and accuracy release. Selective alignment is up to 1.45×
faster, VBEM up to 1.9× faster on large samples, and the fragment-length
distribution is estimated from every concordant pair instead of only the
uniquely mapped ones. No index rebuild is needed: the index format is
unchanged from 2.9.0.

**Results change slightly**, through the fragment-length distribution only.
The speedups below leave results bit-for-bit unchanged.

## Changes to results: fragment-length distribution

salmon estimated the fragment-length distribution (FLD) from uniquely mapped
proper pairs only. That sample is biased: which fragments map uniquely
depends on their length. On a 4M-pair simulation only 29% of fragments were
unique, and they ran long, because short fragments fall more often in exons
shared by several isoforms.

Now every concordant pair adds one length:
- **A uniquely mapped pair:** its only length, as before.
- **A pair with several placements:** the shortest length any of them
  implies. The longer alternatives usually come from isoforms that retain an
  intron or extend an exon between the mates, so the shortest is usually the
  true one. Ties go to the lowest transcript id, so the choice doesn't depend
  on placement order.

The same rule applies everywhere salmon builds the FLD: the mapping pass,
`--rad` re-derivation (`--fldPolicy derive`, piscem RADs), alignment (BAM)
input and genome projection. A distribution baked into a RAD file therefore
still matches one re-derived from it. `-l A` library-type detection still
counts only uniquely mapped pairs, since an ambiguous pair's placements can
disagree on orientation.

On the 4M-pair simulation, where every fragment's true length is known:

| estimator | fragments used | mean | histogram L1 error |
|---|---|---|---|
| truth | | 180.52 | 0 |
| unique only (2.9.1) | 29% | 182.95 | 0.048 |
| shortest implied length (2.10.0) | all concordant | 180.34 | 0.0056 |

On real data the bias runs either way, depending on which fragments map
ambiguously. The muscle sample's mean rises from 159.5 to 164.4, and
GEUVADIS's falls from 172.8 to 170.1. For an ambiguous fragment the shortest
implied length is never longer than the true one, so a rise means the
unique-only estimate was biased short by at least that much.

The effect on abundances is small, and where it is measurable it is an
improvement:
- Transcript-level MARD against truth improves by 0.0003–0.0004 in selective
  alignment on two of three GENCODE v50 simulations.
- Gene-level MARD improves on all six GENCODE runs (three simulations, two
  mapping modes).
- An Ensembl simulation is unchanged within noise.

(#1259)

## Speed

Full quantification against 2.9.1, wall time in seconds. Both builds are
from source with the same profile, on one machine (AMD EPYC 9575F) pinned to
one socket. Each figure is the median of 3 runs, or a single run at one
thread.

| workload | threads | 2.9.1 | 2.10.0 | speedup |
|---|---|---|---|---|
| muscle (SRR13638690, 4.4M × 2×170 bp), GENCODE v50, selective alignment | 1 | 151.0 | 108.1 | 1.40× |
| | 8 | 20.5 | 14.8 | 1.38× |
| | 16 | 11.4 | 8.5 | 1.34× |
| GEUVADIS (ERR188044, first 10M × 2×75 bp), Ensembl cDNA, selective alignment | 1 | 173.5 | 119.8 | 1.45× |
| | 8 | 24.2 | 17.1 | 1.42× |
| | 16 | 13.2 | 9.5 | 1.39× |
| muscle, `--sketch` | 1 / 8 / 16 | 37.7 / 5.5 / 3.4 | 38.2 / 5.5 / 3.4 | unchanged |
| GEUVADIS, `--sketch` | 1 / 8 / 16 | 48.1 / 7.0 / 4.0 | 47.9 / 6.9 / 4.0 | unchanged |

Sketch mode is unchanged, because inference is a small share of a full
sketch run. Peak memory is unchanged in all of these runs.

### Selective alignment

Scoring a candidate mapping is cheaper. Every change here is exact, and the
mappings are byte-for-byte those of 2.9.1.

- **Flanks without a DP.** A read flank whose ungapped alignment reaches
  `match·length − (gap open + gap extend)` can't be beaten by any gapped path,
  so it is scored without the dynamic program. Otherwise, the same bound
  limits how much reference and how wide a band the DP needs.
- **Early stop.** A candidate below the acceptance threshold is rejected
  whatever its exact score, so flanks are aligned last and stop as soon as the
  read can no longer be accepted.
- **No maximum tracking in ksw2.** The flank and gap DPs read only ksw2's
  end-of-query and end-of-target scores, so the overall maximum, which ksw2
  tracked on every anti-diagonal, is no longer computed (ksw2rs 0.2.1's new
  `KSW_EZ_NO_MAX`).
- **Orphans skipped when a pair exists.** Once a fragment has a valid
  concordant pair to a transcript, its orphan alignments would be discarded
  anyway, so they are no longer aligned. This doesn't apply under
  `--recoverOrphans`, where an orphan can become a pair.
- A read's reverse complement is kept between candidates instead of being
  recomputed.

(#1261, which supersedes #1260)

### VBEM

Under VBEM, most transcripts in a large connected component end up with an
abundance so far below the prior that their expectation is a constant, and
their terms are negligible in every class they belong to. They are now taken
out of the iteration. Their abundances are recovered from the class scales
when they're checked, and they come back if a term stops being negligible.
- On GEUVADIS's largest component (131k transcripts), 43% of the work is
  dead by iteration 114.
- VBEM time at one thread drops from 9.8 to 6.5 s on GEUVADIS and from 3.9 to
  2.1 s on a GENCODE simulation.

Results were byte-identical to 2.9.1's VBEM in every check, including
bootstraps, Gibbs sampling and `--perNucleotidePrior`. The one theoretical
exception: leaving out a term that is at most 10⁻²⁰ of a class's total could
change the last bit of a sum. Plain EM (`--useEM`) is unchanged. (#1261)

## Memory

- In `--rad` mode, the collapsed equivalence classes are freed once packed
  when no bias correction will read them again, lowering peak memory by
  5–10%.
- With bootstraps, the VBEM speedup costs some memory: each replicate running
  at once keeps per-transcript state for its largest component. On GEUVADIS
  with 20 bootstraps at 16 threads, peak memory is 1.44 GB against 1.26 GB
  for 2.9.1, and the run is 1.6× faster (13.0 vs 20.3 s).

## Other changes

- **Clearer errors for unreadable gzip input.** A truncated or corrupt
  `.fastq.gz` used to fail with a bare decoder message such as "compressed
  data ended early", without saying which file. The error now names the file
  and the likely cause, for example `trunc_r2.fq.gz: compressed data ended
  early (the file is probably truncated: an incomplete download or copy?)`.
  (piscem-rs 0.12.2, #1258)
- **`--online`** (the deprecated one-pass mode) now says it will be removed
  in a future release, instead of naming a release that has passed. (#1258)

## Checked

- `cargo test --release --workspace --locked`: 388 passed. Clippy and fmt
  are clean.
- Selective alignment: `quant.sf`, sorted `--writeMappings` output and the
  `meta_info.json` counters are identical with and without the alignment
  changes on 500k muscle pairs, with:
  - the default options;
  - `--softclip`, `--softclipOverhangs`, `--recoverOrphans`;
  - `--fullLengthAlignment`, `--minScoreFraction 0.8`, `--hardFilter`;
  - single-end reads.
- `quant.sf` is identical across thread counts.

## Acknowledgements

Thanks to Benjamin Demaille for the selective-alignment and VBEM work in
#1261.
