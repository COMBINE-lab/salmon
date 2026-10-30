# salmon 2.9.0

A performance release. `--sketch` quantification is up to **2.5× faster**,
selective alignment up to **1.6× faster** on isoform-rich references, and
`salmon index` up to **4.9× faster** with less memory, and indexes are now
byte-for-byte reproducible. Two things change for users:

- **Indexes must be rebuilt.** The index format moves to v3 (a new sshash
  dictionary format); salmon refuses older indexes with the command to rebuild
  them.
- **VBEM now decides convergence for each independent group of transcripts
  separately**, and keeps iterating slowly mixing groups of near-identical
  isoforms until they have actually settled. Abundances shift, mostly within
  those families, and are more accurate on realistic simulations. Plain EM
  (`--useEM`) results are unchanged.

## Speed

Full quantification (`salmon quant --sketch`, wall time in seconds, median
of 3 interleaved runs, 2 at one thread; one machine, AMD EPYC 9575F), against
salmon 2.8.0 (built from its release tag) and kallistors 0.6.0 on the same
reads and reference:

| workload | threads | 2.8.0 | 2.9.0 | kallistors 0.6.0 |
|---|---|---|---|---|
| muscle (SRR13638690, 4.4M × 2×170 bp, trimmed as in kallistors' benchmark), GENCODE v50 | 1 | 105.9 | 43.1 | 42.1 |
| | 8 | 14.6 | 6.1 | 9.2 |
| | 16 | 8.9 | 3.7 | 6.7 |
| GEUVADIS (ERR188044, first 10M × 2×75 bp), Ensembl GRCh38 cDNA | 1 | 82.4 | 53.3 | 58.9 |
| | 8 | 10.9 | 7.7 | 11.3 |
| | 16 | 6.3 | 4.4 | 7.3 |

Peak memory is unchanged from 2.8.0 (0.3–1.4 GB on these runs; kallistors
used 1.0–3.7 GB). Selective alignment at 8 threads: muscle 34.2 → 21.2 s,
GEUVADIS 29.2 → 25.0 s. Building the index at 16 threads: GENCODE v50
64.0 → 13.2 s (kallistors 13.8 s), Ensembl cDNA 17.6 → 4.6 s (kallistors
6.0 s).

Where the time went:

- **Inference.** The EM and VBEM now work on the connected components of the
  equivalence-class graph: transcripts are renumbered densely by component,
  single-transcript classes are folded into constants, and each component is
  iterated in cache-sized blocks. This schedule is kallistors' (see
  Acknowledgements). VBEM's expectations `exp(ψ(·))` are computed by a
  vectorized, branch-free routine instead of a library digamma, so a VBEM
  iteration now costs within 15% of a plain-EM one. Inference time (from the
  intermediate mapping file, 1 / 8 / 16 threads): muscle VBEM 46.6 s → 0.88 /
  0.15 / 0.13 s, GEUVADIS VBEM 3.7 / 0.76 / 0.55 s; plain EM costs about
  what kallistors' does (muscle 5.8 / 0.80 / 0.47 s vs 5.4 / 1.0 / 0.64 s;
  GEUVADIS 3.9 / 0.87 / 0.59 s vs 3.9 / 0.82 / 0.78 s).
- **Mapping (`--sketch`).** piscem-rs 0.12 stops tracking a transcript as
  soon as it misses a hit (it could no longer be accepted), instead of keeping
  a hash-map entry for every occurrence of every k-mer hit; the accepted
  mappings are the same, and mapping is 12–20% faster. Reading the intermediate mapping file back
  and building equivalence classes is also cheaper.
- **Selective alignment.** Grouping candidate anchors by transcript is now
  linear, and repeated hits that extend to an already-seen anchor are skipped:
  on the muscle run, full quantification at 8 threads goes from 34.2 s to
  21.2 s (Ensembl/GEUVADIS: 29.2 s → 25.0 s).
- **Index construction.** cf1-rs batches and prefetches its MPHF lookups, uses
  a faster MPHF for in-memory builds, and skips work repeated across the
  isoforms of a gene; piscem-rs and salmon parse and hash their inputs in one
  parallel pass. GENCODE v50 at 16 threads: 64.0 s → 13.2 s; Ensembl cDNA:
  17.6 s → 4.6 s.

## What changes in the results

### VBEM convergence (default optimizer)

VBEM used to iterate every transcript until the *slowest* group of
transcripts converged. It now iterates each connected component — a set of
transcripts linked by shared reads — on its own, and stops it once the usual
criterion holds over the component (same tolerance, same minimum iterations)
**and** the movement still to come, projected from how fast its steps are
shrinking, is below 0.1% of the component's fragments. The second condition
matters for families of near-identical isoforms (e.g. *TTN*, *NEB*): their
steps are tiny long before they settle, so the criterion alone would stop them
far from converged. Components are independent problems, so this is still the
same model.

On six simulations with realistic expression profiles (three profiles,
GENCODE v50, selective alignment and `--sketch`), the total absolute error
against the true counts falls on every sample, by up to 31% (L1 error 0.145 →
0.101 of all fragments; 95% bootstrap intervals exclude zero), concentrated in
those isoform families; transcript-level Spearman and mean absolute relative
difference are unchanged (within ±0.001) and gene-level estimates are
unchanged. `meta_info.json` reports how the components ended
(`em_components`; see the output-format reference): a handful of very slowly
mixing families can reach the iteration cap after meeting the criterion,
which is reported there and in the log, and `em_converged` stays `true`. The
expert flag `--stopRule global` reproduces 2.8.0's VBEM results.

Results remain deterministic and identical at any `-p`. With `--useEM`, and
with `--emAccel squarem` or `daarem`, results are unchanged from 2.8.0 (to
within floating-point rounding, which does not show at the printed precision).

### A note on `--emAccel squarem`

SQUAREM is not a faster route to the default answer. Its stopping test is
applied across each three-step cycle, extrapolation included, so at the
default tolerance it ends far closer to the optimum than unaccelerated
iteration does. On simulated data, converging that far makes
transcript-level rankings slightly *worse* for any algorithm (more truly
low-expressed transcripts are driven to zero), and even when SQUAREM is
stopped at the same distance from the optimum as the default, it is slower.
The default (no acceleration) is the recommended setting.

### Mapping order within RAD records

With `--writeRad` (and in piscem's own RAD output), the mappings of a
single-end or orphaned read are now listed in transcript-id order (forward
before reverse for the same transcript) instead of hash order. The set of
mappings and their contents are unchanged, paired-end records were already in
this order, and quantification is unaffected.

### Single-end reads and orphans

How a single-end read or an orphan (one mate mapped) is weighted by the
fragments that could have produced it was wrong in the default two-pass path,
and depended on how `-l` was spelled:

- The intermediate RAD stored a reverse-strand read's *right* end where the
  quantifier expects its leftmost base (piscem's convention, which salmon now
  follows explicitly), so the "unseen mate must fit" bound and the read's
  5'-end sequence-bias context were a read length off. Reconstructed 5' ends
  on simulated reads: reverse single-end 0% → 98% exact, orphans 60% → 98%.
- Whether that bound applied was decided by `-l` (`-l A` applied it, `-l U`,
  `SF` and `SR` did not); it now follows each read's actual mapping type, so
  `-l A` and an explicit single-end type give the same result.
- **Departures from C++ salmon.** For single-end libraries the bound now uses
  the `--fldMean`/`--fldSD` default prior widened to an SD of at least 0.4×
  the mean, unless you pass `--fldMean`/`--fldSD`, which are then used as
  given; a narrow default that is too long for the library otherwise crushes
  every read near a transcript end. For orphans the bound is mixed 50/50 with
  "the mate lies off the transcript" (the reason it did not map), capping the
  penalty at ln 2.

Against simulated truth (fragment lengths 180±40, 250±25 and 320±60) at
default settings, single-end Spearman improves on every library type, e.g.
0.774 → 0.809 (stranded, `-l SR`) and 0.753 → 0.800 (unstranded, `-l U`),
and is never significantly worse than 2.8.0 with `-l A`; paired-end results
change only through orphans (MARD −0.0006). If you know your library's
fragment-length distribution, passing it improves single-end accuracy
further.

### Everything else

Apart from the changes above, results are byte-identical to 2.8.0:
selective-alignment mapping output, the mappings themselves in `--sketch`
mode (only their order within single-end/orphan records changes), equivalence
classes, and `quant.sf` under `--useEM` (verified on GENCODE v50 and Ensembl,
paired-end without orphans).

## Indexes

- **Format v3.** sshash-lib 0.8 changes the dictionary format (6.0). Indexes
  built by salmon 2.7.0–2.8.0 (format v2) are refused with a rebuild command.
- **Reproducible builds.** Contigs are numbered by their first appearance in
  the reference, so building the same index twice gives identical files.
- **Memory.** Peak RSS of `salmon index` is lower than 2.8.0's: GENCODE v50
  6.5 → 5.0 GB, Ensembl cDNA 3.1 → 2.6 GB (16 threads). cf1-rs now uses
  sux's PHastR minimal perfect hash for both in-memory and external-memory
  builds; genome-scale (decoy) cDBG construction is about 1.5–1.9× faster.

## Fixes

- **FASTA reads.** Uncompressed FASTA read files were refused in reads mode
  as truncated FASTQ (the truncation check assumed four-line records). They
  are accepted again, and give the same results as the equivalent FASTQ.
- **Read and fragment positions in selective alignment** (#1231, Benjamin
  Demaille). Positions are now derived from each read's projected ends rather
  than from its first and last matching k-mer, which cut a fragment short
  whenever a mismatch or indel fell near a read end. On simulated reads with
  known coordinates, reverse-strand mate positions go from 70% to 98% exact,
  fragment lengths from −2.0 bp mean error to +0.02 bp, and positions of
  reads overhanging a transcript end from ~65% to 98.5% exact; dovetailed
  fragments are now detected at the true rate. Quantification accuracy
  against truth improves slightly (MARD −0.001, with or without `--seqBias`).
  Visible effect on real data: the estimated fragment-length mean is
  ~2–3 bp higher (e.g. 170.1 → 172.7 on GEUVADIS). This is the estimator,
  not the positions: salmon learns the distribution from uniquely mapping
  pairs, which are longer on average than all fragments; the old positions'
  downward bias happened to cancel that. `--sketch` is unaffected.
- **Other fixes from #1231 (Benjamin Demaille):**
  - Paired-end sequence-bias contexts under `--online` used mate 1's
    orientation for the whole fragment, putting the reverse model at the
    fragment's leftmost base whenever mate 1 was the reverse read (half of
    `IU` pairs, all of `ISR`).
  - Reverse-strand single-end reads used an off-by-one (exclusive) end
    coordinate, and reverse orphans used their leftmost base as the 5' end.
  - Recovered orphan pairs (`--recoverOrphans`) now get correct positions.
  - The reference sequence store is checked against the index when loaded.
  - Ambiguity counts are accumulated in 64 bits.
  - `--writeMappings` no longer emits zero-length `0M` CIGAR operations.
  - The bias convolutions guard a zero lower fragment-length bound.
  - Ensembl GFF3 transcripts resolve to genes through `Parent`.
  - The Gibbs sampler's prior honours `--perNucleotidePrior`.
  - An EM in which no transcript is above the cutoff counts as converged,
    and the convergence check uses the new iterate, as C++ salmon does.
  - Numeric `quant` arguments are validated, and C++'s negative `--mp` form
    is accepted.
  - Gibbs sampling reuses its buffers between rounds, and the
    fragment-length distribution is updated once per fragment (speed only).
- **Soft-clipping at transcript ends** (#1239, Benjamin Demaille). With
  `--softclip` / `--softclipOverhangs`, a read overhanging the end of a
  transcript is now scored with the overhang clipped, as the option
  promises; default scoring (no soft-clipping) is unchanged.


## Acknowledgements

This release was prompted by, and borrows from,
[kallistors](https://github.com/deminden/kallistors) (Denis Demin), a Rust
reimplementation of kallisto whose 0.6.0 release reported quantification 2–3×
faster than salmon 2.8.0's `--sketch` mode on some datasets. We used its
published benchmark setup (SRR13638690 trimmed with its fastp settings,
GENCODE v50) to reproduce that comparison, and studied its implementation.

- **The EM schedule is kallistors'.** Iterating each connected component in
  cache-sized tiles while reproducing the global stopping iteration, folding
  single-transcript classes into constants, and renumbering active
  transcripts densely — we found this was the whole of its EM speed advantage
  and adopted it in `salmon-infer` (`components` module).

Thanks to Denis Demin for making kallistors and its benchmarks openly
available, and to Benjamin Demaille for the #1231 and #1239 fixes.

## Dependencies

piscem-rs 0.12, sshash-lib 0.8 (sux 0.15 / epserde 0.13), cf1-rs
<!-- TODO(release): version -->, bramble-rs 0.2, noodles 0.119
(noodles-bam 0.98, noodles-sam 0.93, noodles-bgzf 0.52), paraseq-temp
0.6.2-pre.1, plus semver-compatible updates. Each noodles crate now appears
once in the dependency graph.

## For library users

- `salmon-infer`: the default EM/VBEM path runs on a new component engine
  (`components` module); the accelerated paths keep the sharded kernels.
- `salmon-eqclass`: the class builder is a sharded hash table with lookups by
  borrowed label.
- piscem-rs 0.12: `MappingCache::hit_map` is a `HitMap` (see the piscem-rs
  changelog).
