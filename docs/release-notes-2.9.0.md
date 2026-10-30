# salmon 2.9.0

A performance release. `--sketch` quantification is up to **2.4× faster**,
selective alignment up to **1.4× faster** on isoform-rich references, and
`salmon index` up to **3.7× faster**, with indexes that are now byte-for-byte
reproducible. Two things change for users:

- **Indexes must be rebuilt.** The index format moves to v3 (a new sshash
  dictionary format); salmon refuses older indexes with the command to rebuild
  them.
- **VBEM now converges each independent group of transcripts separately.**
  Abundances shift slightly, mostly between near-identical isoforms, and are
  more accurate on simulated data. Plain EM (`--useEM`) results are unchanged.

<!-- TODO(release): confirm version number (index format break; 2.7.0 set the
precedent of a minor bump for one). -->

## Speed

Full quantification (`salmon quant --sketch`, wall time, median of 3,
interleaved runs on one machine), against salmon 2.8.0 and kallistors 0.6.0
on the same reads and reference:

| workload | threads | 2.8.0 | 2.9.0 | kallistors 0.6.0 |
|---|---|---|---|---|
| muscle (SRR13638690, 4.4M × 2×170 bp), GENCODE v50 | 1 | 106.4 s | 43.3 s | 41.0 s |
| | 8 | 14.8 s | 6.1 s | 9.2 s |
| | 16 | 9.2 s | 3.9 s | 6.6 s |
| GEUVADIS (ERR188044, 10M × 2×75 bp), Ensembl cDNA | 1 | 84.3 s | 58.0 s | 63.2 s |
| | 8 | 11.0 s | 8.1 s | 11.5 s |
| | 16 | 6.7 s | 5.0 s | 7.5 s |

Peak memory is unchanged (0.3–1.3 GB on these runs; kallistors used 1.0–3.7 GB).

<!-- TODO(release): the 2.8.0 column is develop + the dependency upgrades
(piscem-rs 0.11.1 / sshash-lib 0.8), which leaves results byte-identical to
2.8.0; re-run against the released 2.8.0 binary for the final table. -->

Where the time went:

- **Inference.** The EM and VBEM now work on the connected components of the
  equivalence-class graph: transcripts are renumbered densely by component,
  single-transcript classes are folded into constants, and each component is
  iterated in cache-sized blocks. This schedule is kallistors' (see
  Acknowledgements). On the muscle run, VBEM falls from 46.6 s to 0.6 s
  single-threaded, and plain EM costs the same per iteration as kallistors'.
- **Mapping (`--sketch`).** piscem-rs 0.12 tracks only the transcripts a read
  can still map to — kallisto's running intersection of compatible
  transcripts, applied to piscem's position-aware hits — instead of a hash map
  entry for every occurrence of every k-mer hit; 12–20% faster mapping. Reading the intermediate mapping file back
  and building equivalence classes is also cheaper.
- **Selective alignment.** Grouping candidate anchors by transcript is now
  linear, and repeated hits that extend to an already-seen anchor are skipped:
  on the muscle run, full quantification at 8 threads goes from 31.8 s to
  22.9 s (Ensembl/GEUVADIS: 29.7 s → 25.7 s).
- **Index construction.** cf1-rs batches and prefetches its MPHF lookups, uses
  a faster MPHF for in-memory builds, and skips work repeated across the
  isoforms of a gene; piscem-rs and salmon parse and hash their inputs in one
  parallel pass. GENCODE v50 at 16 threads: 47 s → 12.8 s; Ensembl cDNA:
  13.0 s → 4.7 s.

## What changes in the results

### VBEM convergence (default optimizer)

VBEM used to iterate every transcript until the *slowest* group of
transcripts converged. It now stops each connected component — a set of
transcripts linked by shared reads — when that component meets the usual
criterion (same tolerance, same minimum iterations). Components are
independent problems, so this is still the same model; it simply stops doing
extra iterations on groups that have already converged.

On two polyester simulations (76 bp and 50 bp), accuracy against the true
counts improves: Spearman +0.010 / +0.008, mean absolute relative difference
−0.0022 / −0.0018 (95% bootstrap intervals exclude zero). Against 2.8.0's
output on real data, transcript-level Spearman is 0.93 (GENCODE v50, muscle)
and 0.985 (Ensembl, GEUVADIS); the largest differences are between
near-identical isoforms of very large genes (e.g. *NEB*, *TTN*), where the
reads cannot tell them apart. Gene-level estimates agree closely (Spearman
0.998 on muscle).

Results remain deterministic and identical at any `-p`. With `--useEM`, and
with `--emAccel squarem` or `daarem`, results are unchanged from 2.8.0 (to
within floating-point rounding, which does not show at the printed precision).

<!-- TODO(release): add the SQUAREM tolerance finding once measured. -->

### Mapping order within RAD records

With `--writeRad` (and in piscem's own RAD output), the mappings of a
single-end or orphaned read are now listed in transcript-id order (forward
before reverse for the same transcript) instead of hash order. The set of
mappings and their contents are unchanged, paired-end records were already in
this order, and quantification is unaffected.

### Everything else

Apart from the two changes above, results are byte-identical to 2.8.0:
selective-alignment mapping output, the mappings themselves in `--sketch`
mode (only their order within single-end/orphan records changes), equivalence
classes, and `quant.sf` under `--useEM` (verified on GENCODE v50 and Ensembl,
paired- and single-end).

## Indexes

- **Format v3.** sshash-lib 0.8 changes the dictionary format (6.0). Indexes
  built by salmon 2.7.0–2.8.0 (format v2) are refused with a rebuild command.
- **Reproducible builds.** Contigs are numbered by their first appearance in
  the reference, so building the same index twice gives identical files.
- **Memory.** Peak RSS of `salmon index` falls on GENCODE v50
  (6.8 GB → 4.3–4.8 GB) but rises on smaller references, e.g. Ensembl cDNA
  1.9 GB → 2.5 GB, from the faster in-memory MPHF construction.
  <!-- TODO(release): update after the cf1-rs sux-PHastR evaluation. -->

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
GENCODE v50) to reproduce that comparison, and studied its implementation:

- **The EM schedule is kallistors'.** Iterating each connected component in
  cache-sized tiles while reproducing the global stopping iteration, folding
  single-transcript classes into constants, and renumbering active
  transcripts densely — we found this was the whole of its EM speed advantage
  and adopted it in `salmon-infer` (`components` module).
- **The sketch-mode collector applies kallisto's principle** (Bray et al.,
  *Nature Biotechnology* 2016): a read is compatible only with the
  transcripts common to all of its k-mer hits, so only those need tracking.
  kallistors' fast implementation of this intersection is what showed us
  how much it saves.

Thanks also to Benjamin Demaille for the #1231 and #1239 fixes.

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
