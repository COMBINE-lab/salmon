# salmon 2.0.2 (draft — in progress)

A correctness-focused patch that closes the remaining selective-alignment
mapping/quantification gaps against C++ salmon (1.12.1), plus the uni-MEM default
seeding change. **No output-format changes** — `quant.sf` and inferential
replicates are produced as before. Indices built with 2.0.0/2.0.1 still load,
but rebuilding is recommended to pick up the uni-MEM seeding.

These changes were found and validated by a head-to-head parity study against
C++ salmon on simulated human data (polyester, 193,759 transcripts) with matched
indices and poly-A clipping; see `docs/mapping-parity-differences.md`. With them,
C++↔Rust per-read target-set concordance is **99.98%**, and the remaining
differences are genuine co-optimal-paralog ties at the margin.

## Library-type auto-detection is now *applied* (`-l A`)

Previously, in automatic library-type mode the detector sampled the read prefix
and the inferred type was written to the output JSON — but it was **never used to
filter mappings**: the expected format stayed unset for the whole run, so the
strand-compatibility filter was inactive and orientation/strand-incompatible
mappings were counted during quantification. (Explicit `-l` types were filtered
correctly; only `-l A` was affected.)

Now, matching C++ salmon, the detector locks in the most likely library type once
it has seen enough of the prefix (50,000 confidently mapped fragments) and that
type is **enforced as a strand-compatibility filter for the rest of the sample**;
the inferred type continues to be reported in `lib_format_counts.json` /
`meta_info.json`. This is faithful for single-end and paired-end reads, across
orientation (inward/outward/matching) and strandedness (sense/antisense/
unstranded) — the inference and compatibility logic are ports of salmon's
`LibraryTypeDetector` and `compatibleHit`. For an unstranded library (the common
`-l A` result) the filter is a no-op; for a stranded library wrong-strand
mappings are now correctly excluded.

## Concordant-pairing correctness (matches salmon's `joinReadsAndFilter`)

A fragment that is a proper inward (FR) pair on its true transcript can map to a
paralog whose homologous segments are inverted or reordered, so the same pair
reads as **outward** or **same-strand** there. salmon rejects such pairs as
non-concordant; the Rust port had been keeping them as spurious extra
equivalence-class members. Three fixes bring concordancy in line with salmon,
independent of library type:

- **Pair over all per-target chains.** Pairing no longer collapses each mate to a
  single best-coverage chain before pairing, so a valid concordant placement at a
  repeat locus is not lost to a higher-coverage non-pairing one.
- **Emit all coverage-tied pairs.** When repeat loci tie on seed coverage, all are
  forwarded so alignment+finalize keep the locus that actually scores best
  (coverage does not always predict alignment score).
- **Require proper inward, opposite-strand orientation.** Outward pairs (forward
  mate downstream of the reverse mate) are now treated as dovetailed and dropped
  under the default no-dovetail policy, and same-strand pairs are rejected
  (salmon's `satisfiesOri`). Both are recoverable with `--allowDovetail` /
  explicit library types where appropriate.

## Anchored-alignment & chaining correctness

- **Divergent-paralog / hidden-indel scoring.** The anchored (PuffAligner-style)
  scorer no longer absorbs a read-overlap + reference-gap as a free overlap; the
  implied indel is DP-aligned and penalized, so divergent partial-paralog
  placements are correctly rejected (as `--fullLengthAlignment` and C++ do).
- **Contained-anchor chaining.** An anchor fully contained in its predecessor
  (a repeat sub-anchor on a different diagonal) is no longer chained, preventing a
  spurious overlap+gap that tanked otherwise-perfect alignments.

## uni-MEM seeding by default

Seeding now extends k-mer hits to unitig-local maximal exact matches (uni-MEMs)
by default — faster than, and at least as accurate as, the previous sparse k-mer
seeding. `--sparseSeeds` restores the old behavior. Output is equivalent; this is
a performance/representation change.

## Related C++ fixes (salmon 1.12.1)

Two of these were also backported to the final C++ line so the two implementations
agree (pinned pufferfish `17e1ccf`): an **inclusive concordant fragment-length
bound** (a fragment of exactly `maxFragmentLength` is valid) and **orphan-rescue**
when one mate of a pair passes the per-mate score threshold.
