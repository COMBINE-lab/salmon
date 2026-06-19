# salmon 2.0.2 (draft — in progress)

A correctness-focused patch that closes the remaining selective-alignment
mapping/quantification gaps against C++ salmon (1.12.1), plus the uni-MEM default
seeding change. **No output-format changes** — `quant.sf` and inferential
replicates are produced as before. Indices built with 2.0.0/2.0.1 still load,
but rebuilding is recommended to pick up the uni-MEM seeding.

These changes were found and validated by a head-to-head parity study against
C++ salmon on simulated human data (polyester, 193,759 transcripts); see
`docs/mapping-parity-differences.md`. With them, C++↔Rust per-read target-set
agreement is **~99.98%**. The small remaining residual has two parts: a poly-A
clipping mismatch in the older Rust index (resolved by the current binary), and a
marginal seeding-sensitivity tail where Rust retains a few near-optimal (1–2
mismatch-worse) secondary placements that C++'s seeding does not surface — the
same class of boundary difference documented for real data, not co-optimal ties.

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

## Decoy-aware indices: correct decoy boundary, no abundance-phase stall (#1019)

On a decoy-aware index (e.g. GRCh38 transcriptome + genome decoy) two related
defects are fixed:

- **Decoy boundary off-by-N.** The index builder (cf1-rs/piscem) relocates
  sub-`k` "short" transcripts — which carry no k-mers and never tile into the de
  Bruijn graph — to the *end* of the reference list, **after** the decoy block,
  so the references are ordered `[transcripts][decoys][short transcripts]`. The
  recorded `first_decoy_index` had been computed in input order, which on GRCh38
  pointed 87 references *into* the decoy block. As a result ~87 genome decoys
  (including a 250 Mb chromosome) were treated as transcripts. This (a) inflated
  the mapping rate, because reads best-explained by those decoys were counted as
  mapped instead of discarded, and (b) caused the **abundance phase to stall** —
  the expected-GC bias model swept a whole chromosome as if it were a transcript
  (one core pinned for minutes). The decoy span is now located by scanning the
  built index's reference names (`first_decoy_index` + a new `num_decoys`), so
  decoy classification, the bias-model bounds, and decoy-fragment accounting are
  correct. On SRR1039508 (3 M-read subset, `--seqBias --gcBias`) the abundance
  phase now completes in seconds instead of stalling, and the reported mapping
  rate drops from an inflated 98.8 % to **93.7 %**, matching C++ salmon's 94.1 %.
  Indices built before 2.0.2 **must be rebuilt** to pick up the corrected
  metadata (older indices fall back to the previous suffix interpretation).

- **Short transcripts in `quant.sf`.** The sub-`k` short transcripts are recorded
  by name and length and reported in `quant.sf` with 0 reads / 0 TPM (rather than
  silently dropped), so every input transcript appears in the output;
  `meta_info.json` reports the true `num_decoy_targets`.

This was found and fixed using the #1019 reproducer (GRCh38 Ensembl r115 gentrome
+ SRR1039508); the related expected-bias-model parallelization/decoy-exclusion of
PR #1020 is also folded in.

### `--allowDecoyOrphans` (new flag)

By default a fragment whose best transcript alignment is outscored by a decoy
(genome) alignment is discarded as decoy-dominated (matching salmon). On the
SRR1039508 subset, C++ salmon retains slightly more such fragments as transcript
mappings than the Rust port; the new `--allowDecoyOrphans` keeps a fragment's
transcript placement even when a decoy outscores it. On the 1 M-read subset
(no bias) this closes about half the residual gap to C++:

| run | mapping rate | decoy fragments |
| --- | --- | --- |
| C++ salmon (1.12.1) | 94.15 % | 14,128 |
| Rust (default) | 93.71 % | 51,525 |
| Rust `--allowDecoyOrphans` | 93.93 % | 49,342 |

We verified the residual is **not** a Rust seeding/chaining failure: of the 4,664
fragments (per 1 M) that C++ maps to a transcript and Rust does not, **all 4,664
produced chained candidates in Rust** — none are missed seeds — and exactly **one**
is a low-multiplicity, high-score alignment dropped for a non-decoy reason. The
clean low-multiplicity cases (1,391) are decoy-dominated, and **1,383 of them are
mapped by C++ as a single-mate orphan** (one mate exonic, the other mate unmapped
on the transcript) while Rust maps **both mates concordantly to the genome** and
discards the fragment — i.e. these are genomic / pre-mRNA fragments whose intronic
mate the decoy is designed to catch; the Rust port is arguably the more correct of
the two here. The remaining fragments Rust dropped but C++ kept were ~39 % whose
candidates were all filtered — **77 % of those map to more than 200 transcript
loci** (gene families / repeats, where one shared contig is annotated by hundreds
of transcripts). Investigating these uncovered — and we **fixed** — a real bug in
C++'s `maxReadOcc` cap (see below); they are no longer a divergence. Across the
whole set there were essentially **no clean, high-score, low-multiplicity
transcript alignments lost** (one, on the 1 M-read subset).
`--allowDecoyOrphans` is off by default; turn it on to bias the reported rate
toward C++ by retaining the transcript orphan for the genome-concordant fragments.

### `maxReadOcc` cap fix + default 250 + orphan policy (`--orphansRequireUnmappedMate`)

Three related changes, applied to **both** the Rust port and C++ salmon 1.12.1 so
they stay consistent:

- **`maxReadOcc` was a no-op in C++ selective-alignment mode — now fixed.** C++
  computed `tooManyHits = jointHits.size() > maxReadOccs` correctly but then
  `clearAlignments()` on an *already-empty* alignment group, while the read was
  aligned from the untouched `jointHits` with no cap guard — so over-occurring
  fragments were never actually dropped (verified: pre-fix `--maxReadOcc 1` and
  `200` gave identical results). The C++ fix applies the cap **after alignment**,
  on the number of distinct *aligned* non-decoy targets — matching the semantics of
  the Rust port's `maps.len()` cap. (Capping the pre-alignment `jointHits` count
  instead would over-drop gene-family reads that *seed* hundreds of transcripts but
  only align to a few dozen; on byte-identical sim references that mistake inflated
  the residual to 222/43, whereas the post-alignment cap leaves it at 147/50.)
  Verified the cap is now live: C++ `--maxReadOcc 1` drops to 32.9 % mapped
  (uniquely-*aligned* reads only). The cap counts the **total distinct aligned
  mapping set** (concordant + orphan union), which determines fragment ambiguity /
  equivalence-class size.
- **Default raised 200 → 250** in both, a mild relaxation now that the cap is
  actually enforced.
- **New `--orphansRequireUnmappedMate`** (both tools, default off). When a fragment
  has no concordant mapping, orphans are by default reported for **both** mates
  (their union — established salmon behavior). With this flag, a single-mate
  (orphan) mapping is emitted only when the read's mate is **entirely unmapped**, so
  a read that mapped only alongside a mate mapping to a *disjoint* reference set is
  not reported as an orphan. On the 1 M-read subset this removes ~1.2 % of fragments
  (those are genuinely both-mates-map-disjoint pairs).

After these fixes, on the SRR1039508 1 M-read subset at the shared default
(`maxReadOcc 250`): C++ 94.09 %, Rust 93.80 %, and Rust `--allowDecoyOrphans`
94.02 % — i.e. with the decoy-orphan flag the two agree to **0.07 %**, the
remaining difference being the benign seeding tail. With
`--orphansRequireUnmappedMate` the rate drops ~1.2 % in both. On **byte-identical**
simulated references (same clipped FASTA, matched thresholds) the per-read
target-set residual is 147 Rust-superset / 50 C++-superset, unchanged by the cap
fix now that both cap on the aligned-mapping count.

## Related C++ fixes (salmon 1.12.1)

These were also applied to the final C++ line so the two implementations agree:
an **inclusive concordant fragment-length bound** (a fragment of exactly
`maxFragmentLength` is valid) and **orphan-rescue** when one mate of a pair passes
the per-mate score threshold; the **`maxReadOcc` no-op fix** (clear `jointHits` so
the cap actually drops over-occurring fragments), the **default 200 → 250**, and
the new **`--orphansRequireUnmappedMate`** flag described above. The cap /
`jointHits` clearing and the orphan-union gating live in pufferfish
(`SalmonMappingUtils.hpp`, `Util.cpp`, `Util.hpp`); pin bumped accordingly.
