# salmon 2.1.0 (draft — in progress)

A correctness-focused release that closes the remaining selective-alignment
mapping/quantification gaps against C++ salmon (1.12.1), adds N-aware decoy
indexing, and introduces an explicit salmon index *format* version. Released as a
**minor** (not patch) version because of the index-format bump and the breadth of
the decoy / short-transcript / duplicate-symmetry correctness changes.

**Index rebuild required.** salmon 2.1.0 writes (and requires) index format
**v1**, recorded as `index_version` in `info.json`. Indices built by salmon
2.0.0/2.0.1 carry no format version and are **rejected on load** with a rebuild
message — they predate the decoy / short-transcript contiguity guarantee below
and could mis-classify references. `quant.sf` and the inferential-replicate
formats are otherwise unchanged.

## Index format version (`index_version`)

salmon now records its own on-disk index **format version** in `info.json`,
independent of the software-release string (`salmon_version`) and of the
underlying piscem/cf1-rs versions. `SalmonIndex::load` refuses to open an index
older than the minimum it supports (currently v1) and prints an actionable
`salmon index …` rebuild command instead of risking a silent mis-load. The field
is bumped only when a layout/semantics change makes an older index unsafe to
read, so version going forward is explicit rather than inferred.

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
  Indices built before 2.1.0 **must be rebuilt** (they are now rejected on load by
  the index-format-version check above, rather than silently mis-interpreted).

- **Short transcripts in `quant.sf`.** The sub-`k` short transcripts are recorded
  by name and length and reported in `quant.sf` with 0 reads / 0 TPM (rather than
  silently dropped), so every input transcript appears in the output;
  `meta_info.json` reports the true `num_decoy_targets`.

This was found and fixed using the #1019 reproducer (GRCh38 Ensembl r115 gentrome
+ SRR1039508); the related expected-bias-model parallelization/decoy-exclusion of
PR #1020 is also folded in.

- **Sub-`k` decoy references are dropped at index time.** A decoy whose sequence
  is `<= k` after cleaning/poly-A clipping carries no k-mers, so it can never be
  seeded or catch a read. Previously the builder relocated such a reference into
  the trailing "short" block, scattering the decoy region; the contiguous-decoy
  assumption then broke and a transcript adjacent to the decoys could be
  mis-classified as a decoy — its reads dropped and its row omitted from
  `quant.sf` (silent data loss). These useless sub-`k` decoys are now dropped at
  build with a warning, which keeps the decoy block contiguous so the per-fragment
  decoy test stays a single O(1) range check. Sub-`k` *transcripts* are still kept
  and reported at 0 reads. A build-time guard verifies decoy contiguity and aborts
  with a clear error rather than risk a silent mis-classification.

- **Deterministic, input-order index build.** The cDBG builder (cf1-rs) now emits
  its reference tiling in input order via a bounded reorder buffer
  (`synchronize_output`), so the built reference numbering is deterministic and
  decoys — always last in the input — stay one contiguous block by construction
  (the prior task-completion order could scatter a small decoy among transcripts).
  This makes the O(1) decoy range check and the contiguity guard correct by design.

## N-aware decoy indexing (cf1-rs 0.5)

Decoy sequences now **retain their ambiguous (`N`) bases** instead of having them
replaced with pseudo-random ACGT. cf1-rs splits the de Bruijn graph on `N` runs
natively (recording the gaps in the tiling), which avoids seeding spurious k-mers
across assembly gaps and yields a less tangled, smaller, faster-to-build graph;
the reference store keeps the raw bytes and the aligner encodes `N` as a mismatch
(dna5 code 4). Transcripts are still `N`-replaced (matching salmon `FixFasta`). On
the full GRCh38 gentrome (≈5 % N) this removes ~151 M spurious k-mers (−5.7 %),
shrinks the index ~1.3 %, and trims build time/peak-RSS a couple percent; the
savings scale with N content. (`cf1-rs` ≥ 0.5 also adds `--poly-N-stretch` gating:
without it an `N`-containing input now fails loudly rather than corrupting.)

## `duplicate_clusters.tsv` under `--keepDuplicates`

salmon detects exact-sequence-duplicate transcripts and lists them in
`duplicate_clusters.tsv` even when `--keepDuplicates` retains them (downstream
tooling relies on the file). The Rust port had gated both the detection and the
file emission behind *not* keeping duplicates, so a `--keepDuplicates` index had no
`duplicate_clusters.tsv` at all. Duplicates are now detected in both modes (the
cluster list is populated regardless), only *collapsed* when not keeping
duplicates, and the file is always written. (#1015)

## Decoy handling in `--sketch` (pseudoalignment) mode

Decoys are now handled in sketch mode, where they previously **leaked**. Sketch
mappings are built directly from piscem's accepted hits and bypassed the
selective-alignment finalize where all decoy logic lived, so on a decoy-aware
index decoy references entered the equivalence classes as if they were
transcripts: decoy-only fragments were counted as mapped, decoys stole EM mass,
and `num_decoy_fragments` was never recorded. On SRR1039508 (3 M-read subset,
GRCh38 + 194 decoys) this inflated the sketch mapping rate to **97.5 %** with
**0** decoy fragments reported.

Sketch mode now applies the same decoy policy as selective alignment: decoy-only
fragments are dropped and counted as decoy, and decoy targets are removed from the
equivalence class. The corrected rate is **92.6 %** with **147,190** decoy
fragments — in line with selective alignment (93.7 %) — and removing the decoy
mass also recovers real-transcript abundances (nonzero transcripts 49.5 k → 51.1 k,
toward SA's 51.8 k).

- **`--decoyThreshold` is a no-op in sketch mode** (warned if set): pseudoalignment
  returns only equally-best mappings, so the `bestTxp < decoyThreshold * bestDecoy`
  comparison never triggers. A fragment is decoy-dominated only when it maps to
  decoys *and no transcript*.
- **`--allowDecoyOrphans` works in sketch mode.** Because sketch pairing is
  same-tid only, a fragment with one mate on a transcript and the other on a decoy
  forms no concordant pair and would otherwise be dropped as unmapped. With the
  flag, the transcript mate is recovered as an orphan (only when the other mate's
  hits are entirely decoys). On the subset above this recovers +8,270 fragments
  (92.6 % → 92.9 %), matching SA mode's `--allowDecoyOrphans` effect.

## Mapper allocation/perf

The alignment/chaining hot path reuses per-thread ksw2 scratch buffers (a reusable
aligner plus DNA5-encoded query/target buffers) instead of allocating per call, and
`chain_mems` gains exact single-MEM and two-MEM fast paths for the common small
cases. These are score/result-preserving; the two-MEM fast path carries the same
contained-anchor guard as the general DP. (#1015)

## Duplicate-transcript symmetry: consistent fragment-length reads

For sets of exact-duplicate transcripts the per-member read split is statistically
unidentifiable, so it must be handled symmetrically (as C++ salmon does). The Rust
port had been reading the fragment-length distribution's `pmf` directly from the
live, concurrently-updated histogram while building equivalence classes, so two
lookups of the *same* length (for two duplicate transcripts in one fragment) could
return slightly different values. The VBEM Dirichlet prior (`α<1`, sparsity-
inducing) then amplified that tiny asymmetry into a winner-take-all split,
concentrating a duplicate group's mass on one member and reporting fewer expressed
transcripts than C++. The FLD probability is now read from a per-fragment snapshot
(refreshed at mini-batch boundaries, mirroring salmon's cached CMF), so duplicate
transcripts receive identical weights and stay symmetric; quantification is also
now deterministic across thread counts on this path. On GEUVADIS ERR188044 vs C++
(byte-identical references) this raised NumReads `log`-Pearson 0.974→0.992 and
brought the nonzero-transcript counts into agreement (90.9k→95.3k vs C++ 95.3k).

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

## Decoy-orphan handling in selective alignment (`--allowDecoyOrphans`)

Two related fixes to how a fragment that pairs concordantly on the genome decoy
but only orphans onto a transcript (one mate exonic, the other intronic/genomic —
very common with a decoy-aware index) is handled:

- **A concordant decoy pair no longer suppresses a transcript orphan.** The
  orphan-fallback rule discarded *all* orphans whenever any concordant pair
  existed — including a concordant pair to a *decoy*. That destroyed the
  transcript orphan, leaving only the decoy pair, so the fragment was dropped as
  decoy-dominated with no surviving non-decoy mapping — and `--allowDecoyOrphans`
  could not even rescue it (it only acts when a transcript mapping survives). Now
  only a concordant *transcript* pair suppresses orphans; a decoy pair leaves the
  transcript orphan for the decoy-domination logic to adjudicate.
- **`--allowDecoyOrphans` now works as intended.** With the above fixed, the flag
  recovers the transcript orphan when the other mate maps to the genome decoy
  (default still drops it). On SRR1039508 (full) this raises the `--allowDecoyOrphans`
  rate 93.92 % → 94.81 %; the **default rate is unchanged** (byte-identical), and
  the recovered fragments match what C++ keeps as orphans. This is the
  selective-alignment mirror of the sketch-mode decoy-orphan rescue above.

## Equivalence-class weights: log-space normalization (no lost mapped mass)

Per-fragment equivalence-class weights are now normalized in **log** space
(`exp(auxProb − auxDenom)`, as C++ salmon does) rather than linearly (`w/Σw`). The
linear form, guarded by `Σw > 0`, silently produced all-zero weights for a
fragment whose implied lengths all have ~0 fragment-length-distribution
probability (every `w·exp(logFragProb)` underflows to 0); the VBEM then dropped
that equivalence class's count, **losing mapped mass**. The log-space form is
mathematically identical for the normal case (per-class scaling is EM-invariant)
but stays well-defined under total underflow. On SRR1039508 (full) the
mapped-mass loss (sum of `quant.sf` `NumReads` vs `num_mapped`) drops from **190
fragments to 0.1** — matching C++.

## `num_dovetail_fragments` is now reported

Rust always dropped dovetailed concordant pairs under the default no-dovetail
policy (matching C++) but reported `num_dovetail_fragments = 0`, because it only
inspected surviving pairs (never dovetailed after filtering). The counter is now
wired to the pairing stage and reports fragments whose only concordant pairing was
a dovetail. Diagnostic only — no change to mapping or quantification.

## Run-to-run reproducibility

Quantification is **byte-identical run-to-run when single-threaded** (`-p 1`); the
duplicate-transcript symmetry fix (per-fragment FLD snapshot + uniform init) means
exact-duplicate groups converge deterministically. Multi-threaded runs have a
small residual wobble (~0.26 % of assigned reads on SRR1039508, nonzero-transcript
set unchanged) from the **stochastic FLD training** under nondeterministic
fragment→thread scheduling — the same mechanism present in C++ salmon. Measured
head-to-head (`-p 16`, two runs each), Rust is in fact **more reproducible than
C++**: about half the average and total per-transcript variation, and far fewer
transcripts shifting by >2 % (Rust 2.2 % of expressed transcripts vs C++ 5.1 %).
Full parallel determinism (per-fragment-seeded FLD acceptance + order-independent
accumulation) is tracked as a follow-up; it is an improvement *over* C++, not a
correctness gap.

## Related C++ fixes (salmon 1.12.1)

These were also applied to the final C++ line so the two implementations agree:
an **inclusive concordant fragment-length bound** (a fragment of exactly
`maxFragmentLength` is valid) and **orphan-rescue** when one mate of a pair passes
the per-mate score threshold; the **`maxReadOcc` no-op fix** (clear `jointHits` so
the cap actually drops over-occurring fragments), the **default 200 → 250**, and
the new **`--orphansRequireUnmappedMate`** flag described above. The cap /
`jointHits` clearing and the orphan-union gating live in pufferfish
(`SalmonMappingUtils.hpp`, `Util.cpp`, `Util.hpp`); pin bumped accordingly.
