# salmon 1.12.1

A reproducibility, correctness, and decoy-mapping patch on top of 1.12.0. It
resolves the lopsided splitting of duplicate (and other low-abundance,
non-identifiable) transcripts reported in issues #1008 and #1011, fixes the
`--writeMappings` + `--gcBias` crash (#1010), corrects several selective-alignment
and decoy-attribution behaviors (with two new flags, `--allowDecoyOrphans` and
`--orphansRequireUnmappedMate`), and enforces `maxReadOcc` in selective-alignment
mode (default raised to 250). All bug fixes track the Rust salmon (2.x) port,
which several of them are backported from.

## What changed

### Deterministic offline optimizer accumulation
The parallel VBEM/EM M-step (`CollapsedEMOptimizer`) previously accumulated each
transcript's abundance with atomic compare-and-swap adds across threads. Because
floating-point addition is not associative, the per-iteration result depended on
the (nondeterministic) order in which threads landed their updates. It now
accumulates into per-thread private buffers that are reduced into the output in a
fixed (shard-index) order — making the M-step deterministic for a given input and
eliminating the inter-thread CAS contention (a measurable per-iteration speedup
on large references; ~33% faster per iteration on a 193k-target human index at 16
threads).

### Uniform offline initialization by default
The offline optimizer is now seeded **uniformly** by default rather than from the
online (streaming) abundance estimates. The online estimates carry hog-wild,
multithreaded run-to-run noise; when fed as the optimizer seed, the VBEM update
amplifies that noise into arbitrary winner-take-all splits of low-abundance,
non-identifiable transcripts (exact duplicates being the extreme case). On
ground-truth simulations the uniform seed is accuracy-neutral versus the online
seed (no measurable change in Spearman/MARD, including the low-abundance stratum),
while removing the dominant source of run-to-run variability.

`--useOnlineSeed` restores the previous online-seeded initialization (e.g. for
faster convergence). On that path a **signature-equalization pass** is applied to
the seed: transcripts that are mutually non-identifiable (appear in exactly the
same equivalence classes with the same combined weights — exact duplicates and
the like) are assigned an identical seed, so the optimizer does not split them
arbitrarily. This pass is applied **only** under `--useOnlineSeed`; the default
uniform seed is already symmetric and needs no equalization.

### Selective-alignment mapping fixes (pufferfish)

Mapping-side correctness fixes backported from the Rust salmon port. The pufferfish
dependency pin is advanced to `1c78859` over the course of this release; the fixes
below are part of it.

- **Inclusive concordant fragment-length bound.** `joinReadsAndFilter` filtered
  concordant pairs with a strict `fragmentLen < maxFragmentLength`, which dropped
  a fragment whose length is *exactly* `maxFragmentLength`. Such a fragment is at
  — not beyond — the maximum, so it is a valid pair; the strict bound silently
  lost otherwise co-optimal placements (e.g. a read pairing to a paralog whose
  exon spacing puts the fragment right at the bound). The bound is now inclusive
  (`<=`), matching the Rust port. Effect is small and one-directional (recovers
  exactly-at-bound pairs); on a 1M-read human simulation it raised C++↔Rust
  target-set concordance from 99.963% to 99.966% with no regressions.

- **Orphan-rescue when one mate fails alignment.** When a concordant pair has
  exactly one mate passing the per-mate score threshold,
  `PuffAligner::calculateAlignments` now emits that mate as an orphan
  (`PAIRED_END_LEFT`/`RIGHT`) instead of discarding the whole fragment, so a
  strong mate is not lost because its partner is error-laden or mis-oriented —
  matching the Rust port's `m1`/`m2` orphan emission.

- **Decoy/transcript equal-footing seed chaining.** `fillMemCollection` counted
  MEM hits over non-decoy references only, so a mate that maps *only* to a decoy
  (e.g. the intronic mate of a genomic/pre-mRNA fragment) got no chains
  (`findOptChain` early-returns when `maxHits == 0`) and could never form the
  concordant decoy pair its exonic mate belongs to — leaking the genomic fragment
  as a spurious transcript orphan. MEM hits are now counted over **all** references
  (transcript and decoy). On the SRR1039508 3M subset against a GRCh38 decoy index
  this brings C++ decoy attribution in line with the Rust port (93.74% / 154,800
  decoy fragments vs Rust 93.75% / 154,284; previously 94.02% / 42,044).

### Decoy-orphan recovery: `--allowDecoyOrphans`

When a fragment maps concordantly to a decoy but one mate also maps to a
transcript, salmon discards it as decoy-dominated by default. The new
**`--allowDecoyOrphans`** flag (default **off**) keeps the transcript placement
instead, matching the Rust port's flag of the same name. Decoy domination is
enforced at several stages, each relaxed under the flag: `updateRefMappings` keeps
tracking a non-decoy hit below the decoy cutoff; `filterAndCollectAlignments` gates
non-decoy hits on the best non-decoy score (rather than `decoyThresh * bestDecoyScore`)
so a single-mate transcript orphan is not filtered out (the `estAlnProb`/`minAlnProb`
relative filter still prunes them); and the consumer keeps the fragment when a
transcript placement exists. On the full SRR1039508 GRCh38-decoy run this maps
94.36% with the flag vs 93.69% default, closely tracking Rust (94.81% / 93.69%).
**The default path is unchanged.**

### `maxReadOcc` is now enforced in selective-alignment mode (default 250); `--orphansRequireUnmappedMate`

`maxReadOcc` was effectively a **no-op** in selective-alignment mode: the
too-many-hits test was computed on the *pre-alignment* candidate set, but the
discard cleared an already-empty alignment group while the read stayed aligned
from the untouched hit list, so no cap was applied. It is now applied **after
alignment**, on the number of distinct aligned non-decoy targets (matching the Rust
port's `maps.len()` semantics). Capping the pre-alignment candidate count instead
would over-drop gene-family reads that *seed* hundreds of transcripts but *align*
to only a few dozen. The default is raised **200 → 250** (a mild relaxation now
that the cap actually bites). Verified live: `--maxReadOcc 1` now reduces to
uniquely-aligned reads.

Also new: **`--orphansRequireUnmappedMate`** (default **off**) emits a single-mate
(orphan) mapping only when the read's mate is *entirely unmapped*, rather than the
union of both mates' orphans when a pair has no concordant mapping.

### `--writeMappings` + `--gcBias` out-of-bounds crash (#1010)

Issue #1010 reported a segfault when `--writeMappings` and `--gcBias` were used
together (stack trace in `Transcript::gcDesc`). The root cause is a data-flow
hazard between the two features: the SAM writer clamped a right-overhanging
fragment's length by writing the clamped value back into the shared
`QuasiAlignment::fragLen`, and the GC-bias collector then reused that mutated
field to compute its `[start, stop]` window — so a wrapped/corrupted `stop` could
slip past the old `start >= 0 && stop < RefLength` guard and index `GCCount_[]`
out of bounds inside `gcDesc`.

This release closes the hazard at every link:

- **pufferfish `SAMWriter`** no longer mutates `qa.fragLen`; the TLEN clamp is
  computed into a local, so bias collection always sees the original value.
- **GC/positional-bias collection** (`SalmonQuantify.cpp`) computes its window
  from the clamped `fragLengthPedantic(RefLength)` recompute rather than the raw
  `fragLen` field, and tightens the guard to also require `start < RefLength` and
  `stop >= start`.
- **`Transcript::gcDesc`** gains an entry bounds-guard that returns an invalid GC
  descriptor (rather than reading out of bounds) for any out-of-range or
  ill-formed window — the defensive backstop at the faulting access itself.

These changes are defensive in depth and produce identical bias models and
`quant.sf`/SAM output on non-pathological inputs; together they make the
`gcDesc` out-of-bounds access unreachable.

The fix was verified directly: driving the pre-fix `gcDesc` with the exact
malformed windows the hazard produces flags an out-of-bounds `GCCount_` read
(bytes before/after the allocation) under valgrind, while the guarded version
returns an invalid descriptor with no out-of-bounds access. Because the bad index
is usually only a few elements past the array, the read typically lands on a mapped
heap page and returns garbage (silently corrupting the GC-bias model) — it only
**segfaults** when the index happens to fall on an unmapped page, which is why the
crash is data- and layout-dependent. The entry guard costs ~0.3 ns per call
(unmeasurable end-to-end), so it is kept as a permanent backstop at the faulting
site.

## Known limitation: residual run-to-run variability (the FLD-feedback path)

The seed change removes the *dominant, amplified* run-to-run variability, but a
small residual remains, and is **expected**. Its root cause is the
fragment-length distribution (FLD):

- The FLD is trained online over the burn-in window, so *which* fragments train
  it depends on read-arrival/thread-scheduling order.
- The FLD is consumed *during* the online pass to compute the fragment-length
  probability term of the (range-factorized) equivalence-class weights, so those
  weights bake in the FLD's training trajectory, which is also order-dependent.

Because the FLD sits inside the online feedback loop, this residual cannot be
removed by changing how fragment-length observations are weighted (we verified
that uniform, abundance-aware, and deterministic score-confidence weightings all
leave it essentially unchanged). Closing it fully would require a structural
change — deterministic read-delivery order plus freezing the FLD before it feeds
the equivalence-class weights (or a two-pass design) — which is out of scope for
this patch. In practice the residual is a fraction of a percent of expressed
transcripts and is smaller than the run-to-run variability of the pre-1.11
(non-SSHash) line; equivalence-class label sets and counts are unaffected.

## Notes

This is a patch on the final C++ salmon line. The Rust rewrite (salmon 2.0)
already seeds the optimizer uniformly with a deterministic single-writer M-step
by construction, so it is not subject to the duplicate-splitting behavior; the
same FLD-feedback residual applies to it and is likewise documented.
