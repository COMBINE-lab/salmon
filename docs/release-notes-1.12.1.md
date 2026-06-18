# salmon 1.12.1

A focused reproducibility / correctness patch on top of 1.12.0 addressing the
lopsided splitting of duplicate (and other low-abundance, non-identifiable)
transcripts reported in issues #1008 and #1011.

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
