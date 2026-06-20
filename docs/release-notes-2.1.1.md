# salmon 2.1.1

A focused reporting-correctness patch on top of 2.1.0: on a **stranded** library,
salmon over-reported `num_mapped` / `percent_mapped`, breaking the invariant
`Σ NumReads (quant.sf) == num_mapped`. Reported in issue
[#1025](https://github.com/COMBINE-lab/salmon/issues/1025) by
[@BenjaminDEMAILLE](https://github.com/BenjaminDEMAILLE), with the reads-mode fix
in PR [#1026](https://github.com/COMBINE-lab/salmon/pull/1026).

## `num_mapped` is now counted after the strand-compatibility filter (#1025)

A fragment is counted as *mapped* only once it has at least one **strand-compatible**
placement that is actually quantified — matching C++ salmon, and restoring
`Σ NumReads == num_mapped`.

Previously, in reads (selective-alignment / sketch) mode, `num_mapped` was
incremented as soon as a fragment had *any* mapping, **before** the
strand-compatibility filter. On a stranded library a fragment whose every mapping
is strand-incompatible is then dropped (contributes no mass), but it had already
been counted as mapped — inflating `num_mapped`, `percent_mapped`,
`num_compatible_fragments`, and `num_assigned_fragments`. Unstranded libraries
were unaffected (the filter never drops a whole fragment).

On one real paired-end sample against a GRCh38 decoy-aware index, the auto-detected
stranded run reported `num_mapped = 311,967` but `Σ NumReads = 295,597` (5.25%
gap); after the fix the two agree, matching C++ salmon and the forced-unstranded
run.

**Alignment (`-a`) mode had the same bug, larger.** It set `num_mapped` to the
count of *every aligned fragment in the BAM*, before the per-fragment
strand-compatibility filter applied during assignment. On a stranded library
applied to FR alignments this over-reported by ~49% (e.g. `num_mapped = 280,047`
but `Σ NumReads ≈ 143,305`) and likewise broke the invariant; even unstranded was
off by ~0.09% from orphan/empty drops. Alignment mode now counts a fragment as
mapped only when it has a surviving strand-compatible placement, and its console
summary / `salmon_quant.log` distinguish **aligned fragments** (`num_processed`)
from **strand-compatible and quantified** (`num_mapped`) so a (necessarily) <100%
rate on a stranded BAM reads as strand-incompatibility, not lost alignments. The
`meta_info.json` field names are unchanged for downstream-tool compatibility.

## What is *not* affected

- **The point-estimate `quant.sf` for the compatible fragments is unchanged.**
  This was always a counting/reporting bug; the abundances come from the EM over
  the (unchanged) equivalence classes of the strand-compatible fragments.
- **Read counts already summed to `num_mapped` by construction.** salmon's
  Rust M-step distributes exactly `count` fragments per equivalence class (the
  VBEM prior only reweights and adds no mass), so `Σ NumReads == num_mapped`
  holds for the point estimate in every mode without a separate "scaled counts"
  renormalization. (C++ salmon needs its `useScaledCounts` step because its VBEM
  alphas carry the prior pseudo-count; we do not, and so do not replicate its
  mode-dependent enabling/disabling of that step.) Bootstrap/Gibbs replicates
  (`--numBootstraps` / `--numGibbsSamples`, off by default) are likewise made to
  sum to `num_mapped`; with this fix that target is now the correct mapped-fragment
  count on stranded data.

## Acknowledgments

Thanks to **[@BenjaminDEMAILLE](https://github.com/BenjaminDEMAILLE)** for the
clear, well-quantified report in [#1025](https://github.com/COMBINE-lab/salmon/issues/1025)
and the reads-mode fix + regression test in
[#1026](https://github.com/COMBINE-lab/salmon/pull/1026). The alignment-mode
extension and reporting changes were added on top of that PR.
