# salmon 2.1.1

A correctness-and-completeness patch on top of 2.1.0, focused on fragment-count
reporting and the inferential-replicate (`--numBootstraps` / `--numGibbsSamples`)
machinery:

1. **`num_mapped` reporting fix on stranded libraries** (issue
   [#1025](https://github.com/COMBINE-lab/salmon/issues/1025) /
   PR [#1026](https://github.com/COMBINE-lab/salmon/pull/1026), thanks
   [@BenjaminDEMAILLE](https://github.com/BenjaminDEMAILLE)).
2. **Negligible-abundance handling reworked** — the post-EM `min-alpha`
   truncation now *redistributes* the trimmed mass instead of rescaling, and the
   trimmed total is reported as `inference_truncated_mass`.
3. **Inferential-replicate output now aligns with `quant.sf`** — decoy columns are
   excluded and short (sub-`k`) transcripts are included, so `bootstraps.gz` /
   `names.tsv.gz` line up row-for-row with `quant.sf` (matching C++ salmon, fixing
   positional readers such as tximport on decoy-aware indices).
4. **Alignment mode (`-a`) gains bootstrap / Gibbs support** — previously only
   mapping-based quant produced inferential replicates.

---

## 1. `num_mapped` is now counted after the strand-compatibility filter (#1025)

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

## 2. Negligible abundances are redistributed, not rescaled; `inference_truncated_mass` reported

After the EM/VBEM converges, salmon zeroes transcripts whose estimated count falls
below `min-alpha` (1e-8) so vanishing point masses do not appear in `quant.sf`.
Previously the Rust port zeroed those transcripts and then **rescaled** the
survivors back up to the original total — which silently moved a tiny amount of
mass onto transcripts that did not earn it.

2.1.1 instead runs one final masked M-step: the truncated transcripts are pinned
inactive and each equivalence class redistributes its fragment mass only among its
*surviving* members. Mass therefore flows to the genuine co-mapped transcripts of
each dropped one, with no rescale-up. Under VBEM an explicit inactive mask prevents
the prior pseudo-count from reviving a truncated transcript. In the rare case where
**every** transcript in a class is truncated, that class's mass cannot be
reassigned; it is summed into a new `inference_truncated_mass` field in
`meta_info.json` (and a one-line `WARN` is logged) rather than being rescaled away.
On all tested data (easy/hard simulation and the real decoy-aware sample) this
value is `0.0`.

This makes the point estimate's behaviour explicit and mass-faithful, and applies
the same finalize step to every bootstrap replicate. (C++ salmon's analogous step
is its `useScaledCounts` rescale, which it needs because its VBEM alphas carry the
prior pseudo-count; the Rust M-step is mass-conserving by construction — the prior
only reweights and adds no mass — so we do not replicate C++'s mode-dependent
enabling/disabling of that step.)

For Gibbs, each posterior sample is now extracted as *truncate-then-normalize*: the
per-transcript rate is truncated at `min-alpha` first, then survivors are
normalized to the mapped-fragment total. (Previously the order was reversed, which
could lose a sliver of mass.) Gibbs replicates remain salmon's smoothed counts.

## 3. Inferential-replicate output aligns row-for-row with `quant.sf`

`aux_info/bootstrap/{names.tsv.gz, bootstraps.gz}` are now emitted over exactly the
same reference set, in the same order, as the `quant.sf` rows: **decoy references
are excluded** and sub-`k` **short transcripts are included** (always 0). Both
files therefore align positionally *and* by name with `quant.sf`.

Previously the per-replicate vectors were written over the full internal reference
range, which on a decoy-aware index appended the decoy columns (e.g. 194 extra
`f64` per replicate on a GRCh38 gentrome). A positional reader — tximport and
friends read `bootstraps.gz` against the `quant.sf` rows — would then misalign
every replicate by the decoy offset. C++ salmon excludes decoys here, and the Rust
port now matches: `len(bootstraps.gz) == n_replicates × (quant.sf rows) × 8`.

## 4. Alignment mode (`-a`) supports bootstrap and Gibbs

Inferential replicates are produced from the packed equivalence classes, which are
identical in structure whether they come from mapping or from an input BAM. That
path existed only in reads mode; `--numBootstraps` / `--numGibbsSamples` /
`--thinningFactor` now work in alignment mode too, writing the same
`aux_info/bootstrap/` layout (a transcriptome BAM has no decoy or short references,
so every `quant.sf` row is emitted) and recording the replicate count in
`meta_info.json`.

## Validation

Mass conservation holds exactly for the point estimate and every bootstrap
replicate (`Σ == num_mapped`, `inference_truncated_mass = 0`) across easy/hard
simulation and the real decoy-aware sample. Correlations are computed name-aligned.

**Hard simulation vs. ground truth** (193,759 transcripts; 100 bootstrap / 100
Gibbs; `Σ truth ≈ 4.50M`). On expressed transcripts (truth > 0, n = 101,200) the
posterior summaries improve on the point estimate, as expected:

| estimator | Spearman | log-Pearson |
|---|---|---|
| point estimate | 0.849 | 0.822 |
| bootstrap mean | 0.886 | 0.872 |
| Gibbs median | **0.916** | **0.930** |

(Over all 193,759 transcripts the three agree to ~0.02 Spearman: 0.867 / 0.864 /
0.850.) The Gibbs *median* does not sum to the mapped total — medians are not
additive — whereas the point estimate and bootstrap mean conserve mass exactly.

**Real data** (SRR1039508 subset vs. GRCh38 decoy-aware index, 532,617 quant.sf
rows incl. 65 short transcripts; 100 bootstrap / 100 Gibbs). Bootstrap and Gibbs
means vs. the point estimate, expressed transcripts:

| comparison (expressed, point > 1, n = 46,893) | Spearman | log-Pearson |
|---|---|---|
| point vs bootstrap mean | 0.978 | 0.984 |
| point vs Gibbs mean | 0.830 | 0.873 |

Every replicate — bootstrap and Gibbs — sums to exactly `num_mapped`
(2,812,472), as does each posterior mean.

These track C++ salmon: on the same sample the C++ point-vs-Gibbs correlation is
essentially identical (expressed Spearman 0.825 C++ vs 0.823 Rust), and both tools
fall to ~0.42 Spearman over the full ~90%-zero transcript set — an intrinsic
property of rank-correlating Gibbs-smoothed zeros, not a port artifact.

## What is *not* affected

- **The point-estimate `quant.sf` for the compatible fragments is unchanged** by
  the #1025 fix — it was always a counting/reporting bug; the abundances come from
  the EM over the (unchanged) equivalence classes of the strand-compatible
  fragments. The redistribute change (item 2) only alters how sub-`min-alpha`
  transcripts are handled, leaving all non-trivial abundances intact.
- **Read counts sum to `num_mapped` by construction.** salmon's Rust M-step
  distributes exactly `count` fragments per equivalence class (the VBEM prior only
  reweights and adds no mass), so `Σ NumReads == num_mapped` holds for the point
  estimate in every mode.

## Acknowledgments

Thanks to **[@BenjaminDEMAILLE](https://github.com/BenjaminDEMAILLE)** for the
clear, well-quantified report in [#1025](https://github.com/COMBINE-lab/salmon/issues/1025)
and the reads-mode fix + regression test in
[#1026](https://github.com/COMBINE-lab/salmon/pull/1026). The alignment-mode
extension and reporting changes were added on top of that PR.
