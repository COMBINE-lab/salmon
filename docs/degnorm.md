# Degradation normalization (`--degCoverage` + `salmon degnorm`)

RNA in a sequencing library is not intact, and it does not degrade evenly.
Longer transcripts degrade faster than short ones, and how much any given sample
lost depends on how it was handled. A degraded transcript yields fewer fragments
than its abundance warrants, so its count is an underestimate — and because the
shortfall differs per transcript *and* per sample, no global size factor removes
it. A differential-expression test run on those counts reads degradation
differences as expression differences.

[DegNorm](https://nustatbioinfo.github.io/DegNorm/) (Xiong *et al.*, *Genome
Biology* 2019) measures the effect from the shape of each transcript's coverage
curve and corrects the counts by it. salmon implements that model over its own
data, in two stages.

## The workflow

```sh
# 1. Quantify every sample as usual, plus a coverage dump.
for s in sample1 sample2 sample3; do
  salmon quant -i index -l A -1 ${s}_1.fq.gz -2 ${s}_2.fq.gz \
               --degCoverage -o quants/$s
done

# 2. Fit the model across the cohort.
salmon degnorm --quants quants/sample1 quants/sample2 quants/sample3 \
               -o degnorm
```

Stage 1 writes `aux_info/coverage.gz` in each quantification directory. Stage 2
writes, into its output directory:

| File | Contents |
|------|----------|
| `degradation_index.tsv` | `Name` × sample matrix of degradation indices; `NA` for transcripts that were filtered out rather than fitted |
| `counts_raw.tsv` | the `NumReads` column of each sample's `quant.sf`, as a matrix |
| `counts_adjusted.tsv` | `raw / (1 − DI)`, the degradation-corrected counts |
| `degnorm_summary.json` | parameters, per-sample mapped-fragment counts and depth factors, mean index per sample, how many transcripts were fitted |

`--noCounts` reports the indices only, without reading `quant.sf`.

## Why two stages

Degradation is only visible by comparison. One sample's coverage curve is
whatever it is; it takes *other* samples to establish what the curve should have
looked like, and therefore what is missing. So the expensive, streaming half —
turning every fragment placement into coverage — happens during `salmon quant`,
where the placement is already in hand and the reads will not be visited again;
and the cross-sample fit happens afterwards, over a handful of small dumps.

`salmon degnorm` refuses a single sample rather than reporting zeros for it.

## The model

For one transcript, stack the samples' coverage curves into a matrix `K`
(`m` samples × `n` position bins). Absent degradation every sample shows the
same shape scaled by its own abundance, i.e. `K` is rank one. Degradation only
removes coverage, so the fit is the closest rank-one matrix that stays *above*
the data:

```
minimise ‖K̃ − a bᵀ‖_F   subject to   K̃ ≥ K,  a ≥ 0,  b ≥ 0
```

and the degradation index of sample `i` is the fraction of the fitted envelope's
area that its actual curve is missing:

```
DIᵢ = 1 − (Σ_l K_il) / (Σ_l aᵢ b_l)
```

`DI = 0.5` means half of the fragments that transcript should have produced are
gone, so its count is corrected by `1 / (1 − 0.5)`.

Solved by block coordinate descent: fit the best rank-one approximation of `K̃`
by alternating least squares, then project back with `K̃ ← max(K, a bᵀ)`. Entries
that have fallen below `--maskFrac` of the current fit are held out of the least
squares as damaged (never out of the index) — without that, damage in a region
where only some samples carry evidence is absorbed into a sagging envelope, and
samples that lost nothing pick up an index they have not earned. The
`salmon-degnorm` crate documentation works through the example.

## How this differs from DegNorm

Both fit the same model. They do not see the same data, and the differences are
not cosmetic:

* **Transcripts, not genes.** DegNorm reads genome BAMs and works per gene, on
  exon-union coverage with exons shared between genes removed. salmon has no
  genome: it works per transcript, against the same targets as `quant.sf`. For
  gene-level results, aggregate the adjusted counts with your usual
  transcript-to-gene map.
* **Posterior coverage, not primary alignments.** A genome aligner picks one
  primary alignment per read. salmon does not: a fragment compatible with
  several transcripts contributes to each in proportion to its posterior
  probability, the same weighting the bias models use. Coverage curves for
  transcripts that share sequence are therefore smoother than a
  primary-alignment pileup would be.
* **Binned, not per base.** Coverage is accumulated into `--degCovBins`
  equal-width bins per transcript (default 100) rather than per position. Bins
  are a fixed *fraction* of each transcript, which is what makes them comparable
  across samples without interpolation.
* **Baseline handling.** DegNorm has a baseline-selection step for heavily
  degraded genes. salmon holds individual damaged *entries* out of the envelope
  fit instead. Same intent, different procedure, different numbers.
* **No smoothing.** DegNorm smooths coverage curves; salmon fits the binned
  curves directly, relying on the bin width for smoothing.

Do not expect indices from the two tools to agree numerically on the same data.

## Reading the numbers

Three properties of the model are worth knowing before acting on its output.

**The index has a positive floor.** Coverage curves are noisy, and an envelope
constrained to sit above the data sits above the noise too. On simulated intact
libraries the floor runs around 0.1 at 20× coverage and 0.25 at 5×. DI is
therefore a *comparative* measure — between samples, for a given transcript —
not an absolute percentage of degraded molecules. `--minCoverage` (default 0.1)
and `--minLength` (default 200) drop the transcripts where the floor dominates;
raising `--minCoverage` is the first thing to try if the indices look uniformly
inflated.

**Uniform loss is invisible.** A sample that lost the same fraction of coverage
at every position of a transcript is indistinguishable from a sample that simply
has less of it. Only uneven loss is identifiable — which is what degradation
actually looks like, but it does mean the index measures the *shape* of the
damage, not its total.

**Every sample must come from the same index.** Bin `b` of transcript `t` only
means the same thing in two samples if `t` has the same length in both.
`salmon degnorm` checks names and lengths and fails with the first target that
disagrees rather than fitting something meaningless.

## Cost

The accumulator is `(bins + 1) × transcripts × 8` bytes, shared across threads —
about 160 MB for a human transcriptome at the default 100 bins. It is allocated
only when `--degCoverage` is given, and released as soon as the dump is written.
Per fragment the hot path adds four relaxed atomic `fetch_add`s per compatible
mapping (a difference-array encoding of the fragment's span, so the cost does
not grow with how many bins the fragment covers).

The cohort fit is parallel across transcripts and takes seconds for a human
transcriptome.

## Limitations

* Reads mode only. Alignment mode (`-a`) and RAD input do not produce the
  per-fragment posteriors coverage is weighted by; `--degCoverage` warns and is
  ignored there.
* Not available under `--deterministic`, whose mapping pass deliberately skips
  the per-fragment posterior computation.
* Coverage weights come from the abundance-aware online posterior inside the
  model-training window and from score-only weights after it. This affects only
  how multi-mapping fragments are split between compatible transcripts.
