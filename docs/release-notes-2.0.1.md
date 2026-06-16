# salmon 2.0.1

A patch release with two internal accuracy improvements to selective-alignment
quantification. **No breaking changes, no new options, and no output-format
changes** — indices built with 2.0.0 work unchanged (no rebuild needed), and
`quant.sf` / inferential-replicate outputs are produced exactly as before.

## What changed

Both changes refine how fragments are weighted across **multimapping**
transcripts (paralogs and alternative isoforms), bringing the estimates closer to
C++ salmon's:

- **Fragment-length probability in the equivalence-class weights.** The
  selective-alignment equivalence-class weights now include the fragment-length
  term (salmon's `logFragProb`), which had been omitted. Without it the
  conditional weights were flattened for transcripts whose implied insert size
  differs, slightly degrading how shared reads are split.
- **Abundance-aware fragment-length-distribution training.** The fragment-length
  distribution is now trained from fragments sampled by their abundance-aware
  posterior (matching salmon), so reads shared between near-duplicate transcripts
  contribute the dominant transcript's implied length. This concentrates the
  learned distribution as C++ salmon does.

## Impact

On benchmark data the per-transcript agreement with C++ salmon improves, and the
remaining difference sits at the run-to-run noise floor:

- Real data (GEUVADIS `ERR188044`, GRCh38 cDNA): per-transcript Spearman vs C++
  salmon 0.9705 → **0.9710**; NumReads Pearson **0.9996**; mapping rate
  unchanged.
- Simulated ground truth (human, polyester): the per-transcript Spearman gap to
  C++ salmon shrinks to ≤ 0.002 across difficulty levels and both the default
  (VBEM) and `--useEM` optimizers — at/near the noise floor.

Performance and memory are unchanged, and runs remain reproducible (run-to-run
variation is unaffected). Existing 2.0.0 results are still valid; re-quantifying
picks up the slightly improved estimates.

## Installation

Unchanged from 2.0.0:

```sh
curl --proto '=https' --tlsv1.2 -LsSf \
  https://github.com/COMBINE-lab/salmon/releases/latest/download/salmon-cli-installer.sh | sh
# or: cargo install salmon-cli   |   conda install -c bioconda -c conda-forge salmon
```

## Documentation

Full docs at **<https://combine-lab.github.io/salmon>**.
