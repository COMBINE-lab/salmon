# Bug fix: `--seqBias` observed model was never trained (degenerate correction)

## Summary

Sequence-specific bias correction (`--seqBias`) was effectively degenerate: the
**observed** ("foreground") read-bias model was never populated, so it stayed at
its uniform prior. The correction therefore reduced to `1 / P_expected(context)`
— it penalized transcripts with rare transcriptome contexts instead of using the
actual sequence bias of the reads, and it inflated effective lengths (often
beyond the transcript length).

## Impact

With the bug, `--seqBias`:
- collected **zero** observed-bias samples (independent of dataset size);
- produced effective lengths that frequently exceeded the transcript length
  (on the yeast `ERR458493` set: 4,015 of 6,571 transcripts had
  `effLen > length`; mean effLen 1,084 → 1,859);
- applied a correction uncorrelated with the true sequence bias.

After the fix, the observed model trains normally (≈430k samples on that set),
`effLen > length` drops to 117/6,571, and the bias-corrected effective lengths
match an independent re-implementation to within Pearson **r = 0.99943**
(bias-induced Δ correlation 0.06 → **0.97**).

## Root cause

In the `quant` command setup:

1. A local `int32_t numBiasSamples{0}` is bound to the `--numBiasSamples`
   option (default 2,000,000), but the default is only applied to that local by
   `po::notify(vm)`.
2. `QuantPipelineContext` is constructed **before** `po::store`/`po::notify`,
   capturing `numBiasSamples` **by value** — i.e. the init value `0`.
3. `processQuantOptions()` then does `sopt.numBiasSamples.store(0)`, which also
   overwrites the struct default of 1,000,000.

So `sopt.numBiasSamples` is `0` at run time, and the bias-sampling gate
(`needBiasSample && numBiasSamples > 0 && …`) never fires. The observed
`SBModel` receives no `addSequence` calls and remains uniform.

Because the gate depends on the *budget* being `> 0` (not on how many reads are
available), this affected runs of **any** size — a multi-million-read library
collected zero samples just like a small one.

## Fix

Re-sync the parsed value into the pipeline context after the options are
notified, before `processQuantOptions` consumes it:

```cpp
po::notify(vm);
// numBiasSamples is only applied to the bound local by po::notify(); the
// pipeline context captured it by value (0) before parsing. Re-sync it.
pipelineCtx.numBiasSamples = numBiasSamples;
```

(`--numBiasSamples` now also takes effect, which it previously did not.)

## Validation

- Observed model trains: ~430k samples on yeast `ERR458493`; the dumped
  `obs5_seq` is no longer uniform (pos-0 maxprob 0.250 → 0.295).
- `effLen > length` count: 4,015/6,571 → 117/6,571.
- Corrected effective lengths match the independent Rust port:
  EffLen Pearson 0.91 → **0.99943**, bias-Δ correlation **0.97**.
- Clean (non-bias) runs are unaffected (`numBiasSamples` is only used for bias
  sampling); `sample_data` mapping unchanged at 100%.
