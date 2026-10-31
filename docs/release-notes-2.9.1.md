# salmon 2.9.1

A small performance patch, released so that piscem-infer 0.8.0 can build on
it. **Results are bit-for-bit identical to 2.9.0**, and no index rebuild is
needed.

## Speed

The EM and VBEM M-step that `salmon-infer` runs for each connected group of
transcripts is a little cheaper:

- **Plain EM** (`--useEM`, and piscem-infer's default): the step no longer
  bounds-checks every transcript it reads and writes, or tests every product
  for NaN.
  - Every transcript id in a component is checked once, when the components
    are laid out.
  - Because a NaN product always leaves a NaN in the step's output, that
    output is scanned once per step. In the (so far never seen) case that it
    holds one, the step is recomputed exactly as before.
- **VBEM** (the default): with a prior of at least 0.002 (the default is
  0.01), the step uses the same unchecked loop.

On the GENCODE v50 muscle sample (SRR13638690), time spent in the EM only:

| | 2.9.0 | 2.9.1 |
|---|---|---|
| plain EM, 1 thread (1420 iterations) | 6.62 s | 6.18 s |
| plain EM, 8 threads | 0.98 s | 0.93 s |
| VBEM, 1 thread | 1.013 s | 0.968 s |

In a full `salmon quant` run, the EM is a small share of the time.

## Checked

- `quant.sf` is byte-identical to 2.9.0 on that sample with `--sketch -p 8`,
  both with the default VBEM and with `--useEM`, and is identical across
  thread counts.
- `cargo test --release --workspace --locked`: 383 passed. Clippy
  (`-D warnings`) and fmt are clean.

## For library users

`salmon-infer`'s public API is unchanged.
