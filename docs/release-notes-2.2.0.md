# salmon 2.2.0

A feature release that **decouples mapping from quantification** through the RAD
format, adds a **deterministic** quantification mode, and makes the **bias
models deterministic and faster**. Released as a **minor** (not patch) version
because of the new `--writeRad` / `--rad` / `--deterministic` modes and the new
`libradicl` dependency, even though existing workflows and outputs are
unchanged.

**No index rebuild required.** The on-disk index format is unchanged from 2.1.0
(`index_version` v1); `quant.sf` and the inferential-replicate formats are
unchanged. New behavior is opt-in via the flags below — a plain
`salmon quant -i … -1 … -2 …` run produces the same output as 2.1.2 (the one
exception is the bias models; see *Bias determinism* below).

Depends on the published **`libradicl` 0.14.0** (chunk compression +
backpatchable header tags).

## Decoupled mapping and quantification (`--writeRad`, `--rad`)

salmon can now separate the map and quant phases through the RAD format:

- **`--writeRad <PATH>`** writes per-fragment mappings to a RAD file (sketch or
  selective-alignment profile chosen automatically from the mapping mode).
  Quantification still runs; add **`--skipQuant`** to map only. The file is
  **piscem `map-bulk` compatible**.
- **`--rad <PATH>`** quantifies a RAD file directly, in parallel — whether it was
  written by salmon **or** produced by piscem `map-bulk`/sketch mode. No `-i` is
  needed; reference names travel in the RAD header.

This lets you map once and quantify many times, quantify a shared piscem RAD with
salmon's EM, or split the two phases across machines.

## Deterministic quantification (`--deterministic`)

RAD-mode quantification is **byte-identical across thread counts and runs**, and
a new **`--deterministic`** reads mode brings that guarantee to FASTQ input: it
maps the reads once to an intermediate RAD, then quantifies from it with a fixed
fragment-length distribution (no second mapping pass; the intermediate is removed
on success unless `--keepRad` or `--writeRad <PATH>` is given).

The determinism comes from making the **computation itself** order-independent
rather than sorting records:

- equivalence-class weights accumulate in **fixed-point `u128`** integers
  (integer addition is associative, so thread-count- and order-independent),
  materialized back to `f64` once at the end;
- the fragment-length distribution is built from **integer count histograms** in
  a fixed length order, then frozen;
- salmon-written RAD **bakes** the order-independent FLD, initial abundances, and
  resolved library format into the header (via `libradicl` 0.14 backpatchable
  tags), so it requantifies in **one pass**; a piscem RAD (nothing baked) derives
  a unique-fragment FLD in a first pass, then quantifies — **two passes**. A
  `baked_flags` marker drives the dispatch, and the baked library format is used
  for `-l A` concordance filtering without re-inference.

`-l A` library-type detection in `--deterministic` mode is itself made
order-independent so the baked format is stable.

## Bias-aware RAD, and deterministic + faster bias models

Sequence, GC, and positional bias correction now work on RAD input
(`--rad … --seqBias/--gcBias/--posBias`), sharing the same bias-correction tail
as the alignment (BAM) path. `--deterministic` bias no longer requires `-t`: the
reference sequences are sourced from the index for the second pass.

**Bias determinism.** Bias observed models (sequence/GC/positional) now
accumulate per-fragment mass in **fixed-point integer** buffers and materialize
once, so the trained models — and therefore the bias-corrected effective lengths
— are **byte-identical across thread counts**. (Previously the per-thread
fragment partition made the f64 accumulation order-dependent, giving tiny
thread-count-dependent wobble.) This is the one behavioral change versus 2.1.2
for bias runs; absolute accuracy versus ground truth is unchanged.

**Bias performance.** The effective-length correction was profiled and sped up
along the part that actually dominates:

- the sequence-bias per-position factor build folds the observed/expected models
  into a single precomputed `obs − exp` log-bias table (one context encode + one
  table sweep instead of two), giving **~1.9× on the seqBias correction sweep**;
- the length-sweep convolution for any non-GC bias (sequence and/or positional)
  is computed once as an FFT cross-correlation, `O(L log L)` instead of
  `O(L · n_len)`.

These are reassociation-level changes (the FFT and the fused table sum in a
different but equally valid order); they are deterministic and do not change
accuracy versus truth.

## RAD chunk compression (`--radCompress`)

RAD output can be compressed per chunk:

- **`--radCompress=lz4|zstd|none`** (default **`lz4`**), and **`--noCompressRad`**
  to force uncompressed;
- on by default for `--writeRad` output and the `--deterministic` intermediate.

Compression is **transparent and lossless** — chunks are decompressed in
`libradicl`'s reader thread, so every consumer (salmon, alevin-fry, piscem-infer)
is unchanged, and a RAD with no codec tag (every pre-2.2.0 / piscem-produced file)
reads as uncompressed automatically. On a 36 M-read human RAD, lz4 is ≈ 1.25×
smaller and zstd ≈ 1.9×; lz4 write is neutral-to-faster (parallel per-worker
compression plus less write I/O), while the requantified `quant.sf` is
byte-identical regardless of codec.

## Supporting improvements

- **Default minimizer length `m = 19`** (matches piscem) for new indices.
- **Lazy reference-sequence load** — sketch-mode runs without sequence/GC bias no
  longer load the multi-GB decoy `refseq` (~3 GB saved on a human + genome-decoy
  index).
- **Reusable per-thread mapping scratch** to cut per-read allocation.
- Diagnostic env vars promoted to **hidden flags** (`--biasSeedEMIters`,
  `--dumpBiasModels`); `--keepFixedFasta` added.
- Depends on the published **`libradicl` 0.14.0** (no local-path patch).

## Attribution

The RAD direction was opened by **@BenjaminDEMAILLE** in PRs
[#1033](https://github.com/COMBINE-lab/salmon/pull/1033),
[#1034](https://github.com/COMBINE-lab/salmon/pull/1034), and
[#1037](https://github.com/COMBINE-lab/salmon/pull/1037), which this release
**supersedes**. It adopts his core design — decouple map/quant through RAD, a
`--writeRad` producer and `--rad` consumer, piscem bulk/sketch RAD as a
first-class input, and deterministic RAD quantification. It differs in the
determinism mechanism (order-independent computation + a fixed, baked FLD rather
than an external merge sort of a RAD mapping store) and adds the
information-adaptive 1-vs-2-pass dispatch, bias-aware RAD, and chunk compression.

## Validation

- **Determinism:** RAD quant and `--deterministic` produce byte-identical
  `quant.sf` across thread counts (`-p 1` ≡ `-p 16`), with and without bias.
- **Compression:** requant of `none`/`lz4`/`zstd` RAD is byte-identical (sims and
  a 36 M-read human dataset); the full reader-decompression path and the
  codec-tag-presence behavior are covered by tests.
- **Backward compatibility:** real piscem `map-bulk` RAD (no codec tag)
  quantifies with mass conserved (Σ NumReads = fragment count).
- **Accuracy / concordance:** sim accuracy versus truth and 36 M cross-method
  concordance match the direct reads-mode path; piscem→RAD tracks salmon sketch.
- Reads, alignment (BAM/SAM), and RAD paths all exercised end-to-end across
  core/bias/bootstrap/Gibbs/`--skipQuant`; `cargo test --workspace`, `cargo fmt`,
  and `clippy` clean.
