# salmon 2.2.1

The recommended 2.2.x release: it carries the full 2.2.0 feature set **plus** a
security fix for a `lz4_flex` advisory pulled in by the new RAD compression. If
you are on 2.2.0, upgrade. No index rebuild is required and `quant.sf` is
unchanged.

## Security fix (the reason for the patch)

- **`lz4_flex` advisory** ([GHSA-vvp9-7p8x-rfvv](https://github.com/advisories/GHSA-vvp9-7p8x-rfvv))
  — the RAD chunk de/compression introduced in 2.2.0 transitively pulled
  `lz4_flex` 0.10.0 (via `libradicl`), which carries a **HIGH** advisory: LZ4
  decompression could leak uninitialized memory or reuse the output buffer.
  salmon now depends on **`libradicl` 0.14.1**, which bumps `lz4_flex` to 0.13.1
  (the fix landed in 0.11.0). Drop-in: the RAD on-disk format is unchanged (LZ4's
  block format is stable, so RAD files written by 2.2.0 still decompress) and the
  compressed-RAD round-trip stays byte-identical.

---

The rest of this note summarizes the 2.2 feature set (introduced in 2.2.0 and
shipped here), since 2.2.1 is the release most users will land on.

## Decoupled mapping and quantification (`--writeRad`, `--rad`)

salmon can now separate mapping from quantification through the **RAD** format:

- **`--writeRad <PATH>`** writes per-fragment mappings (sketch or
  selective-alignment profile, chosen automatically). Quantification still runs;
  add **`--skipQuant`** to map only. The file is **piscem `map-bulk`-compatible**.
- **`--rad <PATH>`** quantifies a RAD file directly and in parallel — salmon-
  written **or** piscem `map-bulk`/sketch output. No `-i` is needed; reference
  names travel in the RAD header.

This lets you map once and quantify many times, quantify a shared piscem RAD with
salmon's EM, or split the two phases across machines. Quantifying a
piscem-produced RAD is fast and lean — **16.6 s and 0.52 GB RSS** for the quant
half of a 36 M-read human sample (mapping already done by piscem).

## Deterministic quantification (`--deterministic`) — reproducible *and* faster

RAD-mode quantification is **byte-identical across thread counts and runs**, and
`--deterministic` brings that guarantee to FASTQ input: it maps the reads once to
an intermediate RAD, then quantifies from it with a fixed fragment-length
distribution. Determinism comes from making the computation order-independent
(fixed-point `u128` equivalence-class weight accumulation; integer-count
fragment-length histograms built once in a fixed order) rather than from sorting
records.

**It is reproducibility at *negative* cost.** Because the single mapping pass
feeds a fast deterministic streaming requant instead of the heavier online
dual-phase inference, `--deterministic` is consistently **faster and lighter**
than the default online path. Measured on a 36 M-read GEUVADIS sample
(ERR188044, paired-end) against a decoy-aware human index, 16 threads:

| run | wall | peak RSS |
| --- | --- | --- |
| C++ salmon 1.9 / 2.1.2 (selective alignment) | 97.9 s | 9.9 GB |
| salmon 2.2 direct (selective alignment) | 94.7 s | 10.0 GB |
| **salmon 2.2 `--deterministic` (selective alignment)** | **90.7 s** | **9.8 GB** |
| salmon 2.2 direct (sketch) | 53.1 s | 9.8 GB |
| **salmon 2.2 `--deterministic` (sketch)** | **45.2 s** | **9.7 GB** |

So `--deterministic` runs **~4% faster than the default in selective-alignment
mode and ~15% faster in sketch mode** (and ~7% faster than C++ salmon in SA),
while also using slightly less memory — and unlike either non-deterministic path,
its `quant.sf` is bit-for-bit identical regardless of `-p`. Reproducibility is no
longer a speed/accuracy tradeoff.

`--deterministic` works with bias correction and does not require `-t` (the
reference sequences for the second pass are taken from the index). The
intermediate RAD is removed on success unless `--keepRad` (or `--writeRad PATH`).

## Bias-aware RAD, and deterministic + faster bias models

Sequence, GC, and positional bias correction now work on RAD input
(`--rad … --seqBias/--gcBias/--posBias`), sharing the same correction tail as the
alignment (BAM) path.

- **Deterministic bias models.** Bias observed models accumulate per-fragment
  mass in fixed-point integers, so the trained models — and the bias-corrected
  effective lengths — are **byte-identical across thread counts**.
- **Faster bias correction.** The sequence-bias per-position factor build folds
  the observed/expected models into a single precomputed `obs − exp` log-bias
  table (one context encode + one table sweep instead of two), giving **~1.9× on
  the seqBias correction sweep**; the length-sweep convolution for any non-GC
  bias is computed once as an FFT cross-correlation (`O(L log L)`). Both are
  reassociation-level and accuracy-vs-truth is unchanged.

## RAD chunk compression (`--radCompress`)

RAD output (and the `--deterministic` intermediate) can be compressed per chunk:
**`--radCompress=lz4|zstd|none`** (default `lz4`), `--noCompressRad` to force
uncompressed. Transparent and lossless — chunks are decompressed in the reader,
so all consumers (salmon, alevin-fry, piscem-infer) are unchanged, and a RAD with
no codec tag (any pre-2.2 / piscem file) reads as uncompressed automatically. On a
36 M-read human RAD, **lz4 is ≈ 1.25× smaller and zstd ≈ 1.9×**, with the
requantified `quant.sf` byte-identical regardless of codec.

## Compatibility

- **No index rebuild.** On-disk index format is unchanged from 2.1.x; `quant.sf`
  and the inferential-replicate formats are unchanged.
- A plain `salmon quant -i … -1 … -2 …` run is unaffected, except that the bias
  models are now thread-count-deterministic (accuracy vs truth unchanged).

## Attribution

The RAD direction was opened by **@BenjaminDEMAILLE** in PRs
[#1033](https://github.com/COMBINE-lab/salmon/pull/1033),
[#1034](https://github.com/COMBINE-lab/salmon/pull/1034), and
[#1037](https://github.com/COMBINE-lab/salmon/pull/1037), which the 2.2 series
supersedes. It adopts his core design — decouple map/quant through RAD, a
`--writeRad` producer and `--rad` consumer, and piscem bulk/sketch RAD as a
first-class input — and differs in the determinism mechanism (order-independent
computation + a fixed, baked FLD rather than an external merge sort) plus the
information-adaptive 1-vs-2-pass dispatch, bias-aware RAD, and chunk compression.
