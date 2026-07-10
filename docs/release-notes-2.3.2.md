# salmon 2.3.2

A performance-and-polish patch on 2.3.1. **No index rebuild is required and
`quant.sf` is unchanged.** Mapping output is byte-identical to 2.3.1 (verified on
real data), and the `--deterministic` guarantee — byte-identical results across
thread counts — is preserved and unchanged from 2.3.1.

## Faster mapping

- **Length-dispatched radix sort for uni-MEM grouping** (#1044). The per-read
  grouping of uni-MEMs by `(transcript, orientation)` now uses a radix sort sized
  to the key range instead of a comparison sort, shaving ~2% off mapping time on
  typical bulk runs. The change is purely a hot-path optimization: results are
  **output-neutral** (byte-identical `quant.sf` versus 2.3.1 at single thread, and
  byte-identical `--deterministic` output across thread counts).

## Richer `meta_info.json` and input diagnostics

- **More run metadata in `aux_info/meta_info.json`** (#1044). The file gains
  several additive fields — `detected_library_type`, `num_orphan`,
  `em_converged`, `num_em_iterations`, `range_factorization_bins`,
  `pos_bias_correct`, `peak_rss_kb`, `total_time_seconds` — and a structured
  `diagnostics` array. All existing keys are unchanged, so downstream parsers
  (e.g. tximport, pyroe) keep working.
- **Always-on library-type-mismatch detection** (#1044). salmon now always
  evaluates whether the observed fragment orientations match the requested
  library type and surfaces it via `detected_library_type` / `diagnostics`. It is
  a pure observer — it does not change quantification.
- The richer metadata and zero-cost bad-input diagnostics now also cover the
  alignment (`-a`) and `--rad` input paths.

## Configurable index-build scratch (`--sshashTmpDir`, `--ramLimit`)

- **`salmon index --sshashTmpDir <DIR>`** (#1048) controls where the SSHash
  minimizer sort writes its on-disk scratch. It now defaults to a `sshash_tmp`
  subdirectory of `--tmpdir` (or the index output directory when `--tmpdir` is
  unset) and is cleaned up after the build — fixing the previous behavior of
  leaving an empty `sshash_tmp/` directory in the current working directory. Use
  the flag to place the sort scratch on a separate or faster disk.
- **`salmon index --ramLimit <GiB>`** (#1048) caps the RAM the external minimizer
  sort may use (default `8`); a smaller value uses less memory but spills to disk
  sooner. This is a cohesive replacement for the fix in #1045 (thanks
  [@BenjaminDEMAILLE](https://github.com/BenjaminDEMAILLE), whose report motivated
  the change).

## Other changes

- **Clarified `salmon quant --threads` help text** to note the auxiliary I/O
  thread (#1047, #993) — thanks [@BenjaminDEMAILLE](https://github.com/BenjaminDEMAILLE).
- Internal: salmon is now available in **homebrew-core** (`brew install salmon`),
  so the separate brewsci auto-bump workflow was removed.

## Upgrading

Drop-in from 2.3.1: no index rebuild, and quantification results are unchanged.
The new `--sshashTmpDir` / `--ramLimit` options only affect `salmon index`
resource placement, not index content.
