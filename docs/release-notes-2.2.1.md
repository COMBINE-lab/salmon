# salmon 2.2.1

A security patch on top of 2.2.0.

1. **`lz4_flex` advisory fix** ([GHSA-vvp9-7p8x-rfvv](https://github.com/advisories/GHSA-vvp9-7p8x-rfvv))
   — the RAD chunk de/compression added in 2.2.0 pulled `lz4_flex` 0.10.0 (via
   `libradicl`), which carries a HIGH advisory: LZ4 decompression could leak
   uninitialized memory or reuse the output buffer. salmon now depends on
   `libradicl` 0.14.1, which bumps `lz4_flex` to 0.13.1 (the fix landed in
   0.11.0).

No functional changes. The RAD on-disk format is unchanged — LZ4's block format
is stable, so RAD files written by 2.2.0 still decompress — and the
compressed-RAD round-trip remains byte-identical (`quant.sf` unaffected). No
index rebuild required.
