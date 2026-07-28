---
title: Output format specification
description: Every file salmon 2.0 writes, including exact byte layouts and downstream compatibility.
---

This page specifies **every file** salmon 2.0 writes to a `quant` output
directory, including exact byte layouts for the binary files, so that downstream
tooling (tximport, fishpond/swish, custom scripts) can rely on a stable contract.

Files fall into two compatibility classes:

- **C++-compatible** — byte/loadable-compatible with C++ salmon (≤ 1.12.0) so
  existing downstream tools keep working unchanged. These formats are stable and
  changes to them are breaking.
- **Documented Rust format** — diagnostic dumps not consumed by standard
  downstream tools. The layout is documented here and is stable within the 2.x
  series, but is not promised to match C++ salmon byte-for-byte.

All multi-byte binary values are **little-endian**. All `*.gz` files are standard
gzip streams; "uncompressed payload" describes the bytes *after* gunzip.

```
<output_dir>/
├── quant.sf                       (C++-compatible)   TSV
├── cmd_info.json                  (C++-compatible)   JSON
├── lib_format_counts.json         (C++-compatible)   JSON
├── logs/
│   └── salmon_quant.log           (Rust format)      text
├── libParams/
│   └── flenDist.txt               (C++-compatible)   TSV (one line)
└── aux_info/
    ├── meta_info.json             (C++-compatible)   JSON
    ├── ambig_info.tsv             (C++-compatible)   TSV
    ├── fld.gz                     (Rust format)      gzip → i32[]
    ├── eq_classes.txt.gz          (C++-compatible)   gzip → text   [--dumpEq]
    ├── observed_bias.gz           (Rust format)      gzip → i32[]  [legacy stub]
    ├── observed_bias_3p.gz        (Rust format)      gzip → i32[]  [legacy stub]
    ├── expected_bias.gz           (Rust format)      gzip → f64[]  [legacy stub]
    ├── obs5_seq.gz obs3_seq.gz    (Rust format)      gzip → f64[]  [--seqBias]
    ├── exp5_seq.gz exp3_seq.gz    (Rust format)      gzip → f64[]  [--seqBias]
    ├── obs_gc.gz  exp_gc.gz        (Rust format)      gzip → f64[]  [--gcBias]
    ├── obs5_pos.gz obs3_pos.gz     (Rust format)      gzip → header+f64[]  [--posBias]
    ├── exp5_pos.gz exp3_pos.gz     (Rust format)      gzip → header+f64[]  [--posBias]
    └── bootstrap/                  [--numBootstraps / --numGibbsSamples]
        ├── names.tsv.gz            (C++-compatible)   gzip → text
        └── bootstraps.gz           (C++-compatible)   gzip → f64[]
```

## C++-compatible files

### `quant.sf` — primary abundance table (TSV)

A header line followed by one row per **non-decoy** transcript, in index order
(decoy references, index ≥ first-decoy, are excluded — matching C++ salmon).

```
Name<TAB>Length<TAB>EffectiveLength<TAB>TPM<TAB>NumReads
```

| Column | Type | Notes |
| --- | --- | --- |
| `Name` | string | transcript id, as in the index |
| `Length` | integer | transcript length in nucleotides |
| `EffectiveLength` | float | `--sigDigits` decimals (default 3) |
| `TPM` | float | fixed 6 decimals |
| `NumReads` | float | `--sigDigits` decimals (default 3) |

This is the file tximport reads by default; it is a drop-in replacement for the
C++ output.

### `aux_info/meta_info.json` — run metadata (JSON)

Pretty-printed JSON. tximport keys off several fields (`num_bootstraps`,
`num_valid_targets`, `eq_class_properties`, `index_seq_hash`, …). Fields:

| Field | Type | Meaning |
| --- | --- | --- |
| `salmon_version` | string | e.g. `"2.0.0"` |
| `samp_type` | string | `"bootstrap"`, `"gibbs"`, or `"none"` |
| `opt_type` | string | `"em"` or `"vb"` |
| `num_libraries` | int | currently `1` |
| `library_types` | string[] | detected/declared library type(s) |
| `frag_dist_length` | int | number of FLD length bins |
| `frag_length_mean` / `frag_length_sd` | float | observed fragment length stats |
| `seq_bias_correct` / `gc_bias_correct` / `pos_bias_correct` | bool | bias correction enabled (`--seqBias` / `--gcBias` / `--posBias`) |
| `mapping_type` | string | `"mapping"` (SA) or `"pseudo"` (sketch) |
| `keep_duplicates` | bool | index built with `--keepDuplicates` |
| `index_seq_hash` / `index_name_hash` | string | SHA-256 (hex) of reference seqs / names |
| `index_seq_hash512` / `index_name_hash512` | string | SHA-512 variants |
| `index_decoy_seq_hash` / `index_decoy_name_hash` | string | decoy hashes |
| `num_valid_targets` | int | non-decoy transcript count |
| `num_decoy_targets` | int | decoy count |
| `num_eq_classes` | int | equivalence-class count |
| `eq_class_properties` | string[] | e.g. `["gzipped"]` or `["range_factorized","gzipped"]` |
| `length_classes` | int[] | length-class boundaries (u32) |
| `num_processed` / `num_mapped` | int | fragments observed / mapped |
| `percent_mapped` | float | |
| `num_decoy_fragments` | int | |
| `num_bootstraps` | int | inferential-replicate count (0 if none) |
| `start_time` / `end_time` | string | asctime |

(Plus `quant_errors`, `num_bias_bins`, `serialized_eq_classes`,
`num_dovetail_fragments`, `num_fragments_filtered_vm`,
`num_alignments_below_threshold_for_mapped_fragments_vm`, and `call`.)

#### salmon-rs extensions (beyond the C++ field set)

These fields are added by the Rust implementation; downstream tools that only
know the C++ schema simply ignore them. All are computed from end-of-run
aggregates (no per-fragment cost) and are emitted for every mode — reads,
`-a` alignment, and `--rad`.

| Field | Type | Meaning |
| --- | --- | --- |
| `num_orphan` | int | mapped fragments placed as orphans (one mate only); reads mode |
| `range_factorization_bins` | int | range-factorization bins used (0 = disabled) |
| `num_em_iterations` | int | EM/VBEM iterations actually run |
| `em_converged` | bool | whether the relative-difference criterion was met before the iteration cap |
| `detected_library_type` | string \| null | library format observed by the auto-detector (reads / `--rad`); `null` if not observed (e.g. `-a` transcriptomic BAM) |
| `total_time_seconds` | float | wall-clock seconds for the quantification call |
| `peak_rss_kb` | int | peak resident set size in KiB (Linux `VmHWM`; 0 elsewhere) |
| `diagnostics` | object[] | structured run diagnostics (see below) |

##### `diagnostics` — machine-readable run warnings

An array of `{ "code", "severity", "message" }` objects surfacing likely
bad-input conditions, so pipelines can gate on `code`/`severity` rather than
scraping the log (the same messages are also written to
`logs/salmon_quant.log`). `severity` is `"error"` or `"warning"`. Codes:

| `code` | severity | fires when |
| --- | --- | --- |
| `no_input_fragments` | error | zero fragments were processed (empty/unreadable input) |
| `no_fragments_mapped` | error | zero fragments mapped (reference almost certainly wrong) |
| `very_low_mapping_rate` | warning | < 10% of fragments mapped |
| `low_mapping_rate` | warning | < 30% of fragments mapped |
| `library_type_mismatch` | warning | an explicit `-l` disagrees with the observed library format (reads / `--rad`) |

An empty array means no red flags were detected.

### `aux_info/ambig_info.tsv` — per-transcript read partition (TSV)

Header `UniqueCount<TAB>AmbigCount`, then one row per quantified transcript in
index order (the same `num_valid_targets` set as `quant.sf`). `UniqueCount` =
fragments mapping uniquely to that transcript; `AmbigCount` = fragments mapping
ambiguously.

### `libParams/flenDist.txt` — fragment-length distribution (TSV, one line)

A single line: the normalized fragment-length PMF as tab-separated values in
scientific notation, one value per length bin from 0 to the max fragment length.

### `aux_info/eq_classes.txt.gz` — equivalence classes (gzip text) [--dumpEq]

Only written with `--dumpEq` or `--dumpEqWeights`. Uncompressed payload is
salmon's text format:

```
num_transcripts
num_eq_classes
<transcript name 0>
<transcript name 1>
...                      (num_transcripts names, index order)
<class line>             (num_eq_classes lines)
...
```

- **`--dumpEq`** (collapsed by transcript set): each class line is
  `groupSize` TAB `tid_0` TAB … TAB `tid_{g-1}` TAB `count`.
- **`--dumpEqWeights`**: interleaves the per-transcript combined weights before
  the count: `groupSize` TAB `tid_0 … tid_{g-1}` TAB `w_0 … w_{g-1}` TAB `count`.

Transcript ids index into the name list above.

### `aux_info/bootstrap/` — inferential replicates (gzip) [--numBootstraps / --numGibbsSamples]

Written when bootstrap or Gibbs sampling is requested; byte-compatible with C++
salmon's `GZipWriter::writeBootstrap<double>`.

- **`names.tsv.gz`** — uncompressed payload is the transcript names,
  **tab-separated** on a single line terminated by a newline, in the same order
  and set as the `quant.sf` rows: decoy references are **excluded** and sub-`k`
  "short" transcripts are **included** (always 0). Aligns positionally with both
  `quant.sf` and `bootstraps.gz`.
- **`bootstraps.gz`** — uncompressed payload is raw **`f64` little-endian** values
  with no header: `n_replicates` samples written contiguously, each sample being
  `num_valid_targets` values in `quant.sf` row order (decoys excluded, shorts
  included). Total length = `n_replicates × num_valid_targets × 8` bytes.
  `samp_type` in `meta_info.json` records whether the replicates are `bootstrap`
  or `gibbs`. Available in both mapping-based and **alignment-based** (`-a`)
  quantification.

### `cmd_info.json` / `lib_format_counts.json` (JSON)

`cmd_info.json` records the invocation: `salmon_version`, `index`, `libType`,
`output`, `mates1`, `mates2`, `unmatedReads`, `threads`, `sketch`.

`lib_format_counts.json` records library-format detection: `read_files`,
`expected_format`, `compatible_fragment_ratio`, `num_compatible_fragments`,
`num_assigned_fragments`, `num_frags_with_concordant_consistent_mappings`,
`num_frags_with_inconsistent_or_orphan_mappings`, `strand_mapping_bias`.

## Documented Rust-format files (diagnostic)

These are bias-model and FLD diagnostic dumps, not read by standard downstream
tools. Layouts are documented for completeness.

### `aux_info/fld.gz`

Uncompressed payload: an array of **`i32` little-endian** counts, one per
fragment-length bin. Where C++ salmon draws 10,000 Monte-Carlo samples from the
log-PMF, the port writes the deterministic expected histogram
`round(10000 · pmf[len])` (same type and layout).

### Sequence-bias dumps [--seqBias]

`obs5_seq.gz`, `obs3_seq.gz`, `exp5_seq.gz`, `exp3_seq.gz`: uncompressed payload
is an array of **`f64` little-endian** values — the flattened observed/expected
5′ and 3′ sequence-bias context tables.

### GC-bias dumps [--gcBias]

`obs_gc.gz`, `exp_gc.gz`: uncompressed payload is an array of **`f64`
little-endian** values (observed/expected GC mass bins).

### Positional-bias dumps [--posBias]

`obs5_pos.gz`, `obs3_pos.gz`, `exp5_pos.gz`, `exp3_pos.gz`: uncompressed payload
is a header followed by the model bins:

```
[u32 num_models][u32 bins_per_model] then num_models × bins_per_model f64 LE, row-major
```

Each "model" is one length-class's positional distribution.

### Legacy seq-bias stubs

`observed_bias.gz` (`i32 [0]`), `observed_bias_3p.gz` (`i32 [0]`), and
`expected_bias.gz` (`f64 [1.0]`) are single-element stubs preserved for the
legacy simple-count seq-bias model the port does not implement (it uses the
`SBModel` context model instead). They exist so consumers expecting these
filenames do not error.

### `logs/salmon_quant.log`

A concise human-readable run summary (version, start/end time, library type,
mapping type, observed/mapped fragment counts, mapping rate, equivalence-class
count, fragment-length mean/sd). Downstream tools key off the JSON metadata, not
this log.

## Index directory

The 2.0 index is the piscem-rs format and is **not** compatible with C++ salmon
(pufferfish) indices — they must be rebuilt. Pointing 2.0 at a C++ index (or C++
salmon at a 2.0 index) produces a clear, actionable error. See
[what changed in 2.0](../../migrating/from-cpp/).

## RAD output (`--writeRad`)

`--writeRad <PATH>` (and the `--deterministic` intermediate) write a **RAD** file
— the Reduced Alignment Data format defined by
[libradicl](https://github.com/COMBINE-lab/libradicl) and shared with piscem and
alevin-fry. The base format (prelude, file-level tag definitions, and the
sequence of `[u32 nbytes][u32 nrec][records…]` chunks) is libradicl's contract;
salmon writes a bulk profile that piscem `map-bulk` can also produce and read.
See the [RAD I/O guide](../../guides/rad-and-determinism/) for usage.

Two salmon-relevant details on top of the base format:

- **Baked header tags.** A salmon-written RAD records, as file-level tags, an
  order-independent fragment-length distribution, initial abundances, and the
  resolved library format, plus a `baked_flags` marker. A `--rad` reader consumes
  these to quantify in a single pass and to apply `-l A` concordance filtering
  without re-inference; a piscem RAD carries none of them and is handled with an
  extra derivation pass.
- **Chunk compression (`chunk_codec`).** With compression enabled (the default;
  see `--radCompress`), each chunk's record payload is compressed independently
  and the `nbytes` field is the **compressed** size; the uncompressed size is
  carried inside the payload. A file-level `chunk_codec` byte tag selects the
  codec — `1` = LZ4, `2` = zstd. **A missing `chunk_codec` tag means codec `0`
  (uncompressed)**, so every RAD produced before this feature, and every piscem
  RAD, reads as uncompressed automatically. Decompression happens in the reader,
  so record parsing downstream is identical for compressed and uncompressed RADs.
## Mapping alignment output

`salmon quant --writeMappings <FILE>` writes the retained mappings as unsorted
SAM records. `--writeBam <FILE>` writes semantically equivalent records as
BGZF-compressed BAM. The options are mutually exclusive. Mapping workers run in
parallel, so record order is unspecified; sort the result when a downstream tool
requires coordinate or query-name order.
