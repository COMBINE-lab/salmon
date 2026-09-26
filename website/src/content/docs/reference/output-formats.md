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
| `library_types` | string[] | the library format that was **applied** — the one you gave with `-l`, or the detected one under `-l A`. For what the detector observed, read `detected_library_type` below. |
| `frag_dist_length` | int | number of FLD length bins |
| `frag_length_mean` / `frag_length_sd` | float | fragment length stats |
| `frag_length_source` | string | where the distribution above came from: `reads`, `alignments`, `rad_baked`, `rad_baked_prior`, `rad_derived`, or `prior` |
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

`num_bias_bins` is retained for C++ compatibility and is always `0`. In C++
salmon it reports the size of the legacy 4096-bin k-mer read-bias histogram
behind `observed_bias.gz`; this port replaced that model with the sequence, GC
and positional models, so the quantity it names does not exist here. The shapes
of the models that *do* run are reported in the three salmon-rs fields above.

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
| `seq_bias_context_length` | int | sequence-bias model context length, in nucleotides; omitted when `--seqBias` was not requested |
| `gc_bias_bins` | int | total GC-bias bins, i.e. `--numGCBins` × `--conditionalGCBins`; omitted when `--gcBias` was not requested |
| `pos_bias_bins` | object | positional-bias shape, `{ "length_classes", "bins_per_class" }`; omitted when `--posBias` was not requested |
| `total_time_seconds` | float | wall-clock seconds for the quantification call |
| `peak_rss_kb` | int | peak resident set size in KiB (Linux `VmHWM`; 0 elsewhere) |
| `diagnostics` | object[] | structured run diagnostics (see below) |

##### `diagnostics` — machine-readable run warnings

An array of `{ "code", "severity", "message" }` objects surfacing likely
bad-input conditions, so pipelines can gate on `code`/`severity` rather than
scraping the log. This array is the durable record: the same messages are
written to stderr as the run proceeds, but `logs/salmon_quant.log` carries only
a fixed run summary (version, times, counts, fragment-length statistics) and
not these diagnostics. `severity` is `"error"` or `"warning"`. Codes:

| `code` | severity | fires when |
| --- | --- | --- |
| `no_input_fragments` | error | zero fragments were processed (empty/unreadable input) |
| `no_fragments_mapped` | error | zero fragments mapped (reference almost certainly wrong) |
| `very_low_mapping_rate` | warning | < 10% of fragments mapped |
| `low_mapping_rate` | warning | < 30% of fragments mapped |
| `library_type_mismatch` | warning | an explicit `-l` disagrees with the observed library format (reads / `--rad`) |
| `no_concordant_mappings` | warning | no fragment's observed format agreed with an unstranded expected type — for a paired library, check the reads are properly paired |
| `strand_bias_unstranded` | warning | strand bias > 1% under an unstranded library type (see `strand_mapping_bias` in `lib_format_counts.json`) |
| `high_incompatible_fraction` | warning | more than 5% of judged fragments were incompatible with the expected library type |
| `bam_score_tags_unusable` | warning | deterministic alignment mode scored by `AS`, but the BAM's `AS` tags were missing, partial, or constant — placements were weighted uniformly; consider `--errorModel -t` |

An empty array means no red flags were detected. The library-format codes are
emitted by every mode that writes `lib_format_counts.json`, from the same
conditions the log warnings use, so the log and the report cannot disagree.

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
`num_incompatible_fragments`, `num_assigned_fragments`,
`num_frags_with_concordant_consistent_mappings`,
`num_frags_with_inconsistent_or_orphan_mappings`, `strand_mapping_bias`, and
the raw per-format observed counts (`IU`, `ISF`, `ISR`, `OU`, `OSF`, `OSR`,
`MU`, `MSF`, `MSR`, `U`, `SF`, `SR`).

`num_incompatible_fragments` is an addition to the C++ field set (every other
field keeps its C++ name and meaning). It counts fragments that mapped but had
**no** placement on the strand the declared library type expects:

- `num_compatible_fragments`: fragments with at least one strand-compatible
  mapping.
- `num_incompatible_fragments`: fragments that mapped only on the wrong strand.
  With the default `--incompatPrior 0` these are discarded, so they are **not**
  included in `num_assigned_fragments`; with `--incompatPrior > 0` they are kept
  at reduced weight and are included.
- `compatible_fragment_ratio`: `num_compatible_fragments` over the sum of the
  two, i.e. over the fragments the strand filter actually judged.
- Per-format keys (`ISF`, `ISR`, …): how many mapped fragments had at least one
  placement *observed* as that format, tallied from the raw placements before
  any filtering (a fragment with placements of several formats counts under
  each). `expected_format` is the resolved type the strand filter used — the
  detected type under `-l A`, else the declared one.
- `num_frags_with_concordant_consistent_mappings` /
  `num_frags_with_inconsistent_or_orphan_mappings`: derived from that
  histogram, with the C++ meaning — fragments whose observed format agrees
  with the expected one (exactly the expected format for a stranded type,
  either strandedness of the expected orientation for an unstranded one)
  versus everything else, orphaned mates included.
- `strand_mapping_bias`: of the two strandedness variants of the expected
  orientation (e.g. `ISF`/`ISR` for an inward type), the fraction observed as
  the sense variant. Near 0.5 for a genuinely unstranded library; salmon warns
  when an unstranded run strays more than 1% from it, the classic
  double-stranded-contamination / mis-declared-strandedness check.

On a stranded protocol the wrong-strand fraction is the cheapest available
measure of double-stranded input (genomic DNA carry-over from incomplete DNase
digestion): such fragments map to either strand with equal probability, so an
excess of incompatible fragments implies a comparable amount of contamination
sitting *on* the expected strand, where nothing can distinguish it from RNA.
Comparing the ratio across samples of one cohort flags the affected libraries
without a separate strandedness pass over an alignment.

Both counts cover only fragments compared against a known expected format. Under
`-l A` the fragments consumed before the library type is inferred are in neither
tally, so their sum can be lower than the mapped total; under an unstranded type
every fragment is compatible and the ratio is 1.

The fields are measured in every mode that writes this file: read-based
quantification (selective alignment and `--sketch`, which applies the same
strand filter), alignment-based quantification (`-a`), RAD-input
quantification, and `--deterministic` (whose phase-2 requant rewrites the file
with tallies measured from the RAD's stored orientations). In every mode the
fragment is judged before the filter acts on it, so on a stranded `-l` a
wrong-strand library shows both signals at once: a low
`compatible_fragment_ratio` explaining a collapsed mapping rate.

`-l A` is likewise resolved from the data in every mode: the one-pass reads
modes and online alignment mode sample a prefix with the library-type detector
(detect, lock in, filter the rest of the run), while `--deterministic` and
RAD-input quantification detect from the whole run's records and filter the
quantification pass. `expected_format` always reports the resolved type.

One caveat for alignment input (`-a`): detection treats the reported
alignments as an unfiltered sample of the library, but they are whatever the
aligner chose to report. If the aligner was itself configured to keep only one
orientation or strand, detection mirrors that upstream filter rather than the
library — it cannot conclude `IU` when wrong-strand alignments were already
excluded before salmon ever saw them. For a BAM/SAM filtered this way, pass
the library type explicitly; salmon logs a warning to this effect whenever
`-l A` is used with alignment input.

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

Three salmon-relevant details on top of the base format:

- **Baked header tags.** A salmon-written RAD records, as file-level tags, an
  order-independent fragment-length distribution, initial abundances, and the
  resolved library format, plus a `baked_flags` marker. A `--rad` reader consumes
  these to quantify in a single pass and to apply `-l A` concordance filtering
  without re-inference; a piscem RAD carries none of them and is handled with an
  extra derivation pass.
- **Mapping-pass provenance.** A RAD holds only the fragments that *mapped*, so
  nothing in the records can say how many were *observed*. Without that, a
  requant reports a 100% mapping rate by construction. salmon therefore records
  what its mapping pass saw, as file-level tags:

  | tag | meaning |
  | --- | --- |
  | `num_processed` | fragments observed, mapped or not |
  | `num_dovetail_fragments` | fragments whose mappings were dovetailed |
  | `num_fragments_filtered_vm` | fragments dropped by the score filter |
  | `num_alignments_below_threshold_vm` | below-threshold alignments among mapped fragments |
  | `num_decoy_fragments` | fragments whose best mapping was to a decoy |
  | `index_seq_hash`, `index_name_hash`, and the `512` / `decoy` variants | identity of the index the mappings were made against |
  | `keep_duplicates` | how that index was built (omitted if the index predates recording it) |
  | `mapping_type` | `mapping`, `pseudo` or `alignment` |
  | `source_programs` | the source BAM's `@PG` lines, for an alignment-derived RAD |

  The counters are only final at end of pass, so their slots are reserved and
  backpatched at finalize; the `BAKED_MAP_COUNTERS` bit in `baked_flags` says
  they were filled, since a reserved slot and a genuine count of zero are
  otherwise indistinguishable. The index identity is known up front and written
  directly, flagged by `BAKED_INDEX_PROV`. `mapping_type` is recorded rather than
  inferred from the record profile, because a BAM-derived RAD is written in the
  selective-alignment profile although its fragments came from an aligner.

  Quantifying a RAD that lacks these — a piscem RAD, or one written by salmon
  before 2.5.0 — still works. salmon warns, names the affected fields, and
  records the gap in `meta_info.json` (see below).

- **Alignment provenance.** A RAD built from a BAM carries that BAM's `@PG`
  header chain, so a requant can report *how* the alignments were made
  — the aligner, its version and its command line — rather than only that they
  were alignments. It reappears as `alignment_provenance` in `meta_info.json`,
  for both `-a` and a requant of a BAM-derived RAD.

  The chain is stored as one string tag holding the `@PG` lines joined by
  newlines. RAD imposes no restriction on the content — a string tag is an
  explicit length followed by bytes — so the only constraints come from that
  representation, and both are handled rather than assumed away: tabs, newlines
  and backslashes inside a value are escaped so the framing stays unambiguous and
  reversible, and because the length prefix is a `u16`, a chain longer than 64 KiB
  has whole trailing records dropped (with a warning) rather than being allowed to
  wrap the length and corrupt the file.
- **Chunk compression (`chunk_codec`).** With compression enabled (the default;
  see `--radCompress`), each chunk's record payload is compressed independently
  and the `nbytes` field is the **compressed** size; the uncompressed size is
  carried inside the payload. A file-level `chunk_codec` byte tag selects the
  codec — `1` = LZ4, `2` = zstd. **A missing `chunk_codec` tag means codec `0`
  (uncompressed)**, so every RAD produced before this feature, and every piscem
  RAD, reads as uncompressed automatically. Decompression happens in the reader,
  so record parsing downstream is identical for compressed and uncompressed RADs.

### Metadata completeness

Every `meta_info.json` carries `meta_info_complete`. When `false`,
`incomplete_meta_info_fields` lists each field a complete description would have
carried, together with the upstream that could not supply it and why:

```json
"meta_info_complete": false,
"incomplete_meta_info_fields": [
  {
    "field": "keep_duplicates",
    "source": "index",
    "reason": "this index predates recording how it was built; rebuild the index to record it"
  }
]
```

`source` is `index`, `rad` or `bam`. The point is to distinguish "this run
observed nothing" from "nobody could tell this run", and to name what to fix.

`@PG` is optional in SAM and many tools write none, so an alignment input with no
`@PG` chain is an ordinary outcome rather than an error — but it is still recorded
(`alignment_provenance`, source `bam`), as is a chain that names its programs
without carrying any `CL` field
(`alignment_provenance[].command_line`). Both are reported identically whether
the alignments are read from a BAM directly or through a BAM-derived RAD.
Fields that are *not applicable* are not listed: alignment mode has no salmon
index, so its absent index hashes are by design rather than a gap.

Values that are unknown are omitted rather than guessed, so `keep_duplicates`
being absent means exactly that — not `false`.

## Mapping alignment output

`salmon quant --writeMappings <FILE>` (alias `--writeSam`) writes the retained
mappings as unsorted SAM records. `--writeBam <FILE>` writes semantically
equivalent records as BGZF-compressed BAM. The options are mutually exclusive,
and `--writeSam` names the SAM format explicitly so both formats have a
format-named flag. Both share one header and one record builder, so the two
formats cannot drift apart.

Both options apply to the read-mapping path only. In alignment mode (`-a`) and
RAD-input mode (`--rad`) they are accepted but ignored, with a warning.

### What a record contains

The header is unchanged: `@HD`, one `@SQ` per reference with its length, and
`@PG`.

Each record carries:

| Field | Value |
|-------|-------|
| `MAPQ` | Derived from the number of placements, on STAR's scale: 255 for a uniquely placed fragment, 3 for two placements, 1 for three or four, 0 beyond. `-q 255` therefore means "uniquely placed", as it does for STAR output. |
| `CIGAR` | A base-level alignment with real `I` and `D` operations, and `S` where a read overhangs a transcript end. The placement is re-aligned with traceback when the record is written, which is work a run without mapping output never does. |
| `NM` / `MD` | Edit distance and the reference bases at each mismatch, matching what `samtools calmd` recomputes from the input FASTA. Indexing has to replace non-ACGT bases with ACGT, since the k-mer structure cannot hold them, but the build records where it did so and the record reports the `N` that was there rather than the substitute. |
| `NH` / `HI` | How many placements this fragment has, and which one this record is. |
| `AS` | The score of the alignment in this same record. |
| `XT` | `T` for a transcript placement, `D` for a decoy. In practice always `T`: a fragment whose best placement is a decoy is dropped by scoring before any record is written, so a decoy-aware index changes which fragments appear, not what `XT` says about them. |

`QUAL` is still `*`: salmon does not carry base qualities into the record
builder. That is an omission the format allows, unlike a CIGAR that does not
describe the alignment it sits beside.

In sketch mode (`--sketch`) there is no per-placement alignment score:
pseudoalignment maps by k-mer/equivalence-class compatibility without computing
one, so `AS` is reported as `0` and every reported mapping is a co-equal best
mapping for the read as assessed by the algorithm. A meaningful `AS` appears
only under selective alignment (the default), which computes a banded alignment
score per candidate.

### `MAPQ`

SAM's `MAPQ` is a phred-scaled probability that a placement is wrong, and salmon
does not compute one: its design is to keep every plausible placement and let the
EM apportion them. What it does know is how many placements a fragment has, which
is what STAR reports, so `MAPQ` follows STAR's mapping:

| placements (`NH`) | `MAPQ` |
|---|---|
| 1 | 255 |
| 2 | 3 |
| 3 or 4 | 1 |
| 5 or more | 0 |
| none (a record with no placement) | 0 |

Matching STAR matters because the tools downstream of a salmon BAM are tuned to
it: `samtools view -q 255` means "uniquely placed" to all of them. Before this,
every record carried a constant `1`, so that filter discarded the entire file.

### Compression threads

`--writeBam` compresses BGZF blocks on a small pool of threads that run
alongside the `-p` mapping threads. `--bamCompressThreads <N>` sets the pool
size; left unset, salmon derives it from measured throughput.

One worker compresses about 165 MiB/s of BAM records, and one mapping thread
produces at most about 52 MiB/s of them, so roughly **one compression worker per
3 mapping threads** is the point at which compression stops limiting the run.
That is the default, capped at 8 workers — already around three times the
fastest record production measured at any thread count.

The default deliberately errs upward but stops there. The two failure modes are
not symmetric: one worker too few can halve throughput, because record output
pins to a single worker's 165 MiB/s, while surplus workers simply block on an
empty queue and cost nothing measurable. Going beyond the balance point,
however, spends cores that would otherwise be mapping reads. Raise it if you are
writing BAM to unusually fast storage and have cores to spare; there is no
reason to lower it.

### Record order

Mapping workers run in parallel, and the output makes exactly one ordering
guarantee:

- **All records for a fragment are contiguous.** A worker encodes an entire
  fragment before deciding to hand off its buffer, and buffers are written
  whole, so a fragment is never split or interleaved with another worker's
  records. Tools that stream a file and group by read name — the common reason
  to want `samtools collate` — can read the output directly.
- **Nothing else is ordered.** Fragments appear in whatever order workers finish
  their buffers, which varies with thread count and scheduling. Sort the result
  if a downstream tool needs coordinate or query-name order.

The second point is deliberate rather than incidental: imposing a global record
order would require workers to wait for one another. Only the parts that must be
serial (copying bytes into BGZF blocks, and writing compressed blocks in
submission order) are, while record encoding and block compression both run
fully in parallel.
