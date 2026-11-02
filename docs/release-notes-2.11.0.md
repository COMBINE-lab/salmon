# salmon 2.11.0

A memory and input release:
- Selective alignment holds the reference at two bits per base. Peak memory falls by 36–44% on a GENCODE transcriptome and 29% with a GRCh38 genome decoy, at the same or slightly better speed.
- `quant` and `index` read from FIFOs and process substitution, so a run can be quantified while it downloads.
- A truncated FASTQ is caught while it is read, and the error names the file. Plain FASTQ is no longer read twice.

**Results are bit-for-bit identical to 2.10.0**, and no index rebuild is needed: the index format is unchanged.

## Memory: two-bit reference sequences

Selective alignment needs the reference sequences to extend and score candidate mappings, and it held them at one byte per base. They were the largest allocation of a run. They are now held at two bits per base:
- **Other bytes** (`N` and IUPAC codes) are kept as sorted runs. `refseq.bin` is uppercase, so soft-masking in the input FASTA creates no runs. A GRCh38 genome-decoy index has about a thousand.
- **The sequence is packed as it is read**, so one byte per base is never in memory.
- **Extension compares bases a 64-bit word at a time.** Runs are stepped over a byte at a time, so genome decoys, whose chromosomes all contain some, keep the fast path.
- **Alignment decodes only each gap or flank** it scores.

`--seqBias` and `--gcBias` read whole transcripts, so they keep the byte layout and see no saving. Sketch mode without those flags never loads the sequences, as before.

Peak memory of a full selective-alignment run, 2.10.0 → 2.11.0 (median of 3 runs; one run at one thread):

| workload | 1 thread | 8 threads | 16 threads |
|---|---|---|---|
| muscle (SRR13638690, 4.4M × 2×170 bp), GENCODE v50 | 2.39 → 1.34 GB | 2.66 → 1.65 GB | 3.03 → 1.95 GB |
| GEUVADIS (ERR188044, first 10M × 2×75 bp), Ensembl cDNA | 0.72 → 0.42 GB | 0.96 → 0.84 GB | 1.33 → 1.09 GB |
| GEUVADIS, Ensembl cDNA + GRCh38 genome decoy | | 9.89 → 6.95 GB | 10.20 → 7.26 GB |

Wall time is unchanged or slightly lower in every case (see [Speed](#speed)).

(#1265)

## Reading from streams

`-1/-2`, `-r`, `salmon index -t` and `-d` accept a FIFO or process substitution as well as a file:

```sh
salmon quant -i salmon_index -l A -o SRR000001 \
  -1 <(curl -sSf https://ftp.sra.ebi.ac.uk/vol1/fastq/.../SRR000001_1.fastq.gz) \
  -2 <(curl -sSf https://ftp.sra.ebi.ac.uk/vol1/fastq/.../SRR000001_2.fastq.gz)
```

2.10.0 refused this, with "Illegal seek" or a garbled first byte, and `salmon index -t <(...)` failed with "Failed to read the first two bytes".

- **`quant`:** salmon checked every uncompressed FASTQ for truncation in a separate pass before mapping, and on a pipe that pass took the first bytes. The check now happens in the FASTQ reader, as the input is read.
- **`index`:** a build opened each reference up to three times. It now opens a stream once. When any input is a stream, preprocessing feeds the reference digests from the same single read. Regular files keep their separate hashing pass, so file builds are unchanged in time and memory.
- **Unchanged:** indices built from streams are file-for-file identical to those built from the same files, decoy lists included.

A stream is decoded serially, since the parallel gzip decoder needs a file it can seek in. A plain stream cut exactly between two records cannot be told from a complete one. The new [Reading from streams](https://combine-lab.github.io/salmon/reference/cli/#reading-from-streams) section of the CLI reference covers this and other caveats.

(#1269)

## Truncated FASTQ

The FASTQ reader (paraseq-temp 0.6.3-pre.1, through piscem-rs 0.12.3) now refuses a final record that the input cut short:
- one or two lines of a record;
- a header without its newline;
- a quality line cut short.

Before, a cut one or two lines into the record, or inside its header, was dropped without an error, which is why salmon checked plain FASTQ in a separate pass. Format errors also name their input now. For a stream, the name includes the option it was given to:

```
mapping failed: reads_2.fq: The input ends inside a FASTQ record (48 bytes after the last complete record): it is probably truncated
mapping failed: /dev/fd/62 (-2): The input ends inside a FASTQ record (48 bytes after the last complete record): it is probably truncated
```

Consequences:
- **Plain FASTQ is read once.** On 10M uncompressed GEUVADIS pairs, a run takes 9.25 s instead of 10.33 s at 16 threads, and 16.7 s instead of 17.8 s at 8.
- **Paired inputs report the cut** rather than a mismatch in record counts.
- **A FASTQ that ends in a blank line runs.** The old check counted lines and rejected it.
- **Compressed inputs are unchanged:** a truncated gzip, bzip2, xz or zstd stream already failed in the decompressor.

(#1269, piscem-rs#19)

## Speed

Full quantification, 2.10.0 → 2.11.0, wall time in seconds. Both builds are from source with the same profile, on one machine (AMD EPYC 9575F) pinned to one socket. Each figure is the median of 3 runs, or a single run at one thread. `quant.sf` is identical between the two versions in all 38 runs.

| workload | 1 thread | 8 threads | 16 threads |
|---|---|---|---|
| muscle, selective alignment | 107.0 → 106.1 | 14.82 → 14.66 | 8.43 → 8.39 |
| GEUVADIS, selective alignment | 119.2 → 117.5 | 17.06 → 17.03 | 9.50 → 9.48 |
| GEUVADIS + genome decoy, selective alignment | | 26.43 → 26.55 | 15.36 → 15.25 |
| GEUVADIS, uncompressed FASTQ, selective alignment | | 17.77 → 16.73 | 10.33 → 9.25 |
| muscle, `--sketch` | | 5.46 → 5.45 | 3.40 → 3.40 |
| GEUVADIS, `--sketch` | | 6.83 → 6.84 | 4.05 → 4.09 |

Sketch mode is unchanged in time and memory: without `--seqBias` or `--gcBias` it never loads the reference sequences.

## Fixes

- **`--online --sketch --posBias` no longer panics.** The effective-length correction read each transcript's sequence only for its length, and sketch mode does not load the sequences. It now takes the length directly. (#1265)

## Documentation

- **The inferential replicates guide explains why salmon offers both bootstrap and Gibbs replicates.**
  - Each bootstrap replicate re-runs the optimizer from the same starting point. Where a gene's isoforms are nearly identical, replicates can therefore understate how loosely the data constrain the split between them, and the Gibbs sampler shows that ambiguity more fully.
  - A new section covers what replicates cannot measure: the variance between technical replicate libraries.

  (#1266)

## For library users

The published `salmon-*` crates change in these places:
- `salmon_core::RefProvider::ref_seq` is now `ref_view`, which returns a `RefView` over bytes or packed bases.
- `SalmonIndex::load_with_opts` takes `RefseqLoad::{None, Bytes, Packed}` instead of a `bool`.
- `salmon_model::corrected_effective_length_full` takes the transcript length. The sequence now travels in `BiasInputs::seq`, with the sequence-bias tables.

## Dependencies

- piscem-rs 0.12.3 and paraseq-temp 0.6.3-pre.1. The paraseq changes are offered upstream as noamteyssier/paraseq#102.
- Patch updates to `thiserror`, `serde_json`, `xxhash-rust`, `sha2`, `libc`, `realfft`, `ensure_simd` and `jiff`. Indices and `quant.sf` are byte-identical with them. (#1268)
- The Docker image builds with Rust 1.99. (#1262)
- Website: Astro 7.4, Starlight 0.42.8, smol-toml 1.9.1. (#1267, #1248)

## Checked

- `cargo test --release --workspace --locked`: 395 passed. Clippy (`-D warnings`) and fmt are clean.
- `quant.sf` is identical to 2.10.0 on every benchmark run above, with and without the genome decoy, in selective-alignment and sketch mode.
- Two-bit references:
  - `quant.sf` is identical to byte references on muscle and GEUVADIS, with and without the genome decoy, at 8 and 16 threads, and with `--seqBias --gcBias --posBias`.
  - On a synthetic index with `N` runs, IUPAC codes and a decoy, `quant.sf` and sorted `--writeMappings` output are also identical with:
    - the default options;
    - `--softclip --recoverOrphans`;
    - `--fullLengthAlignment`;
    - `--recoverOrphans --softclipOverhangs`.
- Streams: through FIFOs, `quant.sf` is identical to the file runs:
  - paired, single-end and multiple files per mate;
  - plain and gzip;
  - `--sketch` and `--writeMappings`.
- Every truncation tested, in files and streams, single-end and paired, is refused with the input named.

## Acknowledgements

Thanks to Benjamin Demaille for the two-bit reference packing (#1265), stream input (#1269) and the inferential replicates documentation (#1266).
