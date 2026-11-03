# salmon 2.11.0

A memory, input and correctness release:
- Selective alignment holds the reference at two bits per base. Peak memory falls by 36–45% on a GENCODE transcriptome and 29% with a GRCh38 genome decoy, at the same or slightly better speed.
- `quant` and `index` read from FIFOs and process substitution, so a run can be quantified while it downloads.
- A truncated FASTQ is caught while it is read, and the error names the file. Plain FASTQ is no longer read twice.
- Fixes from an audit of 2.10.0 (#1283): wrong counts in a few situations, mostly in alignment mode (`-a`), and inputs and options that used to give a wrong result without an error.

**What changes in results.** Default runs from reads give the same results as 2.10.0, except in three cases:
- **Stranded, paired-end libraries with selective alignment:** a fragment whose mates map, as orphans, to different transcripts is now judged on each mate's own strand. On the muscle sample below, `-l ISR`, 14 more fragments map (4,139,973 instead of 4,139,959), and 1,130 transcripts change NumReads, moving 2,228 of 4.14M reads (0.05%).
- **Sketch mode with mates shorter than k:** these pairs no longer pick up the previous fragment's hits.
- **SAM/BAM output (`--writeMappings`):** orphans from read 2 are written on their own strand, and under `--sketch` they keep their `AS` score.

Options and inputs outside the defaults change in more places, listed under [Other fixes](#other-fixes):
- `--seqBias` and `--posBias`, for transcripts shorter than every fragment length;
- `--recoverOrphans`;
- single-end reads under `-l A`;
- `--ignoreTxVersion`.

Alignment mode changes in more places, listed under [Alignment mode](#alignment-mode-a). Some runs that used to succeed now stop with an error or exit non-zero: see [Behaviour changes](#behaviour-changes).

No index rebuild is needed: the index format is unchanged.

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
| muscle (SRR13638690, 4.4M × 2×170 bp), GENCODE v50 | 2.42 → 1.34 GB | 2.64 → 1.64 GB | 3.03 → 1.95 GB |
| GEUVADIS (ERR188044, first 10M × 2×75 bp), Ensembl cDNA | 0.71 → 0.42 GB | 1.00 → 0.84 GB | 1.30 → 1.10 GB |
| GEUVADIS, Ensembl cDNA + GRCh38 genome decoy | | 9.88 → 6.94 GB | 10.20 → 7.25 GB |

Wall time is the same, within run-to-run noise, or lower (see [Speed](#speed)).

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
- **Plain FASTQ is read once.** On 10M uncompressed GEUVADIS pairs, a run takes 9.36 s instead of 10.30 s at 16 threads, and 16.7 s instead of 17.8 s at 8.
- **Paired inputs report the cut** rather than a mismatch in record counts.
- **A FASTQ that ends in a blank line runs.** The old check counted lines and rejected it.
- **Compressed inputs are unchanged:** a truncated gzip, bzip2, xz or zstd stream already failed in the decompressor.

(#1269, piscem-rs#19)

## Behaviour changes

Each of these used to run to the end and give a result that was wrong, or easy to misread. Each now stops with an error that says what to fix, or, for the exit code, still writes its output.

- **A run in which no fragment is processed or mapped exits non-zero.** These are the `no_input_fragments` and `no_fragments_mapped` diagnostics. The output files are still written. salmon 1.x also exited 1 here, below its `--minAssignedFrags` threshold. A workflow that runs samples which may map nothing, such as negative controls, must now tolerate the failure. (#1277)
- **`salmon index` refuses:**
  - a decoy list naming a sequence that is not in the FASTA. A list made with `grep '^>'` and never cleaned up used to build with no decoys at all, and the genome was quantified as a transcript. (#1275)
  - a reference name that appears more than once, after the `--gencode` cut. `quant.sf` used to get two rows with that name. (#1288)
  - `-k` above 31, or `-m` not below k, before writing anything. k = 33 used to pass the check and then panic, part-way through writing the index. (#1286)
- **`salmon quant` refuses out-of-range values** for `--vbPrior` (must be above 0), `--fldMax` (at least 1), `--forgettingFactor` (in (0.5, 1]), and the fractions `--incompatPrior`, `--minAlnProb`, `--minScoreFraction`, `--consensusSlack`, `--orphanChainSubThresh`, `--preMergeChainSubThresh` and `--postMergeChainSubThresh` (in [0, 1]). `--scoreExp`, `--decoyThreshold` and `--juncMissDiscount` must be finite and at least 0. `--vbPrior NaN`, for instance, used to run and drop about a third of the reads. (#1286)
- **`salmon quantmerge` refuses** two samples with the same name, such as `--quants a/quant b/quant`, an empty name, a row with fewer than five columns, and a target repeated in one file. Each used to give a matrix that looked complete. Pass `--names` to name samples whose directories end alike. (#1287)
- **`-a` refuses a BAM whose mates are not next to each other**, when at least half the fragments are split, and warns from 5%. Only a header saying `SO:coordinate` was refused before, so a position-sorted file with `SO:unknown`, or no `@HD` line, doubled every count. Group such a file with `samtools collate` or `samtools sort -n`. (#1279)
- **`-a` counts unmapped records in `percent_mapped`.** Aligners keep them by default (bowtie2, and STAR with `--outSAMunmapped Within`), and 2.10.0 left them out, so a library where a quarter of the reads aligned reported 100%. Mapping rates from `-a` now drop to their true value, and the `low_mapping_rate` warning can fire. Quantification is unchanged. (#1278)
- **`-p` above the machine's thread count is capped for mapping too.** The warning said the cap applied, but only the EM used it. (#1284)
- **A rerun into an existing output directory first removes the old `quant.sf` and `quant.genes.sf`.** A rerun killed partway used to leave the old `quant.sf` beside the new metadata. A rerun that fails now leaves no `quant.sf`. `quant.sf` is also written under a temporary name and renamed, so a killed write leaves no truncated file. (#1290)
- **`--decoder parallel=N` is accepted**, as its help already said. `--help` no longer lists the decoder values. (#1293)
- **New warnings:** the alignment error-model flags under `--rad`, which has no CIGARs to score (#1299); a `-t` FASTA sequence whose length differs from the alignment's reference, apart from a clipped poly-A tail (#1300); and `--ignoreTxVersion` collisions (#1304).

## Alignment mode (`-a`)

- **A fragment whose two mates align, as orphans, to different transcripts keeps each mate's own side.** A RAD record stores one fragment type, and 2.10.0 applied it to every hit, so one of the two placements was judged on the other mate's strand. With `-l ISR`, 1000 such fragments went 0 / 1000 instead of 500 / 500. The same fix covers the reads path (see the top of these notes). (#1282)
- **A proper pair now drops the fragment's orphans**, as the reads mapper already did. An orphan's read length used to be read back as a fragment length, and could take as much weight as the pair. (#1282)
- **A single-end BAM under `-l A` is detected as single-end, and its strand is inferred.** 2.10.0 treated it as paired: `--discardOrphans` discarded every read, and a stranded library was quantified as unstranded. `--online` was already right. (#1281)
- **`--errorModel` scores a record without bases (`SEQ *`)** from another record of the same read. Such records are common for secondary alignments, minimap2's in particular. They scored 0, better than any real alignment, so reads moved to whichever transcript was reported second. (#1280)
- **`--hardFilter` keeps a fragment's best placement among those the strand filter keeps.** A better-scoring wrong-strand placement used to set the bar, and the fragment was lost whole. (#1285)
- **Outward library types (`-l OU`, `OSF`, `OSR`) map.** A RAD hit does not record which way a pair points, so every pair read as inward and an outward `-l` dropped them all. Under an outward type only the strand is now judged, so inward pairs are kept too. (#1292)
- **Only `/1` and `/2` are stripped from read names** when grouping mates. Reads named `r/a` and `r/b` used to merge into one fragment. (#1294)

## Speed

Full quantification, 2.10.0 → 2.11.0, wall time in seconds. Both builds are from source with the same profile, on one machine (AMD EPYC 9575F) pinned to one socket. Each figure is the median of 3 runs, or a single run at one thread. `quant.sf` is identical between the two versions in 31 of the 38 runs; the other 7 are muscle with selective alignment, which the orphan fix changes (see the top of these notes).

| workload | 1 thread | 8 threads | 16 threads |
|---|---|---|---|
| muscle, selective alignment | 107.2 → 106.3 | 14.83 → 14.64 | 8.44 → 8.25 |
| GEUVADIS, selective alignment | 119.3 → 117.9 | 17.09 → 17.04 | 9.52 → 9.54 |
| GEUVADIS + genome decoy, selective alignment | | 26.47 → 26.59 | 15.29 → 15.29 |
| GEUVADIS, uncompressed FASTQ, selective alignment | | 17.75 → 16.66 | 10.30 → 9.36 |
| muscle, `--sketch` | | 5.48 → 5.46 | 3.40 → 3.41 |
| GEUVADIS, `--sketch` | | 6.87 → 6.83 | 4.01 → 4.06 |

Sketch mode is unchanged in time and memory: without `--seqBias` or `--gcBias` it never loads the reference sequences.

Two small costs are gone from the mapping pass: the fragment counter shared by all threads is updated once per batch instead of once per fragment (#1308), and a single-end pass that does not quantify no longer builds bias models it then discards (#1297).

In alignment mode, unmapped records are now counted, but on the reader thread: they never reach the workers, so a BAM in which half the records are unmapped is read as fast as before.

## Other fixes

- **`--seqBias` and `--posBias` floor a short transcript's effective length as `--gcBias` does.** For a transcript shorter than every fragment length, the faster path used without `--gcBias` skipped the floor, so a 120 nt transcript got 75 there and 45 with `--gcBias`, a 1.67× difference in TPM. Both now give 45, as C++ salmon does. (#1289)
- **Single-end `-l A` decides the library type from every fragment.** The default mode took the first 50,000 that the threads happened to reach, so a file whose start was not representative could be filtered on the wrong strand. (#1295)
- **`--recoverOrphans` applies the pairing rules to a rescued mate.** A rescue could pair mates that normal pairing rejects as dovetailed, or report a fragment longer than the maximum. (#1296)
- **`--ignoreTxVersion` resolves two versions of a transcript in different genes the same way every run.** An exact name match wins; otherwise the lowest identifier does, and salmon warns. Before, the winner depended on hash order and changed from run to run. (#1304)
- **`--online` writes `ambig_info.tsv` row for row with `quant.sf`**, without decoy rows, and reports a failed final write. (#1298)
- **A `.zstd`, `.lzma` or `.bgzf` suffix is recognised** when deciding whether `-g` names a GTF and whether `-a` names a SAM. (#1303)
- **Reports:** a global-rule EM stopped at the iteration cap blames only the components still changing, not all of them (#1305), and the short EM that seeds the bias models no longer warns that it did not converge (#1306).
- **The EM sums a large component's step size the same way on one thread as on several.** The order differed, so the stopping test could land one iteration apart between `-p 1` and `-p N`. (#1307)

## Fixes

- **Sketch mode no longer reuses the previous fragment's hits for a mate shorter than k.**
  - A mate that short is not mapped, and the step it skipped was the one that cleared that mate's hits from the previous fragment. The pair was then merged with those stale hits.
  - Results depended on read order, and most such fragments were dropped as discordant or, sometimes, counted against the wrong transcripts. Now the short mate contributes nothing, and the fragment counts as an orphan of the other mate.
  - The bug affects adapter-trimmed data whose minimum length is below k (31 by default). On 500k muscle pairs with 3% of R1 mates cut to 20 nt, 2.10.0 mapped 463,057 fragments and changed 2,399 transcripts when the pairs were reversed; 2.11.0 maps 476,293 in either order, with identical results.
  - Selective alignment and single-end sketch were not affected. (#1271)
- **An orphan from read 2 is written to SAM/BAM on its own strand.** Every such record was flagged reverse (0x10), with its sequence reverse-complemented, whatever its orientation. Quantification was not affected. (#1272)
- **A failed `salmon index` no longer destroys the index it was replacing.** `refseq.bin` was emptied before the input was opened, so a mistyped FASTA path left an index that loaded and then crashed. The store is now built under a temporary name and renamed at the end, and `info.json`, written last, marks a complete index. Loading also checks the store against the index. (#1274)
- **Decoys no longer appear in `quant.genes.sf`.** With `-g` and a decoy-aware index, each decoy got a gene row, was listed as a transcript missing from the gene map, and lowered the reported match rate. (#1276)
- **Under `--sketch`, a right orphan's SAM/BAM record keeps its `AS` score**, rather than `AS:i:0`. (#1302)
- **`--online --sketch --posBias` no longer panics.** The effective-length correction read each transcript's sequence only for its length, and sketch mode does not load the sequences. It now takes the length directly. (#1265)

## Documentation

- **The inferential replicates guide explains why salmon offers both bootstrap and Gibbs replicates.**
  - Each bootstrap replicate re-runs the optimizer from the same starting point. Where a gene's isoforms are nearly identical, replicates can therefore understate how loosely the data constrain the split between them, and the Gibbs sampler shows that ambiguity more fully.
  - A new section covers what replicates cannot measure: the variance between technical replicate libraries.

  (#1266)
- The output-formats reference says that an `error` diagnostic makes `quant` exit non-zero, and `MIGRATION.md` notes it beside `--minAssignedFrags`. (#1277)

## For library users

The published `salmon-*` crates change in these places:
- `salmon_core::RefProvider::ref_seq` is now `ref_view`, which returns a `RefView` over bytes or packed bases.
- `SalmonIndex::load_with_opts` takes `RefseqLoad::{None, Bytes, Packed}` instead of a `bool`.
- `salmon_model::corrected_effective_length_full` takes the transcript length. The sequence now travels in `BiasInputs::seq`, with the sequence-bias tables.
- `salmon_rad::SalmonBulkRecord::from_placements` builds a record from its placements' statuses. A hit on the orphan side the record's type does not name sets its `mate_fw` bit, which orphan hits otherwise leave unset. Records with one orphan side are byte-identical to before.
- `salmon_align::AlignQuantResult::quant_rows` returns the rows `quant.sf` reports.
- `salmon_index::build` refuses k above 31, `m` not below k, unmatched decoy names and repeated reference names.
- `salmon_core::quantmerge::merge_quants` returns an error for the inputs `quantmerge` now refuses.
- Removed, with no callers: `salmon_model::corrected_effective_length` and `corrected_effective_length_fft` (use `corrected_effective_length_full`), and `salmon_infer::optimize_with_init` (use `optimize_packed_with_init`). (#1291)
- `salmon_core::write_via_rename` writes a file through a temporary name. (#1290)

## Dependencies

- piscem-rs 0.12.3 and paraseq-temp 0.6.3-pre.1. The paraseq changes are offered upstream as noamteyssier/paraseq#102.
- Patch updates to `thiserror`, `serde_json`, `xxhash-rust`, `sha2`, `libc`, `realfft`, `ensure_simd` and `jiff`. Indices and `quant.sf` are byte-identical with them. (#1268)
- The Docker image builds with Rust 1.99. (#1262)
- Website: Astro 7.4, Starlight 0.42.8, smol-toml 1.9.1. (#1267, #1248)

## Checked

- `cargo test --release --workspace --locked`: 435 passed. Clippy (`-D warnings`) and fmt are clean.
- `quant.sf` is identical to 2.10.0 on every benchmark run above, with and without the genome decoy, in selective-alignment and sketch mode, except muscle with selective alignment (the orphan fix). None of those samples has a mate shorter than 31 nt, so the sketch fix does not apply to them.
- The fixes from the second part of the audit (#1289–#1308) leave `quant.sf`, the mapping counts and run time unchanged at 16 threads on muscle and GEUVADIS (selective alignment and sketch), the genome decoy, GEUVADIS with `--seqBias --posBias`, and GEUVADIS read 1 as single-end under `-l A`; and at one thread on muscle and GEUVADIS with selective alignment.
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
- Alignment mode: on SAM and BAM inputs with 2M mapped pairs, alone or with 2M unmapped pairs, `quant.sf` and run time match the development branch without these fixes, with and without `--online`. `percent_mapped` gives 50% where half the reads are unmapped, rather than 100%.
- GENCODE v50 (670,670 records, with and without `--gencode`) and Ensembl GRCh38 cDNA have no repeated names, so the new name check does not affect them.

## Acknowledgements

Thanks to Benjamin Demaille for:
- the two-bit reference packing (#1265);
- stream input (#1269);
- the inferential replicates documentation (#1266);
- the sketch and SAM fixes (#1271, #1272);
- the audit of 2.10.0 (#1283) and its fixes (#1274–#1282, #1284–#1300, #1302–#1308).
