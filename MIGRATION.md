# Migrating from C++ salmon (≤ 1.12.0) to salmon 2.0 (Rust)

salmon 2.0 is a from-scratch Rust rewrite. It keeps the same core workflow
(`salmon index` → `salmon quant` → `quant.sf`) and the same output formats that
downstream tools consume, but it is a new major version and makes some breaking
changes. This guide maps C++ options/behavior to 2.0.

## Breaking change: rebuild your index

**2.0 uses a new index format and cannot read C++ (pufferfish) indices.** Rebuild
with `salmon index` from 2.0; loading a C++ index (or pointing C++ salmon at a
2.0 index) is detected and rejected with a clear error.

`quant.sf` is unchanged (drop-in for tximport/tximeta). Inferential replicates
(`aux_info/bootstrap/{names.tsv.gz,bootstraps.gz}` for both `--numBootstraps` and
`--numGibbsSamples`) are written in the same format C++ salmon used, so
tximport/fishpond/swish keep working. The bias-model diagnostic dumps in
`aux_info/` (`obs/exp*_seq.gz`, `obs/exp_gc.gz`, `obs/exp*_pos.gz`, `fld.gz`) are
in a documented Rust format (see the docs site); they are not consumed by the
standard downstream R packages.

## Removed subcommands

| C++ | 2.0 |
| --- | --- |
| `salmon alevin` | Removed. Use the [alevin-fry](https://github.com/COMBINE-lab/alevin-fry) ecosystem for single-cell. `salmon alevin …` prints this redirect and exits. |

## Removed options (produce an error + this note)

Passing any of these errors out with a pointer to this guide. They are gone
because the underlying feature was removed or has no equivalent in 2.0.

| Option | Subcommand | Why / alternative |
| --- | --- | --- |
| `--features` | index | Index-feature dump not supported by the cf1-rs/piscem-rs index builder. |
| `--mimicBT2`, `--mimicStrictBT2` | quant | bowtie2-mimicking presets removed; 2.0's selective alignment is the single supported mode. |
| `--minAssignedFrags` | quant | The "zero the output below N assigned fragments" guard was removed; 2.0 reports what it quantifies. |
| `--alternativeInitMode`, `--bootstrapReproject`, `--noGammaDraw` | quant | Inference-internal toggles not present in the 2.0 optimizer/Gibbs implementation. |
| `--numBiasSamples` | quant | 2.0 collects bias samples online (abundance-aware dual-phase), so a fixed sample budget no longer applies. |
| `--auxTargetFile`, `--writeOrphanLinks` | quant -a | Removed alignment-mode features. |

## Accepted but ignored (parse + warn)

These still parse so existing scripts run; 2.0 logs a warning and ignores them
(the behavior is either the default now or handled differently).

- global: `--no-version-check` (and the `SALMON_NO_VERSION_CHECK` environment
  variable) — accepted **before or after** the command (e.g.
  `salmon --no-version-check quant …`), matching C++. It is a silent no-op: 2.0
  never contacts the network to check for a newer release.
- index: `--filterSize`
- quant: `--eqclasses`, `--noFragLengthDist`, `--noSingleFragProb`,
  `--mismatchSeedSkip`, `--disableChainingHeuristic`, `--hitFilterPolicy`,
  `--maxRecoverReadOcc`, `--validateMappings` (selective alignment is the default)
- quant -a: `--mappingCacheMemoryLimit`, `-s/--sampleOut`, `-u/--sampleUnaligned`,
  `--writeQualities`

## New in 2.0

- `--sketch` — alignment-free pseudoalignment mode (faster; quantifies directly
  from k-mer/equivalence-class hits).
- `--sketchStrictOrphans` — in `--sketch`, only orphan a pair when the other mate
  had no matching k-mers (the conservative rule). Default is the relaxed rule
  (orphan when the other mate has no consistent target), which tracks selective
  alignment more closely.
- `--allowDovetail` — now honored in `--sketch` as well (admits dovetailed
  short-insert fragments).
- `--ignoreTxVersion` — compare `-g/--geneMap` identifiers without their trailing
  `.N` version suffix, as tximport's option of the same name does. Off by
  default. Needed for an Ensembl cDNA index against an Ensembl GTF, where the
  FASTA carries the version in the identifier and the GTF puts it in a separate
  `transcript_version` attribute.
- `--degCoverage` / `--degCovBins` and the `salmon degnorm` subcommand —
  degradation normalization in the spirit of
  [DegNorm](https://nustatbioinfo.github.io/DegNorm/). `quant --degCoverage`
  records each transcript's coverage curve while the reads are streaming;
  `salmon degnorm` fits them across a cohort and writes a degradation index per
  transcript per sample, plus one corrected `quant.sf` per sample that tximport
  and DESeq2 read unchanged. Reads mode only. See
  [`docs/degnorm.md`](docs/degnorm.md).

## Behavior differences to be aware of

- **Sketch orphan rule** defaults to the relaxed policy (see `--sketchStrictOrphans`).
- **`-g/--geneMap`: unmatched transcripts are reported differently, not
  aggregated differently.** As in C++ salmon's `aggregateEstimatesToGeneLevel`, a
  transcript with no gene-map entry is emitted as its own single-transcript gene,
  named after the transcript, so nothing quantified is dropped on the way to gene
  level and `tximport` can apply its own transcript-to-gene policy to the whole
  file. What changed is the reporting: C++ warned once *per unmatched transcript*
  (half a million lines on a failed join), while 2.x warns once with the count
  and writes the names to `aux_info/genemap_unmatched_txps.json`, under an
  `unmatched_transcripts` key.
  That file is removed when nothing is unmatched, so its presence always refers
  to the current run. 2.x additionally warns when the match rate is at or below
  50%, says how many transcripts would match without the version suffix, and
  names `--ignoreTxVersion`. Note that a `quant.genes.sf` whose rows are mostly
  stand-ins is transcript-level output under a gene-level file name — the run's
  closing line says so explicitly when any row is a stand-in.
- **Selective-alignment chain pruning:** 2.0 currently defaults
  `--orphanChainSubThresh` and `--postMergeChainSubThresh` to `0.0` (off) — it
  aligns every candidate, which is marginally more sensitive than C++ (which uses
  `0.95`/`0.9`). Quantification is essentially unaffected (per-transcript Pearson
  ≈ 0.999); pass `--orphanChainSubThresh 0.95 --postMergeChainSubThresh 0.9` to
  reproduce C++ mapping counts exactly. See `docs/mapping-parity-differences.md`.

## Unchanged

Index/quant basics, `quant.sf`, `cmd_info.json`, `lib_format_counts.json`,
`aux_info/meta_info.json`, `--libType/-l`, `--threads/-p`, `-k/--kmerLen`,
`-m/--minimizerLen`, `-n/--no-clip` (poly-A clipping, on by default),
`--seqBias`, `--gcBias`, `--posBias`, `--numBootstraps`, `--numGibbsSamples`,
`--useEM`, `--meta` (metagenomic preset), `--dumpEq`, decoys, and
`salmon quantmerge`.
