# CLI flag status — C++ salmon vs. the Rust port

This table tracks every command-line flag of C++ salmon (`salmon index`,
`salmon quant` in reads and alignment modes, `salmon quantmerge`) and how the
Rust port handles it. It is the authoritative reference for drop-in
compatibility.

> **salmon 2.0 (Rust):** a user-facing summary of removed/ignored/new options
> lives in [`MIGRATION.md`](../MIGRATION.md). 2.0 adds a `salmon alevin` stub
> that redirects to alevin-fry, and any unrecognized (⛔) flag now errors with a
> pointer to that migration guide rather than a bare clap message. New 2.0-only
> flags: `--sketch`, `--sketchStrictOrphans`, `--allowDovetail` (in sketch).

## Legend

| Status | Meaning |
|--------|---------|
| ✅ **Implemented** | Functionally equivalent to C++ salmon. |
| 🟰 **Default** | The behavior it selects is now always on; the flag is accepted but redundant. |
| ⚠️ **Accept + warn** | Parsed for command-line compatibility, logs a warning, has no effect (no Rust analog, or genuinely a no-op). |
| ⛔ **Rejected** | Not implemented; passing it is an "unknown argument" error. |
| 🔁 **Covered** | Functionality is reachable via a different flag / is the default. |

Anything not ✅/🟰/🔁 does **not** change results; ⚠️ flags let existing salmon
command lines run unchanged, while ⛔ flags fail loudly so nothing silently
misbehaves.

## Global (top-level) flags

| Flag | Status | Notes |
|------|--------|-------|
| `-q/--quiet`, `-h/--help`, `-V/--version` | ✅ | |
| `--no-version-check` | 🔁 | Accepted before or after the command (C++ puts it before: `salmon --no-version-check <cmd>`); silent no-op — 2.0 never checks for updates. `SALMON_NO_VERSION_CHECK` env is likewise irrelevant. |

## `salmon index`

| Flag | Status | Notes |
|------|--------|-------|
| `-t/--transcripts`, `-i/--index`, `-k/--kmerLen`, `-p/--threads` | ✅ | |
| `--keepDuplicates` | ✅ | |
| `--keepFixedFasta` | ✅ | Keeps the post-N-cleaning FASTA. |
| `--gencode` | ✅ | Truncates names at the first `|`. |
| `-d/--decoys` | ✅ | |
| `--tmpdir` | ✅ | Redirects build intermediates; final index still in `-i`. |
| `-m/--minimizerLen` | ✅ | `--keepIntermediate` also supported (keeps cDBG files). |
| `-n/--no-clip` | ✅ | Default clips poly-A tails (≥10 trailing `A`s → trimmed; all-`A` dropped), matching pufferfish `FixFasta`; `--no-clip` disables. Hashes are pre-clip. |
| `--filterSize` | ⚠️ | cf1-rs cDBG builder does not expose the bloom-filter size. |
| `--features` | ⛔ | |

## `salmon quant` — common / inference

| Flag | Status | Notes |
|------|--------|-------|
| `-l/--libType`, `-o/--output`, `-p/--threads`, `-z/--writeMappings`, `--writeBam` | ✅ | `--writeSam` is a visible alias for `--writeMappings`. `--writeMappings`/`--writeSam` and `--writeBam` are mutually exclusive, and are warned-and-ignored in `-a`/`--rad` modes. Records for a fragment are contiguous; no other order is imposed. Records carry real base-level CIGARs (`I`/`D`/`S`) with `NM` and `MD`, base qualities from the input, `MAPQ` derived from the placement count on STAR's scale (255 unique, 3, 1, 0), and `NH`/`HI`/`XT`/`AS`/`ZW`/`MC`/`MQ` tags. The header declares `VN:1.6 SO:unsorted GO:query`, `@SQ` with `M5` and `UR`, and `@PG`. |
| `--spliced` | ✅ | Report mapping output in genome coordinates with `N` at exon junctions, from `--annotation` (exon structures) and `--genome` (chromosome lengths for `@SQ`). Adds `XS` for the transcript strand. Rust-port feature; not in this salmon build. |
| `--rgLine` | ✅ | `@RG` header line plus `RG:Z` on every record; bwa/STAR spelling, escaped tabs accepted. Rust-port feature; not in this salmon build. |
| `--bamCompressThreads` | ✅ | BGZF compression pool for `--writeBam`; defaults to ~1 worker per 3 mapping threads (max 8), balancing measured deflate throughput against record production. Rust-port feature; not in this salmon build. |
| `--useEM` | ✅ | VBEM is the default; `--useEM` selects plain EM. |
| `--meta` | ✅ | Metagenomic preset: plain EM, no range-factorized eq-classes, uniform init (the Rust EM already inits uniformly). Overrides `--useEM`/`--rangeFactorizationBins`. |
| `--useVBOpt` | 🔁 | VBEM is already the default. |
| `--numBootstraps`, `--numGibbsSamples`, `--thinningFactor` | ✅ | |
| `--vbPrior`, `--perTranscriptPrior`, `--perNucleotidePrior` | ✅ | |
| `--rangeFactorizationBins`, `--incompatPrior`, `--sigDigits` | ✅ | |
| `--dumpEq`, `--dumpEqWeights`, `--writeUnmappedNames` | ✅ | The eq dumps come from the pass that builds the classes, so they work in every RAD-based mode (including the 2.6.0 default and `--rad`) as well as `--online`. `--writeUnmappedNames` is reads-mode only: it warns and is ignored under `-a`/`--rad`, where salmon never sees the unmapped reads. |
| `--initUniform` | ✅ | |
| `--skipQuant` | ✅ | Skips EM + Gibbs/bootstrap + quant.sf; still emits eq-classes/metadata. |
| `-g/--geneMap`, `--genes` | ✅ | Gene-level aggregation. Rust also adds `--ignoreTxVersion` (no C++ equivalent) for Ensembl cDNA + GTF identifier mismatches. |
| `-q/--quiet` | ✅ | |
| `--noLengthCorrection` | ✅ | (`--noEffectiveLengthCorrection`). |
| `--minAssignedFrags` | ⛔ | |
| `--alternativeInitMode`, `--bootstrapReproject`, `--noGammaDraw` | ⛔ | Inference niche. |
| `--eqclasses` | ⚠️ | Quantify-from-eq-class-file *input* mode (≠ `--dumpEq`); not yet implemented. |

## `salmon quant` — fragment-length / bias models

| Flag | Status | Notes |
|------|--------|-------|
| `--seqBias`, `--gcBias`, `--posBias` | ✅ | |
| `--fldMean`, `--fldSD`, `--fldMax` | ✅ | Priors. With `--rad` they are ignored unless `--fldPolicy` says otherwise (a salmon RAD bakes its FLD); salmon warns when an explicitly-supplied value is ignored. |
| `--fldPolicy` | ✅ | `baked` (default) / `derive` / `prior`: where a RAD requant's fragment-length distribution comes from. Rust-port feature; not in this salmon build. |
| `-f/--forgettingFactor` | ✅ | Online-path only: since 2.6.0 the deterministic default runs no online phase, so it warns and is ignored without `--online`. |
| `--numAuxModelSamples` | ✅ | Online-phase model-training window; same `--online`-only caveat as `--forgettingFactor`. |
| `--biasSpeedSamp` | ✅ | GC convolution fragment-length stride. |
| `--numGCBins`, `--conditionalGCBins` | ✅ | |
| `--noBiasLengthThreshold` | ✅ | |
| `--reduceGCMemory` | 🟰 | The rank-bitvector GC representation it selects is now the default (faster + ~2× leaner, identical results). |
| `--numBiasSamples` | ⛔ | No separate seq-bias sample cap (collected within the aux window). |
| `--numPreAuxModelSamples` | ✅ | Pre-aux burn-in window (port default 5000; salmon 1,000,000); same `--online`-only caveat as `--forgettingFactor`. |
| `--noFragLengthDist` | ✅ | Drops the fragment-length term from a placement's probability (the proper-pair PMF term and the orphan ambiguous-length term alike). Honoured on every quantification path — reads, `-a`, `--rad`, deterministic and online — except the deprecated `--online -a`, which warns that it lacks it. |
| `--noSingleFragProb` | ✅ | Replaces an orphan's bounded-CMF fragment-length estimate with salmon's flat `LOG_EPSILON` penalty. Honoured on the default deterministic / RAD paths as well as `--online`; it changes point estimates wherever orphans are contested. |

## `salmon quant` — mapping / selective alignment (reads mode)

| Flag | Status | Notes |
|------|--------|-------|
| `--sketch` | ✅ | Pseudoalignment-only mapping (piscem-rs permissive sketch). Rust-port feature; not in this salmon build. |
| `--uniMEMs`, `--refMEMs` | ✅ | |
| `--ma`, `--mp`, `--go`, `--ge` | ✅ | Alignment scoring. |
| `--minScoreFraction`, `--bandwidth`, `--fullLengthAlignment` | ✅ | |
| `--consensusSlack`, `--maxOccsPerHit`, `--maxReadOcc` | ✅ | |
| `--allowDovetail`, `--recoverOrphans`, `--discardOrphansQuasi` | ✅ | |
| `--decoyThreshold`, `--minAlnProb`, `--hardFilter`, `--scoreExp` | ✅ | |
| `--preMergeChainSubThresh` | ✅ | (Rust default 0.8; salmon 0.75.) |
| `--postMergeChainSubThresh`, `--orphanChainSubThresh` | ✅ | Coverage-based concordant-pair / orphan prune. |
| `--softclip`, `--softclipOverhangs` | ✅ | |
| `--mismatchSeedSkip` | ⚠️ | Post-miss seed skip lives in the piscem-rs seed walk; not yet tunable. |
| `--disableChainingHeuristic` | ⚠️ | Rust uses a loss-less ref-distance early-break, not salmon's numRounds heuristic. |
| `--hitFilterPolicy` | ⚠️ | Rust filters hits after chaining (salmon's AFTER default); accepted, not yet tunable. |
| `--maxRecoverReadOcc` | ⚠️ | Orphan-recovery candidate cap; accepted, not yet implemented. |
| `--validateMappings` | ⚠️ | Deprecated in salmon too (no effect); selective alignment is the default. |
| `--noSA` | 🔁 | Selective alignment is always on. |
| `--mimicBT2`, `--mimicStrictBT2` | ⛔ | Bowtie2-mimicking presets. |
| `--puff`, `--type` | 🔁 | Single index type; obsolete. |

## `salmon quant -a` — alignment mode

| Flag | Status | Notes |
|------|--------|-------|
| `-a/--alignments`, `-t/--targets` | ✅ | BAM/SAM input. |
| `--noErrorModel`, `--numErrorBins` | ✅ | |
| `--discardOrphans` | ✅ | Alignment-mode orphan discard; also honoured when requantifying a RAD (`--rad`), so a requant agrees with the run that wrote the file. |
| `--hardFilter`, `--scoreExp` | ✅ | Also honoured on `-a`/`--rad`, not just in reads mode: the requant applies them to the RAD's `AS`-derived scores exactly as it does to selective-alignment ones. |
| `--ont` | ✅ | Redirect (long-read mode is out of scope, as upstream). |
| `--mappingCacheMemoryLimit` | ⚠️ | Rust streams with bounded buffers; no mass-banking cache. |
| `--sampleUnaligned`/`-u` | ✅ | Adds a `FLAG 0x4` record per unmapped fragment to `--writeMappings`/`--writeBam`, sequence and qualities intact, so the output covers the whole library. No longer requires `--sampleOut`. |
| `--writeQualities` | 🔁 | Base qualities are not optional output: every record carries them whenever the input does, so the flag is accepted and has nothing to switch on. |
| `--sampleOut`/`-s` | ⚠️ | A BAM of posterior-sampled alignments; accepted, not yet implemented. |
| `--auxTargetFile`, `--writeOrphanLinks` | ⛔ | |

## `salmon quantmerge`

| Flag | Status | Notes |
|------|--------|-------|
| `--quants`, `--names`, `--column`, `--missing`, `--genes`, `-o/--output` | ✅ | Complete. |

## Hidden / testing flags (C++ `getHiddenOptions` / `getTestingOptions`)

All ⛔ (rejected) — none affect results, and they target C++-specific internals:
`--readBatchSize`, `--adaptiveReadBatch`, `--maxHashResizeThreads`,
`--progressUpdateMs`, `--disableLiveProgress`, `--emitJoinDedupStats`,
`--noRichEqClasses`, `--noFragLenFactor`, `--disableAlignmentCache`,
`--rankEqClasses`, `--noExtrapolateCounts`, `--auxDir`.

## Out of scope (matching upstream)

`salmon alevin` (single-cell) is a stub that prints the alevin-fry redirect, as
in recent C++ salmon.
