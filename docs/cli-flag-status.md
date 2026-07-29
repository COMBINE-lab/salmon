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
| `-l/--libType`, `-o/--output`, `-p/--threads`, `-z/--writeMappings`, `--writeBam` | ✅ | `--writeSam` is a visible alias for `--writeMappings`. `--writeMappings`/`--writeSam` and `--writeBam` are mutually exclusive, and are warned-and-ignored in `-a`/`--rad` modes. Records for a fragment are contiguous; no other order is imposed. |
| `--useEM` | ✅ | VBEM is the default; `--useEM` selects plain EM. |
| `--meta` | ✅ | Metagenomic preset: plain EM, no range-factorized eq-classes, uniform init (the Rust EM already inits uniformly). Overrides `--useEM`/`--rangeFactorizationBins`. |
| `--useVBOpt` | 🔁 | VBEM is already the default. |
| `--numBootstraps`, `--numGibbsSamples`, `--thinningFactor` | ✅ | |
| `--vbPrior`, `--perTranscriptPrior`, `--perNucleotidePrior` | ✅ | |
| `--rangeFactorizationBins`, `--incompatPrior`, `--sigDigits` | ✅ | |
| `--dumpEq`, `--dumpEqWeights`, `--writeUnmappedNames` | ✅ | |
| `--initUniform` | ✅ | |
| `--skipQuant` | ✅ | Skips EM + Gibbs/bootstrap + quant.sf; still emits eq-classes/metadata. |
| `-g/--geneMap`, `--genes` | ✅ | Gene-level aggregation. |
| `-q/--quiet` | ✅ | |
| `--noLengthCorrection` | ✅ | (`--noEffectiveLengthCorrection`). |
| `--minAssignedFrags` | ⛔ | |
| `--alternativeInitMode`, `--bootstrapReproject`, `--noGammaDraw` | ⛔ | Inference niche. |
| `--eqclasses` | ⚠️ | Quantify-from-eq-class-file *input* mode (≠ `--dumpEq`); not yet implemented. |

## `salmon quant` — fragment-length / bias models

| Flag | Status | Notes |
|------|--------|-------|
| `--seqBias`, `--gcBias`, `--posBias` | ✅ | |
| `--fldMean`, `--fldSD`, `--fldMax` | ✅ | |
| `-f/--forgettingFactor` | ✅ | |
| `--numAuxModelSamples` | ✅ | Online-phase model-training window. |
| `--biasSpeedSamp` | ✅ | GC convolution fragment-length stride. |
| `--numGCBins`, `--conditionalGCBins` | ✅ | |
| `--noBiasLengthThreshold` | ✅ | |
| `--reduceGCMemory` | 🟰 | The rank-bitvector GC representation it selects is now the default (faster + ~2× leaner, identical results). |
| `--numBiasSamples` | ⛔ | No separate seq-bias sample cap (collected within the aux window). |
| `--numPreAuxModelSamples` | ✅ | Pre-aux burn-in window (port default 5000; salmon 1,000,000). |
| `--noFragLengthDist`, `--noSingleFragProb` | ⚠️ | FLD-in-fragment-probability toggles (experimental upstream); accepted, not yet implemented. |

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
| `--discardOrphans` | ✅ | Alignment-mode orphan discard. |
| `--ont` | ✅ | Redirect (long-read mode is out of scope, as upstream). |
| `--mappingCacheMemoryLimit` | ⚠️ | Rust streams with bounded buffers; no mass-banking cache. |
| `--sampleOut`/`-s`, `--sampleUnaligned`/`-u`, `--writeQualities` | ⚠️ | Posterior-sampled BAM output; accepted, not yet implemented. |
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
