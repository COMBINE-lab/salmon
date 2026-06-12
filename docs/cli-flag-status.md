# CLI flag status — C++ salmon vs. the Rust port

This table tracks every command-line flag of C++ salmon (`salmon index`,
`salmon quant` in reads and alignment modes, `salmon quantmerge`) and how the
Rust port handles it. It is the authoritative reference for drop-in
compatibility.

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
| `--filterSize` | ⚠️ | cf1-rs cDBG builder does not expose the bloom-filter size. |
| `--features` | ⛔ | |

## `salmon quant` — common / inference

| Flag | Status | Notes |
|------|--------|-------|
| `-l/--libType`, `-o/--output`, `-p/--threads`, `-z/--writeMappings` | ✅ | |
| `--useEM` | ✅ | VBEM is the default; `--useEM` selects plain EM. |
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
| `--numPreAuxModelSamples` | ⛔ | No separate pre-aux burn-in window. |
| `--noFragLengthDist`, `--noSingleFragProb` | ⚠️ | FLD-in-fragment-probability toggles (experimental upstream); accepted, not yet implemented. |

## `salmon quant` — mapping / selective alignment (reads mode)

| Flag | Status | Notes |
|------|--------|-------|
| `--sketch` (alignment-free), `--uniMEMs`, `--refMEMs` | ✅ | |
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
| `--validateMappings`, `--noSA` | 🔁 | Selective alignment is always on. |
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
