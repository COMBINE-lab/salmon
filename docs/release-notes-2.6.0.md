# salmon 2.6.0

The release that makes **deterministic quantification the default**. A plain
`salmon quant` now maps once to an intermediate RAD and quantifies from it,
producing `quant.sf` files that are **byte-identical across runs and thread
counts** — in every input mode: selective alignment, sketch, and alignment
(`-a`). **No index rebuild is required.**

The decision is evidence-backed rather than aspirational (#1140):

- **Speed**: hardware-dependent, and the honest summary is that deterministic
  wins decisively where it matters most and costs a little elsewhere. On a
  modern many-core machine (EPYC, NVMe) it is as fast or faster in every
  measured cell: on 10M real read pairs, sketch is 13% faster wall / 20% less
  CPU at `-p 16` and 11% faster at `-p 4`; selective alignment ties at `-p 16`
  and is ~7% faster at `-p 4` — *including* the cost of writing and re-reading
  the intermediate.

  On an older CPU (Xeon E5-2699 v4, 26M real fragments) the picture is mixed:
  deterministic is ~4% slower at `-p 16` and ~12% slower at `-p 32`, then
  **twice as fast at `-p 64`**, where the online path's lock contention takes
  over. It is also far more *predictable* there — across four repeats at
  `-p 64` the online path ranged 55–105 s (bimodal), while deterministic held
  39–42 s. When sizing a pipeline, that stability is usually worth more than
  the mean.
- **Accuracy**: statistically indistinguishable from the online path against
  ground truth (fourth-decimal differences, direction inconsistent, across two
  simulated samples and two multimapping regimes).
- **Reproducibility**: `quant.sf` byte-identical at every thread count, on real
  data (26M fragments, 238k transcripts, verified at full `--sigDigits 9`
  precision); the online path produces a different file at every thread count.
  Enforced by a test that compares bit-for-bit across rayon pools of 1–16
  threads, over plain EM, VBEM and SQUAREM.

## The deprecation path

- `--deterministic` is **accepted as a no-op** (it is in scripts and CI configs
  already) and will remain so through at least 2.7.
- `--online` selects the pre-2.6 one-pass path,
  **deprecated** with removal planned for 2.7.0. It announces its deprecation
  when used; `--online --deterministic` is an error rather than a guess.
- Online-only knobs — `--forgettingFactor`, `--numAuxModelSamples`,
  `--numPreAuxModelSamples` — **warn and are ignored** when passed without
  `--online`, instead of being silently inert. Under `--online -a`,
  `--noLengthCorrection`/`--noFragLengthDist` likewise warn (the online
  alignment path never supported them; the default path honours them).
- `aux_info/meta_info.json` records which path ran in a new **`inference_path`**
  field (`deterministic` / `online` / `none` for a mapping-only pass), so the
  transition is auditable from existing output.

## What changes in behavior (beyond the flip itself)

- **`--skipQuant`** now stops after the mapping pass and keeps the RAD as the
  run's output (previously, combined with `--deterministic`, it was silently
  ignored and a full `quant.sf` was written). With `--dumpEq` it warns that no
  classes exist to dump and points at `--rad --dumpEq` on the kept RAD.
- **`-l A` accounting**: the deterministic path detects the library type from
  the *whole run* and filters afterward, so `num_compatible_fragments` covers
  every judged fragment; the online path never counted the ~50k-fragment
  detection prefix. `compatible_fragment_ratio` is unchanged in meaning; the
  counts it is built from are more complete. Threshold-based QC on these fields
  should expect slightly higher absolute counts.
- **Alignment mode scores by the BAM `AS` tag by default** (measured at least
  as accurate as the error model at half the wall time); `--errorModel -t`
  opts into the order-independent counting error model. A BAM whose `AS` tags
  are missing, partial, or constant now draws a warning naming the way out —
  recorded in `meta_info.json`'s `diagnostics` as `bam_score_tags_unusable`,
  not just the log. The *online* error model exists only behind `--online -a`
  and leaves with it.

## Disk behavior of the intermediate

The intermediate RAD is lz4-compressed and deleted on success. Measured on
GENCODE v49 with 26.1M real fragments: **~45 MB per 1M fragments in sketch mode
and ~55 MB in selective alignment** (it grows with the multimapping rate).

**Disk speed is not a factor.** Validated on rotational storage (measured
105 MB/s): writing the RAD is entirely hidden behind mapping compute — the
mapping phase takes 30.55 s with the RAD on a spinning disk versus 30.92 s with
it in `/dev/shm` — and reading it back is bound by decompression rather than
I/O. (Uncompressed, the RAD is 5× larger and still reads 2.7× *faster* than
zstd.) Selective alignment is completely insensitive to disk speed. Fast
scratch storage is not needed for the deterministic default.

New in this release cycle:

- **`--radScratchDir <dir>`** places it on node-local or memory-backed storage
  (e.g. `/dev/shm`), pid-suffixed so concurrent runs can share a scratch
  volume. A `--skipQuant` deliverable always goes to the output directory.
- A **startup preflight** warns when the target volume has less free space than
  the total compressed input, and the writer **checks free space during the
  pass**, so a filling disk fails in seconds with actionable advice
  (`--radScratchDir`, `--radCompress zstd` for ~53% more shrink at ~5.5% wall)
  instead of at the end of mapping. On failure, salmon-owned intermediates
  (including `.partial-*` files) are cleaned up; explicit `--writeRad` paths
  are never touched.

## Accounting and correctness fixes in this cycle

- **`lib_format_counts.json` restored to C++ semantics in every mode** (#1131,
  #1136/#1137, #1138): measured `num_compatible_fragments` /
  `num_incompatible_fragments` with the ratio over judged fragments; the
  per-observed-format histogram (`ISF`, `ISR`, …); format-agreement-based
  concordant/inconsistent counts; a measured `strand_mapping_bias` with the
  classic warnings. Sketch mode now applies the strand-compatibility filter
  (#1136) and `-l A` works in every path, including online alignment mode
  (which detects from the reported alignments — with a warning that an
  upstream orientation filter cannot be recovered from).
- **Flags that were accepted and silently dropped under `--deterministic` now
  work**: `--biasSpeedSamp`, `--noBiasLengthThreshold`, `--initUniform`
  (forwarded; note that every reads-mode path already starts the optimizer
  uniformly, so it toggles behavior only in `--online -a`),
  `--noLengthCorrection` and `--noFragLengthDist` (newly implemented in the
  RAD path — tagged-end protocols in particular), and `--dumpEq` /
  `--dumpEqWeights` (previously wrote a well-formed file declaring **zero**
  classes; now dumps the classes the EM consumed, in every RAD-based mode).
- **`total_time_seconds` in `meta_info.json` spans the whole run** (it
  previously covered only the requant pass — reporting ~0.4s for a run whose
  mapping took minutes).
- The deterministic mapping pass has a **live progress spinner**, like the
  one-pass path.

## The pre-release behavior audit (#1140)

Before tagging, every quant flag was traced per path and the tree swept for
computed-but-unused bookkeeping, docs contradicting code, and paths that
diverged without a reason. 111 findings survived verification and were
addressed across #1153–#1158 — most by fixing the code, some by correcting
documentation that described behavior the code never had, and one by
documenting a divergence deliberately left alone (the deprecated `--online -a`
fragment-length model, which is leaving in 2.7.0 and gets no parity work).
Most of it is invisible from outside; these are the parts that are not.

### Output that changes shape

- **Decoy references no longer appear in `quant.sf`** on any RAD-based path —
  which, after the flip, is the default. Phase 1 writes the full reference
  table (decoys included) into the RAD, and the writer had no decoy boundary to
  filter on, so a decoy-aware index produced decoy rows that `--online` never
  emitted. The boundary is now baked into the RAD header. The same fix applies
  to `aux_info/bootstrap/names.tsv.gz`, the bootstrap matrix, `ambig_info.tsv`,
  and `num_valid_targets`. **If you parse `quant.sf` by row index, the row set
  changes** — back to what the docs always promised.
- **`meta_info.json` stops reporting placeholders**: `samp_type` names the
  posterior that actually ran (tximport keys off it together with
  `num_bootstraps`), `serialized_eq_classes` reflects whether classes were
  dumped, `num_decoy_fragments` reads the counter baked in the RAD instead of a
  hardcoded `0`, `num_decoy_targets` reports the real decoy count instead of a
  hardcoded `0` (it contradicted `num_decoy_fragments` in the same file), and
  `frag_length_source` distinguishes an observed distribution from the
  `--fldMean`/`--fldSD` prior a single-end run keeps — and no longer answers
  `rad_baked` for a plain `-a` run, which leaked the internal intermediate RAD
  into the run's own record. Alignment input reports `alignments`; the
  `rad_baked*` spellings are for a RAD you supplied yourself.
- **`library_types` reports the library format that was *applied*, not the one
  the detector observed.** With an explicit `-l` it reported the detected
  format, contradicting `cmd_info.json` and `lib_format_counts.json` in the same
  output directory: `-l ISR` on ISF-looking data reported `["ISF"]`, while the
  run really did filter as ISR and map nothing. The inferred format keeps its
  own field, so metadata now reports **both** — what was applied
  (`library_types`) and what was observed (`detected_library_type`). **If you
  parse `library_types` to learn what salmon detected, read
  `detected_library_type` instead.**
- **New: the bias models' shapes are recorded.** `meta_info.json` gains
  `seq_bias_context_length`, `gc_bias_bins` and `pos_bias_bins`, each present
  only when that correction ran. `num_bias_bins` keeps its C++ meaning — the
  size of a legacy k-mer histogram this port replaced — and so stays `0`;
  nothing previously recorded what the sequence, GC and positional models
  actually looked like.

- **`logs/salmon_quant.log` describes the run you asked for.** The log is
  written by the pass that quantifies, which the default reads path now shares
  with `-a`, so a plain reads run was headed "alignment mode" and reported its
  unmapped fragments as strand-incompatible. Reads runs are labelled and
  counted as reads runs again.
- **`cmd_info.json` records your invocation**, not phase 2's view of it. The
  requant used to overwrite it with the alignment-mode schema, so `index`,
  `mates1`/`mates2`, `threads` and `sketch` vanished from the run's own record
  and `lib_format_counts.json` reported `read_files: []`.
- **`MAPQ` follows STAR's placement-count scale** (255 unique, then 3 / 1 / 0)
  instead of a constant `1`, so `samtools view -q 255` means "uniquely placed"
  as downstream RNA-seq tooling expects. Previously that filter discarded the
  entire file.
- **Right-orphan records carry their real alignment score.** In
  `--writeMappings`/`--writeBam` output, every record for an orphan whose mapped
  mate is read 2 reported `AS:i:0` for an alignment that actually scored ~180–300;
  it now reports the true score, in both SAM and BAM. Left orphans and proper
  pairs were already correct. This is output-only — `quant.sf` is unaffected.
  (In sketch mode `AS` remains `0` for all records by design: pseudoalignment
  computes no per-placement score, so every mapping is a co-equal best mapping.)

### Estimates that change

Each of these was a flag that was accepted, documented, and silently did
nothing on the default path while working under `--online`. Runs that passed
them will now get different numbers — the intended ones.

- **`--perNucleotidePrior`** was ignored by Gibbs sampling on every path (all
  three call sites hardcoded a per-transcript prior) and inside bootstrap
  replicates (effective lengths never reached the replicate EM, which silently
  falls back to the flat prior). Posterior samples change.
- **`--noSingleFragProb`** stopped at the mapping pass, so it was inert on
  every RAD path. Point estimates change where orphans are contested.
- **`--discardOrphans`** was read only by the online alignment path; the
  deterministic writer now drops half-mapped fragments before they reach the
  RAD.
- **`--scoreExp`** was never forwarded to deterministic `-a`, though phase 2
  applies a score exponent to `AS`-scored placements exactly as it does to
  selective-alignment ones.
- **`--recoverOrphans` pairs are strand-judged.** The rescue site left the
  pair's orientation unset, and the one-pass filter accepts an unset format
  unconditionally, so recovered wrong-strand pairs were counted compatible
  under a stranded `-l`. The deterministic path always judged them. Mapped
  counts drop for the affected fragments, which is the correct answer.
- **`--meta`** promises a uniform optimizer start in its own log line; under
  `--online -a` — the one path with a warm start — it never set it.

### Flags that now say when they do nothing

Roughly twenty accepted-but-inert flag×mode combinations now warn by name
rather than being silently ignored: the selective-alignment scoring and seeding
knobs under `--sketch`, `--sketchStrictOrphans` without `--sketch`,
`-i`/`--index` alongside `-a` or `--rad`, `-t`/`--targets` in reads mode (the
index carries the sequences), `--writeUnmappedNames` outside reads mode, and
the RAD-shaping flags on the two paths that write no RAD. `--genome` and
`--juncMissDiscount` without `--annotation` are now errors rather than dead
flags, and `--online` under `--rad` says it has no effect instead of promising
the one-pass path — it also no longer suppresses the inert-knob warnings it
used to mask.

Two flags were repaired rather than warned about: **`--dumpBiasModels`** works
again on every path (it silently stopped writing `bias_models.txt` when
deterministic became the default, because only the one-pass writer had a dump
site), and **`--rad` now forwards `--numGCBins`, `--conditionalGCBins`,
`--biasSpeedSamp` and `--noBiasLengthThreshold`**, which it accepted while
`--gcBias` itself was wired.

Gene-level output is also consistent now: under `--skipQuant`, `-g` no longer
writes a `quant.genes.sf` full of zeros while announcing success, and in
genome-projection mode a `--skipQuant` run keeps its projected RAD as the
deliverable instead of deleting it and producing nothing at all.

### Faster, same answer

Deterministic phase 1 with `--seqBias`/`--gcBias`/`--posBias` was building the
full observed-bias apparatus — including a rank bitvector over the entire
reference concatenation — and per-fragment model updates on the mapping hot
path, none of which any consumer could read: the requant re-derives observation
from the RAD and the reference. Construction now matches its consumer's gate.
Bias-corrected runs use less memory and less mapping-thread time; estimates are
unchanged.

## Also

- The divan-based alignment-kernel benchmark harness (#1128).
- `master`→`develop` housekeeping; `postponed` label for parked work
  (#1041, #1123, #1132, #1129 record their dispositions).
