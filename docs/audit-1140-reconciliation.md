# Audit reconciliation — #1140 (2.6.0 deterministic-default readiness)

Every item raised anywhere in the #1140 thread — Ben's initial five-break review,
the 111-finding readiness sweep (sections A–G), the flip-PR checklist, the disk
story, and the two validation runs — reconciled against `develop` on 2026-08-24
(`b860b74d`). Each row was checked against the code or behaviour on `develop`,
not against the PR descriptions that claimed to fix it.

**The union of the two tables below is exactly every item the thread raised.**
Resolved: 44. Unresolved: 3, plus 2 closed by decision. The unresolved three are
collected in the tracking issue filed for them; none blocks the determinism flip.

Legend for "verified": **code** = confirmed by reading the code path on
`develop`; **behaviour** = confirmed by running the 2.6.0 binary; **PR** =
resolution merged in a reviewed PR and not independently re-run here.

---

## Resolved

### Ben's initial review (five silent breaks)

| # | Item | Fixed in | Verified |
|---|---|---|---|
| IB1 | `--dumpEq` writes an empty file under `--deterministic` (phase 1 dumps before any class exists) | #1147 | PR |
| IB2 | Five non-online flags dropped by `run_deterministic`: `--initUniform`, `--biasSpeedSamp`, `--noBiasLengthThreshold` (forward gaps), `--noLengthCorrection`, `--noFragLengthDist` (feature gaps) | #1144 (first three), #1148 (last two) | code |
| IB3 | `meta_info.json` `total_time_seconds` reports only phase 2 | #1143 | PR |
| IB4 | Alignment mode with no usable `AS` degrades to uniform weights, silently (3× MARD) | #1145 | code |
| IB5 | `-l A` `num_compatible_fragments` accounting changes by the detection prefix | release notes #1159 | code |

### A. Release-blocking output contract

| # | Item | Fixed in | Verified |
|---|---|---|---|
| A1 | Decoy references leak into `quant.sf` / bootstraps / `ambig_info.tsv` / `num_valid_targets` on every `quantify_rad` path | #1153 | code (`first_decoy_index`/`num_decoys` baked into RAD header; `quant_row_indices` excludes in the shared writer) |
| A2 | Phase 2 clobbers `cmd_info.json` (`index`/`mates`/`threads`/`sketch`), `lib_format_counts.json` ends `read_files: []` | #1153 | PR |
| A3 | `samp_type` hardcoded `"none"` even when bootstraps/Gibbs ran | #1153 | PR |
| A4 | `serialized_eq_classes` hardcoded `false` (actively wrong post-#1147) | #1153 | PR |
| A5 | `num_decoy_fragments` hardcoded `0`; genome projection never bakes `MapCounters` (placeholder `num_processed`) | #1153 | code (`genome_project.rs` bakes real `MapCounters`) |
| A6 | `frag_length_source: "reads"` reported for single-end runs whose FLD is the `--fldMean`/`--fldSD` prior | #1153 | code (`is_paired()` → `Reads` else `Prior`) |

### B. Mode-dispatch input handling

| # | Item | Fixed in | Verified |
|---|---|---|---|
| B7 | `-a` + `-1/-2/-r` silently drops the FASTQs (no clap conflict) | #1153 | code (clap `conflicts_with_all`) |
| B8 | `-i` accepted-but-inert with `-a`/`--rad`, no warning | #1155 | code |
| B9 | `-t` silent no-op in both reads modes; help text claims it enables the error model | #1155 | code |

### C. Default `-a` path gaps

| # | Item | Fixed in | Verified |
|---|---|---|---|
| C10 | `--discardOrphans` inert on deterministic `-a` | #1158 | code (`bam_rad.rs:122` consults `opts.discard_orphans`, threaded to all three workers) |
| C11 | `--scoreExp` not forwarded to deterministic `-a` | #1158 | code |
| C12 | `-l A` upstream-filter caveat warns only on the online `-a` path | #1158 | code (shared `warn_l_a_alignment_input_caveat`, called from both writers) |

### D. Inference knobs

| # | Item | Fixed in | Verified |
|---|---|---|---|
| D13 | `--perNucleotidePrior` ignored by Gibbs (hardcoded `per_transcript_prior: true`) and bootstrap replicates (`eff_lens = None`) | #1158 | code (all three Gibbs sites `!opts.em.per_nucleotide_prior`; bootstrap test `bootstrap_replicates_follow_the_per_nucleotide_prior`) |
| D14 | `--noSingleFragProb` inert on default/RAD paths; stale "[not yet implemented]" help on it and `--noFragLengthDist` | #1158 | code (wired on default/RAD; help text accurate) |
| D15 | `--meta`'s uniform-init leg not applied on `--online -a` | #1158 | code (`opts.init_uniform = args.init_uniform \|\| args.meta` before dispatch) |

### E. The warn-and-ignore matrix

| # | Item | Fixed in | Verified |
|---|---|---|---|
| E | ~20 accepted-but-inert flag×mode combinations, driven by a single mode×flag applicability table; includes `--dumpBiasModels` broken on the new default (write site gated `!skip_quant`) and `--genome`/`--juncMissDiscount` without `--annotation` | #1156 | PR + code spot-checks (sketch warns on scoring knobs; `--dumpBiasModels` restored in the shared writer) |

### F. Cross-path divergences

| # | Item | Fixed in | Verified |
|---|---|---|---|
| F16 | `--recoverOrphans` pairs bypass the strand filter in one-pass but are judged in deterministic mode (#1136 class) | #1154 | PR |
| F18 | Deterministic `-a`/genome projection bake `expected.or(detected)` as the RAD library-format tag; reads-det bakes resolved | #1155 | PR |
| F19 | `inference_truncated_mass` warning exists only on the deprecated paths | #1155 | code (moved to shared writer) |
| F20 | `extra_diagnostics` appended only by `quantify_rad` | #1155 | code |

### G. Docs and dead code

| # | Item | Fixed in | Verified |
|---|---|---|---|
| G | ~15 docs-contradict-code corrections and ~11 dead-code removals, incl. observed-bias construction matching its consumer's gate, `TranscriptGroup.valid` constant `true`, both vacuous `SimplePosBias` tests | #1157 | PR |

### Flip-PR checklist (tracked across review comments)

| # | Item | Fixed in | Verified |
|---|---|---|---|
| CK1 | `score_tag_warning` reaches `meta_info.json`'s diagnostics array, not just the log | #1150 flip | code (`run_deterministic_align` pushes `bam_score_tags_unusable` into `extra_diagnostics`) |
| CK2 | `--noLengthCorrection`/`--noFragLengthDist` under `--online -a`: warn-and-ignore | #1148 | code (`main.rs` warn on the deprecated path) |
| CK3 | `--radScratchDir` × `--skipQuant`: promoted RAD goes to the output dir under its stable name | #1149 | PR |
| CK4 | `--deterministic --skipQuant --dumpEq` warns and points at `--rad --dumpEq` | #1147 | PR |
| CK5 | Live progress spinner wired in deterministic mode (was one-pass-only) | #1150 flip | PR |

### Disk story

| # | Item | Fixed in | Verified |
|---|---|---|---|
| DS1 | Startup free-space preflight against a conservative estimate | #1146 | PR |
| DS2 | Delete the salmon-owned partial RAD (incl. `.partial-*`) on failure in either phase | #1146 | PR |
| DS3 | Periodic free-space check so ENOSPC fails in seconds, not at end-of-pass | #1146 | PR |
| DS4 | `--radScratchDir <dir>` flag (also answers in-memory-store for Linux via `/dev/shm`) | #1146 | PR |

### Validation runs

| # | Item | Fixed in | Verified |
|---|---|---|---|
| V1 | **Release blocker**: `quant.sf` not reproducible across thread counts (shard count was `rayon::current_num_threads()`) | #1167 | code + behaviour (one hash across `-p` 1/4/16/64 at `--sigDigits 9`) |
| V2 | EM 2.14× under-parallelization (missing touched-transcript optimization; dense per-shard clears) | #1172 | code + behaviour (independent A/B: 33% faster at `-p 16`, 42% at `-p 64`) |
| V3 | Release-notes speed claim was hardware-universal; wrong on older CPUs | #1159 | code (Broadwell qualification present) |
| V4 | Documented `--radCompress zstd` and RAD-size figures understated the product | #1159 | code (corrected to ~53%/5.5% and ~45/55 MB per 1M) |

### Closed by decision (not code changes)

| # | Item | Disposition |
|---|---|---|
| F17 | Online `-a` FLD divergence (no CMF length-conditioning, flat orphan penalty, FLD trained pre-filter) | Documented at the training site (#1157), not fixed — the path is scheduled for 2.7.0 removal and gets documentation, not parity work. |
| IB-rec | FASTQ parse errors cannot name the offending *record* | Won't-fix: the record identity is discarded upstream by paraseq before the error surfaces; reporting it is not salmon's responsibility. The input *file* is named (#1166). #1162 closed on this basis. |

---

## Unresolved

All three are open and none blocks the determinism flip. They are collected in
the tracking issue filed alongside this document.

| # | Item | State | Why it is still open |
|---|---|---|---|
| U1 | **EM parallel scaling ceiling.** #1172 removed the large per-iteration constant and the anti-scaling regression (EM was slower at `-p 64` than `-p 16`), but the phase still does not scale meaningfully past ~16 threads — it is memory-bandwidth-bound on the per-shard accumulators. | Open, post-tag | Never a correctness bug. Triaged post-tag from the moment it was found ("Post-tag, in my view" — slow-disk validation comment). #1172 was the optimization pass and closed the regression; the residual ceiling is genuinely harder and was scoped as follow-up. Independently confirmed in the #1172 review: bit-identical across `-p`, ~4% slower than online at `-p 16`, ~2× faster at `-p 64`. |
| U2 | **#1141 — mapping-record self-consistency.** `--writeMappings`/`--writeBam` records where the CIGAR is synthesized from the read length and the `AS` describes a different alignment than the CIGAR beside it (2.4% of records self-contradict on an `AS`-vs-`NM` audit). | Open PR #1141 | Only the MAPQ slice of #1141 landed (#1151 derives MAPQ from the placement count; #1171 sets MAPQ 0 not 255 for unplaced records). The rest was split out before the audit on a "demonstrated need" question, re-raised for 2.6.0 inclusion in the thread, and left as an open scoping decision ("2.7.0 or #1141 only") that was never made. Opt-in output only — does not touch the default `quant.sf` contract. |
| U3 | **#1142 — mapping-output completeness.** `QUAL`, `M5`/`UR`, `MC`/`MQ`, `ZW`, `@RG`, unaligned records — fields SAM permits to be absent. Their value is checkability: with qualities and `M5` the file round-trips through `samtools fastq` and names its reference. | Open PR #1142 | Explicitly optional per the thread ("optional and I am not claiming otherwise"). Stacked behind #1141's decision. |
