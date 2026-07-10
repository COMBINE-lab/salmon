# salmon (Rust) pre-release testing checklist

A reusable checklist of the manual/semi-automated tests we run before tagging a
release, in addition to the in-tree `cargo test` suite and the cross-platform
CI. It exists so each release exercises the same surface area and so regressions
in paths the unit tests do not cover (real data, decoy indices, alignment mode,
messy inputs, parity vs C++) are caught before publishing.

Treat every box as "must pass or have an understood, documented reason." The
commands are templates — substitute your own index/reads paths. Where a test
guards a specific historical bug, the issue/fix is noted so the guard is not
silently dropped.

## 0. Datasets used (assemble once per release)

- **`sample_data/`** (in repo): 15 transcripts, 10k paired 50 bp reads with the
  true transcript in each read name. Fast smoke tests + determinism.
- **Simulated human** (polyester over GENCODE/Ensembl cDNA, ~194k transcripts):
  `clean.fa` transcriptome + `sim_R{1,2}` (easy) / `sim_hard_R{1,2}` (hard) +
  `ground_truth_{easy,hard}.tsv` (`name<TAB>true_count<TAB>length`). Accuracy vs
  a known truth, including short (<k) transcripts.
- **Real + decoy**: GRCh38 cDNA+ncRNA gentrome with the primary assembly as
  decoy (`salmon index … -d decoys.txt --keepDuplicates`), and a public
  paired-end RNA-seq run (we use SRR1039508, ~22.9M PE). Decoy handling, the
  bias/abundance path at scale, and real-data parity vs C++.
- **Reference C++ salmon** (the matching `1.x` release / fixed build) for parity.
- **bowtie2 + samtools** (module/`PATH`) for alignment-mode inputs.

A correlation helper is handy: read two `quant.sf` (or a `quant.sf` and a truth
file), report Spearman, **log**-Pearson (raw Pearson is dominated by a few
high-TPM transcripts), MARD, and nonzero-set Jaccard. Use log-scale metrics for
TPM throughout.

## 1. Build / static gates

- [ ] `cargo build --release` clean
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` clean
- [ ] `cargo test --workspace` green
- [ ] `cargo fmt --all --check` clean
- [ ] CI green on all target platforms (linux x86_64/aarch64, macOS arm64/x86_64)

## 2. Index build + format

- [ ] Build a plain transcriptome index; build a decoy-aware index
      (`-d decoys.txt --keepDuplicates`); both load.
- [ ] `info.json` records `index_version`; loading an index built by a release
      older than the current `MIN_READABLE_INDEX_VERSION` is **rejected** with an
      actionable rebuild message (guard for the index-format-version bump).
- [ ] Short (`< k`) transcripts are retained and appear in `quant.sf` with their
      true `Length`, a sane `EffectiveLength` (≥ 1), and 0 reads. No transcript
      from the input FASTA is missing from `quant.sf`.
- [ ] Decoys ≤ k, decoys with N-runs, and short transcripts compose correctly
      (covered by the `build_composes_*` index unit test; re-confirm it runs).

## 3. Quantification modes (reads)

Each should complete, conserve reads, and produce a sane `quant.sf`.

- [ ] Selective alignment, paired-end (`-l A`)
- [ ] Sketch / pseudoalignment, paired-end (`--sketch`)
- [ ] Single-end (`-r`), both SA and sketch
- [ ] Multiple input files / lanes (`-1 a.fq b.fq -2 c.fq d.fq`), plain **and**
      gzipped — mapped count matches the concatenated single-file run
- [ ] Explicit library types `IU` / `ISF` / `ISR` and auto `A` all run; `ISF`+`ISR`
      partition an unstranded set
- [ ] `--useEM` (vs default VBEM): both run; near-identical ranking on clean data
- [ ] `--numBootstraps N` and `--numGibbsSamples N`: inferential replicate files
      written under `aux_info/bootstrap/`
- [ ] `--skipQuant`: equivalence classes emitted, no `quant.sf`

## 4. Bias / fragment-length / length correction

- [ ] `--seqBias`, `--gcBias`, `--posBias`, and all three composed — each
      completes and shifts quant without breaking the expressed set
- [ ] **`--seqBias --gcBias` on a decoy-aware index completes without a
      single-core stall** (regression guard for the #1019 abundance-phase hang)
- [ ] **Mass conservation**: `Σ NumReads` in `quant.sf` ≈ `num_mapped`
      (loss < ~1 fragment) with bias on a decoy index (guard for the log-space
      eq-class normalization fix)
- [ ] Bias-corrected quant vs C++ is at the expected concordance level
- [ ] `--noFragLengthDist` measurably changes quant (the FLD term is actually
      used) in **both** SA and sketch; `--noSingleFragProb`, `--noLengthCorrection`
      run
- [ ] FLD is trained from data (observed `flenDist` mean differs from the prior)
      in both SA and sketch

## 5. Output / auxiliary surfaces

- [ ] `quant.sf`, `cmd_info.json`, `lib_format_counts.json`,
      `aux_info/meta_info.json` present and structurally complete
- [ ] `-g t2g` gene aggregation (`quant.genes.sf`) and `salmon quantmerge`
- [ ] `--dumpEq` / `--dumpEqWeights` (`eq_classes.txt.gz`), `--writeUnmappedNames`
- [ ] `--writeMappings` SAM: realistic positions/flags; valid in SA **and**
      sketch (sketch must not emit `POS=1` for every record); does not crash with
      `--gcBias` (incl. on a decoy index — Rust has no analog of the C++ #1010
      `--writeMappings`+`--gcBias` segfault)
- [ ] bias model dumps (`obs/exp` seq/gc/pos) written under the bias flags;
      `fld.gz` always

## 6. Alignment-based input (`-a`)

Generate alignments with bowtie2 (transcriptome index), reporting many
multimappers (`-k 100`).

- [ ] `-k 100` multimapping BAM: completes, mapped count and quant match C++ `-a`
- [ ] **Messy BAM without `--no-discordant`/`--no-mixed`** (genuine discordant +
      mixed/singleton + secondary records): does **not** panic; mapped count
      matches C++; discordant/mixed mates degrade to orphan placements, nothing
      dropped or assumed-paired
- [ ] **Coordinate-sorted BAM is rejected up front** with an actionable
      "collate by read name" message (allowing `GO:query` / `SO:queryname`)
- [ ] Messy-vs-clean (`--no-*` flags) differ only gracefully (no crash/garbage) —
      confirms those flags are advisory, not required

## 7. Determinism / robustness / edge cases

- [ ] `-p 1` run twice → byte-identical `quant.sf`
- [ ] `-p N` (N>1) run twice → negligible run-to-run TPM drift
- [ ] gzipped input == plain input (identical `quant.sf`)
- [ ] empty input → graceful error (no panic)
- [ ] malformed FASTQ (record boundary broken) → clear error (no panic)
- [ ] mismatched mate counts → graceful error
- [ ] all-unmapped input (foreign reads) → 0 mapped, valid all-zero `quant.sf`,
      no crash

## 8. Accuracy & parity (the headline numbers)

- [ ] **Sim vs ground truth** (easy + hard, SA + sketch): Spearman/MARD in the
      expected band; SA ≥ sketch; enabling the FLD improves both
- [ ] **Real data vs C++** (SA): mapping rate matches C++ to ~0.01%; TPM
      log-Pearson ≈ 0.99 / Spearman ≈ 0.97
- [ ] **Real data sketch**: mapping rate sane; concordance with SA in the
      expected band
- [ ] Runtime / peak RSS in line with the previous release (no large regression)
- [ ] No errors/panics/unexpected warnings in any run log

## 9. Publish gates

- [ ] All release work merged to `master`; `master` CI green
- [ ] Release notes finalized (`docs/release-notes-<version>.md`, no "draft"
      marker), covering every behavioral change
- [ ] External deps (cf1-rs, piscem-rs, ksw2rs) published on crates.io at the
      versions the workspace requires
- [ ] `scripts/bump_and_publish.sh <version> --dry-run` shows the correct version
      bump, tag, and dependency-ordered crate list; leaf-crate
      `cargo publish -p salmon-core --dry-run` packages + verifies cleanly
- [ ] Publish: `scripts/bump_and_publish.sh <version> --publish` (bump → commit →
      tag → push → cargo-dist binaries → publish the 9 `salmon-*` crates)

## Known gaps / candidates to add

Not yet part of the routine sweep; add as bandwidth allows or if a release
touches the relevant code:

- Mate-pair / outward library types (`MSF`/`OSF`/…) end-to-end.
- `--fullLengthAlignment`, `--mimicBT2` / `--mimicStrictBT2` mapping presets.
- `--noEffectiveLengthCorrection` (distinct from `--noLengthCorrection`).
- FLD-prior knobs (`--fldMean`/`--fldSD`/`--fldMax`) and multimapping caps
  (`--maxReadOcc`/`--maxOccsPerHit`).
- Truncated/corrupt gzip and truncated BAM robustness.
- Format-level diff of `eq_classes.txt.gz` and the binary bias dumps vs C++
  (currently presence + shape are checked, not byte content).
- `--gencode` reference-name munging at index build.
