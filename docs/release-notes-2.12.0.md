# salmon 2.12.0

A cleanup, pairing and speed release:
- **The one-pass `--online` path is removed.** It was deprecated in 2.6.0. The removal does not change default results.
- **`--recoverOrphans` is exact and 28–51× faster.** Its output is byte-identical to the slower search's. A run with the flag now takes 1.17–1.65× as long as the same run without it.
- **`-a` reads SAM or BAM from a pipe**: a FIFO, `<(...)` or `/dev/stdin`.
- **Selective alignment pairs small dovetails by default.** Mates of a short insert that read past each other into adapter are paired when the overhang is at most 32 bases. The adapter bases are clipped and score nothing.
- **`--posBias` no longer gives effective lengths longer than the transcript.** On 500k muscle pairs, 23,206 transcripts had one, up to 1.99×.
- **Faster inferential replicates and bias correction.** The Gibbs posterior phase runs in parallel, `bootstraps.gz` is compressed on every core, and `--gcBias` and `--seqBias` are 14% and 6% faster. Peak memory of a paired-end reads run at `-p 16` is up to 13% lower than in 2.11.0.
- **Fixes from an audit of 2.11.0 (#1314)**, mostly in alignment mode (`-a`), `--skipQuant`, Gibbs sampling and output files.

**What changes in results.** Runs with default options change in these places:
- **Paired-end selective alignment:** small dovetails now pair. On the full muscle sample (`-l ISR`), 42,902 more fragments are concordant, 116 more map, and 9,583 reads move. On 10M GEUVADIS pairs, 16,338 more are concordant and 4,549 reads move. A pair that starts before its transcript also gets the right length, which moves 0.8 reads on 500k muscle pairs.
- **Paired-end `--sketch`:** fragments are one base shorter, so the fragment-length distribution (FLD) mean drops by about one base. On full muscle, 4,078.5 reads move.
- **Single-end selective alignment:** on read 1 of 500k muscle pairs, 4.4 reads move across 61 transcripts.
- **Genome projection (`-a <genome BAM> --annotation`):** the projected pair now comes from one genomic mate pair. 30k reads move on ERR188044 and 55k on muscle.
- **Gibbs samples** (`--numGibbsSamples`) differ from 2.11.0 at the same seed. They are statistically equivalent, and they no longer depend on `-p`.

Runs with these options change in more places:
- `--posBias`;
- `--gcBias` on single-end reads;
- `--allowDovetail`;
- `--recoverOrphans` on reads with lowercase bases;
- Gibbs sampling with `--perNucleotidePrior`;
- `-a` with `--errorModel` or bias correction, when `-t` lacks some references or its sequences carry poly-A tails past the BAM's lengths;
- `-a` with `-l IU` or `-l U` and bias correction;
- the mapping rate that `--skipQuant` reports.

The details are under [Changes to results](#changes-to-results). Results that a 2.11.0 `--online` run gave are compared with the default under [Removed: the one-pass `--online` path](#removed-the-one-pass---online-path).

No index rebuild is needed: the index format is unchanged.

## Changes to results

### Selective alignment: small dovetails pair by default

**A pair runs between the mates' 5′ ends** (#1356). For mates that face each other, the fragment is `[start of the forward-strand mate, end of the reverse-strand mate)`. Each read starts at one end of the molecule, so bases that a read has past the other read's 5′ end are adapter read-through, not insert. This holds for every inward library (`ISF`, `ISR`, `IU`). Pairing, the mapper's fragment start and rescued pairs used to measure the union of the two mates. Outward-facing pairs keep the union.

**Small dovetails pair without `--allowDovetail`** (#1357). When the insert is shorter than the reads, each mate runs past the other's 5′ end. Selective alignment rejected every such pair unless `--allowDovetail`, so these fragments were kept as orphans: 2.4% of fragments on a stranded muscle library (2×301, trimmed) and 0.3% on GEUVADIS. Now:
- **A pair whose overhang is at most 32 bases pairs.** The overhang is `max(fw end − rc end, fw start − rc start, 0)`, and the fragment is `[fw start, rc end)`. Orphan rescue uses the same test. The hidden `--dovetailTolerance` sets the limit (see [New features](#new-features)).
- **Only the insert is scored.** In a dovetailed pair, each mate's bases past the other read's 5′ end are soft-clipped and score 0, with no penalty beyond what `--softclip` gives a clipped end. `--minScoreFraction` is measured against the kept bases, so a short fragment is not filtered for its adapter. This applies to pairs admitted by the tolerance or by `--allowDovetail`, in default and `--fullLengthAlignment` scoring, and in orphan rescue. Pairs where neither mate overhangs the other are untouched.
- **`--writeMappings` shows the clipped bases as soft clips** (e.g. `5S80M10S`), and the mate position moves to the fragment start. RAD output already followed the 5′-end rule.

Without the insert-only scoring, pairing these fragments made accuracy worse. Their two mates share the read-through bases, which align against isoform-specific downstream sequence, so adding both mates' scores counted that noise twice.

Accuracy against truth: 2M simulated fragments, 2×100 reads on GENCODE v50, fragment length N(120, 45), so 31% of fragments are shorter than the reads. "Misassigned" is half the total absolute error in per-transcript counts.

| simulation | run | transcript Spearman | misassigned reads | gene-level misassigned | FLD mean (true 122.5) |
|---|---|---|---|---|---|
| untrimmed adapter + trimming residue | before | 0.7249 | 320,620 | 53,909 | 142.87 |
| | tolerance 32, no insert-only scoring | 0.7175 | 339,140 | 54,328 | 131.43 |
| | **2.12.0** | **0.7398** | **299,306** | **39,990** | **128.33** |
| trimming residue only | before | 0.7354 | 304,673 | 18,685 | 141.52 |
| | 2.12.0 | 0.7431 | 294,403 | 16,967 | 122.83 |
| clean (no read-through) | before | 0.7435 | 294,262 | 17,170 | 122.83 |
| | 2.12.0 | 0.7434 | 294,259 | 17,170 | 122.83 |

On real data, before → after, at `-p 16`:

| | GEUVADIS 200k `IU` | muscle 500k `ISR` | full muscle `ISR` | GEUVADIS 10M `IU` |
|---|---|---|---|---|
| concordant fragments | +333 | +4,757 | +42,902 | +16,338 |
| orphans | 12,021 → 11,720 | 14,227 → 9,482 | 126,762 → 83,976 | 622,544 → 607,975 |
| `num_dovetail_fragments` | 494 → 134 | 12,830 → 6,822 | 113,939 → 59,490 | 24,670 → 6,941 |
| mapped | +28 | +12 | +116 | +1,517 |
| FLD mean | 171.78 → 171.57 | 164.04 → 164.97 | 164.43 → 165.38 | 170.06 → 169.86 |
| reads moved | 106 | 1,207 | 9,583 | 4,549 |

`--dovetailTolerance 0` gives the same output as before #1357 on all four data sets.

**`--allowDovetail` admits inward dovetails only** (#1356). A pair whose forward mate starts at or past the reverse mate's end faces outward. `--allowDovetail` used to form these pairs too. The RAD stores only the two strands, so quantification read them back as inward pairs and added their outer span to the FLD. Outward pairs are now never formed, in pairing or in orphan rescue. Without `--allowDovetail` they keep their current treatment: the better mate stays as an orphan. On real data these pairs are rare and look like repeat-copy artefacts, not circular RNAs. Fragments whose only pairing is outward are 98.5 per million in GEUVADIS and 1,115 per million in muscle, 98% of them in NEB and TTN.

With `--allowDovetail`, before → after, at `-p 8`:

| | GEUVADIS 200k `IU` | muscle 500k `ISR` | full muscle `ISR` |
|---|---|---|---|
| reads moved | 24.1 | 146.7 | 1,025.9 |
| FLD mean | 171.66 → 171.62 | 165.09 → 164.95 | 165.54 → 165.38 |

The 5′-end rule alone, with `--allowDovetail`, moves 347 reads on full muscle and lowers the mean fragment length by 0.18 bases. Mapped counts are unchanged.

### Sketch mode: fragment length and start

Sketch mode (`--sketch`) reported each mapped pair one base too long. For the pair r1 = t[50..130], r2 = rc(t[350..430]), it gave 381 bases and a 3′ end of 430, where selective alignment gives 380 and 429. The off-by-one was in piscem-rs's template length and is fixed there, in piscem-rs 0.12.4. A sketch pair also now starts at its forward mate, rather than at the smaller of the two mate positions: the same 5′-end rule as selective alignment. A start below 0 is clipped to 0 and the fragment is shortened by the same amount.

- Sketch FLD means drop by about one base: 200k GEUVADIS pairs go from 172.7111 to 171.7112, and 500k muscle pairs from 167.6435 to 166.6438.
- With the start rule as well, on full muscle 20,946 transcripts change, 919 of them by at least 1 read, and 4,078.5 reads move.
- Mapping counts are unchanged. Selective-alignment output is unchanged.

(#1360)

### Selective alignment: a pair that starts before its transcript

A pair whose projection starts before its transcript was written to the RAD with its start clamped to 0 but its full length. The record ran past the fragment's real end by the clipped amount, so the reverse mate's 5′ end, the fragment's GC window and its FLD weight were all off by that amount. The start is now clipped and the length shortened to match. Such placements are rare: 322 of 8.7M proper-pair placements on 500k muscle pairs (0.004%), and 9,703 of 26.3M (0.04%) on a GC-bias simulation with realistic abundances. Most overhang by 1 base.

| run | transcripts changed (≥ 1 read) | reads moved | FLD mean |
|---|---|---|---|
| muscle 500k | 452 (0) | 0.8 | 164.0417 → 164.0373 |
| muscle 500k, `--seqBias --gcBias --posBias` | 2,946 (20) | 107.8 | same |
| hard simulation | 7,191 (52) | 242.8 | 198.6669 → 198.6657 |
| easy simulation | 10,008 (0) | 9.9 | 248.6006 → 248.5994 |
| GC simulation, realistic abundances, `--gcBias` | 4,361 (0) | 3.9 | 248.6579 → 248.6565 |

Mapping counts are unchanged. On the hard simulation, well-expressed transcripts go from Spearman 0.98770 to 0.98784 and MARD 0.0845 to 0.0844; the other simulations are identical to 4 decimals. `--sketch` and `-a` are byte-identical. (#1364)

### Genome projection: one mate pair per transcript

In genome-projection mode (`-a <genome BAM> --annotation <GTF>`), the projected record took the first read-1 placement and the first read-2 placement that bramble listed on a transcript. When the aligner reports several alignments of one fragment (STAR multimappers, or repeats inside one gene such as TTN), those two need not be mates. The score, position and fragment length could then come from different alignments, and the result depended on bramble's output order.

The reported pair is now the best-scoring read-1/read-2 pair whose genomic records are mates: each points at the other's start and, when both carry `HI`, the values agree. Without such a pair, salmon falls back to the first read-1 and read-2 placements, as before. Orphans are unchanged. The output no longer depends on bramble's order, and it is the same at `-p 1`, 16 and 32.

- **Reads moved:** 30k on ERR188044 and 55k on muscle, mostly between TTN isoforms.
- **Accuracy** on the genome-quantifiable subset, Spearman / MARD with bramble 0.2.1: the multimapping simulation goes from 0.87304 / 0.51491 to 0.87689 / 0.50399, the easy one from 0.87328 / 0.51425 to 0.87349 / 0.51396, and the hard one from 0.85725 / 0.57589 to 0.85705 / 0.57622 (noise).
- **bramble-rs 0.2.2** pairs mates by `HI` when several alignments of a fragment sit at the same positions. On top of the pair selection, it projects 16.7k more fragments on ERR188044 and 560 to 880 more on the simulations. On the multimapping simulation Spearman goes from 0.87689 to 0.87765 and MARD from 0.50399 to 0.50147.
- **Overhead:** none measurable.

(#1359)

### `--posBias` effective lengths

With `--posBias`, many transcripts got a bias-corrected effective length longer than the transcript. On 500k muscle pairs (`-l A`), 23,206 transcripts did, 5,880 of them expressed, up to 1.99×. TTN-201 (82,380 nt) got 163,808. Two defects stacked:
- **The spline knots were in the wrong place.** Fragment positions are binned into 20 uniform bins, but the spline was fitted through non-uniform knots, so each bin's value was drawn at the wrong position. The obs/exp ratio no longer averaged to 1 along a transcript: 1.28 per end in the longest length class on muscle. The knots now sit at the bin centers.
- **The positional effect was counted twice.** The 5′ and 3′ models each measure the whole positional effect, so their product averaged above 1. A fragment's weight is now the geometric mean of its two ratios, which reproduces both observed marginals and averages to 1.

| data | effective length > length, before | after | reads moved |
|---|---|---|---|
| muscle 500k pairs | 23,206 | 0 | 2,639 (0.6%) |
| GEUVADIS 200k pairs | 9,897 | 0 | 525 (0.3%) |

TTN-201's effective length goes from 163,808 to 81,112, and TTN's gene TPM from 2,793 to 4,947. On a simulation with uniform fragment starts, where `--posBias` should change nothing, the old code made the TPM error 2.2% worse than no correction; the fixed code makes it 0.4% worse. On a simulation with 3′ positional bias, `--posBias` went from Spearman / MARD 0.7663 / 0.474 to 0.7701 / 0.468, against 0.7692 / 0.469 with no correction. Under real positional bias the correction still does not move the error much (see [Known limitations](#known-limitations-and-follow-ups)). Runs without `--posBias` are unchanged. (#1352, #1361)

### Gibbs samples

- **Each Gibbs round now runs in parallel, fixed-size blocks** (4096 transcripts, 2048 equivalence classes), each with its own random stream. Samples differ from earlier versions at the same seed, so `bootstraps.gz` is not byte-identical. They are statistically equivalent: on 10M GEUVADIS pairs with 100 samples, per-transcript posterior means differ from the old sampler's by a median 1.31%, where two independent halves of one old-sampler run differ by 1.34%, and the posterior SDs match (median ratio 0.999). Samples no longer depend on the thread count: `-p 1`, 4 and 16 give the same `bootstraps.gz`. `quant.sf` is unchanged. (#1347)
- **Gibbs sampling uses the optimizer's prior on every path.** The default reads path and alignment mode used `max(1, vbPrior)` per nucleotide. With `--perNucleotidePrior --vbPrior 0.01`, a 1500 bp transcript got about 1500 pseudo-reads instead of about 15. Under `--useEM`, `--perNucleotidePrior` also switched the Gibbs prior to per-nucleotide. The rule is now: the per-nucleotide `vbPrior` for VBEM with `--perNucleotidePrior`, `max(vbPrior, 1)` per transcript for VBEM otherwise, and `1e-3` per transcript for EM. On 200,000 GEUVADIS pairs with `--useEM --perNucleotidePrior`, the posterior mass on transcripts with a zero point estimate drops from 95,531 to 2,314 reads, which is what `--useEM` alone gives. Default flags give byte-identical Gibbs output, and `quant.sf` does not change. (#1316)

### Other options and modes

- **`--gcBias` on single-end reads is skipped, with a warning.** Only proper pairs fill the observed fragment-GC table. With single-end reads the table stayed empty, normalized to flat, and the correction still rescaled effective lengths. On 200,000 GEUVADIS reads run as single-end (`-l U`), `--gcBias` moved 26,208 reads and changed 17,021 effective lengths. The run now gives the same `quant.sf` as without `--gcBias`, logs `--gcBias: no paired fragment was observed (single-end reads?), so fragment-GC correction is skipped`, and writes `gc_bias_correct: false`. Paired-end `--gcBias` is unchanged. (#1315)
- **`-a` with a `-t` that lacks some references:**
  - With `--errorModel`, a placement on a reference with no sequence, or a record with no bases, scored 0, which equals a perfect alignment. It now gets the fragment's worst real score.
  - `--gcBias` read past the end of a missing sequence, into the next transcript's.
  - On 200k GEUVADIS pairs with a `-t` missing every 10th transcript, `--errorModel` gave 49% of the reads to the missing 10%; 2.12.0 gives them 7.4%, against 10.8% with the full `-t`. The result is 6,936 reads away from the full-`-t` run instead of 77,778. `--gcBias` with the partial `-t` is 6,975 reads away instead of 11,461. With the full `-t`, the output does not change.
  - Genome projection gives a transcript with an exon it cannot slice (for example, a `chr1` vs `1` contig mismatch) an empty sequence and warns with the count. It used to skip such exons silently.

  (#1319)
- **`-a` and genome projection with `-l IU` or `-l U` and bias correction** seed the bias models under the strand rule given. An explicit unstranded type used to fall back to the detected format in the seed EM. On 500,000 stranded (ISR) muscle pairs aligned to GENCODE v50, `-a -l IU --seqBias --gcBias` moves 2,179 reads (0.46%). `-l ISR` output is unchanged. (#1325)
- **`-a` with a `-t` FASTA longer than the alignment's reference by a poly-A tail** is cut back to the header length (`@SQ LN`), so bias correction no longer runs over the tail. This happens when a `--writeMappings` SAM is requantified with `-a -t` and the original FASTA; GENCODE v50 has 11,659 transcripts with such tails. On 500,000 muscle pairs with `--seqBias --gcBias --posBias`, 525 expressed tailed transcripts change effective length (median 0.44%, max 9.4%), and 304 reads move. Without bias flags the output is byte-identical. Any other length mismatch is untouched and still warns. (#1338)
- **`--recoverOrphans` compares bases case-insensitively.** With a reverse-strand anchor and a lowercase (soft-masked) partner, every base was a mismatch, so the rescue failed. With every base of 500k muscle and 200k GEUVADIS pairs lower-cased, `--recoverOrphans` used to move 3,169 and 718 reads; the lower-case input now gives exactly the upper-case result. (#1332)
- **`--recoverOrphans --allowDovetail` spans a rescued pair over both mates**, as normal pairing does. An anchor at [200,300) and a partner at [220,280) gave 80, shorter than the read. This moves 0.1 reads on 500k muscle pairs and none on 200k GEUVADIS pairs. Without `--allowDovetail` the output is byte-identical. (#1322)
- **Single-end reads apply the consensus filter before keeping one chain per target.** A target's highest-coverage chain could fail the score filter while another of its chains passed, and the target vanished from the read. On read 1 of 500k muscle pairs as single-end, 4.4 reads move across 61 transcripts, and the mapping count is unchanged. 200k GEUVADIS reads, read 1 or read 2 alone, are byte-identical, and so are paired-end runs. (#1333)
- **`--skipQuant` with a stranded `-l` no longer counts wrong-strand fragments as mapped.** On 200k GEUVADIS pairs with `-l ISR`, `--skipQuant` reported 92.27% mapped and the full run 46.53%; both now report 46.53%, and the orphan count matches too. If every fragment is wrong-strand, `--skipQuant` now stops with the "0 fragments mapped" error, as the full run already did. Full runs are unchanged. (#1334)

## Removed: the one-pass `--online` path

`--online` selected the pre-2.6 one-pass quantification path. It was deprecated in 2.6.0, when the deterministic two-phase path became the default (#1140), and it is now removed. Every `salmon quant` maps to an intermediate RAD and then quantifies that RAD with a fixed fragment-length distribution, so `quant.sf` is byte-identical at any `-p`. The removal does not change default results.

What happens to existing command lines:

- **`--online` is an error.** The message says to drop the flag. It names `--errorModel -t <FASTA>` for error-model scoring in alignment mode and `--radScratchDir <dir>` for placing the intermediate RAD (for example on `/dev/shm` or node-local scratch).
- **`-f/--forgettingFactor`, `--numAuxModelSamples` and `--numPreAuxModelSamples`** tuned only the online phase. They are still accepted, so old scripts and C++-era command lines keep running, but they print a warning, do nothing, and are hidden from `--help`.
- **`--deterministic`** (already the default) and **`--initUniform`** (every EM already starts uniform) are still accepted and do nothing.

If you are moving from `--online` to the default, these differences may show up:

- **Alignment mode (`-a`) scores placements by the BAM `AS` tag.**
  - `--online -a -t <FASTA>` trained an alignment error model by default.
  - To get error-model scoring now, add `--errorModel -t <FASTA>`. It trains an order-independent counting model in a second pass over the BAM.
  - Every benchmark with ground truth found `AS` scoring at least as accurate, in a single pass over the BAM.
- **Piped alignment input with `--errorModel` is not supported.**
  - `--errorModel` reads the BAM twice, so the BAM must be a regular file.
  - Write the BAM to disk first, or keep the default `AS` scoring, which needs only one pass.
  - `--online -a` was the only single-pass error-model path, and it is gone.
- **`-l A` is detected from the whole run.**
  - The default detects the library type from every fragment, then filters against the result when it quantifies.
  - `--online` locked the type in after a prefix of about 50k fragments and never judged the fragments in that prefix.
  - On small inputs, or inputs with mixed orientations, a few more fragments can be counted as incompatible. `num_compatible_fragments` covers every fragment that was judged.
- **Single-end reads without `--fldMean`/`--fldSD` use a widened default prior.**
  - The single-end fragment-length term uses the default prior with its SD widened, so it stays robust when the default does not fit the library.
  - If you know the library's fragment-length distribution, pass `--fldMean`/`--fldSD` and it is used exactly as given.
  - On a simulation whose true distribution equals the default prior, the widening costs some accuracy: TPM Spearman 0.868 against 0.873 for `--online`. Passing the true mean and SD matches `--online` (0.873).

Paired-end accuracy against simulated truth, scored on TPM against true molar abundance, is unchanged from `--online`: Spearman agrees to within 0.0002 and MARD on well-expressed transcripts to within 0.003, in no consistent direction. This holds for selective alignment, `--sketch`, a high-multimapping set and `--seqBias --gcBias --posBias`. On simulations with real GC and 3′ positional bias, with the matching correction, Spearman agrees to within 0.002. On the GC simulation with realistic abundances, MARD on well-expressed transcripts is a little higher than `--online`'s (0.295 against 0.289 with `--gcBias`, from 0.369 with no correction), so the default keeps about 92% of the improvement.

The default uses less CPU in every measured cell and less wall time in all but one: with reads-mode bias correction it takes about 5% more wall time at `-p 16`, which is the rough seed EM it runs before baking the RAD. Measured on 4.65M simulated pairs at `-p 16`, wall time in seconds:

| workload | default | `--online` |
|---|---|---|
| selective alignment | 3.20–3.30 | 3.57 |
| `--sketch` | 1.67–1.69 | 2.04–2.12 |
| `-a`, `AS` scoring (1M pairs) | 0.90 | 1.44 |
| `-a`, error model | 2.99 | 3.32 |
| reads, `--seqBias --gcBias --posBias` | 10.30–10.37 | 9.75–9.83 |

Two 2.12.0 fixes to the `--online` path landed before its removal and have no effect now: `-a --online --posBias` learned no observed positional model (#1320), and `--online` called `--writeRad` inert while it wrote the RAD (#1326).

(#1361)

## For library users

The `salmon-*` crates are versioned with salmon, so these changes ship in 2.12.0 of each crate.

**Breaking:**
- **`salmon_quant::quantify` is now `salmon_quant::map_to_rad(opts, rad_path)`.**
  - It runs only the mapping pass: it writes the RAD and the mapping-pass metadata, and estimates **no abundances**.
  - To quantify, call `salmon_align::quantify_rad` on the RAD it wrote. That is what the `salmon` binary does.
- **`QuantResult` is now `MapResult`**, without the abundance, TPM, posterior and bias fields.
- **`QuantOptions` no longer has** `deterministic_fld`, `forgetting_factor`, `num_aux_model_samples`, `num_pre_aux_model_samples` or `init_uniform`. `write_rad`, `skip_quant` and the eq-class dump flags remain; they are read by the caller's two-phase driver, not by `map_to_rad`.
- **Removed with `--online`** (#1361):
  - `salmon_align::quantify_alignments`, which was the one-pass alignment-mode path. Use `write_alignment_rad` + `quantify_rad`.
  - `AlignQuantOptions.{forgetting_factor, init_uniform, num_aux_model_samples, num_pre_aux_model_samples}`.
  - `salmon_infer::OnlineInference`.
  - `FragmentLengthDistribution::{refresh_online, online_snapshot, online_cmf_snapshot}`.
- **Removed, with no callers outside their own tests** (#1328):
  - `salmon_core::math::log_sub` and `log_sum`;
  - `salmon_eqclass::IdentityHasher` and `IdentityBuildHasher`;
  - `salmon_map::score::finalize_mappings` and `finalize_mappings_counted` (use `finalize_mappings_counted_into`);
  - `impl Default for salmon_rad::FragmentChunkBuf`.
- `salmon_map::PairingConfig` has a new public field, `dovetail_tolerance`; `salmon_map::DEFAULT_DOVETAIL_TOLERANCE` is 32. (#1357)

**Added:**
- `salmon_core::asctime_now`. `salmon_align::asctime_now` stays, as a re-export. (#1327)
- `salmon_infer::GibbsOptions::for_em`, which derives the Gibbs prior from the EM options. (#1316)
- `salmon_rad::clip_pair_start(pos, fragment_len)`, which clips a negative start to 0 and shortens the length to match. (#1364)
- `salmon_model::dumps::write_gz_segments`, the parallel single-member gzip writer behind `bootstraps.gz`. (#1343)
- `salmon_map::align::find_in_window`, the exact fast window search behind `--recoverOrphans`. `align_in_window` is unchanged. (#1355)
- `salmon_map::align_chain_clipped_view`. (#1357)

## New features

### `-a` reads from a pipe

`salmon quant -a` reads SAM or BAM from a FIFO, a process substitution or `/dev/stdin`:

```sh
salmon quant -l A -o sample -a <(samtools view -b ...)
```

- **One open per run.** The input used to be opened up to three times: for the header, for the `-l A` pairedness peek, and for the record pass. On a pipe, each later open started mid-stream. The header is now read once and the record pass streams through the same reader. The `-l A` peek keeps what it read and the pass replays it. This covers the default path and genome projection (`--annotation`).
- **Format by content.** The first 64 KiB decide the format: BGZF that decompresses to `BAM\1` is BAM, anything else is SAM (plain or compressed). The extension is no longer used, so a SAM named `.bam`, or a file with no extension, reads correctly.
- **Clearer errors.** A missing input says "opening alignment input X"; an empty one says "alignment input X is empty".
- **`--errorModel` still needs a file**, since it reads the records twice. On a pipe it stops at once with a clear message.

Before, a FIFO hung, and `<(cat x.bam)` or `cat x.bam | salmon quant -a /dev/stdin` failed with "invalid BAM header". The same input through `<(cat ...)`, `/dev/stdin` and a FIFO gives a `quant.sf` byte-identical to the file. (#1339, #1351)

### `--dovetailTolerance`

A hidden expert flag: the largest overhang, in bases, at which inward mates that extend past each other's 5′ end still pair without `--allowDovetail`. The default is 32, and 0 rejects every dovetail as 2.11.0 did. `--allowDovetail` stays the public switch for any inward dovetail. Single-end runs warn that the flag does nothing. (#1357)

## Speed and memory

### Compared with 2.11.0

Release binaries built from the v2.11.0 tag and from 2.12.0, on x86 (AMD EPYC), 16 pinned cores unless marked. Wall times are medians of three runs for the first two rows and the decoy row, and single runs otherwise.

| Run | Wall (s) | Peak memory (GB) | Mapped | Reads moved |
|---|---|---|---|---|
| full muscle, `ISR` | 8.37 → 8.18 | 1.95 → 1.73 | +116 | 9,589 |
| GEUVADIS 10M pairs, `IU` | 9.45 → 9.44 | 1.10 → 0.97 | +1,517 | 4,586 |
| muscle 500k pairs, `-p 1` | 12.90 → 12.69 | 1.33 → 1.27 | +12 | 1,211 |
| muscle 500k pairs, `--recoverOrphans` | 58.67 → 1.46 | 1.92 → 1.78 | +198 | 4,662 |
| GEUVADIS, `--sketch` | 4.02 → 4.01 | 0.76 → 0.71 | 0 | 5,072 |
| GEUVADIS, decoy-aware GRCh38 index | 15.20 → 15.63 | 7.25 → 7.15 | +1,478 | 8,517 |
| GEUVADIS, `--seqBias --gcBias --posBias` | 24.27 → 20.70 | 1.49 → 1.37 | +1,517 | 223,463 |
| GEUVADIS read 1, single-end `-l A` | 6.92 → 6.28 | 0.86 → 0.89 | 0 | 0 |
| GEUVADIS, `--numBootstraps 30` | 22.90 → 21.20 | 1.52 → 1.47 | +1,517 | 4,586 |
| GEUVADIS, `--numGibbsSamples 100` | 40.56 → 17.38 | 1.11 → 0.97 | +1,517 | 4,586 |
| GEUVADIS `-a` (bowtie2 BAM) | see below | 0.45 → 0.43 | 0 | 0 |

- **`-a`** is identical in output. Its wall time is limited by reading the 2.3 GB BAM and varies by up to 2× between identical runs. With the BAM in memory, CPU time is 59–64 s for 2.11.0 and 55–58 s for 2.12.0 over four runs each.
- **The decoy-aware index is 2–3% slower; mapping time is unchanged.** About 0.15 s is the EM: the newly paired dovetails make its slowest connected component take 3,520 iterations instead of 3,248 (see [Costs](#costs)). About 0.09 s is reading the intermediate RAD, which has four times as many chunks now that reads are mapped in batches of 4,096 (#1346).
- **Peak memory** differences under 5% (single-end, bootstraps, decoy) are within run-to-run noise.
- **Reads moved** with all three bias flags come almost all from the `--posBias` effective-length change (#1352).

Accuracy against simulated truth, on TPM against true molar abundance, 2.11.0 → 2.12.0. Spearman and MARD are over transcripts with at least 100 true fragments:

| Simulation | Spearman | MARD |
|---|---|---|
| easy | 0.99262 → 0.99262 | 0.0764 → 0.0764 |
| hard | 0.98770 → 0.98784 | 0.0845 → 0.0844 |
| realistic GC, `--gcBias` | 0.95510 → 0.95510 | 0.2951 → 0.2952 |
| 3′ positional, `--posBias` | 0.87989 → 0.88468 | 0.4735 → 0.4682 |

Mapped counts on all four simulations are unchanged.

**Per change.** Each figure in the sections below compares the code before and after the change, on real reads unless marked. Unless noted, `quant.sf` is byte-identical before and after.

### `--recoverOrphans`

`--recoverOrphans` was 30–100× slower than the same run without it. For each orphan placement, the missing mate was searched for with a scalar DP over its whole window, with no early exit, no seeding and no reuse across a fragment's targets. That was 98% of the run's CPU. The search is now exact but cheap: a per-fragment memo of windows and window parts already searched, 7-mer seeding into diagonal bands, a floor score from the best diagonal, a pruned DP, and a reference-gap scan four columns at a time (SSE4.1 on x86, NEON on aarch64). Scoring schemes that break these arguments fall back to the plain DP.

Wall time with the flag, median of 3:

| data | `-p` | before | first commit of #1355 |
|---|---|---|---|
| muscle 500k pairs | 16 | 51.3 s | 1.81 s |
| muscle 500k pairs | 1 | 742 s | 17.5 s |
| muscle 4.4M pairs | 16 | 418 s | 10.9 s |
| GEUVADIS 10M pairs | 16 | 689 s | 15.7 s |
| GEUVADIS 1M pairs | 1 | 1125 s | 22.1 s |

Two more commits cut the overhead further. With the flag, a run now takes this much longer than the same run without it:

| data | `-p` | overhead |
|---|---|---|
| muscle 500k pairs | 16 | 1.17× |
| muscle 500k pairs | 1 | 1.21× |
| muscle 4.4M pairs | 16 | 1.20× |
| muscle 4.4M pairs | 1 | 1.20× |
| GEUVADIS 10M pairs | 16 | 1.54× |
| GEUVADIS 1M pairs | 1 | 1.65× |

`quant.sf`, `lib_format_counts.json` and the `num_*` counts in `meta_info.json` are byte-identical to the old search's on muscle and GEUVADIS, at `-p 1` and `-p 16`, with and without the flag. Without the flag nothing is slower, and memory is unchanged. (#1355)

### Inferential replicates

- **Gibbs rounds run in parallel** (#1347). On x86, on 10M GEUVADIS pairs with `--numGibbsSamples 100`, the posterior phase goes from 24.4 s to 7.5 s at `-p 16` (wall 40.7 s to 23.5 s), and to 4.6 s at `-p 32`. With no spare thread per chain, a round runs the same blocks in order: at `-p 4` the posterior phase takes 24.8 s, against 24.5 s before.
- **`bootstraps.gz` is compressed on every core** (#1343). Samples are compressed in parallel, in batches of 64, into one gzip member, so any gzip reader sees the same bytes as before. On x86, on 10M GEUVADIS pairs at `-p 16`, the output phase goes from 5.8 s to 0.3 s with `--numBootstraps 100` and from 6.8 s to 0.35 s with `--numGibbsSamples 100`; the whole run is 9% and 16% faster. Alignment mode with 100 bootstraps saves the same 5.5 s. Peak memory rises by at most 19 MB. The decompressed file is byte-identical; only the compressed bytes differ.

### Bias correction

- **`--gcBias`** (#1345): GC fractions of context and fragment windows come from precomputed tables. On x86, on 10M GEUVADIS pairs at `-p 16`, `--gcBias` goes from 20.34 s to 17.41 s wall (−14%), with the `em_bias` phase from 9.72 s to 6.75 s. With all three bias flags the run is 12% faster, and the stranded muscle sample with `--seqBias --gcBias` 13% faster.
- **`--seqBias`** (#1342): each sequence-bias factor is read from a table of all 4^9 contexts instead of being computed at every position. Every factor is bit-identical. On x86, on 10M GEUVADIS pairs at `-p 16`, `--seqBias` goes from 13.72 s to 12.89 s wall (−6%), and the `em_bias` phase from 3.62 s to 2.90 s (−20%). With all three bias flags the run is 2.4% faster, and on muscle with `--seqBias --gcBias` 3.4% faster. Peak memory does not change beyond noise.

### Mapping

- **FASTQ is read in batches of 4096 records** instead of 16384 (#1346). On x86 (AMD EPYC), on GEUVADIS (10M pairs, `-l IU`) and stranded muscle (`-l ISR`): peak memory falls 10–14% at `-p 8` and `-p 16` and 14–19% at `-p 32` (GEUVADIS `-p 32`: 1605 MB to 1303 MB), and does not change at `-p 1`. Wall time is 0.7% to 5.1% lower at every thread count.
- **Repetitive single-end reads** dedupe their candidates with a hash set above 112 candidates (#1330). On x86, on GEUVADIS read 1 as single-end (10M reads, `-l A`, GENCODE index), the run goes from 80.6 s to 72.9 s at `-p 1` (−9.6%) and from 7.22 s to 6.43 s at `-p 16` (−10.9%). Peak memory does not change. Paired-end mapping does not use this step.
- The mapper keeps one reverse-complement buffer per thread instead of allocating one per read (#1345).
- The rough EM that seeds the bias models no longer copies the per-mapping weights only Gibbs reads (#1329). Peak memory on GEUVADIS (10M pairs, `-p 16`) does not move beyond run-to-run noise.

### Alignment mode

- BAM records no longer build a read name and CIGAR op list that nothing reads (#1340). On x86, on a 2.3 GB GEUVADIS BAM (10M pairs, multi-mapping, `AS` tags), the run is 3.1% faster at `-p 1` (40.3 s to 39.0 s), and CPU time falls about 4% at every thread count. At higher thread counts the single BAM reader sets the pace, so wall time does not move.
- A run opens its input once instead of two or three times (#1351). At `-p 16` on regular files: GEUVADIS BAM 11.80 → 10.54 s (median of 6), muscle BAM 4.58 → 4.38 s, SAM 1.79 → 1.81 s.

### Index

`salmon index` no longer writes piscem's "Tiny" dictionary files (`index.tdct`, `index.tct`), which salmon never reads (#1344). piscem wrote them only for references below 100M distinct k-mers, so human transcriptomes (Ensembl GRCh38 cDNA, GENCODE v50) are not affected. On x86, on yeast cDNA, the index goes from 216 MB to 15 MB on disk, the build from 0.58 s to 0.42 s, and peak memory from 691 MB to 423 MB. On the authors' synthetic 78 MB reference (Apple M4 Max), the index goes from 908 MB to 104 MB, the build from 1.711 s to 0.939 s, and peak memory from 1.8 GiB to 0.9 GiB. Every other index file is byte-identical, and an index built by an earlier version still loads.

### Costs

Pairing small dovetails (#1357) costs, at `-p 16`: GEUVADIS 10M pairs 9.29 → 9.42 s (+1.4%, all in the EM, whose slowest component needs 7,650 iterations instead of 6,968; mapping is flat), and full muscle 8.11 → 8.16 s (+0.7%, mapping the 42,902 new pairs). On the decoy-aware GRCh38 index, the EM's slowest component needs 3,520 iterations instead of 3,248, which costs about 0.15 s.

## Fixes

### Mapping

The mapping fixes that change counts are listed under [Changes to results](#changes-to-results): the 5′-end fragment rule and inward-only `--allowDovetail` (#1356), the sketch fragment length and start (#1360), the negative-start clip (#1364), case-insensitive orphan rescue (#1332), the rescued-pair span (#1322), and the single-end filter order (#1333).

### Alignment mode (`-a`) and genome projection

- **A late failure after a `--skipQuant` genome projection keeps `genome_projected.rad`**, the run's only result. The error path deleted it. (#1317)
- **A pipe as `-a` input** no longer hangs (#1339), and is now read (#1351): see [New features](#new-features).
- The fixes for a partial `-t` (#1319), `-l IU`/`-l U` bias seeding (#1325) and poly-A-tailed FASTA (#1338) are under [Changes to results](#other-options-and-modes).

### Quantification and inference

- **`--skipQuant` reports the full run's mapping rate** under a stranded `-l` (#1334).
- **Single-end `--gcBias`** is skipped with a warning (#1315).
- **Gibbs uses the optimizer's prior** in every driver (#1316).

Details are under [Changes to results](#other-options-and-modes).

### CLI and output

- **`quant.genes.sf`, `--writeMappings` and `--writeBam` are written through a temporary name and a rename**, as `quant.sf` and RAD files already were. On an error the BAM writer still wrote the BGZF EOF block, so `samtools quickcheck` passed on a truncated BAM. Now a failed run leaves no BAM, SAM or gene table that looks complete. A target that exists and is not a regular file (a symlink such as `/dev/stdout`, a FIFO or a device) is written in place, so writing mappings to stdout or a pipe works as before. (#1324)
- **Each run's intermediate RAD has a name of its own.** It was `{stem}-{pid}.rad`, and pids repeat across containers and nodes. Two runs sharing a `--radScratchDir` on network storage could write into the same file. The name now adds the start time in nanoseconds. (#1318)
- **`meta_info.json` writes `samp_type: "none"` when no posterior samples are drawn**, for example under `--skipQuant` with `--numBootstraps 30`. It used to say `"bootstrap"` beside `num_bootstraps: 0`. tximport reads these fields. (#1335)
- **`salmon debug-map` refuses a truncated or corrupt FASTQ** and names the file, as `quant` does. It used to print a partial table and exit 0. Intact files give the same output, and CRLF input and a trailing blank line still work. (#1336)
- **`--skipQuant -g` in default reads mode logs that gene-level aggregation is skipped**, as the other drivers do. (#1337)

### Index

- **A failed `salmon index` removes `cleaned_refs.fa`.** On a preprocessing failure only `refseq.bin.partial` was removed, so a cleaned copy of the input (about 3 GB for a genome-decoy index) stayed behind. (#1331)
- **`salmon index -t idx/cleaned_refs.fa -i idx` is refused** before anything is written. It used to empty its own input before reading it. (#1331)

## Dependencies

- bramble-rs 0.2.2 (from 0.2.0): when several alignments of a fragment sit at the same positions, mates are paired by `HI`. (#1359)
- piscem-rs 0.13.0 retains the paired-end template-length and `SimpleHit::pair_start()` fixes from 0.12.4 (#1360), uses an arena for reference names, and improves strict contig walking. The name arena saves about 10.48 MB on the reviewed full GENCODE index; no general speedup is claimed for Salmon's default permissive mapping. Existing indexes remain readable. Rust callers of `RefInfo::names()` now receive an iterator of `&str` instead of `&[String]`.
- sshash-lib and tiny-dict 0.8.2: lazy streaming minimizers, widened bucket starts, and the empty-skew-partition correctness fix, with existing index formats preserved.
- cf1-rs 0.6.1: batched small-reference output, reuse of walked unitigs, four-base packing, and MPHF construction selected using the configured soft memory budget. The budget is not an RSS limit; large inputs may use slower streaming construction to reduce memory. Output formats and unitig IDs are unchanged from 0.6.0.
- ksw2rs 0.3.1 (from 0.2.1): improved exact-score/CIGAR NEON kernels and the fix for the 0.3.0 `KSW_EZ_NO_MAX` slowdown. It also corrects banded traceback and scalar scoring, so affected CIGARs or scalar results can change. Default SIMD score-only results are unchanged. Salmon's usual short score-only extensions are approximately at parity in the upstream measurements; these kernel results are not an end-to-end Salmon speedup claim. Genome projection through bramble-rs 0.2.2 still uses that crate's ksw2rs 0.2.1 dependency.
- Upstream paraseq 0.6.3 replaces paraseq-temp 0.6.3-pre.1. It includes the fork's named input errors and truncated-FASTQ rejection, and retains resizable worker pools. Downstream libraries sharing parser types with piscem must also use upstream paraseq.
- libradicl 0.21.2. (#1312)
- Patch updates to `flate2` 1.1.11, `zstd` 0.14.4, `jiff` 0.2.45 and `sux` 0.15.4. Index files are byte-identical with them. (#1365)
- Website: Starlight 0.42.9, Astro 7.4.7. (#1310)

## Documentation and internals

- `cargo doc` with `-D warnings` runs in CI, after six broken intra-doc links were fixed. (#1350)
- The determinism guide, CLI reference, output formats, introduction and migration guides no longer describe `--online`. (#1361)
- Helpers that had two to four copies across crates now have one, in `salmon-core`. (#1327)
- Dead and test-only code is removed. (#1328)
- A test pins orphan compatibility for every library type, mate and strand. (#1323)

## Known limitations and follow-ups

- **`--posBias` reshapes effective lengths only within a fragment-length class.** Correcting isoform assignment under positional bias needs a per-fragment positional weight, which is a separate change. (#1352)
- **`--errorModel` cannot read `-a` from a pipe.** Supporting it would mean spooling the whole input. (#1351)
- **Under `-l A`, `--skipQuant` still counts every placed fragment as mapped**, since the mapping pass has no expected format yet. (#1334)
- **A `salmon index` failure after preprocessing** (in the dBG or SSHash build) still leaves the cleaned FASTA and the cDBG files. (#1331)
- **The intermediate RAD is not opened with `create_new`**, so a name clash would not fail loudly. Two runs would have to start in the same nanosecond with the same pid. (#1318)
- **Orphan rescue counts `N` against `N` as a match.** (#1332)

## Checked

- `--recoverOrphans`: `quant.sf`, `lib_format_counts.json` and the `num_*` counts in `meta_info.json` are byte-identical to the old search's on muscle and GEUVADIS at `-p 1` and `-p 16`. `--writeMappings` SAM output is identical on 50k pairs. A property test compares the fast window search with the reference DP on about 10k random cases over 7 scoring schemes, and 20,796 muscle and 313,774 GEUVADIS real search inputs were checked the same way. (#1355)
- The `--online` removal leaves the default `quant.sf` byte-identical to the code before it for reads (plain, bias, sketch, single-end, Gibbs including `bootstraps.gz`, and the high-multimapping set) and for `-a` (`AS`, `--errorModel -t`, and bias). `--skipQuant` metadata is identical, and its RAD requantifies to the identical `quant.sf`. (#1361)
- Pipe input: on a 10M-pair GEUVADIS BAM (`-l IU` and `-l A`) and a muscle SAM and BAM, regular-file `quant.sf` is byte-identical to the code before the change, and the same input through `<(cat ...)`, `/dev/stdin` and a FIFO is byte-identical to the file. (#1351)
- Small dovetails: `--writeMappings` output (9.5M records from the simulation, 17.8M from 500k muscle pairs) has no CIGAR/SEQ length mismatches and no mate-field inconsistencies. `samtools quickcheck` passes, and `sort` and `index` succeed. (#1357)
- Negative-start clip: on the hard simulation, the `--writeRad` run equals the default run, and requantifying with `--rad` gives the same `quant.sf`. `--fldPolicy derive` re-derives a byte-identical `fld.gz` and `quant.sf`. (#1364)
- Genome projection: `quant.sf` is byte-identical when bramble's output order is shuffled, on the easy, hard and multimapping simulations, ERR188044 and muscle. (#1359)
- The speed-ups in #1329, #1330, #1340, #1342, #1343, #1344, #1345 and #1346 leave `quant.sf` byte-identical. The refactors in #1327 and #1328 leave `quant.sf`, the bootstrap and Gibbs draws and the FLD byte-identical on GEUVADIS and on the sample BAM in alignment mode.

## Acknowledgements

Thanks to Benjamin Demaille for:
- the audit of 2.11.0 (#1314) and its fixes and speed-ups (#1315–#1320, #1322–#1340);
- the performance and peak-memory pass (#1341) and its changes (#1342–#1347);
- finding the one-base sketch fragment length (#1321); his test is part of #1360.
