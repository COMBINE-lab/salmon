# salmon 2.12.0

> **Draft.** This file collects the 2.12.0 changes as they merge; it is
> finalized (and this note dropped) at release.

## The one-pass `--online` path is removed

`--online` selected the pre-2.6 one-pass quantification path. It was deprecated
in 2.6.0, when the deterministic two-phase path became the default (#1140), and
it is now removed. Every `salmon quant` maps to an intermediate RAD and then
quantifies that RAD with a fixed fragment-length distribution, so `quant.sf` is
byte-identical at any `-p`. The removal does not change default results.

What happens to existing command lines:

- **`--online` is an error.** The message says to drop the flag. It names
  `--errorModel -t <FASTA>` for error-model scoring in alignment mode and
  `--radScratchDir <dir>` for placing the intermediate RAD (for example on
  `/dev/shm` or node-local scratch).
- **`-f/--forgettingFactor`, `--numAuxModelSamples` and
  `--numPreAuxModelSamples`** tuned only the online phase. They are still
  accepted, so old scripts and C++-era command lines keep running, but they
  print a warning, do nothing, and are hidden from `--help`.
- **`--deterministic`** (already the default) and **`--initUniform`** (every EM
  already starts uniform) are still accepted and do nothing.

If you are moving from `--online` to the default, these differences may show up:

- **Alignment mode (`-a`) scores placements by the BAM `AS` tag.**
  - `--online -a -t <FASTA>` trained an alignment error model by default.
  - To get error-model scoring now, add `--errorModel -t <FASTA>`. It trains an
    order-independent counting model in a second pass over the BAM.
  - Every benchmark with ground truth found `AS` scoring at least as accurate,
    in a single pass over the BAM.
- **Piped alignment input with `--errorModel` is not supported.**
  - `--errorModel` reads the BAM twice, so the BAM must be a regular file.
  - Write the BAM to disk first, or keep the default `AS` scoring, which needs
    only one pass.
  - `--online -a` was the only single-pass error-model path, and it is gone.
- **`-l A` is detected from the whole run.**
  - The default detects the library type from every fragment, then filters
    against the result when it quantifies.
  - `--online` locked the type in after a prefix of about 50k fragments and
    never judged the fragments in that prefix.
  - On small inputs, or inputs with mixed orientations, a few more fragments
    can be counted as incompatible. `num_compatible_fragments` covers every
    fragment that was judged.
- **Single-end reads without `--fldMean`/`--fldSD` use a widened default
  prior.**
  - The single-end fragment-length term uses the default prior with its SD
    widened, so it stays robust when the default does not fit the library.
  - If you know the library's fragment-length distribution, pass
    `--fldMean`/`--fldSD` and it is used exactly as given.
  - On a simulation whose true distribution equals the default prior, the
    widening costs a little accuracy: Spearman 0.817 against 0.827 for
    `--online`. Passing the true mean and SD matches `--online` (0.827).

Paired-end accuracy against simulated truth is unchanged from `--online`, to
the fourth decimal and in no consistent direction. This holds for selective
alignment, `--sketch`, a high-multimapping set and `--seqBias --gcBias
--posBias`. It also holds on simulations with real GC and 3′ positional bias:
with the matching correction, the default and `--online` agree to within 0.001
Spearman. The default uses less CPU in every measured cell and less wall
time in all but one: with reads-mode bias correction it takes about 5% more
wall time at `-p 16`, which is the rough seed EM it runs before baking the RAD.
Measured on 4.65M simulated pairs at `-p 16`.

### Library API

- **`salmon_quant::quantify` is now `salmon_quant::map_to_rad(opts, rad_path)`.**
  - It runs only the mapping pass: it writes the RAD and the mapping-pass
    metadata, and estimates **no abundances**.
  - To quantify, call `salmon_align::quantify_rad` on the RAD it wrote. That is
    what the `salmon` binary does.
- **`QuantResult` is now `MapResult`**, without the abundance, TPM, posterior
  and bias fields.
- **`QuantOptions` no longer has** `deterministic_fld`, `forgetting_factor`,
  `num_aux_model_samples`, `num_pre_aux_model_samples` or `init_uniform`.
  `write_rad`, `skip_quant` and the eq-class dump flags remain; they are read
  by the caller's two-phase driver, not by `map_to_rad`.
- **Removed from the other crates:**
  - `salmon_align::quantify_alignments`, which was the one-pass alignment-mode
    path. Use `write_alignment_rad` + `quantify_rad`.
  - `AlignQuantOptions.{forgetting_factor, init_uniform,
    num_aux_model_samples, num_pre_aux_model_samples}`.
  - `salmon_infer::OnlineInference`.
  - `FragmentLengthDistribution::{refresh_online, online_snapshot,
    online_cmf_snapshot}`.
