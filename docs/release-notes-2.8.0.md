# salmon 2.8.0

A compatibility and correctness release. salmon now reads the versioned
**RAD v2** format that current piscem writes, and the opt-in DAAREM EM
accelerator is reworked: it is now correct and fast with plain EM, and it is
**no longer accepted with VBEM** (the default optimizer). Nothing else about
quantification changes, and **no index rebuild is needed**: the index format
is unchanged from 2.7.0.

## Results: unchanged by default

Against the released 2.7.0, on `sample_data` (paired-end, `-l A`), `quant.sf`
is **byte-identical** in every mode checked:

- default (VBEM)
- `--useEM`
- `--sketch`
- `--emAccel squarem`
- `--useEM --emAccel daarem`
- `--numBootstraps 5`

Indexes built by 2.7.0 and 2.8.0 are byte-identical at `-p 1`.

## Reading RAD v2 (piscem ≥ 0.24)

piscem now writes RAD files in the versioned **v2** format: a `RAD_FILE`
magic number, a spec version, a prelude extension block, and a semantic
*role* on every tag. salmon moves to **libradicl 0.21**, which detects the v2
prelude and reads bulk records by their declared roles. Quantifying piscem's
bulk RAD output (`salmon quant` on a RAD input) works with both the legacy
format and v2.

salmon's own RAD writer is unchanged and still writes the legacy format, so
its output is byte-identical to 2.7.0's. Emitting v2 from salmon is future
work.

## `--emAccel daarem`: reworked, and plain EM only

**Behavior change.** `--emAccel daarem` now requires plain EM (`--useEM` or
`--meta`). With VBEM, the default optimizer, salmon exits with an error that
suggests `--emAccel none` or `squarem`.

Earlier releases could return **inaccurate abundances** with `daarem`:

- **Negative intermediate abundances.** The Anderson extrapolation
  overshot below zero and the negative iterate was kept. On a GENCODE run
  (646k transcripts) this produced a false "converged" far from the answer
  (Spearman 0.90 against a fully converged EM), with some reads lost to
  truncation.
- **A convergence cycle under VBEM.** The acceptance test judged steps by
  the residual ‖F(x) − x‖. VBEM's sparse prior drives near-identical
  isoforms to winner-take-all solutions, and the residual *grows* while the
  losing isoform drains. So the test refused VBEM's own path to the optimum:
  one run repeated the same state every 9,750 iterations until `maxIter`.

DAAREM now follows the objective-based variant of the reference method (the
CRAN `daarem` package's `daarem_base_objfn`):

- A step is accepted unless it lowers the log-likelihood by more than 0.01.
- Extrapolated abundances are floored so they stay positive.
- It stops only when both successive iterates and one plain EM step agree to
  within `rel_diff_tol`.
- Its vector work runs in parallel, and results are identical at any thread
  count.

With plain EM on three GENCODE v50 datasets (646k transcripts, 64 threads):

| dataset | DAAREM, default tolerance | DAAREM, tolerance 1e-3 | SQUAREM, tolerance 1e-4 |
|---|---|---|---|
| SRR21186103 | 2.3 s | 5.6 s (0.11 nats from optimum) | 231 s |
| cerebellum (ERR13232583) | 2.0 s | 3.3 s (0.011) | 153 s |
| testis (SRR11517405) | 1.3 s | 4.6 s (0.053) | 313 s |

Under VBEM, no DAAREM variant was an accelerator. It was slower than plain
VBEM at the default tolerance, and did not converge within 100k iterations at
`1e-4` on two of the three datasets. VBEM's objective has many optima, and
the extrapolation kept crossing between them. With VBEM, use no acceleration
(the default) or `--emAccel squarem`, which pays off at tighter tolerances.

## Library changes (salmon-infer)

- **`EmOptions::validate`** is new. It rejects DAAREM with VBEM, and the
  optimizer enforces it too.
- **Target-less classes are skipped.** An equivalence class with no targets
  is now ignored instead of causing an out-of-bounds panic. salmon never
  produces one, but external callers that build `PackedEqClasses` directly
  can (piscem-infer did, for fragments that failed its strand filter).
- **SQUAREM docs corrected.** They no longer promise the same fixed point
  as plain iteration under VBEM.

## Dependencies

- **libradicl 0.18 → 0.21.** libradicl no longer depends on noodles.
- **noodles unified on 0.116:** noodles-bam 0.95, noodles-sam 0.90, and
  bramble-rs 0.1.9. The build now carries exactly one version of every
  noodles crate. This supersedes #1193 and #1194.
- **Other updates:** zstd 0.14, and patch updates to crossbeam, clap and jiff.

## Also in this release

- **Allocator code moved.** The global allocator now lives in its own module
  (`salmon-cli/src/global_alloc.rs`), and the binary logs which allocator it
  was built with at debug level. `docs/allocator-choice.md` records what
  would be needed to justify changing the default (mimalloc). No behavior
  change. Thanks to @BenjaminDEMAILLE.
- **Website.** The documentation site builds again, and its dependencies were
  updated (Astro 7.3, Starlight 0.42).
- **CI.** GitHub Actions dependencies were updated.
