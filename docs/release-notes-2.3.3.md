# salmon 2.3.3

A features-and-fixes release on 2.3.2. **The default `salmon quant` output is
unchanged** — the new inference accelerators are opt-in (default off), the
profiling additions are pure instrumentation, and the one bug fix is
behavior-neutral. No index rebuild is required.

## Optional EM/VBEM acceleration (`--emAccel`)

salmon's abundance estimation is a fixed-point EM (or VBEM under `--useVBOpt`)
that converges only linearly; on high-ambiguity / slow-mixing problems it can
crawl, and every bootstrap replicate pays the cost again. Two off-the-shelf
accelerators are now available (thanks
[@BenjaminDEMAILLE](https://github.com/BenjaminDEMAILLE), #1051, #1052):

- **`--emAccel squarem`** — SQUAREM (Varadhan & Roland 2008): a squared
  extrapolation from the last two EM maps plus a stabilizing step. Near-zero
  overhead; a small, safe speedup on typical runs and a large one on
  slow-mixing problems.
- **`--emAccel daarem`** — DAAREM (Henderson & Varadhan 2019): damped Anderson
  acceleration over a window of residuals (a multi-secant quasi-Newton step),
  with a self-contained eigensolver (no new dependency). Wins big on hard,
  ill-conditioned problems; because of its per-iteration overhead it is best
  reserved for those cases (it can be slower than plain EM when the base EM
  already converges quickly).

Both reach the **same fixed point** as plain EM/VBEM (validated: TPM Pearson
≈ 0.9998–0.99999 vs. the default, total mass conserved, identical accuracy vs.
simulation ground truth), and both accelerate the point estimate **and** every
bootstrap replicate (Gibbs is unaffected). The **default is `none`**, so
existing results and `--deterministic` output are unchanged; the accelerators
are not byte-identical to `none` (differences are at the convergence-tolerance
level), which is why they are opt-in.

## Profiling harness + per-phase timing

- A reproducible profiling harness (#1051, @BenjaminDEMAILLE): a `profiling`
  build profile, a `criterion` EM/VBEM M-step benchmark, and `scripts/profile.sh`
  (samply).
- Per-phase wall-clock timing on the `salmon::timing` tracing target
  (index_load / mapping / eff_length_collapse / em_bias / posterior / output),
  available inline at `info` or in isolation via `RUST_LOG=salmon::timing=info`.
  It now covers **all** quantification paths — reads (online and
  `--deterministic`), alignment (`-a`), and `--rad` — so the EM phase is visible
  regardless of input mode. Pure instrumentation; no effect on output or
  determinism.

## Fixes

- **Radix-sort scratch fill** (#1055): the `(transcript, orientation)` uni-MEM
  radix grouping added in 2.3.2 constructed a zero-length `Mem` as a
  never-read scratch placeholder, which tripped a debug assertion and **panicked
  debug builds** on real reads. Release builds (all shipped binaries) were
  unaffected — the assertion is compiled out and the placeholder was always
  overwritten — so this changes no results (verified byte-identical `quant.sf`).
  The fix restores usable debug builds and closes the invariant hole.

## Other

- Tidier `salmon --help`: the removed `alevin` subcommand is hidden and the
  `(Rust port)` parenthetical is dropped (#1053, @BenjaminDEMAILLE).

## Upgrading

Drop-in from 2.3.2: no index rebuild, and default quantification results are
unchanged. Try `--emAccel squarem` (or `daarem` for slow-mixing / large
bootstrap workloads) if you want faster inference.
