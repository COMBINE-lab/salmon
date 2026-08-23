# The EM/VBEM phase: architecture, optimizations to date, and where to look next

Written as a briefing for someone evaluating further optimization of the
inference phase. It describes what is there now, what has already been tried
and measured, and which specific things are known to still be on the table —
so that effort goes somewhere new rather than re-deriving what is here.

All line references are against `develop` at the time of writing plus
`fix/em-determinism-and-scaling` (#1167); check them, they drift.

---

## 1. Why this phase matters

On a 26.1M-fragment sketch run against GENCODE v49 (238,442 transcripts), the
wall time splits roughly:

| phase | time | share |
|---|---|---|
| phase 1 — mapping (includes writing the intermediate RAD) | 30.5 s | 45% |
| phase 2 — **EM inference** | 31.9 s → **15.2 s after #1167** | 47% → ~31% |
| phase 2 — RAD read | 3.9 s | 6% |
| index load, eff-length collapse, posterior, output | ~1.5 s | 2% |

Before #1167 the EM was as expensive as all of mapping. It is now roughly half
that, and it is still the second-largest cost and the one that scales worst.

Measured on `newton` (Xeon E5-2699 v4, 88 cores); see §7 for the data.

---

## 2. Where the EM sits

Two entry points, both in `crates/salmon-infer/src/lib.rs`:

- `optimize(...)` — builds a `PackedEqClasses` from a `CollapsedEqClasses` and
  runs to convergence.
- `optimize_packed_with_init(p, opts, parallel, init_alphas, eff_lens)` — the
  one production actually calls, with `parallel = true` everywhere (reads mode
  `salmon-quant/src/lib.rs:1039,1041,1249`; RAD/alignment
  `salmon-align/src/rad.rs:2214,2216,2294` and `salmon-align/src/lib.rs:2461,2463,2500`).

Both funnel into **`run_em_counts(p, counts, opts, parallel, min_iter, init_alphas, eff_lens)`**
(`lib.rs:279`), which owns the buffers, builds the shard plan, selects the
M-step kernel, and drives one of three convergence loops.

Note the argument order: `init_alphas` precedes `eff_lens` and both are
`Option<&[f64]>`, so swapping them type-checks silently. That has bitten us
once already (the bootstrap path passed effective lengths as the initial
abundance vector; `uncertainty.rs:136` now has a comment about it).

### Callers with different needs

- **Point estimate** — `parallel = true`, full iteration budget.
- **Bias seed EM** (`salmon-quant/src/lib.rs:1039`) — a deliberately
  under-converged short run (`min_iter = max_iter = bias_seed_em_iters`,
  `min_alpha = 0.0`) whose output seeds the bias model.
- **Bootstrap replicates** (`uncertainty.rs:136`) — `parallel = false`, because
  the *replicates* are the parallel dimension; each replicate runs a sequential
  EM inside a `par_iter` over replicates. Optimizing the sequential kernel
  therefore matters as much as the parallel one for `--numBootstraps` runs.
- **Gibbs** (`uncertainty.rs`) — separate sampler, shares `PackedEqClasses`
  and the `weights` array rather than `combined`.

---

## 3. The data structure: `PackedEqClasses` (CSR)

`crates/salmon-infer/src/packed.rs`. Flat CSR, built once per run:

```
labels:   Vec<u32>   flat transcript ids; class i spans labels[starts[i]..starts[i+1]]
starts:   Vec<u64>   CSR offsets, len = num_classes + 1
combined: Vec<f64>   weight/effLen per (class, transcript) incidence — EM reads this
weights:  Vec<f64>   raw conditional weights, same layout — Gibbs reads this
counts:   Vec<u64>   per-class fragment count
num_txps: usize
```

`starts` is `u64` deliberately: these are cumulative *incidence* counts, and a
`u32` wrapped silently past 4G incidences and handed the optimizer slices
belonging to other classes (#1097). Costs ~3.5% of the structure.

Scale on the reference dataset: **536,311 classes, 238,442 transcripts**, mean
class size ~2.5, so ~1.3M incidences. Each M-step reads every incidence twice
(denominator, then distribution) — call it ~2.6M f64 reads of `combined` plus
the same number of scattered reads of `alpha_in` and scattered writes to an
accumulator.

`weights` is populated only when a posterior method needs it
(`PosteriorMethod::needs_raw_weights`); a point-estimate run does not pay for
it.

---

## 4. The M-step kernels

Four, dispatched by `(use_vbem, parallel)` in the closure at `lib.rs:~300`:

| | sequential | parallel |
|---|---|---|
| EM | `em_step_seq` | `em_step_par` |
| VBEM | `vbem_step_seq` | `vbem_step_par` |

**Plain EM**, per class `ci` with `count` fragments and members `(tid, w)`:

```
denom = Σ alpha_in[tid] * w
if denom > MIN_EQ_CLASS_WEIGHT:
    inv = count / denom
    for each (tid, w): acc[tid] += alpha_in[tid] * w * inv
```

Single-member classes short-circuit to `acc[tid] += count`. The inner product
is recomputed rather than cached in scratch — the multiply was measured cheaper
than the extra memory traffic inside the parallel closure.

**VBEM** substitutes `exp_theta` for `alpha`:
`exp_theta[i] = exp(digamma(alpha_i + prior_i) − digamma(Σ_j alpha_j + prior_j))`,
computed by `fill_exp_theta` **once per M-step, sequentially**, because it
depends on a global sum over all transcripts. That is an `O(num_txps)` serial
prologue on every iteration, including `num_txps` `digamma` calls — worth
measuring, see §6.

Priors come from `prior_alphas_vec` (`lib.rs:241`): flat `vb_prior` per
transcript, or `vb_prior * max(1, effLen)` under `--perNucleotidePrior`.

**VBEM is the production default** (`opt_type: "vb"`; the CLI sets
`use_vbem = !args.use_em && !args.meta`). `EmOptions::default()` has
`use_vbem: false`, so any benchmark or test built on the default is exercising
the *non*-default kernel. This mattered: the original determinism test did
exactly that.

---

## 5. Parallelism: `ShardPlan` (this is what #1167 changed)

### What it was

```rust
let nshards = rayon::current_num_threads().clamp(1, 64);
// per shard, per iteration:
buf.iter_mut().for_each(|x| *x = 0.0);      // dense clear,  O(num_txps)
...
reduce_shards(shards, alpha_out);            // dense reduce, O(nshards * num_txps)
```

Each shard owned a private dense `num_txps` accumulator and processed a
contiguous slice of classes with non-atomic adds — chosen over a shared
`AtomicF64` array, which contended badly on hot transcripts.

**Two defects.**

1. **Determinism.** `nshards` tracked the thread count, so `-p` changed the
   class partition and therefore the grouping of floating-point sums. On real
   data `-p 4` and `-p 32` produced different `quant.sf`, *at default
   precision*. Confirmed by construction: `-p 64` and `-p 80` (both clamped to
   64) were byte-identical; `-p 16` differed.
2. **Scaling.** Clear + reduce cost `O(nshards × num_txps)` per iteration
   regardless of the data — ~30M ops at 64 shards against ~1.3M incidences of
   real work — and grew with thread count while the work stayed fixed.

### What it is now

`ShardPlan` (`packed.rs`, ~line 400):

- **`nshards` derives from the data alone**:
  `nclasses.div_ceil(MIN_CLASSES_PER_SHARD=4096).clamp(1, MAX_SHARDS=64)`.
  On the reference dataset that is 64 shards regardless of `-p`. Rayon
  schedules them over whatever threads exist — a scheduling decision, which
  cannot change arithmetic.
- **`touched[s]`** — sorted, deduplicated transcript ids reachable from shard
  `s`'s class range, computed **once** (labels never change between
  iterations). The per-iteration clear touches only these.
- **`contrib_start` / `contrib`** — CSR inversion of `touched`
  (transcript → contributing shards, ascending). The reduction is
  `par_iter` over transcripts, each an independent fixed-order sum over only
  the shards that contribute. Parallel *and* deterministic.

Per-iteration overhead is now `O(Σ|touched|)` ≈ `O(incidences)` rather than
`O(nshards × num_txps)`.

### Determinism invariants (do not break these)

Any future change must preserve all four, or `quant.sf` stops being
reproducible across `-p`:

1. Shard count and boundaries depend only on the data.
2. Each shard's accumulation order over its classes is fixed.
3. The reduction sums a transcript's contributors in a fixed (ascending shard)
   order.
4. SQUAREM/DAAREM vector math is sequential (see §6) — if parallelized, it must
   use a fixed-order reduction, not a nondeterministic one.

Enforced by `packed.rs::shard_plan_determinism::em_result_is_independent_of_thread_count`,
which compares bit-for-bit across rayon pools of 1/2/3/8/16 threads over **all
four** configurations (plain EM, VBEM, plain+SQUAREM, VBEM+SQUAREM). It is
verified to *fail* on the pre-#1167 tree (1809 of 5000 transcripts differ
between 1 and 2 threads). Weights in that fixture are deliberately not
exactly representable — with tidy powers of two every grouping agrees and the
test proves nothing.

---

## 6. Convergence and acceleration

### The loop

`max_rel_diff(alpha_in, alpha_out, cutoff)` (`lib.rs:~250`) — max over
transcripts of `|out−in|/out`, considering only transcripts with
`alpha_in > alpha_check_cutoff` (without the cutoff a transcript oscillating
between 1e-30 and 2e-30 blocks convergence forever). **This is a sequential
loop over `num_txps`, run every iteration.**

Buffers are swapped, not copied. `min_iter` floors the convergence check;
`max_iter` is the hard stop.

### Acceleration (`EmAccel`)

- **`None`** (default) — plain fixed-point iteration, one M-step per iteration.
  Output unchanged from historical salmon.
- **`Squarem`** (SqS3) — per cycle: two M-steps `x1=F(x0)`, `x2=F(x1)`, form
  `r = x1−x0` and `v = (x2−x1)−r`, extrapolate
  `xp = x0 − 2αr + α²v` with `α = −‖r‖/‖v‖` clamped to `α ≤ −1`, halving back
  toward −1 if any component goes negative; then a stabilizing third M-step
  `x_next = F(xp)`. Three M-steps per cycle, far fewer cycles. Not
  byte-identical to `None` (different iterate sequence, same fixpoint within
  `rel_diff_tol`).
- **`Daarem`** — damped Anderson acceleration with restarts and
  residual-monotonicity control over a window of recent iterates
  (`daarem.rs`, 539 lines). Better than SQUAREM on ill-conditioned,
  high-dimensional problems.

**The vector/norm math in SQUAREM and DAAREM is sequential over `num_txps`,
deliberately**, so results do not depend on thread count. The comment says this
is "cheap next to a per-class M-step" — that was true when the M-step cost
~30 s. After #1167 the M-step is ~2× cheaper, so the ratio has moved and this
should be re-measured before it is assumed still true.

### Finalization

`finalize_truncate_redistribute` (`lib.rs:~265`) — applies `min_alpha`
truncation as a *mass-preserving redistribution* (a masked final M-step over
the surviving members), not a rescale. Mass that cannot be redistributed
(fully-truncated classes) is reported as `inference_truncated_mass`. A
non-positive `min_alpha` makes it a no-op (used by the bias warm-up).

---

## 7. What has been done, with measurements

| change | effect | where |
|---|---|---|
| CSR `PackedEqClasses` replacing per-class `Vec`s | (pre-existing) | `packed.rs` |
| Per-shard dense accumulators replacing a shared `AtomicF64` array | removed CAS contention on hot transcripts | `packed.rs` |
| `u64` CSR offsets | correctness past 4G incidences (#1097) | `packed.rs` |
| SQUAREM, DAAREM | opt-in, far fewer M-steps | `lib.rs`, `daarem.rs` |
| **Data-derived shard count** (#1167) | **fixed cross-`-p` reproducibility** | `packed.rs` |
| **Sparse clear via `touched`** (#1167) | see below | `packed.rs` |
| **Parallel sparse reduce via CSR inversion** (#1167) | see below | `packed.rs` |

EM phase time, newton, same binary pair, back to back, sketch, 26.1M fragments:

| `-p` | before #1167 | after #1167 | |
|---|---|---|---|
| 16 | 32.58 s | **15.19 s** | 2.14× |
| 32 | 23.69 s | **13.95 s** | 1.70× |

Mapping and RAD-read unchanged, so the effect is isolated to this phase.

End-to-end wall on an EPYC 9575F (10M-pair equivalent workload, sketch):
`-p 4` 1.03×, `-p 16` 1.06×, `-p 32` 1.14×, `-p 64` 1.17×.

---

## 8. Known remaining bottlenecks — the actual open questions

**The EM still does not scale.** 15.19 s at `-p 16` vs 13.95 s at `-p 32` is
9% for double the threads. #1167 removed a large constant factor, not the
ceiling. Candidates, roughly in order of my confidence:

1. **Dense per-shard accumulator memory traffic.** 64 shards × 238,442 × 8 B ≈
   **122 MB** of accumulator working set, against ~55 MB of L3 on the Xeon and
   far more on the EPYC — which is consistent with the fix helping more on the
   EPYC than on newton. We now compute `touched[s]` anyway, so each shard's
   accumulator could be **index-compressed** to its own touched set
   (`buf[local_idx]` with a `tid → local_idx` map, or a sorted dense array
   indexed by rank within `touched[s]`). That shrinks the working set to
   `O(Σ|touched|)` ≈ 10 MB and makes the adds cache-local. This is the single
   most promising item and it preserves all four determinism invariants.
2. **Sequential `max_rel_diff` every iteration** — `O(num_txps)` serial, ~238k
   per iteration. Parallelizable with a fixed-order max-reduction (max is
   associative and exactly representable, so a tree reduction is deterministic).
3. **VBEM's `fill_exp_theta` prologue** — `O(num_txps)` serial including
   `num_txps` `digamma` calls, every M-step, on the *production default* path.
   The global sum must be sequential (or a fixed-order tree reduce), but the
   per-element `exp(digamma(...) − log_norm)` map is embarrassingly parallel.
4. **SQUAREM/DAAREM vector math**, sequential over `num_txps` per cycle. Same
   fixed-order-reduction caveat.
5. **Shard count is capped at 64** and sized by classes, not by incidences.
   Classes vary in size, so contiguous equal-count ranges are not equal-*work*
   ranges; a work-balanced partition (equal incidence counts, still data-derived)
   would reduce tail effects. Must remain independent of `-p`.
6. **`em_step_seq` matters more than it looks** — every bootstrap replicate uses
   it (`--numBootstraps` runs `replicates` sequential EMs inside a parallel
   iterator), so improving the sequential kernel improves posterior runs.

Things deliberately **not** worth revisiting, with reasons:

- *Shared atomic accumulator instead of shards* — tried; CAS contention on hot
  transcripts dominated. That is why shards exist.
- *Caching `alpha_in[tid] * w` in scratch inside the parallel M-step* — tried;
  the extra memory traffic cost more than the recomputed multiply.
- *Making the shard count track threads* — this is the #1167 bug. Never again.

---

## 9. Test data provenance (on `newton`)

Everything needed to reproduce the numbers above is staged on **`newton`**
under **`/mnt/scratch7/rob_disktest/`** (spinning disk, ~105 MB/s; 461 GB free
at time of writing). `/mnt/scratch1` is also rotational (~145 MB/s).

```
/mnt/scratch7/rob_disktest/
├── bin/
│   ├── salmon              pre-#1167 binary (baseline)
│   └── salmon_fixed        post-#1167 binary
├── data/
│   ├── gencode_v49/        salmon index, 3.1 GB, 238,442 transcripts
│   ├── SRR21186103_1.fastq.gz   839 MB
│   └── SRR21186103_2.fastq.gz   837 MB
└── *.tsv, *.log            results from the runs described here
```

**Dataset**: SRR21186103, human, paired-end 2×~100 bp, **26,135,185 fragments**,
71.8% mapping rate against GENCODE v49. Copied from
`/scratch1/rob/rshash_testing/human_reads/` on the submit host; the index came
from `/scratch1/rob/flex_bench/salmon_idx/gencode_v49`.

Scale it produces: **536,311 equivalence classes**, ~1.3M incidences,
**2811 EM iterations** to convergence (deterministic path; the online path
varies, 2555–2700).

### Reproducing an EM measurement

```bash
cd /mnt/scratch7/rob_disktest
./bin/salmon_fixed quant -i data/gencode_v49 -l A \
  -1 data/SRR21186103_1.fastq.gz -2 data/SRR21186103_2.fastq.gz \
  -p 32 --sketch -o out/whatever 2> run.log

# per-phase timings (the log carries ANSI escapes — strip them first, this
# tripped up an earlier parser):
sed 's/\x1b\[[0-9;]*m//g' run.log \
  | grep -oE 'phase="[a-z_]+" elapsed_s=[0-9.]+'
```

Phase names: `index_load`, `mapping`, `eff_length_collapse`, `em_bias`,
`posterior`, `output` (phase 1), then `rad_read`, `em_bias`, `posterior`,
`output` (phase 2). **The EM you care about is the *second* `em_bias`.** The
phases sum to wall time, so an unaccounted remainder means a parse error.

Use `--sketch` for EM work: mapping is cheap, so the EM is a large share and
iteration counts are unaffected.

Kernel/loop selection (verified against `salmon quant --help`):

| flag | effect |
|---|---|
| *(none)* | VBEM — **the production default** |
| `--useEM` | plain EM kernel |
| `--meta` | plain EM + no range factorization (overrides `--useEM`) |
| `--emAccel none\|squarem\|daarem` | convergence loop; `none` is the default |
| `--perNucleotidePrior` | VBEM prior becomes `vb_prior * effLen` |
| `--numBootstraps N` | N sequential EMs inside a parallel iterator (exercises `em_step_seq`) |

Benchmark VBEM unless you specifically mean to measure plain EM — the default
`EmOptions` says otherwise and it is an easy mistake.

### Benchmarking caveats learned the hard way

- **Warm the page cache first.** With 503 GB of RAM, the first run over this
  1.7 GB of input is materially slower than later ones. An early A/B here
  showed a "19% penalty" that was entirely cold-cache input reads attributed to
  whichever arm ran first. Run once to warm, then measure.
- **Do not interleave arms with other work.** Wall-time comparisons taken from
  a script that ran other binaries between arms disagreed with back-to-back
  phase measurements by ~20 points.
- **The online path has enormous variance at high `-p`.** Four repeats at
  `-p 64` gave 103.6, 59.8, 104.5, 55.5 s — bimodal, an 88% spread — while the
  deterministic path held 39–42 s. Any single-run comparison against `--online`
  at high thread counts is meaningless.
- **`salmon --version` on a fresh binary first.** newton is Broadwell with no
  AVX-512; the build config floors at `x86-64-v2` with runtime SIMD dispatch,
  so binaries built on the AVX-512 submit host do run — but a
  `target-cpu=native` build would `SIGILL` there.
- No root on newton, so caches cannot be dropped between runs.

### A smaller local fixture

For quick iteration without newton, the determinism test's synthetic generator
(`packed.rs::shard_plan_determinism`) builds 20,000 classes over 5,000
transcripts with non-representable weights — enough to shard several times and
to expose grouping differences, and it runs in under a second.
