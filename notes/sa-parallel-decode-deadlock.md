# Deadlock: broker-driven resizing + paired parallel decode

**Status: open. Blocks enabling the parallel decoder by default in *either*
mapping mode.** Found while measuring the engagement threshold.

> **Correction.** This was first written up as selective-alignment-specific,
> because sketch completed in 36 s with identical decoder settings. It is not.
> A later run of `--sketch -p 64 --decoder parallel` **hung as well**. Sketch was
> winning a race, not avoiding a bug — which the first draft listed as an open
> question and should not have leaned on. Treat the mode as affecting only *how
> often* it hangs, not *whether*.

## Reproduction

```
salmon quant -i <gencode_v49> -l A \
  -1 SRR21186103_1.fastq.gz -2 SRR21186103_2.fastq.gz \
  -p 64 --decoder parallel            # no --sketch
```

26.1 M fragments, 2 gzip inputs. Hangs. Observed twice: once killed at 10 min,
once left for 6 min at **0.0 cores busy measured over three 5 s deltas** with
134 threads alive. It is a deadlock, not a slow run.

## Bisect

Full input, 200 s timeout, `HUNG` = killed at the limit.

| config | broker | paired | result |
|---|---|---|---|
| `-p 64 parallel` | Adaptive | yes | **HUNG** (x2) |
| `-p 64 parallel=1` (2 slots) | Pinned | yes | ok 58 s |
| `-p 64 parallel=2` (4 slots) | Pinned | yes | ok 54 s |
| `-p 64 parallel=4` (**8 slots**) | Pinned | yes | ok 59 s |
| `-p 16 parallel` | Adaptive | yes | ok 129 s |
| `-p 32 parallel` | Adaptive | yes | ok 75 s |
| `-p 64 parallel` single-end | Adaptive | no | ok 39 s |
| `--sketch -p 64 parallel` | Adaptive | yes | **HUNG** |

**The split is not the variable.** `parallel=4` plans exactly the same 56 mapping
/ 8 decode split as plain `parallel` and completes; the difference is that
`workers_per_file: Some(n)` yields `DecodeAllocation::PinnedPerFile`, so
`adaptive()` is false and **the broker never runs**. `None` yields `Adaptive`
and the broker resizes the mapping pool during the run.

Three conditions are jointly required; removing any one makes it complete:

1. **the broker actively resizing** (Adaptive, not Pinned),
2. **paired input** — paraseq's paired `fill` takes R1 and R2 as *separate*
   locks; single-end has one and does not hang,
3. **enough mapping threads to contend** — `-p 16` and `-p 32` survive, `-p 64`
   does not.

Mapping mode changes only the odds: selective alignment hung 2 of 2, sketch 1 of
2. A slower consumer holds the reader lock longer, widening the window a resize
can land in. This is also why piscem drives the same decoder pool without
hanging — its per-fragment work is far cheaper.

## Evidence

`gdb -p <pid> -batch -ex "thread apply all bt"`, 134 threads
(`sa-deadlock-backtrace.txt` in the run directory):

| threads | blocked in |
|---:|---|
| 62 | `parking_lot::RawMutex::lock_slow` <- `paraseq::parallel::paired::PairedReader` |
| 2 | `rapidgzip_core::pool::PoolMember::acquire` |
| 1 | `rapidgzip_core::reader` `Message::recv` |
| 64 | rayon `wait_until_cold` (idle; salmon's own pool, expected) |

The cycle: mapping workers queue on the paired-reader mutex; the holder is
inside `fill` waiting on rapidgzip; rapidgzip's workers are waiting to acquire a
decode-pool permit. Permits are not being released, so nothing advances.

Paired `fill` locks R1 and R2 *separately* (`paraseq/parallel/paired.rs`), so a
pair sustains two concurrent inflate streams against one shared pool — the
likeliest place for permits to be held by a stream whose consumer is itself
blocked behind the other stream's lock.

## Why this is not only a `--decoder parallel` problem, or an SA problem

With the provisional selective-alignment threshold of 29 per stream and F = 2,
`auto` engages from `-p 58` upward. So **the default path reaches this too** at
any budget at or above ~58:

| `-p` | sketch engages | SA engages |
|---:|---|---|
| 32 | yes | no |
| 64 | yes | **yes** |
| 128 | yes | **yes** |

A lower measured threshold would make it *more* exposed, not less.

## What is not yet known

Answered by the bisect: not the slot count, yes it needs enough mapping
threads, and sketch is *not* safe.

Still open:

* Where exactly the cycle closes. The three necessary conditions are known, but
  not which lock is held across which resize. A resize retires a worker at a
  batch boundary (`try_release_slot`, after the order gate advances); the
  suspect window is a worker inside paired `fill` holding R1 and waiting for R2
  while the pool shrinks under it.
* Whether the fix belongs in paraseq (lock both mates as one unit, or refuse to
  retire inside `fill`), in the broker (do not resize while any worker is
  inside a producer call), or in how salmon sizes the shared pool.
* Whether `-p 48` or similar also hangs — the boundary between 32 (ok) and 64
  (hangs) has not been bisected.

**Do not enable the parallel decoder by default in any mode until this is
resolved.** Pinning slots (`--decoder parallel=N`) is a working escape hatch
today precisely because it disables the broker, which is the feature.
