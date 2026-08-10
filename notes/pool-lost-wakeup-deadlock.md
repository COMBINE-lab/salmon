# Deadlock: broker-driven resizing + paired parallel decode

**Status: RESOLVED — rapidgzip-core 0.3.1 / piscem-rs 0.9.1.** Root cause was a
lost wakeup in rapidgzip's `DecoderPool::release`, which notified waiters
without holding the scheduler lock; the final release could land between a
waiter's admission check and its park, leaving free permits, parked decoders,
and no future notify. Fixed upstream (rapidgzip-rust `e513d6d`, with a
regression test that deadlocks the unfixed pool in under two seconds), and
verified here: all six previously-hanging configurations (SA -p 48/56/64 x2,
sketch -p 64/128) complete in 37–62 s with mapped counts identical to their
serial-decoder runs. Neither paraseq nor salmon needed changes; the paired
R1/R2 locking merely raised the slow-path traffic that made the window
reachable. The history below is kept as found.

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
| `--sketch -p 128 parallel` | Adaptive | yes | **HUNG** |
| `-p 40 parallel` | Adaptive | yes | ok 64 s |
| `-p 48 parallel` | Adaptive | yes | **HUNG** |
| `-p 56 parallel` | Adaptive | yes | **HUNG** |

**The split is not the variable.** `parallel=4` plans exactly the same 56 mapping
/ 8 decode split as plain `parallel` and completes; the difference is that
`workers_per_file: Some(n)` yields `DecodeAllocation::PinnedPerFile`, so
`adaptive()` is false and **the broker never runs**. `None` yields `Adaptive`
and the broker resizes the mapping pool during the run.

Three conditions are jointly required; removing any one makes it complete:

1. **the broker actively resizing** (Adaptive, not Pinned),
2. **paired input** — paraseq's paired `fill` takes R1 and R2 as *separate*
   locks; single-end has one and does not hang,
3. **enough mapping threads to contend** — `-p 16`, `-p 32` and `-p 40` survive;
   `-p 48` and above do not. Nothing about the planned split is discontinuous
   across that boundary (35/5 completes, 42/6 hangs), which is what a widening
   race window looks like rather than a threshold being crossed.

Mapping mode changes only the odds: selective alignment hung 2 of 2 at `-p 64`,
sketch 2 of 3 (the one survivor being the 36 s run that produced the cost-share
measurement). A slower consumer holds the reader lock longer, widening the window a resize
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
* ~~Where the 32/64 boundary lies.~~ Bisected: `-p 40` completes (64 s),
  `-p 48` and `-p 56` hang. The opening splits either side are
  35 map / 5 decode (ok) and 42 map / 6 decode (hangs), so nothing about the
  split itself is discontinuous there — consistent with a race whose window
  simply widens with thread count rather than a threshold effect.

~~Do not enable the parallel decoder by default in any mode until this is
resolved.~~ Resolved as above; the default path is safe with
rapidgzip-core >= 0.3.1.

## Fastest route to a regression test

The existing `tests/pool_resize_flush_safety.rs` already churns the pool between
1 and 16 workers, but over a **single-end** collection — the one arrangement the
bisect shows does *not* hang. Rebuilding that test against a
`CollectionType::Paired` collection should reproduce this in seconds instead of
a 200 s timeout, and would be the right thing to have in hand before attempting
a fix, so the fix can be shown to work rather than assumed to.
