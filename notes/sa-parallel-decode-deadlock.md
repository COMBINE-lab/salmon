# Deadlock: selective alignment + parallel decoder

**Status: open. Blocks enabling the parallel decoder by default in
selective-alignment mode.** Found while measuring the engagement threshold.

## Reproduction

```
salmon quant -i <gencode_v49> -l A \
  -1 SRR21186103_1.fastq.gz -2 SRR21186103_2.fastq.gz \
  -p 64 --decoder parallel            # no --sketch
```

26.1 M fragments, 2 gzip inputs. Hangs. Observed twice: once killed at 10 min,
once left for 6 min at **0.0 cores busy measured over three 5 s deltas** with
134 threads alive. It is a deadlock, not a slow run.

`--sketch` with *identical* decoder settings completes in 36 s, and
`--decoder serial` in selective-alignment mode completes in 55 s. Only the
combination hangs, which points at consumption rate rather than configuration.

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

## Why this is not only a `--decoder parallel` problem

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

* Whether it depends on the 8 opening decode slots, or on any slot count.
* Whether it reproduces at smaller `-p`, or needs enough mapping threads to
  saturate the reader mutex.
* Whether sketch is safe or merely fast enough to never lose the race —
  36 s is not proof of absence.
* Whether it is salmon's wiring, a rapidgzip pool-exhaustion bug, or the
  interaction of paraseq's split R1/R2 locking with a shared pool. piscem drives
  the same pool the same way without hanging, but piscem's consumer is much
  cheaper per fragment, which is exactly the variable that differs.

Bisecting decode slots at fixed `-p`, and testing single-end (one lock rather
than two), would separate these.
