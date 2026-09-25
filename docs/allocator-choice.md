# Which global allocator salmon ships, and what would justify changing it

salmon ships with mimalloc as its `#[global_allocator]`. This note records the
measurement behind that choice, so the question does not have to be re-argued
from first principles, and states the evidence a candidate replacement has to
produce before it is worth carrying in the build.

The short version: the alternatives were measured, none of them beat mimalloc
end to end, and so none of them is in the tree. The measurement is kept here
because a negative result nobody wrote down gets re-litigated every year.

## Where the choice lives

`crates/salmon-cli/src/global_alloc.rs` holds the `#[global_allocator]` and the
`NAME` constant reporting which allocator a build actually got. `main.rs` logs
it at debug level, and `examples/cache_bench.rs` pulls the same file in via
`#[path]`, so the microbenchmark and the real binary can never disagree about
what is in play:

```
cargo build --release -p salmon-cli                      # mimalloc (default)
cargo build --release -p salmon-cli --features sysalloc  # system malloc

RUST_LOG=debug ./salmon swim   #  -> "global allocator: mimalloc"
```

`sysalloc` exists as a diagnostic escape hatch (isolating an allocator-specific
bug, or building where mimalloc will not). It is not a supported configuration.

## The gates a replacement has to clear

An allocator is not a free swap: it is native code linked into every build, a
second thing to keep building on every platform salmon ships to, and a variable
that every later performance measurement has to control for. A microbenchmark
win does not pay for that. Before either adding an allocator as a build feature
or changing the default, a proposal has to show:

1. **An end-to-end win on a real quant run.** Either measurably faster wall
   clock without a non-trivial increase in peak RSS, or the same wall clock
   with a considerable decrease in peak RSS. Phase timings and allocator
   microbenchmarks are diagnostics, not evidence.

2. **That the win holds across the input, output, and option surfaces.** At a
   minimum: selective alignment and sketch mode, with and without decoys,
   with and without mapping output, and across a realistic range of `-p`.
   An allocator that wins on one configuration and loses on another has
   found a workload, not an improvement.

Supporting discipline, since these numbers are easy to get wrong (see the traps
below): repetitions with the arms interleaved rather than grouped, the spread
reported alongside the median, `quant.sf` hashed to prove the arms are
computing the same thing, and a warmed page cache.

This is the same bar applied to other performance proposals in this repository:
the hasher swap in #1126 (5 to 8 percent in isolation, 0 percent end to end)
and the GPU backend in #1041 were both held to it.

## The measurement behind the current default

Round of 2026-08, on the four candidates that were built for it (mimalloc,
jemalloc, snmalloc, system malloc):

- Index: GENCODE v50 human transcripts, 654,828 references, k=31, m=19 (~2 GB
  on disk).
- Reads: 10M simulated 2x100 bp pairs, lognormal transcript abundance, ~250 bp
  mean fragment, 0.5% substitution rate. Deterministic, so every allocator sees
  byte-identical input.
- `salmon quant -l A -p 16`, Apple M4 Max (16 cores), macOS.
- 5 reps, allocators interleaved within each rep rather than grouped, so
  thermal drift and background load hit all four equally.

Only the `mapping` phase allocates heavily (per-read scratch buffers allocated
on the parser thread and freed on a worker thread, the cross-thread free
pattern allocators differ most on). `em_bias` is compute-bound and serves as a
control: an allocator that "wins" there is measuring noise.

Per-rep `mapping` seconds, median, from two independent clean rounds:

| round | mimalloc | jemalloc | snmalloc | system |
|-------|----------|----------|----------|--------|
| A (4-way, 5 reps) | 31.0 | 33.2 | 30.8 | 42.7 |
| B (2-way, 5 reps) | 34.1 | n/a | 36.3 | n/a |

Peak RSS, round A: mimalloc 3.84 GiB, jemalloc 3.66, snmalloc 4.22, system 3.71.

Read it as:

- **System malloc is reliably and substantially worse**: +38% on the mapping
  phase, ~+20% on total wall. This is the finding that justifies shipping some
  replacement allocator at all, and it reproduced in every round.
- **jemalloc does not win.** ~7% behind mimalloc on mapping, consistently. Its
  strengths (fragmentation control, RSS on long-lived large heaps) are not what
  a batch quant run is bounded by, and its macOS support is the weakest of the
  three.
- **mimalloc and snmalloc are tied.** Round A puts snmalloc marginally ahead,
  round B puts mimalloc ahead by 6%; the rounds disagree in *direction*, which
  means the difference is smaller than the noise floor. snmalloc also costs
  ~0.4 GiB more peak RSS.

So: keep mimalloc. Nothing here justifies a switch, and snmalloc's +0.4 GiB is
a real cost against a speedup that does not reproduce.

Note what this round does **not** establish, by gate 2's standard: it is one
machine, one reference, one simulated read set, one option configuration. It is
enough to conclude "no candidate showed a win worth carrying", which is what it
concluded. It would not be enough to justify a switch had one of them looked
good, and a future proposal should not treat it as a template for sufficiency.

## Two traps this measurement fell into

Worth recording, because both produced confident-looking wrong answers first.

**A single-threaded allocation microbenchmark predicts nothing.**
`examples/cache_bench.rs` (alloc+free of a `MappingCache` on one thread) ranks
them snmalloc 43 ns/op, system 70, mimalloc 86, jemalloc 89: snmalloc twice as
fast as mimalloc, and system malloc *beating* mimalloc. Neither survives contact
with the real workload, where system malloc is worst by a wide margin.
Same-thread alloc/free is snmalloc's fast path and macOS's nano-zone's fast
path; it is not what salmon does.

**The first round of the real benchmark was contaminated**: cold page cache
over 4 GB of freshly written reads plus residual load. Mapping times ranged
43-68s where the settled machine gives 29-37s, and the ranking came out
snmalloc > jemalloc > mimalloc, the reverse of the clean result. Absolute
numbers also drift between rounds (round B's mimalloc is 34.1 where round A's
is 31.0 for the same binary), so **only within-round comparisons are valid**.
Warm up, watch the spread, and distrust any ranking whose per-rep variance
exceeds the gap it claims to have found.

## Re-running it

The alternate allocators are not features in the tree, so reproducing the table
above means adding them back locally: `tikv-jemallocator` (pinned to 0.5, to
match the `tikv-jemalloc-sys` that piscem-rs -> sux -> rdst already links, since
two versions of a `links = "jemalloc"` native library cannot coexist in one
graph) and `snmalloc-rs` 0.7, each behind its own feature, with the
`#[global_allocator]` statics in `global_alloc.rs` selected by `cfg`.

For the run itself, `scripts/bench_em_rad.sh` is the model to copy: interleaved
arms, `/usr/bin/time -f '%e\t%U\t%S\t%M'` for wall, CPU and peak RSS, results
appended to a TSV, and `quant.sf` hashed per row. That harness replays a RAD
and therefore excludes mapping, which is exactly the phase an allocator moves,
so an allocator round needs the full-quant equivalent rather than that script
as-is.
