# Which global allocator salmon should ship, and what the alternatives measure

salmon ships with mimalloc as its `#[global_allocator]`. This note records the
measurement behind that choice, so the question does not have to be re-argued
from first principles, and documents the build flags that let anyone re-run it
on their own hardware. Numbers here are measured; where a difference is inside
run-to-run noise that is said rather than rounded into a winner.

## The build flags

`crates/salmon-cli/src/global_alloc.rs` selects one of four allocators. Cargo
features are additive, so the `cfg`s there encode a fixed priority (sysalloc >
snmalloc > jemalloc > mimalloc) and any combination resolves to exactly one
`#[global_allocator]`. The binary logs which one it got at debug level:

```
cargo build --release -p salmon-cli                      # mimalloc (default)
cargo build --release -p salmon-cli --features sysalloc  # system malloc
cargo build --release -p salmon-cli --features jemalloc  # jemalloc
cargo build --release -p salmon-cli --features snmalloc  # snmalloc

RUST_LOG=debug ./salmon swim   #  -> "global allocator: mimalloc"
```

The alternates are diagnostics, off by default, and add nothing to a normal
release build. `tikv-jemallocator` is pinned to 0.5 because piscem-rs -> sux ->
rdst already links `tikv-jemalloc-sys` 0.5, and two versions of a `links =
"jemalloc"` native library cannot coexist in one graph.

## The workload

- Index: GENCODE v50 human transcripts, 654,828 references, k=31, m=19 (~2 GB on disk).
- Reads: 10M simulated 2x100 bp pairs, lognormal transcript abundance, ~250 bp
  mean fragment, 0.5% substitution rate. Deterministic, so every allocator sees
  byte-identical input.
- `salmon quant -l A -p 16`, Apple M4 Max (16 cores), macOS.
- 5 reps, allocators interleaved within each rep rather than grouped, so thermal
  drift and background load hit all four equally.

Only the `mapping` phase allocates heavily (per-read scratch buffers allocated
on the parser thread and freed on a worker thread, the cross-thread free
pattern allocators differ most on). `em_bias` is compute-bound and serves as a
control: an allocator that "wins" there is measuring noise.

## Results

Per-rep `mapping` seconds, and the median, from two independent clean rounds:

| round | mimalloc | jemalloc | snmalloc | system |
|-------|----------|----------|----------|--------|
| A (4-way, 5 reps) | 31.0 | 33.2 | 30.8 | 42.7 |
| B (2-way, 5 reps) | 34.1 | n/a | 36.3 | n/a |

Peak RSS, round A: mimalloc 3.84 GiB, jemalloc 3.66, snmalloc 4.22, system 3.71.

Read it as:

- **System malloc is reliably and substantially worse**: +38% on the mapping
  phase, ~+20% on total wall. This is the finding that justifies shipping *some*
  replacement allocator at all, and it reproduced in every round.
- **jemalloc does not win.** ~7% behind mimalloc on mapping, consistently. Its
  strengths (fragmentation control, RSS on long-lived large heaps) are not what
  a batch quant run is bounded by, and its macOS support is the weakest of the
  three.
- **mimalloc and snmalloc are tied.** Round A puts snmalloc marginally ahead,
  round B puts mimalloc ahead by 6%; the rounds disagree in *direction*, which
  means the difference is smaller than the noise floor. snmalloc also costs
  ~0.4 GiB more peak RSS.

So: keep mimalloc. Nothing here justifies a switch, and the +0.4 GiB is a real
cost against a speedup that does not reproduce.

## Two traps this measurement fell into

Worth recording, because both produced confident-looking wrong answers first.

**A single-threaded allocation microbenchmark predicts nothing.**
`examples/cache_bench.rs` (alloc+free of a `MappingCache` on one thread) ranks
them snmalloc 43 ns/op, system 70, mimalloc 86, jemalloc 89: snmalloc twice as
fast as mimalloc, and system malloc *beating* mimalloc. Neither survives contact
with the real workload, where system malloc is worst by a wide margin. Same-thread
alloc/free is snmalloc's fast path and macOS's nano-zone's fast path; it is not
what salmon does.

**The first round of the real benchmark was contaminated**: cold page cache
over 4 GB of freshly written reads plus residual load. Mapping times ranged 43-68s
where the settled machine gives 29-37s, and the ranking came out
snmalloc > jemalloc > mimalloc, the reverse of the clean result. Absolute numbers
also drift between rounds (round B's mimalloc is 34.1 where round A's is 31.0 for
the same binary), so **only within-round comparisons are valid**. Warm up, watch
the spread, and distrust any ranking whose per-rep variance exceeds the gap it
claims to have found.

## Re-running it

The harness is not checked in (it needs a multi-GB index and read set), but it is
four lines: build the four binaries, simulate reads against a real transcriptome,
and loop `rep { for allocator }` under `/usr/bin/time -l`, reading the per-phase
timings off the `salmon::timing` tracing target. `scripts/profile.sh` already
builds the same workload shape for its GIAB tier.
