# The sshash bucket-bounds cache, and what it is worth to salmon

Notes on an upstream `sshash-lib` optimization, why salmon gets it without any
salmon-side change, and how much it actually moves each mapping mode. Numbers
here are measured, not projected; where a figure is inside the measurement noise
that is said rather than rounded away.

## What the change is

Every sshash dictionary *seed* resolves a minimizer to a bucket before it can
scan it. That resolution is an MPHF evaluation (`PartitionedMphf::get`) followed
by `locate_bucket`, which is two Elias-Fano accesses and therefore two
`DArray::select` calls. Together they were ~25% of cycles in the seed path.

Both steps depend only on the *minimizer value*. Consecutive k-mers share a
minimizer for runs averaging `(k-m+2)/2`, so most seeds re-derive bounds they
just computed. Measured on gencode v49 (k=31, m=19): **71.3% of 46.3M seeds had
the same minimizer as the previous k-mer and a found minimizer**, each redoing
the full resolution for an answer already in hand.

The fix caches the resolved bounds on the `StreamingQuery`, keyed on minimizer
value. Only the *bounds* are cached; the per-k-mer `MinimizerInfo`, which carries
the minimizer's offset within the k-mer, still flows into the bucket scan
untouched, so results are unchanged. The mapping from minimizer to bounds is a
pure function of an immutable dictionary, so entries never go stale and need no
invalidation, not even across reads.

It stacks with the existing within-unitig extension rather than overlapping it:
`seed_with_dict` only runs once `remaining_string_bases == 0`, so the 71.3% is
already the population that survives extension.

## Why salmon inherits it for free

`salmon-map` does not talk to `sshash-lib` directly for queries. `collect.rs`,
`sketch.rs` and `extend.rs` all construct
`piscem_rs::mapping::streaming_query::PiscemStreamingQuery`, which wraps
sshash's `StreamingQuery` — exactly where the cache lives. **No salmon code
changes.** It arrives with a dependency bump:

| crate | salmon pins | needs |
|---|---|---|
| `piscem-rs` | 0.6.4 | rebuild against a patched `sshash-lib` |
| `sshash-lib` | 0.6 | >= the release carrying the cache |

Cargo only honours a `[patch]` table from the workspace root, so `piscem-rs`'s
own patch section does **not** apply transitively — salmon has to redirect
`sshash-lib` itself. That is what the temporary `[patch.crates-io]` block on this
branch does.

## Measured effect

Index: gencode v49 pc transcripts, k=31, 238,442 references.
Reads: SRR21186103 (150 bp human bulk RNA-seq, 88.8% mapping rate).
Machine: AMD EPYC 9575F, `-p 32`.

### Upstream, in isolation

`piscem-rs map-bulk --dict sshash`, 32M reads, 3 reps:

| | CPU |
|---|---|
| without cache | 105.26 s |
| with cache | **94.68 s (-10.0%)** |

### In salmon

| mode | measurement | without | with | delta |
|---|---|---|---|---|
| selective alignment | mapping phase (wall, 8M) | 8.935 s | 8.860 s | -0.8% |
| selective alignment | `--skipQuant` CPU (8M) | 256.11 s | 254.33 s | -0.7% |
| full quant (SA + VBEM) | total CPU (8M) | 343.00 s | 338.08 s | -1.4% |
| **sketch** | mapping phase (wall, 8M) | 2.112 s | 1.999 s | **-5.4%** |
| **sketch** | `--skipQuant` CPU (32M) | 222.14 s | 212.00 s | **-4.6%** |
| sketch | `--skipQuant` CPU (8M) | 60.54 s | 59.60 s | -1.6% |

### Reading these

**The mode is what matters, and the reason is dilution.** For the same 8M reads,
mapping-only CPU is 256 s under selective alignment versus 60 s in sketch mode —
selective alignment is roughly three quarters of SA mapping work. The dictionary
seed path the cache accelerates is a small slice of that, so a 10% upstream win
lands as under 1%. Sketch mode forgoes alignment, the seed path is a much larger
share, and the same change is worth ~5%.

**Sketch is the mode to quote, at roughly 4-5%, with wide error bars.** Three
sketch measurements gave -5.4%, -4.6% and -1.6%; the two longer-running ones
agree near 5% and the 8M CPU run is the noisiest (bcache won only 5 of 8 pairs,
with a 7% spread across repetitions of the *same* binary). Do not quote a
sharper figure than "about 4-5%" from this data.

**Sub-1% figures in this table are not evidence.** Separately built binaries
differ in code and data layout, which produces a stable, directional offset of
around 1% that repetition does not average away. The selective-alignment rows
sit inside that; they are reported for completeness and should be read as "no
measurable change", not as a small win.

## Status

Nothing here changes salmon's behaviour or output — the cache is a pure
memoization and mapping results are unchanged. The only salmon-side action is
the dependency bump once `sshash-lib` ships the change, at which point the
`[patch.crates-io]` block added on this branch should be dropped.
