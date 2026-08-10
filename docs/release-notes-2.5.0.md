# salmon 2.5.0

A feature and performance release on 2.4.1: 28 pull requests across three
streams — an adaptive thread scheduler shared with piscem, a mapping- and
quantification-side performance campaign, and a round of correctness and
input-handling fixes.

**Quantification results are unchanged.** With the serial decoder, `quant.sf`
is byte-identical to 2.4.1 at a fixed thread count. Across thread counts,
results agree far beyond the reported precision (measured maximum deviation
4.5×10⁻¹³ against a reporting step of 5×10⁻⁴; see "Determinism", below).
**No index rebuild is required** (though a rebuild is recommended to help 
with propagating the provenance of the duplicate retention strategy used 
during the index construction).

## Adaptive thread scheduling (#1119, #1093)

`-p` now names one execution-slot budget shared between gzip decoding and
mapping, rather than a mapping-thread count with decompression taken on top.
On compressed input a live controller (piscem-rs's `thread-broker`) measures
each side's cost while the run proceeds and moves slots to where they pay.

Whether the parallel decoder is engaged *at all* is decided per mapping mode,
because the answer depends on how much work a fragment costs relative to
inflating its bytes — and that ratio was measured, not assumed:

| mode | engages when `-p` ≥ (per gzip input) |
|---|---|
| selective alignment | 50 |
| selective alignment `--deterministic` | 50 |
| sketch | 10 |
| sketch `--deterministic` | 9 |

For selective alignment the crossover was bracketed empirically on 26 M read
pairs: the serial decoder wins at `-p 64` (by 12 %), the parallel decoder wins
at `-p 128` (by 7 %). Serial decoding is simply the right choice for SA at most
real budgets, and the default encodes that.

Two new flags:

- `--decoder {auto,serial,parallel,parallel=N}` — override the engagement
  decision. The decision is **per input file**: a plain, FIFO, or otherwise
  non-seekable input falls back to serial decoding alone, without costing the
  other inputs their parallel decoder, and is never probed (a byte read from a
  pipe cannot be re-read).
- `--threadPolicy <FILE>` — a JSON override of the engagement threshold, in the
  same format piscem uses. Unknown fields are an error, not a silent no-op.

Counter integrity under thread migration is guaranteed and tested: worker
threads may now retire mid-run when the controller moves a slot, and every
per-thread tally and output buffer (counters, bias, SAM/BAM/RAD buffers,
unmapped names) commits exactly once regardless — verified by a test that
forces ~134 retire/respawn cycles and checks delivery by fragment id.

This work surfaced a race in rapidgzip's shared decoder pool (a lost wakeup in
its permit release path) that could hang a run outright. It is fixed upstream
in `rapidgzip-core` 0.3.1, which this release requires via piscem-rs 0.9.1.

## Performance

A campaign across the mapping and quantification paths (#1084, #1085, #1091
series; #1094, #1095, #1096, #1099, #1100, #1113, #1114, #1117, #1118, #1120 —
several contributed by @BenjaminDEMAILLE):

- per-worker tallies merged once per thread instead of shared atomic
  read-modify-writes per fragment (up to 7 eliminated per fragment);
- scratch reuse through the alignment, chaining, scoring, RAD and positional-
  bias paths, replacing per-fragment allocation (including the `BTreeMap` →
  flat-scratch placement grouping and reusable ksw2 workspaces);
- the packed CSR equivalence-class layout is built once per run and pre-sized;
- `--writeUnmappedNames` streams through a bounded per-thread buffer instead of
  accumulating every name in memory;
- a fragment's placements are ordered once, not once per consumer.

These are result-preserving changes, and the campaign's tests treat them that
way — e.g. the two dedup implementations are asserted to agree *including on
ties*, a distinction that once changed 55 equivalence classes on a 10 M-fragment
run while every sample-data check passed.

## Determinism

`thread_count_invariant` previously asserted bit-identical results across
thread counts. That assertion was stronger than the code's contract: the EM
merges per-worker shards, so float addition order varies with scheduling, and
the test was flaky. Measured on a deliberately hard fixture (800 transcripts,
60 k fragments, 1–4-way ambiguity, 2–32 threads, 50 runs), the maximum
cross-thread deviation is 4.5×10⁻¹³ — nine orders of magnitude below the
5×10⁻⁴ step at which `quant.sf` reports `NumReads`, and it cannot compound
because each EM iteration renormalises. The test now gates at 10⁻⁹ (six orders
stricter than anything observable in the output), and repeated single-threaded
runs are still required to be bit-identical. PR #1056, which restructured the
reduction to chase bit-equality, is superseded.

## Fixes

- **A failed `--writeRad` no longer leaves an unmarked partial file** readers
  cannot distinguish from a complete one (#1105, #1106).
- **Packed CSR offsets widened to `u64`** ahead of inputs whose incidence
  counts exceed `u32::MAX`; the previous ceiling assertion is gone (#1097,
  #1112, #1100).
- **Debug builds no longer abort on decoy indices** from an over-strong
  orphan-status assertion (#1115, #1116).
- **`use_aux` burn-in gating is deterministic** and consistent between code
  paths (#1089, #1107).
- **`--geneMap` accepts gzip-compressed GTF/GFF, fails before the run instead
  of after it** when the file is unusable, and reports unmatched transcripts
  as one warning with counts and an actionable hint (`--ignoreTxVersion`)
  rather than a line per transcript (#1074, #1076, #1077, #1109).
- **Input validation errors name the offending file and flag** (#1110).

## Compression handling (#1102, #1103)

All input paths now detect compression from content, not file names, via one
shared helper: gzip (including BGZF and multi-member), bzip2, xz and zstd.
This closes the previous inconsistency where `salmon index -t` accepted zstd
but other paths did not. Verified in this release: zstd-compressed reads
(`-1/-2`), zstd targets (`quant -a -t`), zstd SAM (`-a`), and gzip content
under a misleading name all work.

## RAD provenance (#1108)

RAD files now carry provenance file tags: total observed fragments, dovetail
and score-filtered counts, decoy-best counts, and the identity of the index the
mappings were made against. The tags are additive and optional on read — files
from older versions parse unchanged, and absent counters are distinguishable
from zero.

## Dependencies

- `piscem-rs` 0.9.1 (thread broker; requires `rapidgzip-core` 0.3.1 with the
  decoder-pool deadlock fix)
- `paraseq` → `paraseq-temp` 0.5.0-pre.1, a republication of paraseq's
  unreleased dev-0.5.0 branch plus the resizable worker pool the broker acts
  through. The library target is still named `paraseq`. `paraseq` remains the
  official crate; salmon will return to it at 0.5.0.
- `libradicl` 0.17, `noodles` 0.115.

## Documentation

Beyond the feature docs above: an exhaustive commenting pass across the
workspace (#1073, ~3,400 lines — what/how/why comments on the mapping,
inference and I/O internals), the packed CSR layout and M-step documentation,
within-transcript bias notes, and a pre-release testing checklist under
`docs/release-testing-checklist.md`, which this release was run against.
