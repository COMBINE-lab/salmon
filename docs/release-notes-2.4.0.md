# salmon 2.4.0

A feature release on 2.3.4. It adds BAM mapping output, an explicit
fragment-length policy for RAD requantification, and a matching `--writeSam`
spelling for the existing SAM output, alongside a full dependency refresh.

**Quantification results are unchanged.** `quant.sf`, index provenance hashes,
`fld.gz`, and bootstrap replicates are byte-identical to 2.3.4, and Gibbs
samples agree to within one floating-point ULP. **No index rebuild is required.**

## `--writeBam`: BAM mapping output (#1065)

`salmon quant --writeBam <FILE>` writes the retained mappings as unsorted,
BGZF-compressed BAM. The records are semantically equivalent to `--writeMappings`
SAM output — paired, orphan and single-end flags, spoofed overhang-aware CIGARs,
sequence orientation, mate coordinates, template length, mapping quality, and
the `NH`/`HI`/`XT`/`AS` tags.

SAM and BAM now share one format-neutral record builder and one header, so the
two cannot drift apart. `--writeMappings` output is byte-identical to 2.3.4.

The two options are mutually exclusive.

This feature was proposed and first implemented by @BenjaminDEMAILLE in #1029.

### Record order

Exactly one ordering guarantee: **all records for a fragment are contiguous.** A
worker encodes a whole fragment before handing off its buffer, and buffers are
written whole, so a fragment is never split or interleaved with another worker's
records — tools that stream and group by read name can consume the output
directly.

Nothing else is ordered. Fragments appear in whatever order workers finish, which
varies with thread count. Sort the result if a downstream tool needs coordinate
or query-name order.

### `--bamCompressThreads`

BGZF compression runs on a small dedicated pool alongside the `-p` mapping
threads. The default is derived from measured throughput — one worker compresses
about 165 MiB/s of BAM records while one mapping thread produces at most about
52 MiB/s, so roughly one compression worker per three mapping threads, capped at
eight. `--bamCompressThreads <N>` overrides it.

The default deliberately errs upward but stops at that balance point: one worker
too few can halve throughput, while surplus workers merely block on an empty
queue — but going further spends cores that would otherwise be mapping reads.

## `--writeSam`: a format-named alias for `--writeMappings` (#1065)

`--writeSam` is a visible alias for `-z/--writeMappings`, so SAM output has a
format-named flag to pair with `--writeBam`. It shares the same argument, so it
inherits the `--writeBam` mutual exclusion.

## `--fldPolicy`: fragment-length distribution control for `--rad` (#1062, #1066)

salmon's RAD writer bakes its fragment-length distribution into the header —
**including under `--skipQuant`** — and the reader preferred it over
`--fldMean`/`--fldSD`/`--fldMax`. That default is right for reproducing the
writing run, but it meant those flags were silently ignored on every
salmon-produced RAD, with no way to override. A fragment-length sensitivity
analysis could run to completion having perturbed nothing.

salmon now **warns** when a baked distribution supersedes explicitly-supplied
FLD flags, naming the flags, the baked mean and SD, whether they were observed,
and the remedy. Runs that inherit the defaults stay quiet.

`--fldPolicy` chooses where the distribution comes from:

| Policy | Behavior |
| --- | --- |
| `baked` *(default)* | Use the RAD's baked distribution. Exact parity with the writing run; behavior unchanged from 2.3.4. |
| `derive` | Ignore it; rebuild from this RAD's own uniquely-mapped proper pairs. |
| `prior` | Ignore both; `--fldMean`/`--fldSD` alone determine it. |

Use `prior` for a sensitivity analysis — it is the only setting under which
varying `--fldMean` changes the result.

`aux_info/meta_info.json` now records `frag_length_source` in every mode:
`reads`, `alignments`, `rad_baked`, `rad_baked_prior`, `rad_derived`, or `prior`.

### Single-end RADs

A single-end run has no fragment lengths to measure, so what it bakes is its own
`--fldMean`/`--fldSD` prior. Reading such a RAD back therefore inherits *the
writing run's* prior. salmon reports this as `rad_baked_prior` to distinguish it
from an observed distribution; pass `--fldPolicy prior` to apply your own values.

Reported by @Hugo-Polloli.

## Docker images include `procps` (#1064)

The published images now ship `procps`, so `ps` is available and
[Nextflow's container requirements](https://www.nextflow.io/docs/latest/container.html)
are met. Without it every task failed with *"Command `ps` required by nextflow
to collect task metrics cannot be found"*.

Note that the image entrypoint is `salmon`, so a workflow manager supplying its
own command needs the entrypoint overridden — `docker.entrypointOverride = true`
for Nextflow. See the installation docs.

Reported by @rich7409.

## Dependency refresh (#1068, #1070)

A full audit of all 42 declared dependencies. Everything is now a single version
in the tree, which required publishing four COMBINE-lab crates first:
`bramble-rs` 0.1.7, `libradicl` 0.15.0, `cf1-rs` 0.5.1 and `sshash-lib` 0.6.2.

Notable moves: `libradicl` 0.15, `noodles-sam`/`-bam`/`-bgzf` 0.87/0.92/0.49,
`needletail` 0.7, `sha2` 0.11, `indicatif` 0.18, `scc` 3, `statrs` 0.19, and the
`rand` family to 0.10 (collapsing three tree versions to one).

Every step was verified output-neutral: index provenance hashes, `fld.gz`,
`quant.sf` and bootstraps byte-identical, Gibbs within one ULP.

`indicatif` 0.18 also drops `number_prefix`, clearing RUSTSEC-2025-0119.

### BGZF thread control

noodles 0.48 deprecated `MultithreadedReader/Writer::with_worker_count` and made
it a no-op — block work now dispatches through `rayon::spawn`. Left alone, BGZF
compression would have landed on salmon's global rayon pool (sized to `-p` and
already mapping), and both the derived worker count and `--bamCompressThreads`
would have controlled nothing. salmon now runs each BGZF reader and writer inside
a dedicated rayon pool, which also keeps compression off the mapping threads.
Reported upstream as zaeleus/noodles#409.

## Upgrading

Drop-in from 2.3.4: **no index rebuild**, and quantification results are
unchanged for every mode. `--writeMappings` output is byte-identical.

The one behavior change is a new warning: if you pass `--fldMean`, `--fldSD` or
`--fldMax` alongside `--rad` on a salmon-produced RAD, salmon now tells you those
values are being superseded by the baked distribution. The result is the same as
2.3.4 — use `--fldPolicy` if you want the flags to take effect.
