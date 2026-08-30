# salmon 2.7.0

The release that adopts the **unified canonical k-mer dictionary** — a port of
C++ SSHash v6.0.0 into `sshash-rs` 0.7 / `piscem-rs` 0.10. This is an **index
format break**: index format moves to **v2**, and every existing salmon index
must be rebuilt with `salmon index`. Quantification results are **unchanged —
provably**: against the released 2.6.0, `quant.sf` and the per-fragment RAD
mapping records are *byte-identical* on every validation dataset, from
`sample_data` (selective alignment and sketch, with and without decoys) up to
GRCh38 cDNA × 36M real read pairs (193,760 transcripts, identical
`num_mapped`). What you get for the rebuild is a smaller, faster index.

**You must rebuild your indexes.** Both gates say so explicitly: a v1 index is
rejected with an actionable message naming the salmon that built it, and (new
in this release) an index built by a *future* salmon is rejected with an
"upgrade salmon" message instead of a cryptic I/O error. `info.json` now also
records `sshash_format_version`, so future dictionary changes are diagnosable
without reading binary headers.

## What changed

SSHash previously had two indexing modalities: *regular* and *canonical*
(the salmon/pufferfish convention, where a k-mer and its reverse complement
share one entry). Canonical answers cost a 4/3 factor in super-k-mer density
because two independent minimizer processes were compared under a third,
unrelated order.

SSHash v6 keeps a single order: the minimizer of a k-mer is the position
minimizing the hash of the **canonical m-mer** at each position. This is
mirror-equivariant — a k-mer and its reverse complement select literally the
same m-mer, so one probe answers both strands and reports orientation — at
the plain forward density. With the centre-closest tie-break (Cologni &
Pibiri, Proposition 24) the scheme is also *forward*, which lets the builder
enforce one-super-k-mer-per-position as a hard invariant. There is no
modality flag anymore, in salmon or below it.

The dictionary header now records the build seed and a hasher guard, so a
hash-library drift fails loudly at load instead of silently collapsing the
mapping rate.

## Numbers (vs 2.6.0, same machine)

- **Index build** (GRCh38 cDNA): ~23% lower peak RSS (4.2 → 3.3 GB), ~3%
  faster wall.
- **Index size**: total dictionary bytes shrink 8–13% across test references.
- **Mapping**: the dictionary's streaming lookups are ~33% faster (single
  minimizer iterator plus v6's same-minimizer memos, which target exactly the
  negative-in-positive streams sequencing errors produce). End to end:
  sketch-mode map-only ~5.5% faster wall / ~5.7% less CPU on 36M real pairs;
  full sketch quant ~3% faster wall.

## Fixes that rode along

- Wide-k (k ≥ 33) streaming lookups produced wrong results in release builds
  (64-bit arithmetic on 128-bit k-mer state); nothing in salmon's default
  k = 31 path was affected, but the class is now fixed and pinned by tests.
- A non-default dictionary build seed used to produce an index that returned
  zero hits (queries hard-coded the default seed).

## Validation

The core invariant — the string set and its numbering are untouched, only the
index layer above changes — was tested exhaustively rather than assumed:
per-k-mer results byte-identical to the previous implementation over ~30M
lookups; byte-identical to C++ SSHash v6 itself at k = 31/47/63; property
tests (mirror-equivariance, forwardness, incremental-vs-reference) against an
independent implementation of the spec; and the end-to-end bit-identity runs
above. Full narrative: `docs/sshash-v6-port.md`.

## Compatibility notes

- `salmon quant` against a pre-2.7 index directory fails with a message
  telling you to rebuild; older salmon against a 2.7 index fails inside the
  dictionary loader with an "incompatible format version" error.
- No salmon command-line options changed. (In the underlying `piscem`/`sshash`
  tools, the removed `--canonical` modality flags are accepted for one release
  and warn.)
