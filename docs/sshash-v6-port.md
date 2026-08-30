# The sshash v6 port (index format v2)

Salmon's k-mer dictionary, `sshash-lib`, was updated to 0.7 — a port of C++
SSHash v6.0.0 (jermp/sshash PR #93). This is an **index-format break**:
indices built by salmon ≤ 2.4.x (index format v1) must be rebuilt with
`salmon index`. Both gates fire with actionable messages: salmon's
`info.json` version check, and sshash's own `.ssi` header check underneath it.

## What changed

SSHash previously had two indexing modalities: *regular* (minimizer = argmin
of the m-mer hash) and *canonical* (pick the smaller value between the best
forward and best reverse m-mer, so a k-mer and its reverse complement share a
bucket). Canonical answers — the only kind mapping needs — cost a 4/3 factor
in super-kmer density because the two minimizer processes were compared under
a third, unrelated order.

v6 keeps a single order: the minimizer of x is the locus i minimizing
`h(κ(i))`, where `κ(i) = min(m-mer_i, rc(m-mer_i))` is the **canonical m-mer
at that locus**. κ is invariant under reverse complementation and its
sequence reverses under it, so x and rc(x) select literally the same m-mer —
same bucket, one probe, orientation reported — at the plain forward density
`2/(k−m+2)`. Ties are broken centre-closest (Cologni & Pibiri, Prop. 24),
which also makes the scheme *forward*: one super-kmer per sampled position,
enforced by a build-time check. There is no modality knob any more.

The `.ssi` header (now format 5.0) also records the build seed and a
hasher-magic guard, so a rapidhash behavior change fails loudly at load
instead of silently zeroing the mapping rate.

## What it means for salmon

- **Results: none.** The SPSS text and unitig numbering are untouched, so
  every `LookupResult` a mapping consumes is identical. Validated end to end:
  `quant.sf` is bit-identical (SA and sketch modes) between the pre-port and
  post-port stacks on `sample_data`, and bit-identical across all 193,760
  transcripts with identical `num_mapped` (33,008,042) on GRCh38 cDNA ×
  ERR188044 (36M read pairs).
- **Index**: total `.ssi` bytes shrink 8–13%; `salmon index` on GRCh38 cDNA
  is ~6% faster at ~17% lower peak RSS.
- **Mapping**: sshash streaming lookups are ~33% faster (single minimizer
  iterator + v6's same-minimizer memos for the negative-in-positive streams
  sequencing errors produce); human sketch quant ~3% faster wall overall.
- Also fixed upstream in sshash-rs 0.7: K ≥ 33 streaming produced wrong
  results in release builds (u64 truncation of u128 k-mer state), and a
  non-default build seed produced an index that returned zero hits.

## Validation trail (sshash-rs repo)

- `tests/minimizer_properties.rs`: mirror-equivariance, forwardness, and
  incremental-vs-brute-force against an independent from-the-spec reference.
- `scripts/oracle_diff.sh`: per-k-mer dumps byte-identical to a 0.6.4 oracle
  over ~30M lookups (point everywhere; streaming for K ≤ 31 — old streaming
  was wrong at K > 31).
- `scripts/cpp_dump.cpp`: per-k-mer dumps byte-identical to C++ SSHash v6
  itself at k = 31, 47, 63.
- piscem-rs: full suite + `tiny_sshash_parity` on a freshly built index.
