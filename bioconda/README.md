# bioconda recipe drafts

These are **drafts** of the bioconda recipes for the salmon 2.0 (Rust) cutover.
They are kept here in the salmon repo so they are ready to submit; the canonical
copies live in
[`bioconda/bioconda-recipes`](https://github.com/bioconda/bioconda-recipes)
under `recipes/`.

## `salmon/` — salmon 2.0 (Rust)

Replaces the current C++ `salmon` recipe at the 2.0 release. Builds the `salmon`
binary from the `salmon-cli` crate with `cargo install --locked` (reproducible
via the committed `Cargo.lock`). Requires the published crates.io versions of
cf1-rs / piscem-rs / ksw2rs (release plan R2/R3).

**Before submitting:** fill in the `sha256` of the `v2.0.0` source tarball.

## `salmon-cpp/` — final C++ salmon (1.12.0)

Preserves the C++ line for reproducibility / emergency fixes. It is the C++
1.12.0 recipe renamed to `salmon-cpp`, with pufferfish pinned to `5dce7f41` (the
orientation fix shipped in 1.12.0). It installs a `salmon` binary, so it
**conflicts with the `salmon` package** — install one or the other.

## Submission order

1. **Now:** the autobump PR
   ([#66276](https://github.com/bioconda/bioconda-recipes/pull/66276)) ships C++
   `salmon` 1.12.0 (with the pufferfish-pin fix already pushed).
2. **At 2.0 release (plan R7/R8):** open a PR that (a) replaces `recipes/salmon`
   with the Rust recipe here, and (b) adds `recipes/salmon-cpp` from here.

The `doc_url` for the 2.0 recipe points at the new Astro docs site
(`https://combine-lab.github.io/salmon`).
