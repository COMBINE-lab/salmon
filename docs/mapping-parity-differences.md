# Mapping parity vs. C++ salmon — an orientation bug the Rust port uncovered

This note records a real mapping-sensitivity difference we found between the
Rust port and C++ salmon on real short-read data, the investigation that traced
it, and its resolution. **The difference turned out to be a genuine bug in C++
salmon (pufferfish's SSHash k-mer lookup), now fixed upstream.** With the fix,
C++ salmon and the Rust port agree.

## TL;DR

On clean/simulated data the two tools always agreed (both map 100% of
`sample_data`, counts Pearson **r = 0.99998**). On **real short reads** the Rust
port mapped ~2% more reads than C++ salmon. That extra ~2% was not a Rust
"sensitivity choice" — it exposed a **k-mer-orientation bug** in pufferfish's
SSHash streaming lookup that caused C++ salmon to mis-place and discard a class
of reads. After fixing pufferfish (commit `5dce7f4`, salmon pin bumped), C++
salmon maps the same reads: **83.48% → 85.55%**, matching the Rust port to
within 1 read, with no change on clean data. The Rust port (built on
piscem-rs, which derives orientation correctly) was right all along.

## Dataset

- Reads: `ERR458493` (Gierliński/Schurch *S. cerevisiae* benchmark),
  1,093,957 single-end 51 bp reads.
- Reference: Ensembl R64-1-1 cDNA, 6,612 transcripts.
- Tools: C++ salmon 1.11.4 vs. this Rust port, both selective-alignment mode,
  `-l U`, no bias correction.

## What we observed

| metric | value |
| --- | --- |
| C++ salmon mapping rate (1.11.4, before fix) | 83.48% (913,271 reads) |
| Rust port mapping rate | 85.55% (935,850 reads) |
| C++ salmon mapping rate (after pufferfish fix) | **85.55% (935,851 reads)** |
| `NumReads` Pearson (rust vs salmon, confidently mapped) | 0.9989 |
| `EffectiveLength` Pearson | 0.9996 |

The well-supported transcripts always quantified near-identically; the only
material difference was the ~22,600 reads salmon dropped that the Rust port
mapped. Spot-checking showed those reads were genuine exact / near-exact matches
to a transcript — not spurious Rust mappings.

## The investigation (including the false starts)

We were honest about two wrong turns before landing on the cause; they're worth
recording because each taught us something and the tooling we built to rule them
out is what finally cracked it.

1. **`mismatchSeedSkip` (partial symptom, not the cause).** Lowering
   pufferfish's `mismatchSeedSkip` from 3 to 1 recovered most of the gap
   (83.48% → 85.00%). This *looked* like a seeding-skip bug, but it was really
   just collecting **more** seeds per read, which diluted the real (orientation)
   error — reads with several seeds were never affected because their
   correctly-oriented seeds dominated chain selection.

2. **Flank scoring (a misdiagnosis).** A synthetic length sweep suggested salmon
   under-scored long read flanks; we even prototyped a change to the flank
   alignment cutoff. It was net-neutral-to-negative and didn't fix the target
   reads — because the reads weren't being mis-*scored*, they were being
   mis-*placed*.

3. **The tool that cracked it: `writeMappings`.** salmon's `--writeMappings`
   (`-z`) was silently dropping SAM records (an unrelated flush bug we fixed —
   see below), so we couldn't see *where* salmon was placing these reads. With
   that fixed, the ground-truth placement was immediate: salmon was aligning the
   read on the **reverse strand at the wrong locus** (e.g. read `ERR458493.850`:
   salmon `REV YGR043C:402`, score 24/102; truth/Rust `FWD YGR043C:421`,
   score 96/102). A direct seed-landscape check confirmed the read's *only*
   31-mer seed was a unique **forward** hit — so the reverse placement had no
   seed support and was simply mis-oriented.

We also ruled out the **seed representation** itself: re-running the Rust mapper
with sparse fixed-`k` anchors, reference-extended MEMs, and true
unitig-constrained uni-MEMs gave identical results (85.55%, `NumReads`
r ≥ 0.99999995) — so granularity was never the issue.

## Root cause (C++ salmon / pufferfish)

pufferfish computes, per k-mer hit, whether the read k-mer lies forward on the
unitig (`hitFW`). Two lookup styles exist:

- **Non-streaming** lookups query the *canonical* k-mer word
  (`dict_.lookup(mer.getCanonicalWord())`); the returned orientation is
  canonical-relative and is converted to the query frame with
  `hitFW = (orientation == forward) == fwIsCanonical`.
- The **streaming** lookup queries the *raw query k-mer string*
  (`sq.lookup(kmer_str)`); its orientation is already query-relative, so the
  correct value is just `hitFW = (orientation == forward)`.

The streaming `getRefPos` overload mistakenly reused the canonical-relative
formula (`… == fwIsCanonical`). For a **non-canonical** query k-mer this flips
`hitFW`, tagging the seed reverse. The aligner then aligns the read's reverse
complement at the wrong position, scores it below threshold, and drops it. A
read survives only if it has *other* seeds whose correct orientation wins chain
selection — so the bug was invisible except for reads whose support reduced to a
**single non-canonical k-mer** (common for short reads with one error near an
end).

This was introduced with the **SSHash index refactor (salmon ≥ 1.11.0)**;
pre-1.11 (BooPHF) pufferfish derives orientation via `isEquivalent` and is not
affected. The Rust port is unaffected because piscem-rs derives k-mer
orientation directly from the query.

## Fix

pufferfish `5dce7f4` (branch `codex/for-salmon`): in the streaming overload,
`hitFW = (res.kmer_orientation == sshash::constants::forward_orientation)` —
i.e. drop the `== fwIsCanonical` correction (a no-op for canonical k-mers, a
fix for non-canonical ones). salmon's pufferfish pin is bumped to this commit.
A standalone release note is in `docs/release-note-orientation-fix.md`.

Validation: yeast 83.48% → 85.55% (matches the Rust port within 1 read);
`sample_data` 100% unchanged; `NumReads` correlation on confidently-mapped
transcripts unchanged.

## Related fix in this tree: `writeMappings` record loss

`--writeMappings`/`-z` wrapped the SAM output stream in an `ostream_sink_mt`
with `force_flush = false`; mapping records (written per batch) only reached
disk when the `ofstream` buffer overflowed, so the trailing buffer was lost when
the stream wasn't explicitly closed on the active quant path. For small inputs
this dropped **all** records (the large `@SQ` header overflowed the buffer and
masked it). Fixed by enabling `force_flush` on the sink
(`src/util/QuantOptionsUtils.cpp`). This fix is what made the placement
diagnosis above possible.

## Reproduction

Minimal: a 51 bp read equal to a transcript window with a single mismatch just
past the first 31 bp leaves exactly one error-free k-mer; if that k-mer is
non-canonical, pre-fix salmon places the read reverse at the wrong locus
(score 24/102) and post-fix places it forward at the correct locus (96/102).
`salmon-cli debug-map -i <index> -r <reads>` (Rust) reports per-read placement,
seed coverage, and alignment score for cross-checking.
