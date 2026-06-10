# Mapping-sensitivity fix: k-mer orientation in the SSHash streaming lookup

## Summary

We fixed a bug in pufferfish's SSHash-backed k-mer lookup that caused a subset
of reads to be mis-placed on the wrong strand and subsequently discarded. The
fix measurably improves mapping sensitivity on real short-read data with **no
change on clean data** and no change to the quantification of confidently-mapped
transcripts.

## Impact

On a standard real-data benchmark (Gierliński/Schurch *S. cerevisiae*,
`ERR458493`, 1,093,957 single-end 51 bp reads; Ensembl R64-1-1 cDNA):

| | mapping rate |
| --- | --- |
| before fix | 83.48% (913,271 reads) |
| after fix | **85.55% (935,851 reads)** |

That is ~22,600 additional correctly-mapped reads (+2.1 percentage points). On
clean simulated data (the bundled `sample_data`) the mapping rate is unchanged
at 100%, and per-transcript `NumReads` for confidently-mapped transcripts are
unchanged — the recovered reads are genuine matches that were previously lost,
not new spurious mappings.

## Symptom

Affected reads were ones whose seed support reduced to a **single k-mer** (e.g.
short reads with a sequencing error close to one end, leaving one error-free
k-mer window). Such a read could be reported as **unmapped** or mapped to the
**wrong strand/locus** with a poor alignment score, even when it was an exact or
near-exact match to a transcript. Reads with several seeds were unaffected,
because their correctly-oriented seeds dominated chain selection — which is why
the bug was easy to miss in aggregate.

## Root cause

pufferfish derives, for each k-mer hit, whether the read k-mer lies **forward**
on the unitig (`hitFW`). There are two lookup styles:

- **Non-streaming** lookups query the *canonical* form of the k-mer
  (`dict_.lookup(mer.getCanonicalWord())`), so the returned orientation is
  relative to the canonical word and must be converted to the query's frame:
  `hitFW = (orientation == forward) == fwIsCanonical`.
- The **streaming** lookup queries the *raw query k-mer string*
  (`sq.lookup(kmer_str)`), so its orientation is **already** relative to the
  query k-mer: `hitFW = (orientation == forward)`.

The streaming `getRefPos` overload mistakenly reused the **canonical-relative**
formula (`… == fwIsCanonical`). For a **non-canonical** query k-mer
(`fwIsCanonical == false`) this **flips** `hitFW`, so the read's seed is tagged
reverse. The aligner then aligns the read's reverse complement at the wrong
position, scores it below threshold, and drops it.

## Fix

In the streaming `getRefPos` overload, take the orientation directly:

```cpp
// streaming lookup orientation is query-relative
bool hitFW = (res.kmer_orientation == sshash::constants::forward_orientation);
```

For canonical query k-mers this is identical to the previous expression
(`x == true`), so only the mis-handled non-canonical case changes. The two
non-streaming overloads (canonical-word lookups) are correct and unchanged.

pufferfish commit: `5dce7f4` on branch `codex/for-salmon`.

## Affected versions

**SSHash-based releases only (salmon ≥ 1.11.0).** The bug lives in the SSHash
streaming-query path, which was introduced with the SSHash index refactor in
1.11.0. Earlier releases (BooPHF-based pufferfish, ≤ 1.10.x) derive orientation
via `CanonicalKmer::isEquivalent`, which is query-relative and correct, and are
**not** affected.

## Validation

- Yeast `ERR458493`: 83.48% → 85.55% mapped; matches an independent mapper to
  within 1 read.
- `sample_data` (clean, simulated): 100% mapped, unchanged.
- `NumReads` correlation on confidently-mapped transcripts unchanged
  (Pearson r ≈ 0.999).
- Minimal reproduction: a 51 bp read equal to a transcript window with a single
  mismatch just past the first 31 bp (so only one error-free k-mer remains, and
  that k-mer is non-canonical) was reported reverse-strand at the wrong locus
  (score 24/102) before the fix and forward at the correct locus
  (score 96/102) after.
