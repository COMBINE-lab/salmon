# salmon 2.3.4

A bug-fix release on 2.3.3. It fixes a single-end handling bug in
**alignment mode** (`salmon quant -a`). Paired-end workflows, mapping-based
(reads) quantification, and single-end runs with `--libType U` are unaffected
and produce identical results. No index rebuild is required.

## Fix: single-end alignment inputs are dropped under stranded library types (#1057)

In alignment mode, salmon classified each lone aligned record using only the
BAM `0x40` (first-in-pair) flag, and never consulted `0x1` ("template has
multiple segments" / paired). A **genuine single-end read** sets none of `0x1`,
`0x40`, or `0x80`, so it was mislabeled as a paired-end *right orphan* rather
than single-end. Under a stranded single-end library type (`--libType SF` or
`SR`) that pseudo-orphan classification can never satisfy the single-end
strandedness check, so **every alignment was discarded** and salmon reported
`0 fragments mapped`. `--libType U` was the only setting that worked, because it
accepts either orientation unconditionally.

The fix ([#1058](https://github.com/COMBINE-lab/salmon/pull/1058)) keys the
single-end vs. orphan decision on the `0x1` flag, per the SAM/BAM
specification: a record with `0x1` unset is treated as genuine single-end, and
its own alignment strand (`0x10`) is what the strand filter uses. A record that
*is* part of a pair but reported alone (a true orphan) is still classified
left/right by `0x40`, unchanged.

After the fix, single-end alignment inputs quantify correctly under all three
single-end library types:

- **`--libType SR`** accepts reverse-strand reads (BAM flag `16`),
- **`--libType SF`** accepts forward-strand reads (flag `0`),
- **`--libType U`** accepts both (unchanged).

This matches the behavior of the mapping-based (reads) path and the internal
genome-projection path, both of which already classified single-end reads
correctly.

## Upgrading

Drop-in from 2.3.3: no index rebuild. Results are unchanged for paired-end
inputs, reads-mode quantification, and single-end `--libType U`. If you feed a
single-end transcriptomic BAM/SAM with a stranded library type (`SF`/`SR`),
salmon now quantifies it correctly instead of reporting `0 fragments mapped`.
