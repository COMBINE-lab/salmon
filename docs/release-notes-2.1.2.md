# salmon 2.1.2

A bug-fix patch on top of 2.1.1.

1. **Fixed a panic during selective-alignment mapping** (issue
   [#1038](https://github.com/COMBINE-lab/salmon/issues/1038)) — on some indexes
   `salmon quant` aborted with
   `range start index N out of range for slice of length M` at
   `salmon-map/src/align.rs`. Root-caused to the seed collector in `piscem-rs`
   and fixed there
   ([piscem-rs v0.6.1](https://github.com/COMBINE-lab/piscem-rs/releases/tag/piscem-rs-v0.6.1));
   salmon now depends on `piscem-rs 0.6.1`.

---

## 1. Mapping panic on certain indexes (#1038)

### Symptom

With some reference indexes (reported on a multi-million-transcript Trinity
assembly), `salmon quant` panicked early in the mapping phase:

```
thread '<unnamed>' panicked at crates/salmon-map/src/align.rs:...:
range start index 259 out of range for slice of length 254
```

### Cause

`salmon`'s sketch/permissive seed collector (in `piscem-rs`) anchors a read
position to a contig with a hash lookup, then *skips ahead* along the contig and
verifies the continuation by direct sequence comparison rather than re-hashing.
The skip distance and the verified position are expressed in the **anchor's**
coordinate frame.

When recording the skipped hit, the collector advanced the *last recorded hit*
instead of the anchor. Those coincide in the common case, but when an earlier
skip/walk detour had already recorded a hit at a read position ahead of the
current one, the last recorded hit had a different (larger) contig offset.
Advancing it by the anchor's skip distance produced a contig position past the
end of the contig (`contig_pos + k > contig_len`) — a k-mer offset that does not
exist in the contig. The selective-alignment scorer then sliced the reference at
that out-of-range position and panicked. (C++ piscem has the same code but masks
the bad position with pufferfish's reference-boundary clamp; the Rust scorer
indexed it directly.)

### Fix

Fixed in `piscem-rs` (v0.6.1): the skip search now projects the hit from the
**open-search anchor**, so the recorded position equals the one actually verified
by sequence comparison and always stays within the contig (the skip distance is
bounded by the anchor's distance to the contig end). There is no change in the
common case; on the reproducing 2.86M-transcript Trinity index the mapping rate
is identical (88.40%) and the panic is gone.

As defense in depth, `align_chain` now `debug_assert`s that every seed lies
within the read and reference, and in release builds skips an out-of-bounds
candidate instead of allowing an out-of-bounds slice — so a future coordinate
bug degrades gracefully rather than aborting a run.
