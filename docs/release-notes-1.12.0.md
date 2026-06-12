# salmon 1.12.0 — release notes

**This is the final release of the C++ implementation of salmon.** Future
development continues in **salmon 2.0**, a from-scratch Rust rewrite that is
faster, easier to build and install, and that fixes the issues below natively.
The 1.12.0 line will remain available for reproducibility and emergency fixes
(see "Continuity" below), but new features and performance work land in 2.0.

1.12.0 is a **correctness-focused** release. Every fix here was discovered while
building and cross-validating the Rust rewrite against C++ salmon, and then
backported to C++ so existing 1.11.x users benefit immediately. Indices built
with salmon 1.11.x are **compatible** — no reindexing is required for 1.12.0.

---

## Correctness fixes

### 1. Selective alignment: k-mer orientation in the SSHash streaming lookup
The SSHash streaming k-mer lookup reused the canonical-relative orientation
formula for a query-relative result, flipping `hitFW` for **non-canonical** query
k-mers. Reads whose seed support reduced to a single k-mer (e.g. short reads with
one sequencing error) could be mis-placed on the wrong strand and discarded.

- **Impact:** yeast `ERR458493` (1.09M 51bp reads) mapping rate **83.48% →
  85.55%** (~22,600 reads recovered); clean `sample_data` unchanged at 100%;
  `NumReads` on confidently-mapped transcripts unchanged. Recovered reads are
  genuine matches, not new spurious mappings.
- **Affected versions:** SSHash-based releases only (salmon ≥ 1.11.0).
- pufferfish pinned to the streaming `getRefPos` orientation fix (`5dce7f4`).

### 2. `--seqBias`: observed sequence-bias model was never trained
`numBiasSamples` was zeroed in the pipeline context before option parsing applied
its default, so the bias-sampling gate never fired. The observed ("foreground")
model stayed uniform and the correction degenerated to `1 / P_expected(context)`,
inflating effective lengths (often beyond transcript length). Fixed for **both**
the reads path and the alignment-based (`-a`) path.

- **Impact:** yeast `ERR458493` — observed-bias samples collected **0 → ~430k**;
  transcripts with `effLen > length` **4,015/6,571 → 117/6,571**; bias-corrected
  effective lengths now match an independent re-implementation at Pearson
  **0.91 → 0.99943**. `--numBiasSamples` now also takes effect (it previously did
  not). Clean non-bias runs are unaffected.

### 3. Alignment mode: pair mates by mate fields, not stream adjacency
`BAMQueue` assumed paired-end records are adjacent in the BAM stream. Aligners
that emit all R1 records then all R2 records within a name-collated group
(e.g. `bowtie2 -k`/`-a`) violated this, causing salmon to pair read1-with-read1
(same-orientation) and silently drop fragments. Mates are now paired via
reciprocal `RNEXT`/`PNEXT` (disambiguated by the `HI` tag) using a small
read-ahead interleave buffer.

- **Impact:** a yeast `bowtie2 -k` BAM recovered **1,918,938 → 2,000,000**
  fragments (~81,062 previously dropped); interleaved sample BAM output
  bit-identical. RSEM shares the adjacency assumption but asserts on violation;
  salmon now handles it.

### 4. Positional bias: additive-smoothed factor + 3′ coordinate fix
Corrected the positional-bias factor to use additive smoothing and fixed the 3′
coordinate computation used when applying the positional model.

### 5. `--writeMappings`: always flush SAM records to disk
The SAM/BAM sink could retain buffered records at shutdown; output is now
force-flushed so all mapping records are written.

### 6. `BAMQueue`: do not pre-allocate a 2M-object pool at startup
Removed an unconditional ~2M-object pool allocation at process start, reducing
startup memory and latency for alignment-mode runs.

---

## Continuity

- **salmon 1.12.0 is the last C++ release.** Development moves to **salmon 2.0**
  (Rust): a single portable binary, `cargo`/conda installation, the same
  `quant.sf` output (drop-in for tximport), compatible inferential-replicate
  output (bootstrap/Gibbs), and a new alignment-free `--sketch` mode.
- The C++ source is preserved on the `cpp` branch; a `salmon-cpp` conda package
  may be provided for users who need the C++ line for reproducibility or
  emergency fixes.
- **salmon 2.0 indices must be rebuilt** (new index format); 1.12.0 indices are
  unchanged from 1.11.x.

## Acknowledgements

These fixes came out of the salmon Rust rewrite effort and the cross-validation
between the two implementations.
