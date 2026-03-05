# Benchmark Harness

This directory defines the benchmark contract for Salmon's modernization work.

- CI should run only small smoke-style benchmark inputs that are safe to commit.
- Maintainer signoff for performance-sensitive changes must use external real bulk RNA-seq datasets that are not stored in this repository.
- Required comparisons for signoff:
  - wall-clock runtime
  - peak RSS
  - output parity for supported bulk workflows

Current acceptance targets:

- no repeatable wall-clock regression greater than 5%
- no repeatable peak RSS regression greater than 10%

Changes that require benchmark reruns:

- FASTA/Q parser changes
- allocator changes
- compression backend changes
- alignment I/O changes
- major dependency upgrades on hot paths
