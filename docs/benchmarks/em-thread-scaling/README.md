# Deterministic EM thread-scaling measurements

Raw tab-separated results for the measurements summarized in
`docs/em-architecture-and-optimization.md`. All runs used the 26,135,185-
fragment SRR21186103 / GENCODE v49 workload and were collected on 2026-08-23.

- `local-final.tsv` and `newton-final.tsv`: five alternating baseline/candidate
  repetitions at 1/2/4/8/16/32/64 threads. Baseline is the measured
  `4c8b0a98` (final rebased `f45ec7a2`); candidate is measured `04c0a93b`
  (final rebased `27c19a71`).
- `*-config-64class-vs-64inc.tsv`: two repetitions at 16/32/64 threads. The
  `baseline` arm is compressed/class-balanced 64; `candidate` is
  compressed/incidence-balanced 64.
- `*-config-128-vs-256.tsv`: the same sweep, with incidence-balanced 128 as
  `baseline` and incidence-balanced 256 as `candidate`.
- `local-bootstrap.tsv`: two alternating repetitions with 32 bootstrap
  replicates at 16 and 32 threads. The candidate is measured `0f9bc742` (final
  rebased `2477cd28`); bootstrap executes the same sequential kernel as the
  selected 128-shard build.
- `perf-counters.tsv`: one full-process `perf stat` run per arm at 16/32/64
  threads on each host, recording generic cycles, instructions, and cache
  misses alongside the isolated `em_bias` phase time.

`quant.sf` files are not committed because each is about 38 MiB. Their SHA-256
values are retained in every final-result row. The candidate hash must be
unique per host across all thread counts; the document records the numerical
comparison and mass/convergence checks.
