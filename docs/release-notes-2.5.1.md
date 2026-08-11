# salmon 2.5.1

A dependency and hardening release on 2.5.0. **Quantification results are
unchanged** — verified on real data (26.1 M read pairs): mapped counts
identical in every mapping mode, and `quant.sf` byte-identical between the
serial and parallel decoders under `--deterministic`. **No index rebuild is
required.**

## The worker pool is now bounded by construction

2.5.0's adaptive threading resizes the mapping pool while a run is in flight.
The pool it shipped on (`paraseq-temp` 0.5.0-pre.1) spawned workers on demand
when running workers observed a shortfall — and upstream review of that design
found the observation racy: on small inputs, the worker count could
stochastically exceed the configured ceiling (nine concurrent workers in an
eight-worker pool, deterministically reproducible at 200 records).

This release moves to `paraseq-temp` 0.5.0-pre.2 via piscem-rs 0.9.2: upstream
paraseq's released 0.5.0 plus the redesigned pool offered upstream as
[noamteyssier/paraseq#78](https://github.com/noamteyssier/paraseq/pull/78).
Every worker the pool can run is spawned once with a stable index; a worker
over the target parks at a batch boundary and is woken by growth. Nothing
spawns after startup, so the thread ceiling cannot be exceeded — the bound is
structural, not enforced by counting. A parked worker drops its record
buffers, so batch memory is proportional to workers actually running.

Salmon's resize-safety test was reworked to match: exactly-once record
delivery is still asserted under continuous resizing, and the test now
requires each shrink to be *observed* taking effect, so it cannot pass
vacuously against a pool that ignores resize requests.

## Also

- piscem-rs 0.9.1 → 0.9.2 (carries the pool floor; 0.9.1's rapidgzip-core
  0.3.1 deadlock fix is retained).
- salmon's own manifest now declares the pre.2 pool directly, rather than
  inheriting the floor transitively.

Measured on the full evaluation set at `-p 64` with the parallel decoder:
selective alignment 58 s, sketch 35 s, both `--deterministic` variants within
the same band — all wall times within the 2.5.0 range.
