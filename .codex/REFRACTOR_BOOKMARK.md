# Refactor Checkpoint (2026-03-08)

Branch: `codex/develop-refactor`

## Completed in this checkpoint

- Continued module ownership cleanup under `src/{io,model,util,tests}`.
- Kept `AlignmentIO` as the sole htslib boundary layer with more wrapper logic moved out of headers.
- Continued FASTQ/FASTA modernization around FQFeeder:
  - switched call sites to `salmon_fqfeeder::ParserConfig`.
  - retained lightweight compatibility aliases in `FastxReader.hpp` required by `SAMWriter`/`SalmonMappingUtils`.
- Simplified dependency resolution toward modern CMake package-first behavior.
- Added dependency floor checks for:
  - `spdlog/fmt`
  - `parallel_hashmap`
  - `nlohmann/json`
  - `Eigen >= 3.4`

## Verified at checkpoint

- Configure: `cmake -S . -B build/codex-check`
- Build: `cmake --build build/codex-check -j8 --target salmon unitTests`
- Tests: `ctest --test-dir build/codex-check --output-on-failure -R "unit_tests|salmon_read_test_quasi|salmon_alevin_stub"`

## Resume From Here

1. Remove remaining Staden-era compatibility artifacts after confirming no remaining call paths.
2. Keep shrinking heavy headers (`ReadExperiment.hpp`, `AlignmentLibrary.hpp`) by moving non-template logic into `.cpp`/`.inl`.
3. Continue top-level include shim removal only where include-policy and unit/integration tests remain green.
4. Quant hot-path optimization pass (read-based mode):
   - parser chunking + queue contention
   - mapping-stage allocations and temporary buffers
   - log/progress overhead in worker loops
   - merge/reduce behavior across mini-batches
