# Tests

This repository now targets bulk RNA-seq only.

Current automated coverage lives in:

- `unitTests` for C++ unit coverage
- `salmon_read_test_quasi` CTest integration smoke test for index+quant
- `salmon_alevin_stub` CTest migration test (`salmon alevin` must fail and print breadcrumbs)
- `tests/test_quant.nf` for bulk quantification integration coverage

Single-cell `alevin` test assets are no longer part of the supported test surface.
