# Header Ownership Map

This file tracks migration of direct `include/*` headers into clear module ownership under `include/salmon/internal/*`.

## Policy

- New project headers must not be added directly under `include/`.
- Every project header has one owning module.
- Keep compatibility shims at old paths during migration:
  - old header path remains,
  - old header only includes new canonical header,
  - add a short `TODO` with planned shim removal milestone.
- Third-party headers remain outside module ownership and are tracked separately.

## Top-Level Project Headers

| Current header | Owner module | Target canonical path | Status |
|---|---|---|---|
| `EquivalenceClassBuilder.hpp` | `quant` | `include/salmon/internal/quant/EquivalenceClassBuilder.hpp` | migrated (top-level shim removed) |
| `GenomicFeature.hpp` | `model` | `include/salmon/internal/model/GenomicFeature.hpp` | migrated (top-level shim removed) |
| `TranscriptCluster.hpp` | `quant` | `include/salmon/internal/quant/TranscriptCluster.hpp` | migrated (top-level shim removed) |
| `SalmonOpts.hpp` | `config` | `include/salmon/internal/config/SalmonOpts.hpp` | migrated (top-level shim removed) |
| `SalmonDefaults.hpp` | `config` | `include/salmon/internal/config/SalmonDefaults.hpp` | migrated (top-level shim removed) |
| `SalmonConfig.hpp` | `config` | `include/salmon/internal/config/SalmonConfig.hpp` | migrated (top-level shim removed) |
| `SalmonIndexVersionInfo.hpp` | `index` | `include/salmon/internal/index/SalmonIndexVersionInfo.hpp` | migrated (top-level shim removed) |
| `IndexVersionInfo.hpp` | `index` | `include/salmon/internal/index/IndexVersionInfo.hpp` | migrated (top-level shim removed) |
| `SufficientStatisticsQueue.hpp` | `quant` | `include/salmon/internal/quant/SufficientStatisticsQueue.hpp` | migrated (top-level shim removed) |
| `ReadProducer.hpp` | `io` | `include/salmon/internal/io/ReadProducer.hpp` | migrated (top-level shim removed) |
| `FragmentList.hpp` | `io` | `include/salmon/internal/io/FragmentList.hpp` | migrated (top-level shim removed) |
| `PCA.hpp` | `model` | `include/salmon/internal/model/PCA.hpp` | migrated (top-level shim removed) |
| `MicroOpts.hpp` | `cli` | `include/salmon/internal/cli/MicroOpts.hpp` | migrated (top-level shim removed) |
| `CommonTypes.hpp` | `util` | `include/salmon/internal/util/CommonTypes.hpp` | migrated (top-level shim removed) |
| `BinaryLiteral.hpp` | `util` | `include/salmon/internal/util/BinaryLiteral.hpp` | migrated (top-level shim removed) |
| `backtrace.hpp` | `util` | `include/salmon/internal/util/backtrace.hpp` | migrated (top-level shim removed) |
| `BAMQueue.tpp` | `alignment` | `include/salmon/internal/alignment/BAMQueue.tpp` | migrated (top-level shim removed) |
| `zstr.hpp` | `output` | `include/salmon/internal/output/zstr.hpp` | migrated (top-level shim removed) |

## Third-Party / Vendored Headers at Top-Level

These remain non-module-owned for now and should be tracked under vendoring policy:

- `edlib.h`, `xxhash.h`, `json.hpp`, `httplib.hpp`, `pcg_random.hpp`, `pcg_extras.hpp`, `kseq*.h*`, `format.h`, `spline.h`, `peglib.h`, `strict_fstream.hpp`, `atomicops.h`, `concurrentqueue.h`, `readerwriterqueue.h`, `blockingconcurrentqueue.h`, `lightweightsemaphore.h`, `cuckoohash_*`, `bucket_container.hh`, `default_hasher.hh`, `dbg.hpp`, `make_unique.hpp`, `fastapprox.h`, `posix.h`, `ezETAProgressBar.hpp`.

## Execution Order

1. `config + index` (high include fanout, biggest consistency win).
2. `quant + alignment` core types.
3. `model + io + util` cleanup.
4. Remove obsolete shims after all internal includes are rewritten and one release cycle passes.

## Mechanical Checklist Per Header

1. Add canonical header at target path.
2. Convert old top-level header into forwarding shim.
3. Rewrite includes in `include/`, `src/`, and `tests/` to canonical path.
4. Build and run tests.
5. Mark status here (`migrated` or `shim-only`).
