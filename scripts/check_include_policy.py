#!/usr/bin/env python3

import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

# We only guard implementation-facing code; compatibility shims are allowed to
# include legacy top-level paths while migration is in progress.
SCAN_DIRS = [
    REPO_ROOT / "src",
    REPO_ROOT / "tests",
    REPO_ROOT / "include" / "salmon" / "internal",
]

# Top-level vendor/provider headers that should no longer be included directly.
# Implementation code must include salmon/vendor/* wrappers so provider
# ownership and update policy are centralized.
FORBIDDEN_TOPLEVEL_VENDOR_HEADERS = {
    "concurrentqueue.h",
    "edlib.h",
    "kseq++.hpp",
    "kseq.h",
    "xxhash.h",
    "atomicops.h",
    "blockingconcurrentqueue.h",
    "bucket_container.hh",
    "cuckoohash_config.hh",
    "cuckoohash_map.hh",
    "cuckoohash_util.hh",
    "dbg.hpp",
    "default_hasher.hh",
    "ezETAProgressBar.hpp",
    "fastapprox.h",
    "format.h",
    "httplib.hpp",
    "json.hpp",
    "lightweightsemaphore.h",
    "make_unique.hpp",
    "pcg_extras.hpp",
    "pcg_random.hpp",
    "peglib.h",
    "posix.h",
    "readerwriterqueue.h",
    "spline.h",
    "strict_fstream.hpp",
}

# Legacy Salmon internal headers that were removed from top-level include/.
# Implementation code must include canonical salmon/internal/* paths instead.
FORBIDDEN_LEGACY_INTERNAL_HEADERS = {
    "BinaryLiteral.hpp",
    "CommonTypes.hpp",
    "EquivalenceClassBuilder.hpp",
    "FragmentList.hpp",
    "GenomicFeature.hpp",
    "IndexVersionInfo.hpp",
    "MicroOpts.hpp",
    "PCA.hpp",
    "ReadProducer.hpp",
    "SalmonConfig.hpp",
    "SalmonDefaults.hpp",
    "SalmonIndexVersionInfo.hpp",
    "SalmonOpts.hpp",
    "SufficientStatisticsQueue.hpp",
    "TranscriptCluster.hpp",
    "backtrace.hpp",
    "zstr.hpp",
    "BAMQueue.tpp",
    "BarcodeGroup.hpp",
}

INCLUDE_RE = re.compile(r'^\s*#\s*include\s*[<"]([^">]+)[">]')


def iter_source_files():
    suffixes = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx", ".tpp"}
    for scan_dir in SCAN_DIRS:
        if not scan_dir.exists():
            continue
        for p in scan_dir.rglob("*"):
            if p.is_file() and p.suffix in suffixes:
                yield p


def main() -> int:
    violations = []
    for path in iter_source_files():
        try:
            lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
        except OSError:
            continue
        for lineno, line in enumerate(lines, start=1):
            m = INCLUDE_RE.match(line)
            if not m:
                continue
            include_target = m.group(1).strip()
            if "/" in include_target:
                continue
            if (
                include_target in FORBIDDEN_TOPLEVEL_VENDOR_HEADERS
                or include_target in FORBIDDEN_LEGACY_INTERNAL_HEADERS
            ):
                violations.append((path, lineno, include_target))

    if not violations:
        print(
            "include-policy: OK (no forbidden top-level vendor or legacy internal includes found)"
        )
        return 0

    print("include-policy: found forbidden top-level vendor includes:")
    for path, lineno, include_target in violations:
        print(
            f"  {path.relative_to(REPO_ROOT)}:{lineno}: "
            f'#include "{include_target}"'
        )
    print(
        "\nUse salmon/vendor/* wrappers for vendored headers and "
        "salmon/internal/* paths for Salmon internal headers."
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
