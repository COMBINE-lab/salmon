#!/usr/bin/env python3

import re
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]

CODE_SCAN_DIRS = [
    REPO_ROOT / "src",
    REPO_ROOT / "tests",
    REPO_ROOT / "include" / "salmon" / "internal",
]

CMAKE_SCAN_FILES = [
    REPO_ROOT / "CMakeLists.txt",
    REPO_ROOT / "src" / "CMakeLists.txt",
]
CMAKE_SCAN_FILES.extend((REPO_ROOT / "cmake").rglob("*.cmake"))
CMAKE_SCAN_FILES.extend(REPO_ROOT.rglob("CMakeLists.txt"))

SKIP_PATH_PARTS = {
    "build",
    "external",
    ".git",
    ".codex",
}

ALLOWED_ALIGNMENT_ADAPTER_FILES = {
    REPO_ROOT / "include" / "salmon" / "internal" / "io" / "AlignmentIO.hpp",
    REPO_ROOT / "src" / "io" / "AlignmentIO.cpp",
}

LEGACY_CMAKE_KNOBS = (
    "FETCH_BOOST",
    "BOOST_RECONFIGURE",
    "BOOST_WILL_RECONFIGURE",
    "TBB_RECONFIGURE",
    "TBB_WILL_RECONFIGURE",
    "TBB_TARGET_EXISTED",
    "DO_QUIET_MAKE",
    "NO_RTM",
    "NO_IPO",
    "CONDA_BUILD",
)

FASTXPARSER_INCLUDE_RE = re.compile(r'^\s*#\s*include\s*[<"]FastxParser\.hpp[">]')
FASTXPARSER_PATH_INCLUDE_RE = re.compile(r'^\s*#\s*include\s*[<"]include/FastxParser\.hpp[">]')
SCRAM_SYMBOL_RE = re.compile(r"\bscram_[A-Za-z0-9_]*\b|\bSAM_hdr\b")
RAW_ALIGNMENT_FIELD_RE = re.compile(r"->raw\b")
RAW_ALIGNMENT_HEADER_TEXT_RE = re.compile(r"\bsam_hdr_str\s*\(")
HTSLIB_INCLUDE_RE = re.compile(r'^\s*#\s*include\s*[<"]htslib/')


def iter_source_files():
    suffixes = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx", ".tpp", ".inl"}
    for scan_dir in CODE_SCAN_DIRS:
        if not scan_dir.exists():
            continue
        for p in scan_dir.rglob("*"):
            if p.is_file() and p.suffix in suffixes:
                yield p


def scan_legacy_cmake_knobs():
    violations = []
    patterns = [re.compile(rf"\b{re.escape(knob)}\b") for knob in LEGACY_CMAKE_KNOBS]
    for path in CMAKE_SCAN_FILES:
        if not path.exists() or not path.is_file():
            continue
        if any(part in SKIP_PATH_PARTS for part in path.parts):
            continue
        try:
            lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
        except OSError:
            continue
        for lineno, line in enumerate(lines, start=1):
            for knob, pat in zip(LEGACY_CMAKE_KNOBS, patterns):
                if pat.search(line):
                    violations.append((path, lineno, f"legacy CMake knob `{knob}`"))
                    break
    return violations


def scan_forbidden_includes_and_symbols():
    violations = []
    fastx_adapter = REPO_ROOT / "include" / "salmon" / "internal" / "io" / "FastxReader.hpp"
    for path in iter_source_files():
        try:
            lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
        except OSError:
            continue

        for lineno, line in enumerate(lines, start=1):
            if FASTXPARSER_INCLUDE_RE.match(line) and path != fastx_adapter:
                violations.append((path, lineno, "direct include of FastxParser.hpp"))
            if FASTXPARSER_PATH_INCLUDE_RE.match(line) and path != fastx_adapter:
                violations.append((path, lineno, "direct include of include/FastxParser.hpp"))

            if path not in ALLOWED_ALIGNMENT_ADAPTER_FILES and SCRAM_SYMBOL_RE.search(line):
                violations.append((path, lineno, "raw alignment legacy symbol (scram_/SAM_hdr) outside adapter"))
            if path not in ALLOWED_ALIGNMENT_ADAPTER_FILES and RAW_ALIGNMENT_FIELD_RE.search(line):
                violations.append((path, lineno, "direct AlignmentIO raw-field access outside adapter"))
            if path not in ALLOWED_ALIGNMENT_ADAPTER_FILES and RAW_ALIGNMENT_HEADER_TEXT_RE.search(line):
                violations.append((path, lineno, "direct sam_hdr_str usage outside adapter"))
            if path not in ALLOWED_ALIGNMENT_ADAPTER_FILES and HTSLIB_INCLUDE_RE.match(line):
                violations.append((path, lineno, "direct htslib include outside adapter"))

    return violations


def main() -> int:
    violations = []
    violations.extend(scan_legacy_cmake_knobs())
    violations.extend(scan_forbidden_includes_and_symbols())

    if not violations:
        print("refactor-policy: OK")
        return 0

    print("refactor-policy: violations found:")
    for path, lineno, reason in violations:
        print(f"  {path.relative_to(REPO_ROOT)}:{lineno}: {reason}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
