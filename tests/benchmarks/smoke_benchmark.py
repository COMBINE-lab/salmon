#!/usr/bin/env python3

import argparse
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description="Run Salmon smoke benchmark on committed fixtures.")
    p.add_argument("--salmon", required=True, help="Path to salmon executable")
    p.add_argument("--repo-root", required=True, help="Path to salmon repository root")
    p.add_argument("--output", required=True, help="JSON output path")
    return p.parse_args()


def ensure_fixtures(repo_root: Path):
    sample_dir = repo_root / "sample_data"
    if sample_dir.exists():
        return sample_dir
    tarball = repo_root / "sample_data.tgz"
    if not tarball.exists():
        raise RuntimeError(f"Missing fixture archive: {tarball}")
    subprocess.run(["tar", "xzvf", str(tarball.name)], cwd=repo_root, check=True)
    if not sample_dir.exists():
        raise RuntimeError(f"Failed to materialize fixture directory: {sample_dir}")
    return sample_dir


def run_step(name: str, cmd, cwd: Path):
    ru_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    ru_after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if proc.returncode != 0:
        raise RuntimeError(
            f"{name} failed ({proc.returncode})\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )
    return {
        "name": name,
        "command": cmd,
        "wall_seconds": elapsed,
        "ru_maxrss_before": ru_before,
        "ru_maxrss_after": ru_after,
        "ru_maxrss_delta": max(0, ru_after - ru_before),
    }


def main():
    args = parse_args()
    salmon = Path(args.salmon).resolve()
    repo_root = Path(args.repo_root).resolve()
    output = Path(args.output).resolve()

    if not salmon.exists():
        raise RuntimeError(f"Missing salmon executable: {salmon}")

    fixture_dir = ensure_fixtures(repo_root)
    index_dir = fixture_dir / "sample_benchmark_index"
    quant_dir = fixture_dir / "sample_benchmark_quant"
    if index_dir.exists():
        shutil.rmtree(index_dir)
    if quant_dir.exists():
        shutil.rmtree(quant_dir)

    steps = []
    steps.append(
        run_step(
            "index",
            [str(salmon), "index", "-t", "transcripts.fasta", "-i", str(index_dir.name)],
            fixture_dir,
        )
    )
    steps.append(
        run_step(
            "quant",
            [
                str(salmon),
                "quant",
                "-l",
                "A",
                "-i",
                str(index_dir.name),
                "-1",
                "reads_1.fastq",
                "-2",
                "reads_2.fastq",
                "-o",
                str(quant_dir.name),
            ],
            fixture_dir,
        )
    )

    result = {
        "host": {
            "platform": platform.platform(),
            "python": sys.version,
            "cwd": str(repo_root),
        },
        "acceptance_targets": {
            "wall_clock_regression_max": 0.05,
            "peak_rss_regression_max": 0.10,
        },
        "notes": {
            "ru_maxrss_units": "bytes on macOS, KiB on Linux",
            "scope": "CI smoke benchmark only; maintainer signoff requires external bulk datasets",
        },
        "steps": steps,
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2))
    print(f"smoke benchmark complete: {output}")
    for s in steps:
        print(
            f"{s['name']}: wall={s['wall_seconds']:.3f}s "
            f"ru_maxrss_after={s['ru_maxrss_after']} "
            f"ru_delta={s['ru_maxrss_delta']}"
        )


if __name__ == "__main__":
    main()
