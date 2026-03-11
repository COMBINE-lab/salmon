#!/usr/bin/env python3
import argparse
import json
import re
import subprocess
import time
from pathlib import Path


ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
PROCESSED_RE = re.compile(r"processed\s+([0-9,]+)\s+fragments", re.IGNORECASE)


def strip_ansi(text: str) -> str:
    return ANSI_RE.sub("", text)


def run_window(args: argparse.Namespace) -> dict:
    out_dir = Path(args.out_dir)
    out_dir.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(args.salmon),
        "quant",
        "-lA",
        "-i",
        str(args.index),
        "-1",
        str(args.mates1),
        "-2",
        str(args.mates2),
        "--threads",
        str(args.threads),
        "-o",
        str(out_dir),
        "--adaptiveReadBatch",
    ]

    started = time.monotonic()
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        cwd=str(args.workdir) if args.workdir else None,
    )

    captured = ""
    timed_out = False
    try:
        out, _ = proc.communicate(timeout=args.seconds)
        captured = out or ""
    except subprocess.TimeoutExpired:
        timed_out = True
        proc.terminate()
        try:
            out, _ = proc.communicate(timeout=5)
            captured = out or ""
        except subprocess.TimeoutExpired:
            proc.kill()
            out, _ = proc.communicate()
            captured = out or ""
    elapsed = time.monotonic() - started

    plain = strip_ansi(captured)
    matches = PROCESSED_RE.findall(plain)
    processed = int(matches[-1].replace(",", "")) if matches else 0
    hits_per_frag = None
    hpf_match = re.findall(r"hits per frag:\s*([0-9.]+)", plain)
    if hpf_match:
        hits_per_frag = float(hpf_match[-1])

    return {
        "command": cmd,
        "timed_out": timed_out,
        "elapsed_seconds": elapsed,
        "processed_fragments": processed,
        "fragments_per_second": (processed / elapsed) if elapsed > 0 else 0.0,
        "hits_per_frag": hits_per_frag,
        "exit_code": proc.returncode,
    }


def main() -> int:
    p = argparse.ArgumentParser(description="Run a fixed-duration Salmon quant window and report throughput.")
    p.add_argument("--salmon", required=True, type=Path)
    p.add_argument("--index", required=True, type=Path)
    p.add_argument("--mates1", required=True, type=Path)
    p.add_argument("--mates2", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--seconds", type=int, default=30)
    p.add_argument("--workdir", type=Path, default=None)
    p.add_argument("--json-out", type=Path, default=None)
    args = p.parse_args()

    result = run_window(args)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
