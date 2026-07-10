#!/usr/bin/env bash
#
# profile.sh — reproducible salmon quant profiling harness (Phase 1).
#
# Builds the `profiling` cargo profile (release + line-table debug info, unstripped
# — see [profile.profiling] in the root Cargo.toml; the shipped release/dist binary
# is unaffected) and runs `salmon quant` under a sampling profiler on two workloads:
#   1. mapping-dominated : plain quant
#   2. inference-dominated: quant --numBootstraps 100
# Per-phase wall-clock (index_load / mapping / eff_length_collapse / em_bias /
# posterior / output) is emitted on the `salmon::timing` tracing target and shown
# inline; isolate it with RUST_LOG=salmon::timing=info.
#
# Two data tiers:
#   - smoke (default): the repo's sample_data.tgz — tiny, good for *relative* phase
#     attribution and correctness, too small for stable wall-clock deltas.
#   - giab           : a realistic human RNA-seq sample + full transcriptome index.
#     See "GIAB tier" below to fetch the inputs, then pass them in.
#
# Usage:
#   scripts/profile.sh                              # smoke tier, both workloads
#   scripts/profile.sh -p 8                         # pin thread count
#   scripts/profile.sh -i <index_dir> -1 R1.fq.gz -2 R2.fq.gz   # giab tier
#   PROFILER=none scripts/profile.sh               # just time it, no profiler
#
# Profiler: samply (default; `cargo install samply`) — no sudo, opens the Firefox
# profiler UI. On Linux you can set PROFILER=flamegraph (needs perf). PROFILER=none
# skips sampling and only reports the per-phase timings.
#
# ── GIAB tier: fetching realistic inputs ─────────────────────────────────────────
# Reads: a GIAB reference individual's RNA-seq. NA12878 (= HG001) LCL RNA-seq from
#   the GEUVADIS project is the most canonical/downloadable. Resolve a paired-end
#   run on the ENA browser (https://www.ebi.ac.uk/ena/browser) under study PRJEB3366
#   / ArrayExpress E-GEUV-1, filtering sample = NA12878, library_strategy = RNA-Seq,
#   then download the run's *_1.fastq.gz / *_2.fastq.gz (typically 2x75 bp, tens of
#   millions of pairs). VERIFY the exact ERR run accession on ENA at download time —
#   do not hard-code one. (Alternative: HG002/NA24385 RNA-seq via an SRA search.)
# Transcriptome: GENCODE human transcripts, e.g.
#   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/  ->  gencode.vXX.transcripts.fa.gz
#   (or Ensembl cDNA, ~194k transcripts). Build the index once:
#     salmon index -t gencode.vXX.transcripts.fa.gz -i target/gencode_index -p 8
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)"
INDEX_DIR=""
R1=""
R2=""
BOOTSTRAPS=100
PROFILER="${PROFILER:-samply}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -p) THREADS="$2"; shift 2 ;;
    -i) INDEX_DIR="$2"; shift 2 ;;
    -1) R1="$2"; shift 2 ;;
    -2) R2="$2"; shift 2 ;;
    --bootstraps) BOOTSTRAPS="$2"; shift 2 ;;
    -h|--help) sed -n '2,50p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done

echo "==> building profiling profile (release + debug symbols)"
cargo build --profile profiling -p salmon-cli
SALMON="$REPO_ROOT/target/profiling/salmon"

OUT="$REPO_ROOT/target/profile-run"
mkdir -p "$OUT"

# Resolve the data tier: explicit -i/-1/-2 => GIAB tier; otherwise unpack sample_data.
if [[ -z "$INDEX_DIR" ]]; then
  echo "==> smoke tier: unpacking sample_data.tgz and building index"
  tar -xzf sample_data.tgz -C "$OUT"
  INDEX_DIR="$OUT/sample_index"
  R1="$OUT/sample_data/reads_1.fastq"
  R2="$OUT/sample_data/reads_2.fastq"
  "$SALMON" index -t "$OUT/sample_data/transcripts.fasta" -i "$INDEX_DIR" >/dev/null
fi

# Wrap the quant invocation in the chosen profiler.
run() {
  local tag="$1"; shift
  echo "==> [$tag] threads=$THREADS  profiler=$PROFILER"
  case "$PROFILER" in
    samply)     samply record -o "$OUT/$tag.profile.json" -- "$@" ;;
    flamegraph) cargo flamegraph --profile profiling -o "$OUT/$tag.svg" -- "${@:2}" ;;
    none)       "$@" ;;
    *) echo "unknown PROFILER: $PROFILER" >&2; exit 2 ;;
  esac
}

COMMON=(-l A -i "$INDEX_DIR" -1 "$R1" -2 "$R2" -p "$THREADS")

echo "### workload 1: mapping-dominated"
run map "$SALMON" quant "${COMMON[@]}" -o "$OUT/out_map"

echo "### workload 2: inference-dominated (--numBootstraps $BOOTSTRAPS)"
run boot "$SALMON" quant "${COMMON[@]}" --numBootstraps "$BOOTSTRAPS" -o "$OUT/out_boot"

echo "==> done. Profiles + quant outputs under $OUT"
echo "    Per-phase timings above are on the 'salmon::timing' target"
echo "    (RUST_LOG=salmon::timing=info to isolate them)."
