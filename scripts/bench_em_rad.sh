#!/usr/bin/env bash
# Reproducible A/B harness for the RAD-input inference phase.
#
# The input RAD is created once and reused, so every arm sees identical packed
# equivalence classes and the measurement excludes mapping. Results are appended
# to results.tsv; quant.sf is retained for bitwise and numerical comparisons.

set -euo pipefail

usage() {
  sed -n '2,34p' "$0"
  exit "${1:-0}"
}

BASELINE=""
CANDIDATE=""
RAD=""
OUT=""
THREADS="1,2,4,8,16,32,64"
REPEATS=5
CPU_ORDER=""
EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --baseline) BASELINE="$2"; shift 2 ;;
    --candidate) CANDIDATE="$2"; shift 2 ;;
    --rad) RAD="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --repeats) REPEATS="$2"; shift 2 ;;
    --cpu-order) CPU_ORDER="$2"; shift 2 ;;
    -h|--help) usage 0 ;;
    --) shift; EXTRA=("$@"); break ;;
    *) echo "unknown argument: $1" >&2; usage 2 ;;
  esac
done

[[ -x "$BASELINE" ]] || { echo "baseline is not executable: $BASELINE" >&2; exit 2; }
[[ -x "$CANDIDATE" ]] || { echo "candidate is not executable: $CANDIDATE" >&2; exit 2; }
[[ -f "$RAD" ]] || { echo "RAD input does not exist: $RAD" >&2; exit 2; }
[[ -n "$OUT" ]] || { echo "--out is required" >&2; exit 2; }
[[ ! -e "$OUT" ]] || { echo "refusing to reuse output directory: $OUT" >&2; exit 2; }

mkdir -p "$OUT"
"$BASELINE" --version > "$OUT/baseline-version.txt"
"$CANDIDATE" --version > "$OUT/candidate-version.txt"
lscpu > "$OUT/lscpu.txt"
{
  echo "rad=$RAD"
  echo "threads=$THREADS"
  echo "repeats=$REPEATS"
  echo "cpu_order=$CPU_ORDER"
  printf 'extra='; printf ' %q' "${EXTRA[@]}"; echo
} > "$OUT/config.txt"

printf 'arm\tthreads\trepeat\torder\tem_s\trad_read_s\twall_s\tuser_s\tsys_s\tmax_rss_kb\tem_iters\tconverged\tquant_sha256\n' > "$OUT/results.tsv"

affinity_for() {
  local n="$1"
  [[ -n "$CPU_ORDER" ]] || return 0
  local cpus selected
  IFS=',' read -r -a cpus <<< "$CPU_ORDER"
  (( ${#cpus[@]} >= n )) || { echo "--cpu-order has fewer than $n CPUs" >&2; exit 2; }
  selected="${cpus[0]}"
  for ((i=1; i<n; i++)); do selected+=",${cpus[i]}"; done
  printf '%s' "$selected"
}

phase_seconds() {
  local phase="$1" log="$2"
  sed 's/\x1b\[[0-9;]*m//g' "$log" \
    | grep -oE "phase=\"${phase}\" elapsed_s=[0-9.]+" \
    | tail -n 1 \
    | sed 's/.*elapsed_s=//'
}

json_scalar() {
  local key="$1" file="$2"
  sed -nE "s/.*\"${key}\"[[:space:]]*:[[:space:]]*([^,}]+).*/\1/p" "$file" | tail -n 1
}

run_one() {
  local arm="$1" binary="$2" threads="$3" repeat="$4" order="$5"
  local tag="t${threads}-r${repeat}-o${order}-${arm}"
  local run_out="$OUT/$tag" log="$OUT/$tag.log" timing="$OUT/$tag.time"
  local affinity
  affinity="$(affinity_for "$threads")"
  local prefix=()
  [[ -z "$affinity" ]] || prefix=(taskset -c "$affinity")

  /usr/bin/time -f '%e\t%U\t%S\t%M' -o "$timing" \
    env RUST_LOG=salmon::timing=info "${prefix[@]}" \
    "$binary" quant --rad "$RAD" -l A -p "$threads" --sigDigits 9 \
    -o "$run_out" "${EXTRA[@]}" 2> "$log"

  local wall user sys rss em rad_read iters converged hash
  IFS=$'\t' read -r wall user sys rss < "$timing"
  em="$(phase_seconds em_bias "$log")"
  rad_read="$(phase_seconds rad_read "$log")"
  iters="$(json_scalar num_em_iterations "$run_out/aux_info/meta_info.json")"
  converged="$(json_scalar em_converged "$run_out/aux_info/meta_info.json")"
  hash="$(sha256sum "$run_out/quant.sf" | awk '{print $1}')"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$arm" "$threads" "$repeat" "$order" "$em" "$rad_read" "$wall" "$user" \
    "$sys" "$rss" "$iters" "$converged" "$hash" >> "$OUT/results.tsv"
}

IFS=',' read -r -a thread_counts <<< "$THREADS"

# Warm both instruction paths and the RAD page cache; warmups are deliberately
# outside results.tsv.
warm_threads="${thread_counts[0]}"
run_one baseline "$BASELINE" "$warm_threads" 0 0
run_one candidate "$CANDIDATE" "$warm_threads" 0 1
sed -i '/\t0\t[01]\t/d' "$OUT/results.tsv"

for threads in "${thread_counts[@]}"; do
  for ((repeat=1; repeat<=REPEATS; repeat++)); do
    if (( repeat % 2 == 1 )); then
      run_one baseline "$BASELINE" "$threads" "$repeat" 1
      run_one candidate "$CANDIDATE" "$threads" "$repeat" 2
    else
      run_one candidate "$CANDIDATE" "$threads" "$repeat" 1
      run_one baseline "$BASELINE" "$threads" "$repeat" 2
    fi
  done
done

echo "wrote $OUT/results.tsv"
