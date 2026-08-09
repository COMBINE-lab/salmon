#!/usr/bin/env bash
#
# decode_crossover.sh — find where the parallel gzip decoder starts paying, in
# each of salmon's two mapping modes.
#
# WHY THIS EXISTS
#
# `-p` is one budget shared between inflating the input and mapping it. The
# serial decoder inflates inline on the mapping threads and is therefore
# work-conserving: a thread that just inflated a batch maps it and is never
# idle. It gives F concurrent inflate streams for free (F = gzip inputs). The
# parallel decoder gives d streams but *takes* d mapping threads, so it pays
# exactly when d > F.
#
# The broker converges to d* = N · P/(P+C) (producer/consumer busy time), so
#
#     d* > F   <=>   N > F · (1 + C/P)
#
# i.e. the engagement threshold is 1 + C/P and scales linearly with per-fragment
# consumer cost. Selective alignment does far more work per fragment than sketch
# at identical decode cost, so it must engage later. This script measures *how
# much* later instead of guessing, because the arithmetic says a 4x heavier
# consumer moves the threshold from 8 to ~29 — far outside the 12-16 range one
# would pick by eye.
#
# WHAT IT DOES
#
# For each mode x thread count, runs the same input twice — serial decoder and
# parallel decoder — and reports parallel/serial wall-time ratio. The crossover
# is the smallest -p at which the ratio drops below 1. Dividing that by F gives
# min_threads_per_stream for that mode.
#
# METHOD NOTES (learned the hard way in piscem)
#
#   * Runs are counterbalanced (serial,parallel then parallel,serial) and the
#     better of each pair is used: this machine showed up to ~40% run-to-run
#     spread on short inputs.
#   * Any run under MIN_SECONDS is reported but excluded from the verdict. The
#     whole point of adaptive scheduling is long runs; a 20-second run is
#     dominated by startup and converges too late to matter.
#   * CPU/wall is checked against -p on every run. Three separate invalid
#     measurements in the piscem work would have been caught by that one ratio,
#     so a run that overspends its budget is flagged, not quietly averaged in.
#
# USAGE
#   scripts/decode_crossover.sh -i <index> -1 R1.fq.gz -2 R2.fq.gz [-o outdir]
#   THREADS="8 16 32 64" scripts/decode_crossover.sh ...
set -uo pipefail

SALMON=${SALMON:-target/release/salmon}
THREADS=${THREADS:-"8 16 32 64"}
MODES=${MODES:-"sketch sa"}
MIN_SECONDS=${MIN_SECONDS:-60}
OUT=${OUT:-crossover-out}
IDX=""; R1=""; R2=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) IDX="$2"; shift 2;;
    -1) R1="$2"; shift 2;;
    -2) R2="$2"; shift 2;;
    -o) OUT="$2"; shift 2;;
    *) echo "unknown argument: $1" >&2; exit 2;;
  esac
done
[[ -n "$IDX" && -n "$R1" && -n "$R2" ]] || { sed -n '/^# USAGE/,/^set /p' "$0"; exit 2; }
[[ -x "$SALMON" ]] || { echo "no salmon binary at $SALMON (build --release first)" >&2; exit 2; }

mkdir -p "$OUT"
TSV="$OUT/raw.tsv"
: > "$TSV"
printf 'mode\tthreads\tdecoder\trep\twall\tcpu\tavg_threads\tover_budget\n' >> "$TSV"

# F: gzip inputs. Paired input is two files, both of which must be inflated.
F=2

run_one() { # mode threads decoder rep
  local mode=$1 t=$2 dec=$3 rep=$4
  local o="$OUT/$mode-t$t-$dec-r$rep"
  rm -rf "$o"; mkdir -p "$o"
  local extra=(); [[ "$mode" == "sketch" ]] && extra+=(--sketch)
  /usr/bin/time -f "%e %U" -o "$o/.time" \
    "$SALMON" quant -i "$IDX" -l A -1 "$R1" -2 "$R2" -p "$t" \
      --decoder "$dec" "${extra[@]}" -o "$o" >"$o/out.log" 2>"$o/err.log"
  local rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "  FAILED mode=$mode -p=$t decoder=$dec rep=$rep (see $o/err.log)" >&2
    printf '%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\n' "$mode" "$t" "$dec" "$rep" >> "$TSV"
    return 1
  fi
  read -r wall cpu < "$o/.time"
  local avg over
  avg=$(python3 -c "print(f'{$cpu/max($wall,1e-9):.1f}')")
  over=$(python3 -c "print('YES' if $cpu/max($wall,1e-9) > $t*1.05 else 'no')")
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$mode" "$t" "$dec" "$rep" "$wall" "$cpu" "$avg" "$over" >> "$TSV"
  [[ "$over" == "YES" ]] && echo "  WARNING: mode=$mode -p=$t $dec used $avg threads against a budget of $t" >&2
  echo "  $mode -p=$t $dec rep$rep: wall=${wall}s cpu=${cpu}s (${avg} avg threads)"
  return 0
}

for mode in $MODES; do
  for t in $THREADS; do
    echo "== mode=$mode -p=$t =="
    # Counterbalanced: order within each pair is reversed on the second rep so a
    # warming page cache or a thermal ramp cannot favour one arm systematically.
    run_one "$mode" "$t" serial   1
    run_one "$mode" "$t" parallel 1
    run_one "$mode" "$t" parallel 2
    run_one "$mode" "$t" serial   2
  done
done

echo
echo "=== crossover summary (ratio = best parallel wall / best serial wall) ==="
python3 - "$TSV" "$MIN_SECONDS" "$F" <<'PY'
import sys, csv, collections
tsv, min_s, F = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
best = collections.defaultdict(dict)
with open(tsv) as fh:
    for r in csv.DictReader(fh, delimiter='\t'):
        if r['wall'] == 'NA':
            continue
        k, w = (r['mode'], int(r['threads'])), float(r['wall'])
        best[k][r['decoder']] = min(best[k].get(r['decoder'], 1e18), w)

print(f"{'mode':<8}{'-p':>5}{'serial':>10}{'parallel':>10}{'ratio':>8}  verdict")
print('-'*56)
crossover = {}
for (mode, t) in sorted(best, key=lambda x: (x[0], x[1])):
    d = best[(mode, t)]
    if 'serial' not in d or 'parallel' not in d:
        continue
    s, p = d['serial'], d['parallel']
    ratio = p / s
    short = ' (SHORT - excluded)' if max(s, p) < min_s else ''
    verdict = 'parallel wins' if ratio < 1 else 'serial wins'
    print(f"{mode:<8}{t:>5}{s:>10.1f}{p:>10.1f}{ratio:>8.3f}  {verdict}{short}")
    if not short and ratio < 1 and mode not in crossover:
        crossover[mode] = t

print()
for mode in ('sketch', 'sa'):
    if mode in crossover:
        t = crossover[mode]
        print(f"  {mode}: parallel first wins at -p {t}  ->  min_threads_per_stream = {t}/{F} = {t//F}")
    else:
        print(f"  {mode}: no crossover in the swept range - threshold is above max(-p)")
print()
print("  Set these in decode.rs (ENGAGEMENT_SKETCH / ENGAGEMENT_SELECTIVE_ALIGNMENT).")
print("  A constant fitted to the EDGE of the sweep is unproven: if the crossover")
print("  lands on the largest -p tested, widen the range before believing it.")
PY
