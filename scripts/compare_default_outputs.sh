#!/usr/bin/env bash
# Compare the outputs of two salmon binaries, byte for byte, on the bundled
# sample data, over the default index build and a set of quant configurations.
#
# Used to show that a change which adds opt-in behaviour leaves every default
# output untouched (e.g. `salmon index --unspliced`: without the option the
# index and all quant outputs must be identical to the base release).
#
#   scripts/compare_default_outputs.sh <base salmon binary> <new salmon binary>
#
# Indices are built with one thread: a multi-threaded build is not byte
# reproducible even between two runs of the same binary (the piscem tables'
# layout depends on thread scheduling), while quant results are.
#
# Both binaries run in sibling directories with identical relative paths, so
# the paths recorded in cmd_info.json / meta_info.json match. Only fields that
# measure the run itself are masked before comparing: wall-clock timestamps,
# durations, peak RSS, and the logs/ directory. Everything else, including
# every binary index file, must match exactly.
set -euo pipefail

BASE=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
NEW=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

tar -xzf "$REPO/sample_data.tgz" -C "$WORK"
S=../sample_data
LOG=$WORK/run.log
# A decoy-aware "gentrome": the transcripts, then a pseudo-chromosome made of
# every transcript reverse-complemented and joined by N runs, as the decoy.
python3 - "$WORK/sample_data" <<'PY'
import sys
d = sys.argv[1]
recs, name = [], None
for line in open(f"{d}/transcripts.fasta"):
    line = line.strip()
    if line.startswith(">"):
        recs.append([line[1:].split()[0], ""])
    elif line:
        recs[-1][1] += line
comp = str.maketrans("ACGTacgt", "TGCAtgca")
decoy = ("N" * 50).join(s.translate(comp)[::-1] for _, s in recs)
with open(f"{d}/gentrome.fa", "w") as f:
    for n, s in recs:
        f.write(f">{n}\n{s}\n")
    f.write(f">chrDecoy\n{decoy}\n")
open(f"{d}/decoys.txt", "w").write("chrDecoy\n")
PY

run_all() {
  local bin=$1 dir=$2
  mkdir -p "$dir" && cd "$dir"
  "$bin" index -t $S/transcripts.fasta -i idx -p 1 >>"$LOG" 2>&1
  "$bin" index -t $S/transcripts.fasta -i idx_k23_dups -k 23 --keepDuplicates -p 1 >>"$LOG" 2>&1
  "$bin" index -t $S/gentrome.fa -d $S/decoys.txt -i idx_decoy -p 1 >>"$LOG" 2>&1
  local q=(quant -i idx -l A -1 $S/reads_1.fastq -2 $S/reads_2.fastq -p 4)
  "$bin" "${q[@]}" -o q_default --dumpEq >>"$LOG" 2>&1
  "$bin" "${q[@]}" -o q_bias --seqBias --gcBias --posBias >>"$LOG" 2>&1
  "$bin" "${q[@]}" -o q_boot --numBootstraps 4 >>"$LOG" 2>&1
  "$bin" "${q[@]}" -o q_online --online >>"$LOG" 2>&1
  "$bin" "${q[@]}" -o q_sketch --sketch >>"$LOG" 2>&1
  "$bin" quant -i idx_k23_dups -l A -r $S/reads_1.fastq -p 4 -o q_single >>"$LOG" 2>&1
  "$bin" quant -i idx_decoy -l A -1 $S/reads_1.fastq -2 $S/reads_2.fastq -p 4 \
    -o q_decoy --gcBias >>"$LOG" 2>&1
  "$bin" quant -t $S/transcripts.fasta -l A -a $S/sample_alignments.bam -p 4 \
    -o q_align >>"$LOG" 2>&1
  cd - >/dev/null
}

run_all "$BASE" "$WORK/base"
run_all "$NEW" "$WORK/new"

# Mask the fields that measure the run rather than describe its result.
normalize() {
  python3 - "$1" <<'PY'
import json, sys
p = sys.argv[1]
d = json.load(open(p))
for k in ("start_time", "end_time", "total_time_seconds", "peak_rss_kb"):
    if k in d:
        d[k] = "<masked>"
print(json.dumps(d, indent=2, sort_keys=False))
PY
}

status=0
n=0
while IFS= read -r f; do
  rel=${f#"$WORK/base/"}
  case "$rel" in */logs/*) continue ;; esac
  n=$((n + 1))
  other="$WORK/new/$rel"
  if [[ ! -e "$other" ]]; then
    echo "MISSING in new: $rel"; status=1; continue
  fi
  if [[ "$rel" == *meta_info.json ]]; then
    if ! diff -q <(normalize "$f") <(normalize "$other") >/dev/null; then
      echo "DIFFERS: $rel"; diff <(normalize "$f") <(normalize "$other") | head -20; status=1
    fi
  elif ! cmp -s "$f" "$other"; then
    echo "DIFFERS: $rel"; status=1
  fi
done < <(find "$WORK/base" -type f | sort)

# A default run must not grow new files either.
while IFS= read -r f; do
  rel=${f#"$WORK/new/"}
  case "$rel" in */logs/*) continue ;; esac
  [[ -e "$WORK/base/$rel" ]] || { echo "EXTRA in new: $rel"; status=1; }
done < <(find "$WORK/new" -type f | sort)

if [[ $status -eq 0 ]]; then
  echo "identical: $n files compared (index x3, quant x8), logs and run timings masked"
fi
exit $status
