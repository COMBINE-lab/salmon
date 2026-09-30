#!/usr/bin/env bash
# Public benchmark of `salmon index --unspliced` on total (rRNA-depleted) RNA.
#
# Builds three decoy-aware GENCODE indices from the same references:
#   gentrome   transcripts + the whole genome as decoys (the usual setup)
#   intron     gentrome + --unspliced intron (flank = read length - 1)
#   premrna    gentrome + --unspliced premrna
# quantifies every library against each, and records wall time and peak
# memory of every step. `bench_unspliced.py` then summarizes the outputs.
#
#   scripts/bench_unspliced.sh --salmon BIN --ref DIR --fastq DIR --out DIR \
#       [--threads N] [--read-length L] [--modes "gentrome intron premrna"]
#
# --ref must hold GENCODE's GRCh38.primary_assembly.genome.fa.gz,
# gencode.vNN.transcripts.fa.gz and gencode.vNN.primary_assembly.annotation.gtf.gz.
# --fastq holds <run>_1.fastq.gz / <run>_2.fastq.gz pairs; every pair is used.
# Extra quant options go after `--` (default: -l A --gcBias --seqBias).
#
# The run used in the PR: GENCODE v50, GSE100127 (ABRF rRNA-depletion study,
# Universal Human Reference RNA, technical duplicate libraries, 2x76):
#   SRR5693352 SRR5693353  (Ribo-Zero Gold, site 2)
#   SRR5693399 SRR5693400  (Ribo-Zero Gold, site 4)
set -euo pipefail

SALMON="" REF="" FASTQ="" OUT="" THREADS=16 READ_LEN=76
MODES="gentrome intron premrna"
QUANT_OPTS=(-l A --gcBias --seqBias)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --salmon) SALMON=$(cd "$(dirname "$2")" && pwd)/$(basename "$2"); shift 2 ;;
    --ref) REF=$(cd "$2" && pwd); shift 2 ;;
    --fastq) FASTQ=$(cd "$2" && pwd); shift 2 ;;
    --out) mkdir -p "$2"; OUT=$(cd "$2" && pwd); shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --read-length) READ_LEN="$2"; shift 2 ;;
    --modes) MODES="$2"; shift 2 ;;
    --) shift; QUANT_OPTS=("$@"); break ;;
    -h|--help) sed -n '2,24p' "$0"; exit 0 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
  esac
done
[[ -n "$SALMON" && -n "$REF" && -n "$FASTQ" && -n "$OUT" ]] || { sed -n '2,24p' "$0"; exit 2; }

GENOME=$(ls "$REF"/*primary_assembly.genome.fa.gz)
TXOME=$(ls "$REF"/gencode.v*.transcripts.fa.gz)
GTF=$(ls "$REF"/gencode.v*.primary_assembly.annotation.gtf.gz)
TIMES="$OUT/times.tsv"
[[ -s "$TIMES" ]] || printf "step\twall_s\tmax_rss_bytes\n" > "$TIMES"

# Run a command under the platform's time(1), appending wall time and peak
# RSS (bytes) to times.tsv.
timed() {
  local step=$1; shift
  local log="$OUT/logs/$step.log"
  mkdir -p "$OUT/logs"
  local start end rss
  start=$(date +%s)
  if [[ "$(uname)" == Darwin ]]; then
    /usr/bin/time -l "$@" > "$log" 2>&1
    rss=$(awk '/maximum resident set size/ {print $1}' "$log")
  else
    /usr/bin/time -v "$@" > "$log" 2>&1
    rss=$(awk -F: '/Maximum resident set size/ {print $2 * 1024}' "$log")
  fi
  end=$(date +%s)
  printf "%s\t%s\t%s\n" "$step" "$((end - start))" "$rss" >> "$TIMES"
}

# gentrome FASTA and decoy list (decoys last, as salmon requires)
if [[ ! -s "$OUT/gentrome.fa.gz" ]]; then
  gzip -dc "$GENOME" | awk '/^>/ {print substr($1, 2)}' > "$OUT/decoys.txt"
  cat "$TXOME" "$GENOME" > "$OUT/gentrome.fa.gz"
fi

for mode in $MODES; do
  idx="$OUT/index_$mode"
  [[ -s "$idx/info.json" ]] && continue
  extra=()
  case "$mode" in
    gentrome) ;;
    intron) extra=(--unspliced intron --genome "$GENOME" --gtf "$GTF" --readLength "$READ_LEN") ;;
    premrna) extra=(--unspliced premrna --genome "$GENOME" --gtf "$GTF") ;;
    *) echo "unknown mode $mode" >&2; exit 2 ;;
  esac
  timed "index_$mode" "$SALMON" index -t "$OUT/gentrome.fa.gz" -d "$OUT/decoys.txt" \
    --gencode -p "$THREADS" -i "$idx" ${extra[@]+"${extra[@]}"}
done

for r1 in "$FASTQ"/*_1.fastq.gz; do
  run=$(basename "$r1" _1.fastq.gz)
  r2="$FASTQ/${run}_2.fastq.gz"
  for mode in $MODES; do
    q="$OUT/quant/$mode/$run"
    [[ -s "$q/quant.sf" ]] && continue
    timed "quant_${mode}_$run" "$SALMON" quant -i "$OUT/index_$mode" -1 "$r1" -2 "$r2" \
      -p "$THREADS" -o "$q" "${QUANT_OPTS[@]}"
  done
done
echo "done: $OUT"
