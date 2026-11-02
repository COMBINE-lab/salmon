#!/usr/bin/env bash
#
# quant_ena_stream.sh — quantify a public study from ENA without keeping its FASTQ.
#
# Each run's FASTQ are read straight from ENA through process substitution, so
# they never land on disk: only the output directories do (plus, while a run is
# being quantified, salmon's intermediate mappings inside its output directory).
#
# A run counts as done only when salmon exited cleanly AND it processed exactly
# the number of reads ENA announces. A stream can end early without anything
# noticing (a plain FASTQ has no envelope to say it was cut), so the count is
# what makes a result trustworthy. Done runs get a `.ok` marker and are skipped
# when the script is run again, which is also how to retry the ones that failed.
#
# USAGE
#
#   quant_ena_stream.sh -i INDEX -o OUTDIR (-a ACCESSION | -t REPORT.tsv)
#                       [-j JOBS] [-p THREADS] [-l LIBTYPE] [-- SALMON_QUANT_ARGS...]
#
#   -a  study or run accession ENA knows (PRJNA..., SRP..., SRR...); GEO series
#       (GSE...) must be translated to one of these first
#   -t  an ENA filereport saved earlier, with the columns run_accession,
#       read_count and fastq_ftp, in any order
#   -j  runs quantified at once (default 4)
#   -p  threads per run (default 3)
#   -l  library type (default A, detected from the reads)
#   --  anything else for `salmon quant`, e.g. -- --seqBias --gcBias --posBias
#
# The salmon binary is $SALMON (default: salmon). It must be able to read FIFOs,
# which salmon 2.10.0 and earlier cannot ("Illegal seek").
#
# A run that is not plainly single or paired (ENA lists X.fastq.gz, X_1 and X_2
# when it has orphans) is quantified from X_1 and X_2 only, and then fails the
# count check, since read_count includes the orphans. That is on purpose: it is
# for a person to decide what that run should mean.

set -u

usage() {
    sed -n '/^# USAGE/,/^# The salmon binary/p' "$0" | sed -e 's/^# \{0,1\}//' -e '$d' >&2
    exit 2
}

INDEX="" OUTDIR="" ACCESSION="" REPORT="" JOBS=4 THREADS=3 LIBTYPE=A
while getopts "i:o:a:t:j:p:l:h" opt; do
    case "$opt" in
    i) INDEX=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    a) ACCESSION=$OPTARG ;;
    t) REPORT=$OPTARG ;;
    j) JOBS=$OPTARG ;;
    p) THREADS=$OPTARG ;;
    l) LIBTYPE=$OPTARG ;;
    *) usage ;;
    esac
done
shift $((OPTIND - 1))
[ -n "$INDEX" ] && [ -n "$OUTDIR" ] || usage
if { [ -n "$ACCESSION" ] && [ -n "$REPORT" ]; } || { [ -z "$ACCESSION" ] && [ -z "$REPORT" ]; }; then
    usage
fi

SALMON=${SALMON:-salmon}
EXTRA_ARGS=""
[ "$#" -gt 0 ] && EXTRA_ARGS=$(printf '%q ' "$@")
mkdir -p "$OUTDIR" || exit 1

if [ -n "$ACCESSION" ]; then
    REPORT=$OUTDIR/ena_report.tsv
    curl -sSf --retry 5 --retry-all-errors -o "$REPORT" \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$ACCESSION&result=read_run&fields=run_accession,read_count,fastq_ftp&format=tsv" ||
        {
            echo "could not fetch the ENA report for $ACCESSION" >&2
            exit 1
        }
fi

# No `--retry` on the streams: curl cannot take back what it already wrote to a
# pipe, so a retry would feed the same reads twice. A stalled or dropped transfer
# fails the run instead, and running the script again resumes.
stream() {
    case "$1" in
    *://*) curl -sSf --speed-limit 1024 --speed-time 120 "$1" ;;
    *) curl -sSf --speed-limit 1024 --speed-time 120 "https://$1" ;;
    esac
}

one() {
    local run=$1 expected=$2 url1=$3 url2=${4:-} dir="$OUTDIR/$1" log="$OUTDIR/$1.log" extra processed
    [ -f "$dir/.ok" ] && return 0
    eval "extra=($EXTRA_ARGS)"
    case "$expected" in
    '' | *[!0-9]*)
        echo "FAILED $run: ENA announces no read count, nothing to check against"
        return 1
        ;;
    esac
    # Salmon reuses an output directory; a partial earlier attempt must not count.
    rm -f "$dir/aux_info/meta_info.json"
    : >"$log"
    if [ -n "$url2" ]; then
        "$SALMON" quant -i "$INDEX" -l "$LIBTYPE" -p "$THREADS" ${extra[@]+"${extra[@]}"} -o "$dir" \
            -1 <(stream "$url1" 2>>"$log") -2 <(stream "$url2" 2>>"$log") >>"$log" 2>&1
    else
        "$SALMON" quant -i "$INDEX" -l "$LIBTYPE" -p "$THREADS" ${extra[@]+"${extra[@]}"} -o "$dir" \
            -r <(stream "$url1" 2>>"$log") >>"$log" 2>&1
    fi
    local status=$?
    if [ "$status" -ne 0 ]; then
        echo "FAILED $run: salmon exited with $status, see $log"
        return 1
    fi
    processed=$(sed -n 's/.*"num_processed": *\([0-9]*\).*/\1/p' "$dir/aux_info/meta_info.json" 2>/dev/null)
    if [ "$processed" != "$expected" ]; then
        echo "FAILED $run: read ${processed:-0} reads, ENA announces $expected"
        return 1
    fi
    touch "$dir/.ok"
    echo "ok $run $processed"
}
export -f one stream
export INDEX OUTDIR JOBS THREADS LIBTYPE SALMON EXTRA_ARGS

# One "run count url1 url2" line per run, with the columns found by name. An empty
# count is kept as "-" so that it fails on its own line rather than shifting the
# fields; `one` then refuses it as not being a number.
awk -F'\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) col[$i] = i
        if (!col["run_accession"] || !col["read_count"] || !col["fastq_ftp"]) {
            print "the report needs run_accession, read_count and fastq_ftp columns" > "/dev/stderr"
            exit 2
        }
        next
    }
    $col["fastq_ftp"] == "" {
        print "FAILED " $col["run_accession"] ": ENA lists no FASTQ" > "/dev/stderr"
        failed = 1
        next
    }
    {
        n = split($col["fastq_ftp"], urls, ";")
        r1 = urls[1]
        r2 = ""
        if (n == 2) r2 = urls[2]
        if (n > 2) for (i = 1; i <= n; i++) {
            if (urls[i] ~ /_1\.f(ast)?q\.gz$/) r1 = urls[i]
            if (urls[i] ~ /_2\.f(ast)?q\.gz$/) r2 = urls[i]
        }
        count = $col["read_count"] == "" ? "-" : $col["read_count"]
        # No trailing blank for a single-end run: xargs -L would read it as a line
        # continuation and glue the next run onto this one.
        print $col["run_accession"], count, r1 (r2 == "" ? "" : " " r2)
    }
    END { exit failed }
' "$REPORT" >"$OUTDIR/.runs" || report_status=$?

xargs -P "$JOBS" -L 1 bash -c 'one "$@"' _ <"$OUTDIR/.runs"
jobs_status=$?
rm -f "$OUTDIR/.runs"

done_runs=0
for marker in "$OUTDIR"/*/.ok; do
    [ -f "$marker" ] && done_runs=$((done_runs + 1))
done
total=$(($(wc -l <"$REPORT") - 1))
echo "$done_runs of $total runs done in $OUTDIR"
[ "$jobs_status" -eq 0 ] && [ "${report_status:-0}" -eq 0 ]
