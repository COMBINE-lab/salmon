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
#   quant_ena_stream.sh -i INDEX -o OUTDIR (-a IDS [-a IDS...] | -t REPORT.tsv)
#                       [-n] [-j JOBS] [-p THREADS] [-l LIBTYPE] [-- SALMON_QUANT_ARGS...]
#
#   -a  what to quantify, comma-separated and repeatable, in any mix of:
#         GEO series or samples     GSE312440, GSM9345836
#         ENA/SRA/DDBJ projects     PRJNA..., SRP..., ERP..., DRP...
#         samples and experiments   SAMN..., SRS..., SRX...
#         runs                      SRR..., ERR..., DRR...
#       Every ID has to resolve, and every sample of a GEO series has to have a
#       run at ENA, or nothing is quantified: a study quietly missing samples is
#       worse than an error. The list is saved as OUTDIR/ena_report.tsv.
#   -t  an ENA filereport saved earlier, with the columns run_accession,
#       read_count and fastq_ftp, in any order
#   -n  stop after resolving: list the runs, quantify nothing
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

INDEX="" OUTDIR="" IDS="" REPORT="" DRY_RUN="" JOBS=4 THREADS=3 LIBTYPE=A
while getopts "i:o:a:t:j:p:l:nh" opt; do
    case "$opt" in
    i) INDEX=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    a) IDS="$IDS $(printf '%s' "$OPTARG" | tr ',' ' ')" ;;
    n) DRY_RUN=1 ;;
    t) REPORT=$OPTARG ;;
    j) JOBS=$OPTARG ;;
    p) THREADS=$OPTARG ;;
    l) LIBTYPE=$OPTARG ;;
    *) usage ;;
    esac
done
shift $((OPTIND - 1))
[ -n "$INDEX" ] && [ -n "$OUTDIR" ] || usage
IDS=$(printf '%s' "$IDS" | tr 'a-z' 'A-Z' | xargs)
if { [ -n "$IDS" ] && [ -n "$REPORT" ]; } || { [ -z "$IDS" ] && [ -z "$REPORT" ]; }; then
    usage
fi

SALMON=${SALMON:-salmon}
EXTRA_ARGS=""
[ "$#" -gt 0 ] && EXTRA_ARGS=$(printf '%q ' "$@")
mkdir -p "$OUTDIR" || exit 1

ENA_API=https://www.ebi.ac.uk/ena/portal/api
ENA_FIELDS=run_accession,read_count,fastq_ftp,sample_alias

# Fetch a list, whole or not at all. Captured before anything reads it, since in
# `curl | ...` a transfer dropped midway would pass for a shorter list, and so
# for a study with fewer samples.
fetch() {
    local body
    body=$(curl -sSf --retry 5 "$@") || return 1
    printf '%s\n' "$body"
}

# The runs of one accession ENA indexes directly, without the header line. An
# accession ENA does not know is not an error to its API, just an empty report.
ena_runs() {
    local body
    body=$(fetch "$ENA_API/filereport?accession=$1&result=read_run&fields=$ENA_FIELDS&format=tsv") || return 1
    printf '%s\n' "$body" | tail -n +2
}

# The runs of the GEO samples on stdin. ENA has no idea what a GSM is as an
# accession, but keeps it as the sample alias. Batched because a large series
# would not fit in one query.
ena_runs_of_gsm() {
    xargs -n 100 | while read -r batch; do
        query=$(printf '%s\n' $batch | sed 's/.*/sample_alias="&"/' | sed '$!s/$/ OR /' | tr -d '\n')
        body=$(fetch "$ENA_API/search" -d result=read_run -d limit=0 \
            -d "fields=$ENA_FIELDS" -d format=tsv --data-urlencode "query=$query") || exit 1
        printf '%s\n' "$body" | tail -n +2
    done
}

# The samples of a GEO series. Asking GEO rather than following its BioProject
# link gives the list to check ENA's answer against, and covers a SuperSeries,
# whose samples sit under the BioProjects of its SubSeries.
geo_samples() {
    local body
    body=$(fetch "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=$1&targ=self&form=text&view=brief") || return 1
    printf '%s\n' "$body" | sed -n 's/^!Series_sample_id = //p' | tr -d '\r'
}

# Fail, naming them, if any of the GSMs in $1 has no run in the report in $2.
check_gsm_found() {
    awk -F'\t' -v id="$3" '
        NR == FNR { found[$4]; next }
        !($0 in found) { missing = missing " " $0; n++ }
        END {
            if (n) print "no run at ENA for " n " sample(s) of " id ":" missing > "/dev/stderr"
            exit n > 0
        }' "$2" "$1"
}

resolve_ids() {
    local id runs samples
    runs=$OUTDIR/.resolve.runs
    samples=$OUTDIR/.resolve.samples
    printf 'run_accession\tread_count\tfastq_ftp\tsample_alias\n'
    for id in $IDS; do
        : >"$runs"
        case "$id" in
        *[!A-Z0-9_.-]*)
            echo "not an accession: $id" >&2
            return 1
            ;;
        GSE*)
            geo_samples "$id" >"$samples" || {
                echo "could not reach GEO for $id" >&2
                return 1
            }
            [ -s "$samples" ] || {
                echo "GEO has no samples for $id" >&2
                return 1
            }
            ena_runs_of_gsm <"$samples" >"$runs" || {
                echo "could not reach ENA for $id" >&2
                return 1
            }
            check_gsm_found "$samples" "$runs" "$id" || return 1
            ;;
        GSM*)
            printf '%s\n' "$id" | ena_runs_of_gsm >"$runs" || {
                echo "could not reach ENA for $id" >&2
                return 1
            }
            ;;
        *)
            ena_runs "$id" >"$runs" || {
                echo "ENA refused $id, or could not be reached" >&2
                return 1
            }
            ;;
        esac
        [ -s "$runs" ] || {
            echo "ENA has no runs for $id" >&2
            return 1
        }
        cat "$runs"
    done
}

if [ -n "$IDS" ]; then
    REPORT=$OUTDIR/ena_report.tsv
    # IDs can overlap (a study and one of its runs): a run is quantified once.
    resolve_ids >"$REPORT.tmp"
    resolved=$?
    rm -f "$OUTDIR/.resolve.runs" "$OUTDIR/.resolve.samples"
    if [ "$resolved" -ne 0 ]; then
        rm -f "$REPORT.tmp"
        exit 1
    fi
    awk -F'\t' '!seen[$1]++' "$REPORT.tmp" >"$REPORT"
    rm -f "$REPORT.tmp"
fi

if [ -n "$DRY_RUN" ]; then
    echo "$(($(wc -l <"$REPORT") - 1)) runs in $REPORT"
    cut -f1,2 "$REPORT"
    exit 0
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
