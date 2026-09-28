#!/usr/bin/env bash
# Show the SLURM state of the jobs from the latest (or a given) submission.
#
# usage: pipeline_status.sh <outdir> [jobs_<timestamp>.tsv]
set -euo pipefail

OUTDIR=${1:?usage: pipeline_status.sh <outdir> [jobs_<timestamp>.tsv]}
MANIFEST=${2:-$(ls -t "$OUTDIR"/pipeline_info/jobs_*.tsv 2> /dev/null | head -n 1)}
[[ -f $MANIFEST ]] || { echo "no job manifest found in $OUTDIR/pipeline_info" >&2; exit 1; }

echo "Jobs from $MANIFEST"
ids=$(awk -F'\t' 'NR>1 {print $1}' "$MANIFEST" | paste -sd, -)
[[ -n $ids ]] || { echo "manifest is empty"; exit 0; }

sacct -X -n -P -j "$ids" --format=JobID,State,Elapsed,MaxRSS 2> /dev/null > "${TMPDIR:-/tmp}/sacct.$$" || true
awk -F'\t' -v OFS='\t' '
    NR==FNR { split($0, a, "|"); state[a[1]]=a[2]; el[a[1]]=a[3]; next }
    FNR==1  { print "job_id", "step", "id", "state", "elapsed"; next }
    { s = ($1 in state) ? state[$1] : "UNKNOWN"; print $1, $2, $3, s, (($1 in el) ? el[$1] : "-"); count[s]++ }
    END { printf "\n"; for (s in count) printf "%-12s %d\n", s, count[s] > "/dev/stderr" }
' FS='|' "${TMPDIR:-/tmp}/sacct.$$" FS='\t' "$MANIFEST" | column -t -s $'\t'
rm -f "${TMPDIR:-/tmp}/sacct.$$"

echo
echo "Failed job logs:"
awk -F'\t' 'NR>1 {print $1"\t"$5}' "$MANIFEST" | while IFS=$'\t' read -r id logf; do
    st=$(sacct -X -n -P -j "$id" --format=State 2> /dev/null | head -n 1 || true)
    [[ $st == FAILED* || $st == TIMEOUT* || $st == OUT_OF_ME* || $st == CANCELLED* ]] && echo "  $st  $logf"
done || true
