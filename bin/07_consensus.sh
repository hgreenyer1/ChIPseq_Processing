#!/usr/bin/env bash
# Step 7 (per antibody): nf-core style consensus peak set. Main (blacklist
# filtered) peaks from every IP library of the antibody are merged, and reads
# are counted in the merged peaks with featureCounts for downstream
# differential binding (e.g. DESeq2).
#
# usage: 07_consensus.sh <antibody> <narrow|broad> <all_paired 0|1> <ip_library> [ip_library ...]
source "${CHIPSEQ_PIPELINE_DIR}/bin/common.sh"

ANTIBODY=$1
PEAK_TYPE=$2
ALL_PAIRED=$3
shift 3
LIBS=("$@")
init_step consensus "$ANTIBODY"

OUT=$OUTDIR/macs2/$PEAK_TYPE/consensus/$ANTIBODY
mkdir -p "$OUT"
prefix=$OUT/$ANTIBODY.consensus_peaks

bams=()
for lib in "${LIBS[@]}"; do
    p=$(main_peaks "$lib" "$PEAK_TYPE")
    [[ -f $p ]] || die "missing peaks $p"
    awk -v l="$lib" 'BEGIN{OFS="\t"} {print $1,$2,$3,l}' "$p"
done | sort -k1,1 -k2,2n -S "$(job_mem_75)" -T "$STEP_TMP" > "$STEP_TMP/all_peaks.bed"
for lib in "${LIBS[@]}"; do bams+=("$(lib_bam "$lib")"); done

# chr start end libraries n_libraries
bedtools merge -i "$STEP_TMP/all_peaks.bed" -c 4,4 -o distinct,count_distinct \
    | awk -v m="$CONSENSUS_MIN_LIBS" 'BEGIN{OFS="\t"} $5>=m' > "$prefix.bed"
n_cons=$(wc -l < "$prefix.bed")
log "$n_cons consensus peaks from ${#LIBS[@]} libraries"

# Presence/absence matrix (for UpSet plots etc.)
libs_csv=$(IFS=,; echo "${LIBS[*]}")
awk -v libs="$libs_csv" 'BEGIN{
        OFS="\t"; n=split(libs, L, ",");
        printf "interval\tchr\tstart\tend\tn_libraries"; for (i=1;i<=n;i++) printf "\t%s", L[i]; print ""
    }
    {
        delete has; split($4, s, ","); for (k in s) has[s[k]]=1;
        printf "%s:%d-%d\t%s\t%d\t%d\t%d", $1, $2+1, $3, $1, $2, $3, $5
        for (i=1;i<=n;i++) printf "\t%d", (L[i] in has);
        print ""
    }' "$prefix.bed" > "$prefix.boolean.tsv"

# featureCounts on the consensus peaks
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"} {print $1":"$2+1"-"$3,$1,$2+1,$3,"+"}' \
    "$prefix.bed" > "$prefix.saf"
pe_args=()
if [[ $ALL_PAIRED == 1 ]]; then
    # subread >= 2.0.2 needs --countReadPairs to count fragments
    fc_help=$(featureCounts 2>&1 || true)
    if [[ $fc_help == *--countReadPairs* ]]; then
        pe_args=(-p --countReadPairs)
    else
        pe_args=(-p)
    fi
fi
featureCounts -F SAF -O --fracOverlap 0.2 -T "$CPUS" "${pe_args[@]}" --tmpDir "$STEP_TMP" \
    -a "$prefix.saf" -o "$prefix.featureCounts.txt" "${bams[@]}"

finish_step
