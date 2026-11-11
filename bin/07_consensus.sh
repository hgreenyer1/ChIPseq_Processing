#!/usr/bin/env bash
# Step 7 (per antibody): nf-core style consensus peak set. Main (blacklist
# filtered) peaks from every IP library of the antibody are merged, and reads
# are counted in the merged peaks with featureCounts for downstream
# differential binding (e.g. DESeq2).
#
# usage: 07_consensus.sh <antibody> <narrow|broad> <all_paired 0|1> <ip_library> [ip_library ...]
#   all_paired=0 means at least one library is single-end; read types are then
#   detected per library
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
# subread >= 2.0.2 needs --countReadPairs to count fragments
pe_args=(-p)
fc_help=$(featureCounts 2>&1 || true)
[[ $fc_help == *--countReadPairs* ]] && pe_args+=(--countReadPairs)

# featureCounts cannot mix single-end and paired-end BAMs in one run: count
# each kind separately and join the columns back in library order
is_paired_bam() {
    local fs=$OUTDIR/bowtie2/merged_library/qc/$1.filt.dedup.flagstat n
    if [[ -s $fs ]]; then
        n=$(awk '/paired in sequencing/ {print $1}' "$fs")
    else
        n=$(samtools view -c -f 1 "$(lib_bam "$1")")
    fi
    (( ${n:-0} > 0 ))
}
pe_bams=() se_bams=()
for lib in "${LIBS[@]}"; do
    if [[ $ALL_PAIRED == 1 ]] || is_paired_bam "$lib"; then
        pe_bams+=("$(lib_bam "$lib")")
    else
        se_bams+=("$(lib_bam "$lib")")
    fi
done
fc_common=(-F SAF -O --fracOverlap 0.2 -T "$CPUS" --tmpDir "$STEP_TMP" -a "$prefix.saf")

if (( ${#se_bams[@]} == 0 )); then
    featureCounts "${fc_common[@]}" "${pe_args[@]}" -o "$prefix.featureCounts.txt" "${pe_bams[@]}"
elif (( ${#pe_bams[@]} == 0 )); then
    featureCounts "${fc_common[@]}" -o "$prefix.featureCounts.txt" "${se_bams[@]}"
else
    log "Counting ${#pe_bams[@]} paired-end and ${#se_bams[@]} single-end libraries separately"
    featureCounts "${fc_common[@]}" "${pe_args[@]}" -o "$prefix.pe.featureCounts.txt" "${pe_bams[@]}"
    featureCounts "${fc_common[@]}" -o "$prefix.se.featureCounts.txt" "${se_bams[@]}"
    order=$(IFS=,; echo "${bams[*]}")
    paste "$prefix.pe.featureCounts.txt" "$prefix.se.featureCounts.txt" | awk -F'\t' -v order="$order" '
        BEGIN { OFS="\t"; n=split(order, o, ",") }
        NR==1 { print "# featureCounts: paired-end libraries counted as fragments, single-end as reads"; next }
        NR==2 { for (i=1; i<=NF; i++) if (!($i in col)) col[$i]=i }
        { out=$1; for (i=2; i<=6; i++) out=out OFS $i
          for (k=1; k<=n; k++) out=out OFS $(col[o[k]])
          print out }' > "$prefix.featureCounts.txt"
    rm -f "$prefix.pe.featureCounts.txt" "$prefix.se.featureCounts.txt"
fi

finish_step
