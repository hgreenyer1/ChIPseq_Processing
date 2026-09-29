#!/usr/bin/env bash
# Step 3 (per library): merge the runs (lanes / re-sequencing) of a library and
# apply ENCODE filtering:
#   - drop unmapped, secondary, QC-fail reads and MAPQ < MAPQ_THRESHOLD
#   - paired-end: keep properly paired reads only, fix mate info, drop orphans
#   - mark duplicates (Picard), compute library complexity (NRF, PBC1, PBC2)
#   - remove duplicates (unless KEEP_DUPLICATES=true)
#
# usage: 03_filter.sh <library_id> <paired 0|1> <run_id> [run_id ...]
source "${CHIPSEQ_PIPELINE_DIR}/bin/common.sh"

LIB_ID=$1
PAIRED=$2
shift 2
RUN_IDS=("$@")
init_step filter "$LIB_ID"

BAM_WORK=$WORKDIR/bowtie2
OUT=$OUTDIR/bowtie2/merged_library
QC=$OUT/qc
mkdir -p "$OUT" "$QC"
T=$STEP_TMP
SORT_MEM=$(sort_mem_per_thread)

# ---- merge runs -------------------------------------------------------------
run_bams=()
for r in "${RUN_IDS[@]}"; do
    [[ -s $BAM_WORK/$r.sorted.bam ]] || die "missing aligned BAM for run $r"
    run_bams+=("$BAM_WORK/$r.sorted.bam")
done
merged=$T/$LIB_ID.merged.bam
if (( ${#run_bams[@]} > 1 )); then
    log "Merging ${#run_bams[@]} runs"
    samtools merge -@ "$CPUS" -f "$merged" "${run_bams[@]}"
else
    cp "${run_bams[0]}" "$merged"
fi
samtools index "$merged"
samtools flagstat -@ "$CPUS" "$merged" > "$QC/$LIB_ID.merged.flagstat"
total_reads=$(samtools view -c -F 0x900 "$merged")
mapped_reads=$(samtools view -c -F 0x904 "$merged")

# Optional chromosome removal (e.g. chrM) as a region list of the kept contigs
region_args=()
if [[ -n ${FILTER_CHRS// /} ]]; then
    # shellcheck disable=SC2086
    printf '%s\n' $FILTER_CHRS > "$T/drop_chrs.txt"
    samtools idxstats "$merged" | awk 'NR==FNR{drop[$1]=1; next} !($1 in drop) && $1!="*" {print $1"\t0\t"$2}' \
        "$T/drop_chrs.txt" - > "$T/keep.bed"
    region_args=(-L "$T/keep.bed")
fi

# ---- ENCODE filtering ------------------------------------------------------
filt=$T/$LIB_ID.filt.bam
if [[ $PAIRED == 1 ]]; then
    log "Filtering (paired-end): -F 1804 -f 2 -q $MAPQ_THRESHOLD, fixmate, drop orphans"
    samtools view -@ "$CPUS" -F 1804 -f 2 -q "$MAPQ_THRESHOLD" "${region_args[@]}" -u "$merged" \
        | samtools sort -n -@ "$CPUS" -m "$SORT_MEM" -T "$T/nsort" -o "$T/filt.nsort.bam" -
    samtools fixmate -@ "$CPUS" -r "$T/filt.nsort.bam" "$T/fixmate.bam"
    rm -f "$T/filt.nsort.bam"
    samtools view -@ "$CPUS" -F 1804 -f 2 -u "$T/fixmate.bam" \
        | samtools sort -@ "$CPUS" -m "$SORT_MEM" -T "$T/csort" -o "$filt" -
    rm -f "$T/fixmate.bam"
else
    log "Filtering (single-end): -F 1804 -q $MAPQ_THRESHOLD"
    samtools view -@ "$CPUS" -F 1804 -q "$MAPQ_THRESHOLD" "${region_args[@]}" -b -o "$filt" "$merged"
fi
rm -f "$merged" "$merged.bai"

# ---- mark duplicates -------------------------------------------------------
log "Picard MarkDuplicates"
dupmark=$T/$LIB_ID.dupmark.bam
# shellcheck disable=SC2086
$PICARD_CMD -Xmx"$(job_mem_75)" MarkDuplicates \
    INPUT="$filt" OUTPUT="$dupmark" METRICS_FILE="$QC/$LIB_ID.MarkDuplicates.metrics.txt" \
    VALIDATION_STRINGENCY=LENIENT ASSUME_SORT_ORDER=coordinate REMOVE_DUPLICATES=false \
    TMP_DIR="$T"
rm -f "$filt"
dup_pct=$(awk -F'\t' '/^## METRICS CLASS/ {getline; for (i=1;i<=NF;i++) if ($i=="PERCENT_DUPLICATION") c=i; getline; print $c; exit}' \
    "$QC/$LIB_ID.MarkDuplicates.metrics.txt")

# ---- library complexity (ENCODE: on filtered reads, before de-duplication) --
log "Library complexity"
# Input: count of reads per unique position. Output: TotalReadPairs DistinctReadPairs
# OneReadPair TwoReadPairs NRF PBC1 PBC2
pbc_awk='BEGIN{mt=0;m0=0;m1=0;m2=0}
    ($1==1){m1++} ($1==2){m2++} {m0++; mt+=$1}
    END{
        nrf  = (mt>0) ? sprintf("%.6f", m0/mt) : "NA";
        pbc1 = (m0>0) ? sprintf("%.6f", m1/m0) : "NA";
        pbc2 = (m2>0) ? sprintf("%.6f", m1/m2) : "inf";
        printf "%d\t%d\t%d\t%d\t%s\t%s\t%s\n", mt, m0, m1, m2, nrf, pbc1, pbc2
    }'
pbc_file=$QC/$LIB_ID.pbc.qc
if [[ $PAIRED == 1 ]]; then
    samtools sort -n -@ "$CPUS" -m "$SORT_MEM" -T "$T/pbcsort" -o "$T/dupmark.nsort.bam" "$dupmark"
    bedtools bamtobed -bedpe -i "$T/dupmark.nsort.bam" 2> /dev/null \
        | awk 'BEGIN{OFS="\t"} $1!="chrM" {print $1,$2,$4,$6,$9,$10}' \
        | sort -S "$(job_mem_75)" -T "$T" | uniq -c | awk "$pbc_awk" > "$pbc_file"
    rm -f "$T/dupmark.nsort.bam"
else
    bedtools bamtobed -i "$dupmark" \
        | awk 'BEGIN{OFS="\t"} $1!="chrM" {print $1,$2,$3,$6}' \
        | sort -S "$(job_mem_75)" -T "$T" | uniq -c | awk "$pbc_awk" > "$pbc_file"
fi
sed -i '1i TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF\tPBC1\tPBC2' "$pbc_file"
read -r _ _ _ _ NRF PBC1 PBC2 < <(tail -n 1 "$pbc_file")

# ---- final BAM ---------------------------------------------------------------
final=$(lib_bam "$LIB_ID")
flags=1804
is_true "$KEEP_DUPLICATES" && flags=780
if [[ $PAIRED == 1 ]]; then
    samtools view -@ "$CPUS" -F "$flags" -f 2 -b -o "$final" "$dupmark"
else
    samtools view -@ "$CPUS" -F "$flags" -b -o "$final" "$dupmark"
fi
rm -f "$dupmark"
samtools index "$final"
samtools flagstat -@ "$CPUS" "$final" > "$QC/$LIB_ID.filt.dedup.flagstat"
samtools idxstats "$final" > "$QC/$LIB_ID.filt.dedup.idxstats"
samtools stats -@ "$CPUS" "$final" > "$QC/$LIB_ID.filt.dedup.stats"
final_reads=$(samtools view -c "$final")

if is_true "$CLEAN_INTERMEDIATES"; then
    for b in "${run_bams[@]}"; do rm -f "$b" "$b.bai"; done
fi

mapped_pct=$(awk -v m="$mapped_reads" -v t="$total_reads" 'BEGIN{printf "%.2f", (t>0)?100*m/t:0}')
write_metrics "$METRICS_DIR/libraries/$LIB_ID/1_filter.tsv" \
    paired "$PAIRED" \
    total_reads "$total_reads" \
    mapped_reads "$mapped_reads" \
    mapped_pct "$mapped_pct" \
    dup_fraction "$dup_pct" \
    final_reads "$final_reads" \
    NRF "$NRF" PBC1 "$PBC1" PBC2 "$PBC2"

finish_step
