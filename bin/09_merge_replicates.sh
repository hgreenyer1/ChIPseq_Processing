#!/usr/bin/env bash
# Merged replicates (per sample with >= 2 replicates, IP or control):
# pools the filtered, de-duplicated BAMs of all replicates of a sample (ENCODE
# style pooling, no second duplicate removal: reads from different libraries
# are not PCR duplicates of each other) and makes a coverage bigWig.
# Peaks on the merged BAM are called by 05_call_peaks.sh ... replicate.
#
# usage: 09_merge_replicates.sh <sample> <paired 0|1> <library> <library> [library ...]
source "${CHIPSEQ_PIPELINE_DIR}/bin/common.sh"

SAMPLE=$1
PAIRED=$2
shift 2
LIBS=("$@")
init_step merge "$SAMPLE"

OUT=$OUTDIR/bowtie2/merged_replicate
QC=$OUT/qc
BW_OUT=$OUTDIR/bigwig/merged_replicate
mkdir -p "$QC" "$BW_OUT"

bams=() fl_sum=0 fl_n=0
for lib in "${LIBS[@]}"; do
    b=$(lib_bam "$lib")
    [[ -s $b ]] || die "missing filtered BAM $b"
    bams+=("$b")
    fl=$(lib_fraglen "$lib")
    if [[ -n $fl ]]; then fl_sum=$(( fl_sum + fl )); fl_n=$(( fl_n + 1 )); fi
done

log "Merging ${#LIBS[@]} replicates: ${LIBS[*]}"
bam=$(merged_bam "$SAMPLE")
samtools merge -@ "$CPUS" -f "$bam" "${bams[@]}"
samtools index "$bam"
samtools flagstat -@ "$CPUS" "$bam" > "$QC/$SAMPLE.mRp.flagstat"
samtools idxstats "$bam" > "$QC/$SAMPLE.mRp.idxstats"
samtools stats -@ "$CPUS" "$bam" > "$QC/$SAMPLE.mRp.stats"
reads=$(samtools view -c "$bam")

# Fragment length of the pool: mean of the replicate SPP estimates (ENCODE)
fraglen=""
rm -f "$QC/$SAMPLE.fraglen.txt"
if (( fl_n > 0 )); then
    fraglen=$(( fl_sum / fl_n ))
    echo "$fraglen" > "$QC/$SAMPLE.fraglen.txt"
fi

if [[ $PAIRED == 1 ]]; then
    extend=(--extendReads)
else
    extend=(--extendReads "${fraglen:-$DEFAULT_FRAGLEN}")
fi
log "bamCoverage ($BIGWIG_NORM)"
bamCoverage -b "$bam" -o "$BW_OUT/$SAMPLE.$BIGWIG_NORM.bigWig" \
    --binSize "$BIGWIG_BINSIZE" --normalizeUsing "$BIGWIG_NORM" \
    --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" "${extend[@]}" \
    --numberOfProcessors "$CPUS"

write_metrics "$METRICS_DIR/merged/$SAMPLE/1_merge.tsv" \
    replicates "$(IFS=,; echo "${LIBS[*]}")" \
    paired "$PAIRED" \
    reads "$reads" \
    fraglen "${fraglen:-NA}"

finish_step
