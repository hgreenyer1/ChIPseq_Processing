#!/usr/bin/env bash
# Step 5 (per IP library, and per IP sample for merged replicates):
#   - main MACS2 call (nf-core): narrowPeak (q < MACS_QVAL) or broadPeak
#     (--broad-cutoff MACS_BROAD_CUTOFF), blacklist filtered, FRiP
#   - relaxed MACS2 call (ENCODE, p < RELAXED_PVAL, top RELAXED_NPEAKS) used for
#     replicate reproducibility, plus fold-enrichment and -log10(p) bigWigs
#   - deepTools plotFingerprint of IP vs control
#
# usage: 05_call_peaks.sh <id> <control> <paired 0|1> <narrow|broad> [library|replicate]
#   library   (default): <id> is a library, <control> a control library or "none"
#   replicate          : <id> is a sample whose replicates were merged by
#                        09_merge_replicates.sh; <control> is "none",
#                        "merged:<control_sample>" or a comma separated list of
#                        control libraries (pooled here)
source "${CHIPSEQ_PIPELINE_DIR}/bin/common.sh"

ID=$1
CTRL_SPEC=$2
PAIRED=$3
PEAK_TYPE=$4
LEVEL=${5:-library}

if [[ $LEVEL == "replicate" ]]; then
    init_step mpeaks "$ID"
    IP_BAM=$(merged_bam "$ID")
    FRAGLEN=$(merged_fraglen "$ID")
    BASE=$OUTDIR/macs2/merged_replicate
    FP_OUT=$OUTDIR/deeptools/plotFingerprint/merged_replicate
    METRICS_FILE=$METRICS_DIR/merged/$ID/2_peaks.tsv
else
    init_step peaks "$ID"
    IP_BAM=$(lib_bam "$ID")
    FRAGLEN=$(lib_fraglen "$ID")
    BASE=$OUTDIR/macs2
    FP_OUT=$OUTDIR/deeptools/plotFingerprint
    METRICS_FILE=$METRICS_DIR/libraries/$ID/3_peaks.tsv
fi
[[ -s $IP_BAM ]] || die "missing IP BAM $IP_BAM"

# ---- control ---------------------------------------------------------------------
CTRL_BAM="" CTRL_LABEL="none"
case $CTRL_SPEC in
    none) ;;
    merged:*)
        CTRL_LABEL=${CTRL_SPEC#merged:}
        CTRL_BAM=$(merged_bam "$CTRL_LABEL")
        ;;
    *,*)
        # controls of the replicates span several samples: pool them here
        IFS=, read -ra ctrl_libs <<< "$CTRL_SPEC"
        ctrl_bams=()
        for c in "${ctrl_libs[@]}"; do ctrl_bams+=("$(lib_bam "$c")"); done
        log "Pooling controls ${ctrl_libs[*]}"
        CTRL_BAM=$STEP_TMP/pooled_control.bam
        samtools merge -@ "$CPUS" -f "$CTRL_BAM" "${ctrl_bams[@]}"
        samtools index "$CTRL_BAM"
        CTRL_LABEL="pooled_control"
        ;;
    *)
        CTRL_LABEL=$CTRL_SPEC
        CTRL_BAM=$(lib_bam "$CTRL_SPEC")
        ;;
esac
[[ -z $CTRL_BAM || -s $CTRL_BAM ]] || die "missing control BAM $CTRL_BAM"

# MACS inputs (handles an IP and control of different read types)
macs_inputs "$IP_BAM" "$CTRL_BAM" "$PAIRED" "$FRAGLEN"
ctrl_args=()
[[ -n $MACS_CTRL ]] && ctrl_args=(-c "$MACS_CTRL")
log "MACS2 input: ${MACS_FMT[*]}, control: $CTRL_LABEL"

# ---- main peak call ------------------------------------------------------------
MAIN_OUT=$BASE/$PEAK_TYPE
mkdir -p "$MAIN_OUT/qc"
if [[ $PEAK_TYPE == "broad" ]]; then
    type_args=(--broad --broad-cutoff "$MACS_BROAD_CUTOFF")
else
    type_args=(-q "$MACS_QVAL")
fi
log "MACS2 main ${PEAK_TYPE} peak call"
# shellcheck disable=SC2086
"$MACS_CMD" callpeak -t "$MACS_TREAT" "${ctrl_args[@]}" "${MACS_FMT[@]}" -g "$MACS_GSIZE" \
    -n "$ID" --outdir "$MAIN_OUT" --tempdir "$STEP_TMP" --keep-dup all \
    "${type_args[@]}" $MACS_EXTRA_ARGS
raw_peaks=$MAIN_OUT/${ID}_peaks.${PEAK_TYPE}Peak
peaks=$MAIN_OUT/${ID}_peaks.bfilt.${PEAK_TYPE}Peak
blacklist_filter "$raw_peaks" "$peaks"
n_peaks=$(wc -l < "$peaks")
log "$n_peaks ${PEAK_TYPE} peaks after blacklist filtering"

# FRiP: fraction of (filtered, de-duplicated) reads falling in peaks
cut -f1-3 "$peaks" > "$STEP_TMP/peaks.bed"
total=$(samtools view -c -F 4 "$IP_BAM")
if (( n_peaks > 0 )); then
    in_peaks=$(samtools view -c -F 4 -L "$STEP_TMP/peaks.bed" "$IP_BAM")
else
    in_peaks=0
fi
frip=$(awk -v a="$in_peaks" -v b="$total" 'BEGIN{printf "%.4f", (b>0)?a/b:0}')
printf 'id\treads\treads_in_peaks\tFRiP\n%s\t%s\t%s\t%s\n' "$ID" "$total" "$in_peaks" "$frip" \
    > "$MAIN_OUT/qc/$ID.FRiP.txt"

# ---- relaxed peaks + signal tracks (ENCODE) --------------------------------------
REL_OUT=$BASE/relaxed/$ID
log "MACS2 relaxed call (p < $RELAXED_PVAL)"
call_relaxed_peaks "$IP_BAM" "$CTRL_BAM" "$ID" "$REL_OUT" "$PAIRED" "$FRAGLEN" bdg
n_relaxed=$(wc -l < "$REL_OUT/$ID.relaxed.narrowPeak")
# p < 0.01 always gives many peaks on real data; none means broken input
(( n_relaxed > 0 )) || die "MACS called no peaks even at p < $RELAXED_PVAL for $ID vs $CTRL_LABEL. \
Check that the IP and control are the right files (see the MACS read counts above)."

SIG_OUT=$BASE/signal
mkdir -p "$SIG_OUT"
treat_bdg=$REL_OUT/${ID}_treat_pileup.bdg
ctrl_bdg=$REL_OUT/${ID}_control_lambda.bdg
log "Fold-enrichment and p-value signal tracks"
"$MACS_CMD" bdgcmp -t "$treat_bdg" -c "$ctrl_bdg" -m FE --outdir "$STEP_TMP" --o-prefix "$ID"
bdg_to_bigwig "$STEP_TMP/${ID}_FE.bdg" "$SIG_OUT/$ID.fc_signal.bigWig"
rm -f "$STEP_TMP/${ID}_FE.bdg"
# ENCODE: scale the p-value track by the smaller library (in millions of reads)
# (counted on the reads MACS used, i.e. read 1 only when read types are mixed)
treat_total=$total
[[ $MACS_TREAT != "$IP_BAM" ]] && treat_total=$(samtools view -c -F 4 "$MACS_TREAT")
ctrl_total=$treat_total
[[ -n $MACS_CTRL ]] && ctrl_total=$(samtools view -c -F 4 "$MACS_CTRL")
sval=$(awk -v a="$treat_total" -v b="$ctrl_total" 'BEGIN{m=(a<b)?a:b; printf "%.6f", m/1000000}')
"$MACS_CMD" bdgcmp -t "$treat_bdg" -c "$ctrl_bdg" -m ppois -S "$sval" --outdir "$STEP_TMP" --o-prefix "$ID"
bdg_to_bigwig "$STEP_TMP/${ID}_ppois.bdg" "$SIG_OUT/$ID.pval_signal.bigWig"
rm -f "$STEP_TMP/${ID}_ppois.bdg" "$treat_bdg" "$ctrl_bdg"

# ---- fingerprint ----------------------------------------------------------------
if is_true "$RUN_FINGERPRINT" && [[ -n $CTRL_BAM ]]; then
    mkdir -p "$FP_OUT"
    # mates define the fragment for paired reads; single reads get this length
    extend=(--extendReads "${FRAGLEN:-$DEFAULT_FRAGLEN}")
    log "plotFingerprint"
    plotFingerprint -b "$IP_BAM" "$CTRL_BAM" --labels "$ID" "$CTRL_LABEL" \
        --plotFile "$FP_OUT/$ID.plotFingerprint.pdf" \
        --outRawCounts "$FP_OUT/$ID.plotFingerprint.raw.txt" \
        --outQualityMetrics "$FP_OUT/$ID.plotFingerprint.qcmetrics.txt" \
        --skipZeros --numberOfSamples "$FINGERPRINT_SAMPLES" "${extend[@]}" --numberOfProcessors "$CPUS"
fi

write_metrics "$METRICS_FILE" \
    control "$CTRL_LABEL" \
    peak_type "$PEAK_TYPE" \
    peaks "$n_peaks" \
    FRiP "$frip" \
    relaxed_peaks "$n_relaxed"

finish_step
