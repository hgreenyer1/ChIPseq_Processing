#!/usr/bin/env bash
# Step 5 (per IP library):
#   - main MACS2 call (nf-core): narrowPeak (q < MACS_QVAL) or broadPeak
#     (--broad-cutoff MACS_BROAD_CUTOFF), blacklist filtered, FRiP
#   - relaxed MACS2 call (ENCODE, p < RELAXED_PVAL, top RELAXED_NPEAKS) used for
#     replicate reproducibility, plus fold-enrichment and -log10(p) bigWigs
#   - deepTools plotFingerprint of IP vs control
#
# usage: 05_call_peaks.sh <ip_library_id> <control_library_id|none> <paired 0|1> <narrow|broad>
source "${CHIPSEQ_PIPELINE_DIR}/bin/common.sh"

LIB_ID=$1
CTRL_ID=$2
PAIRED=$3
PEAK_TYPE=$4
init_step peaks "$LIB_ID"

IP_BAM=$(lib_bam "$LIB_ID")
[[ -s $IP_BAM ]] || die "missing IP BAM $IP_BAM"
CTRL_BAM=""
ctrl_args=()
if [[ $CTRL_ID != "none" ]]; then
    CTRL_BAM=$(lib_bam "$CTRL_ID")
    [[ -s $CTRL_BAM ]] || die "missing control BAM $CTRL_BAM"
    ctrl_args=(-c "$CTRL_BAM")
fi

FRAGLEN=$(lib_fraglen "$LIB_ID")
MACS_IN=$(macs_input_args "$PAIRED" "$FRAGLEN")
log "MACS2 input: $MACS_IN"

# ---- main peak call ------------------------------------------------------------
MAIN_OUT=$OUTDIR/macs2/$PEAK_TYPE
mkdir -p "$MAIN_OUT/qc"
if [[ $PEAK_TYPE == "broad" ]]; then
    type_args=(--broad --broad-cutoff "$MACS_BROAD_CUTOFF")
else
    type_args=(-q "$MACS_QVAL")
fi
log "MACS2 main ${PEAK_TYPE} peak call"
# shellcheck disable=SC2086
"$MACS_CMD" callpeak -t "$IP_BAM" "${ctrl_args[@]}" $MACS_IN -g "$MACS_GSIZE" \
    -n "$LIB_ID" --outdir "$MAIN_OUT" --tempdir "$STEP_TMP" --keep-dup all \
    "${type_args[@]}" $MACS_EXTRA_ARGS
raw_peaks=$MAIN_OUT/${LIB_ID}_peaks.${PEAK_TYPE}Peak
peaks=$(main_peaks "$LIB_ID" "$PEAK_TYPE")
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
printf 'library\treads\treads_in_peaks\tFRiP\n%s\t%s\t%s\t%s\n' "$LIB_ID" "$total" "$in_peaks" "$frip" \
    > "$MAIN_OUT/qc/$LIB_ID.FRiP.txt"

# ---- relaxed peaks + signal tracks (ENCODE) --------------------------------------
REL_OUT=$OUTDIR/macs2/relaxed/$LIB_ID
log "MACS2 relaxed call (p < $RELAXED_PVAL)"
call_relaxed_peaks "$IP_BAM" "$CTRL_BAM" "$LIB_ID" "$REL_OUT" "$PAIRED" "$FRAGLEN" bdg
n_relaxed=$(wc -l < "$(relaxed_peaks "$LIB_ID")")

SIG_OUT=$OUTDIR/macs2/signal
mkdir -p "$SIG_OUT"
treat_bdg=$REL_OUT/${LIB_ID}_treat_pileup.bdg
ctrl_bdg=$REL_OUT/${LIB_ID}_control_lambda.bdg
log "Fold-enrichment and p-value signal tracks"
"$MACS_CMD" bdgcmp -t "$treat_bdg" -c "$ctrl_bdg" -m FE --outdir "$STEP_TMP" --o-prefix "$LIB_ID"
bdg_to_bigwig "$STEP_TMP/${LIB_ID}_FE.bdg" "$SIG_OUT/$LIB_ID.fc_signal.bigWig"
rm -f "$STEP_TMP/${LIB_ID}_FE.bdg"
# ENCODE: scale the p-value track by the smaller library (in millions of reads)
ctrl_total=$total
[[ -n $CTRL_BAM ]] && ctrl_total=$(samtools view -c -F 4 "$CTRL_BAM")
sval=$(awk -v a="$total" -v b="$ctrl_total" 'BEGIN{m=(a<b)?a:b; printf "%.6f", m/1000000}')
"$MACS_CMD" bdgcmp -t "$treat_bdg" -c "$ctrl_bdg" -m ppois -S "$sval" --outdir "$STEP_TMP" --o-prefix "$LIB_ID"
bdg_to_bigwig "$STEP_TMP/${LIB_ID}_ppois.bdg" "$SIG_OUT/$LIB_ID.pval_signal.bigWig"
rm -f "$STEP_TMP/${LIB_ID}_ppois.bdg" "$treat_bdg" "$ctrl_bdg"

# ---- fingerprint ----------------------------------------------------------------
if is_true "$RUN_FINGERPRINT" && [[ -n $CTRL_BAM ]]; then
    FP_OUT=$OUTDIR/deeptools/plotFingerprint
    mkdir -p "$FP_OUT"
    if [[ $PAIRED == 1 ]]; then
        extend=(--extendReads)
    else
        extend=(--extendReads "${FRAGLEN:-$DEFAULT_FRAGLEN}")
    fi
    log "plotFingerprint"
    plotFingerprint -b "$IP_BAM" "$CTRL_BAM" --labels "$LIB_ID" "$CTRL_ID" \
        --plotFile "$FP_OUT/$LIB_ID.plotFingerprint.pdf" \
        --outRawCounts "$FP_OUT/$LIB_ID.plotFingerprint.raw.txt" \
        --outQualityMetrics "$FP_OUT/$LIB_ID.plotFingerprint.qcmetrics.txt" \
        --skipZeros --numberOfSamples "$FINGERPRINT_SAMPLES" "${extend[@]}" --numberOfProcessors "$CPUS"
fi

write_metrics "$METRICS_DIR/libraries/$LIB_ID/3_peaks.tsv" \
    control "$CTRL_ID" \
    peak_type "$PEAK_TYPE" \
    peaks "$n_peaks" \
    FRiP "$frip" \
    relaxed_peaks "$n_relaxed"

finish_step
