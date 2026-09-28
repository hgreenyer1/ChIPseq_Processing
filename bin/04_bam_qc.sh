#!/usr/bin/env bash
# Step 4 (per library, IP and control): strand cross-correlation with
# phantompeakqualtools (NSC, RSC, fragment length; ENCODE) and a normalised
# coverage bigWig with deepTools bamCoverage.
#
# usage: 04_bam_qc.sh <library_id> <paired 0|1>
source "${CHIPSEQ_PIPELINE_DIR}/bin/common.sh"

LIB_ID=$1
PAIRED=$2
init_step bamqc "$LIB_ID"

BAM=$(lib_bam "$LIB_ID")
[[ -s $BAM ]] || die "missing filtered BAM $BAM"
SPP_OUT=$OUTDIR/phantompeakqualtools
BW_OUT=$OUTDIR/bigwig
mkdir -p "$SPP_OUT" "$BW_OUT"

fraglen="" nsc="NA" rsc="NA"
rm -f "$SPP_OUT/$LIB_ID.fraglen.txt"

# ---- cross-correlation ------------------------------------------------------
if is_true "$RUN_SPP"; then
    log "phantompeakqualtools on up to $SPP_NREADS reads"
    ta=$STEP_TMP/$LIB_ID.subsample.tagAlign.gz
    # ENCODE: paired-end libraries use read 1 only
    if [[ $PAIRED == 1 ]]; then
        samtools view -b -f 64 "$BAM" | bedtools bamtobed -i stdin
    else
        bedtools bamtobed -i "$BAM"
    fi | awk 'BEGIN{OFS="\t"} $1!="chrM" {$4="N"; $5="1000"; print}' \
       | shuf_seeded "$PSEUDOREP_SEED" -n "$SPP_NREADS" | gzip -nc > "$ta"

    if run_spp.R -c="$ta" -p="$CPUS" -filtchr="chrM" -rf -tmpdir="$STEP_TMP" \
            -savp="$SPP_OUT/$LIB_ID.spp.pdf" -out="$SPP_OUT/$LIB_ID.spp.out" >&2; then
        # col3 = comma separated fragment length candidates (first is best)
        fraglen=$(cut -f3 "$SPP_OUT/$LIB_ID.spp.out" | cut -d, -f1)
        nsc=$(cut -f9 "$SPP_OUT/$LIB_ID.spp.out")
        rsc=$(cut -f10 "$SPP_OUT/$LIB_ID.spp.out")
        if [[ $fraglen =~ ^[0-9]+$ ]] && (( fraglen >= MIN_FRAGLEN )); then
            echo "$fraglen" > "$SPP_OUT/$LIB_ID.fraglen.txt"
        else
            warn "SPP fragment length '$fraglen' < MIN_FRAGLEN, it will not be used"
        fi
    else
        warn "run_spp.R failed for $LIB_ID (often low depth); continuing without NSC/RSC"
    fi
fi

# ---- coverage track ----------------------------------------------------------
if [[ $PAIRED == 1 ]]; then
    extend=(--extendReads)
else
    fl=$(lib_fraglen "$LIB_ID")
    extend=(--extendReads "${fl:-$DEFAULT_FRAGLEN}")
fi
log "bamCoverage ($BIGWIG_NORM)"
bamCoverage -b "$BAM" -o "$BW_OUT/$LIB_ID.$BIGWIG_NORM.bigWig" \
    --binSize "$BIGWIG_BINSIZE" --normalizeUsing "$BIGWIG_NORM" \
    --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" "${extend[@]}" \
    --numberOfProcessors "$CPUS"

write_metrics "$METRICS_DIR/libraries/$LIB_ID/2_bamqc.tsv" \
    fraglen "${fraglen:-NA}" NSC "$nsc" RSC "$rsc"

finish_step
