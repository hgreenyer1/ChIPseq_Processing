#!/usr/bin/env bash
# Step 2 (per sequencing run): align trimmed reads with bowtie2 and
# coordinate-sort.
#
# usage: 02_align.sh <run_id> <library_id> <paired 0|1>
source "${CHIPSEQ_PIPELINE_DIR}/bin/common.sh"

RUN_ID=$1
LIB_ID=$2
PAIRED=$3
init_step align "$RUN_ID"

TRIM_WORK=$WORKDIR/trimgalore
BAM_WORK=$WORKDIR/bowtie2
LOG_DIR=$OUTDIR/bowtie2/logs
mkdir -p "$BAM_WORK" "$LOG_DIR"

if [[ $PAIRED == 1 ]]; then
    trimmed=("$TRIM_WORK/${RUN_ID}_1_val_1.fq.gz" "$TRIM_WORK/${RUN_ID}_2_val_2.fq.gz")
    reads_args=(-1 "${trimmed[0]}" -2 "${trimmed[1]}" -X "$MAX_FRAGMENT_LENGTH")
else
    trimmed=("$TRIM_WORK/${RUN_ID}_trimmed.fq.gz")
    reads_args=(-U "${trimmed[0]}")
fi
for f in "${trimmed[@]}"; do [[ -s $f ]] || die "missing trimmed reads $f"; done

# Leave a couple of cores to samtools sort
bt2_threads=$(( CPUS > 2 ? CPUS - 2 : 1 ))
log "bowtie2 with ${bt2_threads} threads"
# shellcheck disable=SC2086
bowtie2 -p "$bt2_threads" --mm -x "$BOWTIE2_INDEX" "${reads_args[@]}" \
        --rg-id "$RUN_ID" --rg "SM:$LIB_ID" --rg "PL:ILLUMINA" --rg "LB:$LIB_ID" \
        $BOWTIE2_ARGS 2> "$LOG_DIR/${RUN_ID}.bowtie2.log" \
    | samtools sort -@ 2 -m "$(sort_mem_per_thread 2)" -T "$STEP_TMP/${RUN_ID}" \
        -o "$BAM_WORK/${RUN_ID}.sorted.bam" -
samtools index "$BAM_WORK/${RUN_ID}.sorted.bam"
sed 's/^/    /' "$LOG_DIR/${RUN_ID}.bowtie2.log" >&2

if is_true "$CLEAN_INTERMEDIATES"; then
    rm -f "${trimmed[@]}"
fi

finish_step
