#!/usr/bin/env bash
# Step 1 (per sequencing run): FastQC on raw reads, adapter/quality trimming
# with Trim Galore (which also runs FastQC on the trimmed reads).
#
# usage: 01_trim.sh <run_id> <fastq_1> [fastq_2]
source "${CHIPSEQ_PIPELINE_DIR}/bin/common.sh"

RUN_ID=$1
FQ1=$2
FQ2=${3:-}
init_step trim "$RUN_ID"

RAW_QC=$OUTDIR/fastqc/raw
TRIM_OUT=$OUTDIR/trimgalore
TRIM_WORK=$WORKDIR/trimgalore
mkdir -p "$RAW_QC" "$TRIM_OUT/fastqc" "$TRIM_WORK"

# Symlink inputs under the run id so every report is labelled consistently
if [[ -n $FQ2 ]]; then
    ln -sf "$FQ1" "$STEP_TMP/${RUN_ID}_1.fastq.gz"
    ln -sf "$FQ2" "$STEP_TMP/${RUN_ID}_2.fastq.gz"
    reads=("$STEP_TMP/${RUN_ID}_1.fastq.gz" "$STEP_TMP/${RUN_ID}_2.fastq.gz")
    paired_arg=(--paired)
else
    ln -sf "$FQ1" "$STEP_TMP/${RUN_ID}.fastq.gz"
    reads=("$STEP_TMP/${RUN_ID}.fastq.gz")
    paired_arg=()
fi

log "FastQC (raw)"
fastqc --quiet --threads "$CPUS" --outdir "$RAW_QC" "${reads[@]}"

# Trim Galore scales poorly beyond ~4 cores
cores=$(( CPUS > 4 ? 4 : CPUS ))
log "Trim Galore"
# shellcheck disable=SC2086
trim_galore --cores "$cores" --gzip "${paired_arg[@]}" \
    --fastqc --fastqc_args "--quiet --outdir $TRIM_OUT/fastqc" \
    --output_dir "$TRIM_WORK" $TRIM_GALORE_ARGS "${reads[@]}"
mv "$TRIM_WORK/${RUN_ID}"*_trimming_report.txt "$TRIM_OUT/"

finish_step
