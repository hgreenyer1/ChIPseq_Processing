#!/usr/bin/env bash
# Step 8 (once, after everything else, even if some jobs failed): collect the
# per-library / per-sample metrics into summary tables and run MultiQC.
#
# usage: 08_multiqc.sh
source "${CHIPSEQ_PIPELINE_DIR}/bin/common.sh"
init_step multiqc all

QC=$OUTDIR/qc
MQC_CUSTOM=$STEP_TMP/custom_content
mkdir -p "$QC" "$MQC_CUSTOM"

# ---- library table: one row per library, one column per metric ---------------
lib_table=$QC/library_qc_summary.tsv
shopt -s nullglob
lib_dirs=("$METRICS_DIR"/libraries/*/)
if (( ${#lib_dirs[@]} > 0 )); then
    for d in "${lib_dirs[@]}"; do
        lib=$(basename "$d")
        for f in "$d"/*.tsv; do awk -v l="$lib" 'BEGIN{OFS="\t"} {print l,$1,$2}' "$f"; done
    done | awk -F'\t' '
        !($2 in seen_k) { seen_k[$2]=1; keys[++nk]=$2 }
        !($1 in seen_l) { seen_l[$1]=1; libs[++nl]=$1 }
        { v[$1 SUBSEP $2]=$3 }
        END {
            printf "library"; for (k=1;k<=nk;k++) printf "\t%s", keys[k]; print ""
            for (l=1;l<=nl;l++) {
                printf "%s", libs[l]
                for (k=1;k<=nk;k++) printf "\t%s", ((libs[l] SUBSEP keys[k]) in v) ? v[libs[l] SUBSEP keys[k]] : "NA"
                print ""
            }
        }' > "$lib_table"
    log "Wrote $lib_table"
    {
        echo "# id: 'chipseq_library_qc'"
        echo "# section_name: 'ChIP-seq library QC (ENCODE metrics)'"
        echo "# description: 'Read counts after ENCODE filtering, library complexity (NRF, PBC1, PBC2), strand cross-correlation (NSC, RSC) and FRiP. ENCODE targets: NRF > 0.9, PBC1 > 0.9, PBC2 > 10, NSC > 1.05, RSC > 0.8, FRiP > 0.01.'"
        echo "# plot_type: 'table'"
        cat "$lib_table"
    } > "$MQC_CUSTOM/chipseq_library_qc_mqc.tsv"
fi

# ---- reproducibility table -----------------------------------------------------
sample_files=("$METRICS_DIR"/samples/*.tsv)
if (( ${#sample_files[@]} > 0 )); then
    repro_table=$QC/reproducibility_summary.tsv
    awk 'FNR==1 && NR!=1 {next} {print}' "${sample_files[@]}" > "$repro_table"
    log "Wrote $repro_table"
    {
        echo "# id: 'chipseq_reproducibility'"
        echo "# section_name: 'Replicate reproducibility (ENCODE)'"
        echo "# description: 'IDR (narrow) or naive overlap (broad) between true replicates (Nt) and pooled pseudoreplicates (Np). Rescue and self-consistency ratios <= 2 pass.'"
        echo "# plot_type: 'table'"
        cat "$repro_table"
    } > "$MQC_CUSTOM/chipseq_reproducibility_mqc.tsv"
fi
shopt -u nullglob

# ---- MultiQC ------------------------------------------------------------------------
log "MultiQC"
multiqc --force --outdir "$OUTDIR/multiqc" --filename multiqc_report.html \
    --config "$CHIPSEQ_PIPELINE_DIR/assets/multiqc_config.yaml" \
    --ignore "$WORKDIR" --ignore "$OUTDIR/logs" --ignore "$OUTDIR/pipeline_info" \
    "$OUTDIR" "$MQC_CUSTOM"

finish_step
