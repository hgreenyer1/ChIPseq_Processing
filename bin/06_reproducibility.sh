#!/usr/bin/env bash
# Step 6 (per IP sample, i.e. all replicates of one sample): ENCODE replicate
# reproducibility on the relaxed peak sets.
#
#   true replicates   : every pair of replicates           -> Nt (best pair)
#   pooled pseudoreps : pooled reads split in two halves   -> Np
#   self pseudoreps   : each replicate split in two halves -> N1, N2, ...
#
# Comparisons use IDR (narrow samples) or naive overlap (broad samples; a
# pooled peak must overlap a peak in both sets by >= 50% of either peak).
#   conservative set = best true-replicate set
#   optimal set      = the larger of the true-replicate and pooled-pseudorep sets
#   rescue ratio     = max(Np,Nt) / min(Np,Nt)
#   self-consistency = max(Ni) / min(Ni)
# Reproducibility passes if both ratios are <= 2, is borderline if one is, and
# fails otherwise. Unreplicated samples fall back to the self-pseudorep set.
#
# usage: 06_reproducibility.sh <sample> <narrow|broad> <design.tsv>
#   design.tsv columns: ip_library  control_library|none  paired(0|1)
source "${CHIPSEQ_PIPELINE_DIR}/bin/common.sh"

SAMPLE=$1
PEAK_TYPE=$2
DESIGN=$3
init_step repro "$SAMPLE"

METHOD=$REPRO_METHOD
if [[ $METHOD == "auto" ]]; then
    [[ $PEAK_TYPE == "broad" ]] && METHOD="overlap" || METHOD="idr"
fi
OUT=$OUTDIR/reproducibility/$SAMPLE
CMP=$OUT/comparisons
PEAKS_WORK=$STEP_TMP/peaks
mkdir -p "$OUT" "$CMP" "$PEAKS_WORK"
log "Reproducibility for $SAMPLE using $METHOD"

reps=() ctrls=() paired=()
while IFS=$'\t' read -r lib ctrl pe; do
    [[ -z $lib ]] && continue
    reps+=("$lib"); ctrls+=("$ctrl"); paired+=("$pe")
done < "$DESIGN"
n=${#reps[@]}
(( n > 0 )) || die "no replicates in $DESIGN"

ctrl_bam_of() { [[ $1 == "none" ]] && echo "" || lib_bam "$1"; }

# split_bam <in.bam> <out1.bam> <out2.bam>: random half of the read names (mates stay together)
split_bam() {
    local in=$1 o1=$2 o2=$3 names=$STEP_TMP/names.txt nhalf
    samtools view "$in" | cut -f1 | sort -u -S "$(job_mem_75)" -T "$STEP_TMP" \
        | shuf_seeded "$PSEUDOREP_SEED" > "$names"
    nhalf=$(( $(wc -l < "$names") / 2 ))
    head -n "$nhalf" "$names" > "$names.1"
    samtools view -@ "$CPUS" -b -N "$names.1" -o "$o1" -U "$o2" "$in"
    samtools index "$o1"; samtools index "$o2"
    rm -f "$names" "$names.1"
}

# compare <peaks_a> <peaks_b> <pooled_peaks> <out_prefix>
# Writes <out_prefix>.<method>.narrowPeak (blacklist filtered), echoes its count.
compare() {
    local a=$1 b=$2 pooled=$3 prefix=$4 res="$4.$METHOD.narrowPeak"
    if [[ $METHOD == "idr" ]]; then
        local thr
        thr=$(awk -v p="$IDR_THRESHOLD" 'BEGIN{print -log(p)/log(10)}')
        if run_idr --samples "$a" "$b" --peak-list "$pooled" --input-file-type narrowPeak \
               --output-file "$prefix.idr.txt" --rank "$IDR_RANK" \
               --soft-idr-threshold "$IDR_THRESHOLD" --plot --use-best-multisummit-IDR \
               --log-output-file "$prefix.idr.log" >&2; then
            awk -v t="$thr" 'BEGIN{OFS="\t"} $12>=t {if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' \
                "$prefix.idr.txt" | sort -u | sort -k7,7gr > "$prefix.unfilt"
            gzip -nf "$prefix.idr.txt"
        else
            # IDR cannot fit its model on very small peak sets; anything else is a real error
            local na nb
            na=$(wc -l < "$a"); nb=$(wc -l < "$b")
            (( na < 100 || nb < 100 )) \
                || die "IDR failed for $(basename "$prefix") ($na and $nb input peaks), see log above"
            warn "IDR failed for $(basename "$prefix") with only $na / $nb input peaks; reporting 0 peaks"
            : > "$prefix.unfilt"
        fi
    else
        # narrowPeak (10 cols) x narrowPeak: b coords in $12,$13, overlap length in $21
        local keep='BEGIN{FS=OFS="\t"} {s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) print}'
        bedtools intersect -wo -a "$pooled" -b "$a" | awk "$keep" | cut -f1-10 | sort -u \
            | bedtools intersect -wo -a stdin -b "$b" | awk "$keep" | cut -f1-10 | sort -u \
            | sort -k7,7gr > "$prefix.unfilt"
    fi
    blacklist_filter "$prefix.unfilt" "$res"
    rm -f "$prefix.unfilt"
    wc -l < "$res"
}

# ---- self pseudoreplicates -----------------------------------------------------
declare -a self_n
for i in "${!reps[@]}"; do
    lib=${reps[$i]}
    rel=$(relaxed_peaks "$lib")
    [[ -s $rel ]] || die "missing relaxed peaks $rel"
    log "Self pseudoreplicates of $lib"
    split_bam "$(lib_bam "$lib")" "$STEP_TMP/$lib.pr1.bam" "$STEP_TMP/$lib.pr2.bam"
    fl=$(lib_fraglen "$lib")
    for pr in pr1 pr2; do
        call_relaxed_peaks "$STEP_TMP/$lib.$pr.bam" "$(ctrl_bam_of "${ctrls[$i]}")" "$lib.$pr" \
            "$PEAKS_WORK" "${paired[$i]}" "$fl"
        rm -f "$STEP_TMP/$lib.$pr.bam" "$STEP_TMP/$lib.$pr.bam.bai"
    done
    self_n[i]=$(compare "$PEAKS_WORK/$lib.pr1.relaxed.narrowPeak" "$PEAKS_WORK/$lib.pr2.relaxed.narrowPeak" \
        "$rel" "$CMP/$lib.pr1_vs_pr2")
    log "  $lib self-pseudorep peaks: ${self_n[i]}"
done

rescue="NA" self_cons="NA" flag="NA" nt="NA" np="NA" best_pair="NA"
if (( n >= 2 )); then
    # ---- pooled replicates and controls ------------------------------------------
    log "Pooling replicates"
    rep_bams=() all_pe=1 fl_sum=0 fl_n=0
    for i in "${!reps[@]}"; do
        rep_bams+=("$(lib_bam "${reps[$i]}")")
        [[ ${paired[$i]} == 1 ]] || all_pe=0
        fl=$(lib_fraglen "${reps[$i]}")
        [[ -n $fl ]] && { fl_sum=$(( fl_sum + fl )); fl_n=$(( fl_n + 1 )); }
    done
    pooled_fl=""
    (( fl_n > 0 )) && pooled_fl=$(( fl_sum / fl_n ))
    samtools merge -@ "$CPUS" -f "$STEP_TMP/pooled.bam" "${rep_bams[@]}"
    samtools index "$STEP_TMP/pooled.bam"

    mapfile -t uniq_ctrls < <(printf '%s\n' "${ctrls[@]}" | grep -vx none | sort -u)
    pooled_ctrl=""
    if (( ${#uniq_ctrls[@]} == 1 )); then
        pooled_ctrl=$(lib_bam "${uniq_ctrls[0]}")
    elif (( ${#uniq_ctrls[@]} > 1 )); then
        ctrl_bams=()
        for c in "${uniq_ctrls[@]}"; do ctrl_bams+=("$(lib_bam "$c")"); done
        samtools merge -@ "$CPUS" -f "$STEP_TMP/pooled_ctrl.bam" "${ctrl_bams[@]}"
        samtools index "$STEP_TMP/pooled_ctrl.bam"
        pooled_ctrl=$STEP_TMP/pooled_ctrl.bam
    fi

    log "Relaxed peaks on pooled replicates"
    call_relaxed_peaks "$STEP_TMP/pooled.bam" "$pooled_ctrl" "${SAMPLE}.pooled" "$PEAKS_WORK" "$all_pe" "$pooled_fl"
    pooled_peaks=$PEAKS_WORK/${SAMPLE}.pooled.relaxed.narrowPeak
    cp "$pooled_peaks" "$OUT/"

    log "Pooled pseudoreplicates"
    split_bam "$STEP_TMP/pooled.bam" "$STEP_TMP/pooled.pr1.bam" "$STEP_TMP/pooled.pr2.bam"
    rm -f "$STEP_TMP/pooled.bam" "$STEP_TMP/pooled.bam.bai"
    for pr in pr1 pr2; do
        call_relaxed_peaks "$STEP_TMP/pooled.$pr.bam" "$pooled_ctrl" "${SAMPLE}.pooled.$pr" \
            "$PEAKS_WORK" "$all_pe" "$pooled_fl"
        rm -f "$STEP_TMP/pooled.$pr.bam" "$STEP_TMP/pooled.$pr.bam.bai"
    done
    np=$(compare "$PEAKS_WORK/${SAMPLE}.pooled.pr1.relaxed.narrowPeak" \
        "$PEAKS_WORK/${SAMPLE}.pooled.pr2.relaxed.narrowPeak" "$pooled_peaks" "$CMP/${SAMPLE}.pooled_pr1_vs_pr2")
    log "  pooled pseudorep peaks (Np): $np"

    # ---- true replicates -----------------------------------------------------
    nt=-1
    for (( i = 0; i < n; i++ )); do
        for (( j = i + 1; j < n; j++ )); do
            pair="${reps[$i]}_vs_${reps[$j]}"
            c=$(compare "$(relaxed_peaks "${reps[$i]}")" "$(relaxed_peaks "${reps[$j]}")" \
                "$pooled_peaks" "$CMP/$pair")
            log "  $pair peaks: $c"
            if (( c > nt )); then nt=$c; best_pair=$pair; fi
        done
    done
    log "  best true-replicate pair (Nt): $best_pair with $nt peaks"

    cp "$CMP/$best_pair.$METHOD.narrowPeak" "$OUT/$SAMPLE.conservative.$METHOD.narrowPeak"
    if (( nt >= np )); then
        cp "$CMP/$best_pair.$METHOD.narrowPeak" "$OUT/$SAMPLE.optimal.$METHOD.narrowPeak"
        optimal_src=$best_pair
    else
        cp "$CMP/${SAMPLE}.pooled_pr1_vs_pr2.$METHOD.narrowPeak" "$OUT/$SAMPLE.optimal.$METHOD.narrowPeak"
        optimal_src="pooled_pseudoreplicates"
    fi

    ratio() { awk -v a="$1" -v b="$2" 'BEGIN{mx=(a>b)?a:b; mn=(a<b)?a:b; if (mn>0) printf "%.3f", mx/mn; else print "inf"}'; }
    rescue=$(ratio "$np" "$nt")
    smin=$(printf '%s\n' "${self_n[@]}" | awk 'NR==1 || $1<m {m=$1} END{print m}')
    smax=$(printf '%s\n' "${self_n[@]}" | awk 'NR==1 || $1>m {m=$1} END{print m}')
    self_cons=$(ratio "$smax" "$smin")
    flag=$(awk -v r="$rescue" -v s="$self_cons" 'BEGIN{
        bad = (r=="inf" || r+0>2) + (s=="inf" || s+0>2);
        print (bad==0) ? "pass" : (bad==1) ? "borderline" : "fail" }')
else
    warn "$SAMPLE has a single replicate; using its self-pseudoreplicate peaks"
    lib=${reps[0]}
    cp "$CMP/$lib.pr1_vs_pr2.$METHOD.narrowPeak" "$OUT/$SAMPLE.conservative.$METHOD.narrowPeak"
    cp "$CMP/$lib.pr1_vs_pr2.$METHOD.narrowPeak" "$OUT/$SAMPLE.optimal.$METHOD.narrowPeak"
    optimal_src="self_pseudoreplicates"
fi

n_cons=$(wc -l < "$OUT/$SAMPLE.conservative.$METHOD.narrowPeak")
n_opt=$(wc -l < "$OUT/$SAMPLE.optimal.$METHOD.narrowPeak")
self_str=$(IFS=,; echo "${self_n[*]}")
{
    printf 'sample\tmethod\treplicates\tNt\tNp\tself_pseudorep_peaks\tbest_pair\tconservative_peaks\toptimal_peaks\toptimal_source\trescue_ratio\tself_consistency_ratio\treproducibility\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SAMPLE" "$METHOD" "$n" "$nt" "$np" "$self_str" \
        "$best_pair" "$n_cons" "$n_opt" "$optimal_src" "$rescue" "$self_cons" "$flag"
} > "$OUT/$SAMPLE.reproducibility.qc.tsv"
mkdir -p "$METRICS_DIR/samples"
cp "$OUT/$SAMPLE.reproducibility.qc.tsv" "$METRICS_DIR/samples/$SAMPLE.tsv"
log "conservative: $n_cons, optimal: $n_opt, reproducibility: $flag"

finish_step
