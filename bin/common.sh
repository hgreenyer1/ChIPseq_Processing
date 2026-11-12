#!/usr/bin/env bash
# =============================================================================
# Helpers shared by every step script in bin/.
#
# Step scripts are run by SLURM from a spool copy, so they locate this file via
# CHIPSEQ_PIPELINE_DIR, which submit_chipseq.sh exports together with
# CHIPSEQ_CONFIG and CHIPSEQ_OUTDIR.
# =============================================================================

set -euo pipefail

log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2; }
warn() { log "WARNING: $*"; }
die()  { log "ERROR: $*"; exit 1; }

# init_step <step> <id>
# Loads the config, the software environment and sets the common paths.
init_step() {
    STEP_NAME=$1
    STEP_ID=$2
    : "${CHIPSEQ_PIPELINE_DIR:?not set - submit through submit_chipseq.sh}"
    : "${CHIPSEQ_CONFIG:?not set - submit through submit_chipseq.sh}"
    : "${CHIPSEQ_OUTDIR:?not set - submit through submit_chipseq.sh}"

    # shellcheck source=/dev/null
    source "$CHIPSEQ_CONFIG"
    # conda/module activation scripts often reference unset variables
    if declare -F setup_env > /dev/null; then
        set +eu
        setup_env
        set -eu
    fi

    OUTDIR=$CHIPSEQ_OUTDIR
    WORKDIR=$OUTDIR/work
    DONE_DIR=$OUTDIR/pipeline_info/done
    METRICS_DIR=$OUTDIR/qc/metrics
    CPUS=${SLURM_CPUS_PER_TASK:-1}
    MEM_MB=${SLURM_MEM_PER_NODE:-4096}
    STEP_TMP=$WORKDIR/tmp/${STEP_NAME}.${STEP_ID}.${SLURM_JOB_ID:-$$}
    mkdir -p "$WORKDIR" "$DONE_DIR" "$METRICS_DIR" "$STEP_TMP"
    export TMPDIR=$STEP_TMP
    export LC_ALL=C

    rm -f "$DONE_DIR/${STEP_NAME}.${STEP_ID}.done"
    trap 'log "FAILED: ${STEP_NAME} ${STEP_ID} (line ${LINENO}, exit $?)"' ERR
    trap 'rm -rf "$STEP_TMP"' EXIT
    log "Starting ${STEP_NAME} for ${STEP_ID} on $(hostname) (job ${SLURM_JOB_ID:-local}, ${CPUS} cpus, ${MEM_MB} MB)"
}

# finish_step: mark the step as complete so --resume can skip it
finish_step() {
    touch "$DONE_DIR/${STEP_NAME}.${STEP_ID}.done"
    log "Finished ${STEP_NAME} for ${STEP_ID}"
}

is_true() { [[ ${1,,} == "true" || ${1,,} == "yes" || $1 == "1" ]]; }

# Memory per thread for `samtools sort -m` (~60% of the job split over threads)
sort_mem_per_thread() {
    local threads=${1:-$CPUS} mb
    mb=$(( MEM_MB * 6 / 10 / threads ))
    (( mb < 256 )) && mb=256
    echo "${mb}M"
}

# Memory for GNU sort -S / java -Xmx (~75% of the job)
job_mem_75() { echo "$(( MEM_MB * 3 / 4 ))M"; }

# write_metrics <file> <key> <value> [<key> <value> ...]
write_metrics() {
    local file=$1; shift
    mkdir -p "$(dirname "$file")"
    : > "$file"
    while (( $# >= 2 )); do
        printf '%s\t%s\n' "$1" "$2" >> "$file"
        shift 2
    done
}

# Reproducible randomness for shuf: shuf_seeded <seed> [shuf args...]
shuf_seeded() {
    local seed=$1; shift
    if command -v openssl > /dev/null; then
        shuf --random-source=<(openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt < /dev/zero 2> /dev/null) "$@"
    else
        shuf "$@"
    fi
}

# ---------------------------------------------------------------------------
# Standard output locations
# ---------------------------------------------------------------------------
lib_bam()      { echo "$OUTDIR/bowtie2/merged_library/$1.filt.dedup.bam"; }
lib_fraglen()  { local f="$OUTDIR/phantompeakqualtools/$1.fraglen.txt"; [[ -s $f ]] && cat "$f" || echo ""; }
relaxed_peaks() { echo "$OUTDIR/macs2/relaxed/$1/$1.relaxed.narrowPeak"; }
main_peaks()   { echo "$OUTDIR/macs2/$2/$1_peaks.bfilt.$2Peak"; }   # <lib> <narrow|broad>
# merged replicates (per sample)
merged_bam()     { echo "$OUTDIR/bowtie2/merged_replicate/$1.mRp.bam"; }
merged_fraglen() { local f="$OUTDIR/bowtie2/merged_replicate/qc/$1.fraglen.txt"; [[ -s $f ]] && cat "$f" || echo ""; }

# ---------------------------------------------------------------------------
# MACS2
# ---------------------------------------------------------------------------

# macs_input_args <paired 0|1> <fraglen>
# Paired-end: fragments come from the mates (BAMPE). Single-end: reads are
# shifted by the SPP fragment length estimate (ENCODE), or DEFAULT_FRAGLEN.
macs_input_args() {
    local paired=$1 fraglen=$2
    if [[ $paired == 1 ]]; then
        echo "-f BAMPE"
    elif [[ $fraglen =~ ^[0-9]+$ ]] && (( fraglen >= MIN_FRAGLEN && fraglen <= ${MAX_FRAGLEN:-500} )); then
        echo "-f BAM --nomodel --shift 0 --extsize $fraglen"
    else
        echo "-f BAM --nomodel --shift 0 --extsize $DEFAULT_FRAGLEN"
    fi
}

# bam_is_paired <bam>: true if the BAM contains paired-end reads (stops at the first one)
bam_is_paired() {
    local n
    n=$({ samtools view -f 1 "$1" 2> /dev/null || true; } | head -n 1 | wc -l)
    (( n > 0 ))
}

# read1_bam <bam>: read-1-only copy of a paired-end BAM in STEP_TMP (cached), echoes its path
read1_bam() {
    local out
    out=$STEP_TMP/$(basename "$1" .bam).read1.bam
    if [[ ! -s $out ]]; then
        samtools view -@ "$CPUS" -b -f 64 -o "$out" "$1"
        samtools index "$out"
    fi
    echo "$out"
}

# mean_insert_size <paired bam>: mean fragment length of the first 200k proper pairs
mean_insert_size() {
    { samtools view -f 66 "$1" 2> /dev/null || true; } | head -n 200000 \
        | awk '{t = ($9 < 0) ? -$9 : $9; if (t > 0) {s += t; n++}} END {if (n) printf "%d", s / n}'
}

# macs_inputs <treat_bam> <ctrl_bam|""> <paired 0|1> <fraglen>
# Sets MACS_TREAT, MACS_CTRL and MACS_FMT (array of format/model args).
# MACS reads treatment and control in a single format, and reading a
# single-end BAM as BAMPE turns every read into a 0 bp fragment (no peaks at
# all). So when exactly one side is paired-end, both are used as single reads:
# the paired-end BAM is reduced to read 1 (one read per fragment) and reads are
# extended by the fragment length.
macs_inputs() {
    local treat=$1 ctrl=$2 paired=$3 fraglen=$4 ctrl_paired=$3
    MACS_TREAT=$treat
    MACS_CTRL=$ctrl
    if [[ -n $ctrl ]]; then
        if bam_is_paired "$ctrl"; then ctrl_paired=1; else ctrl_paired=0; fi
    fi
    if [[ $paired != "$ctrl_paired" ]]; then
        if [[ $paired == 1 ]]; then
            # the real fragment length is known from the pairs (more reliable
            # than the cross-correlation estimate, especially for broad marks)
            fraglen=$(mean_insert_size "$treat")
            log "Treatment is paired-end but control is single-end: using read 1 of the treatment, extended to the mean insert size (${fraglen:-unknown} bp)"
            MACS_TREAT=$(read1_bam "$treat")
        else
            log "Treatment is single-end but control is paired-end: using read 1 of the control"
            MACS_CTRL=$(read1_bam "$ctrl")
        fi
        paired=0
    fi
    read -ra MACS_FMT <<< "$(macs_input_args "$paired" "$fraglen")"
}

# call_relaxed_peaks <treat_bam> <ctrl_bam|""> <name> <outdir> <paired> <fraglen> [bdg]
# ENCODE-style relaxed call (p < RELAXED_PVAL, top RELAXED_NPEAKS by p-value)
# used for IDR / overlap. Writes <outdir>/<name>.relaxed.narrowPeak.
# Pass "bdg" as the 7th argument to also keep MACS2 bedGraphs (-B --SPMR).
call_relaxed_peaks() {
    local treat=$1 ctrl=$2 name=$3 outdir=$4 paired=$5 fraglen=$6 bdg=${7:-}
    local ctrl_args=() bdg_args=()
    macs_inputs "$treat" "$ctrl" "$paired" "$fraglen"
    [[ -n $MACS_CTRL ]] && ctrl_args=(-c "$MACS_CTRL")
    [[ $bdg == "bdg" ]] && bdg_args=(-B --SPMR)
    mkdir -p "$outdir"
    # shellcheck disable=SC2086
    "$MACS_CMD" callpeak -t "$MACS_TREAT" "${ctrl_args[@]}" "${MACS_FMT[@]}" \
        -g "$MACS_GSIZE" -n "$name" --outdir "$outdir" --tempdir "$STEP_TMP" \
        -p "$RELAXED_PVAL" --keep-dup all --call-summits "${bdg_args[@]}" $MACS_EXTRA_ARGS
    sort -k8,8gr -S "$(job_mem_75)" -T "$STEP_TMP" "$outdir/${name}_peaks.narrowPeak" \
        | awk -v n="$RELAXED_NPEAKS" 'BEGIN{OFS="\t"} NR<=n { $4="Peak_"NR; if ($2<0) $2=0; print }' \
        > "$outdir/$name.relaxed.narrowPeak"
    rm -f "$outdir/${name}_peaks.narrowPeak" "$outdir/${name}_summits.bed" "$outdir/${name}_peaks.xls"
}

# run_idr [idr args...]
# idr 2.0.4.2 uses numpy aliases (numpy.int, ...) removed in numpy >= 1.24,
# which macs3/deepTools require. Unless IDR_CMD is set, run idr's main() with
# the aliases restored, using the interpreter from idr's own shebang.
run_idr() {
    if [[ -n ${IDR_CMD:-} ]]; then
        # shellcheck disable=SC2086
        $IDR_CMD "$@"
        return
    fi
    local idr_bin interp
    idr_bin=$(command -v idr) || die "idr not found"
    read -r -a interp < <(head -n 1 "$idr_bin" | sed 's/^#!//')
    [[ ${#interp[@]} -gt 0 ]] || interp=(python3)
    "${interp[@]}" -c '
import sys, numpy
for name, alias in (("int", int), ("float", float), ("bool", bool), ("object", object)):
    if not hasattr(numpy, name):
        setattr(numpy, name, alias)
from idr.idr import main
sys.argv[0] = "idr"
main()
' "$@"
}

# blacklist_filter <in.bed> <out.bed>
blacklist_filter() {
    if [[ -n ${BLACKLIST:-} && -s ${BLACKLIST:-} ]]; then
        bedtools intersect -v -a "$1" -b "$BLACKLIST" > "$2"
    else
        cp "$1" "$2"
    fi
}

# bdg_to_bigwig <in.bdg> <out.bigWig>: clip to chromosome ends, sort, convert
bdg_to_bigwig() {
    local bdg=$1 bw=$2
    bedtools slop -i "$bdg" -g "$CHROM_SIZES" -b 0 | bedClip stdin "$CHROM_SIZES" "$bdg.clip"
    sort -k1,1 -k2,2n -S "$(job_mem_75)" -T "$STEP_TMP" "$bdg.clip" > "$bdg.sorted"
    bedGraphToBigWig "$bdg.sorted" "$CHROM_SIZES" "$bw"
    rm -f "$bdg.clip" "$bdg.sorted"
}
