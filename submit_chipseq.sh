#!/usr/bin/env bash
# =============================================================================
# submit_chipseq.sh
#
# Submit the whole ChIP-seq pipeline for every sample in a samplesheet to
# SLURM as a graph of dependent jobs (sbatch --dependency=afterok).
#
#   per run      : trim (FastQC + Trim Galore) -> align (bowtie2)
#   per library  : filter (merge runs, ENCODE filtering, MarkDuplicates,
#                  library complexity) -> bamqc (SPP cross-correlation, bigWig)
#   per IP lib   : peaks (MACS2 main + relaxed calls, FRiP, signal tracks,
#                  fingerprint)                  [waits for its control too]
#   per sample   : repro (IDR / naive overlap across replicates, pseudoreps)
#   per antibody : consensus (merged peaks + featureCounts)
#   once         : multiqc (runs after everything, even if some jobs failed)
#
# Run with -h for usage. See README.md for details.
# =============================================================================
set -euo pipefail

PIPELINE_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

usage() {
    cat << EOF
Usage: $(basename "$0") -i samplesheet.csv -o outdir [options]

Required:
  -i, --input FILE     samplesheet (CSV, nf-core/chipseq format)
  -o, --outdir DIR     output directory

Options:
  -c, --config FILE    config file (default: $PIPELINE_DIR/conf/chipseq.config)
  -n, --dry-run        validate the samplesheet and print the jobs, submit nothing
  -r, --resume         skip steps that already finished in a previous run
                       (marked in <outdir>/pipeline_info/done/)
  -h, --help           show this help

Samplesheet columns (header required, order free):
  sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate[,peak_type]
  - one row per FASTQ (pair); rows with the same sample+replicate are merged
  - leave fastq_2 empty for single-end data
  - control rows (inputs) leave antibody, control, control_replicate empty
  - peak_type (optional) is narrow or broad; default from PEAK_TYPE_DEFAULT
EOF
}

die()  { echo "ERROR: $*" >&2; exit 1; }
warn() { echo "WARNING: $*" >&2; }
is_true() { [[ ${1,,} == "true" || ${1,,} == "yes" || $1 == "1" ]]; }

# ---- arguments ----------------------------------------------------------------
SAMPLESHEET="" OUTDIR="" CONFIG="$PIPELINE_DIR/conf/chipseq.config"
DRY_RUN=0 RESUME=0
while (( $# )); do
    case $1 in
        -i|--input)   SAMPLESHEET=${2:?}; shift 2 ;;
        -o|--outdir)  OUTDIR=${2:?}; shift 2 ;;
        -c|--config)  CONFIG=${2:?}; shift 2 ;;
        -n|--dry-run) DRY_RUN=1; shift ;;
        -r|--resume)  RESUME=1; shift ;;
        -h|--help)    usage; exit 0 ;;
        *) usage >&2; die "unknown option $1" ;;
    esac
done
[[ -n $SAMPLESHEET && -n $OUTDIR ]] || { usage >&2; die "--input and --outdir are required"; }
[[ -f $SAMPLESHEET ]] || die "samplesheet not found: $SAMPLESHEET"
[[ -f $CONFIG ]] || die "config not found: $CONFIG"
SAMPLESHEET=$(realpath "$SAMPLESHEET")
CONFIG=$(realpath "$CONFIG")
OUTDIR=$(realpath -m "$OUTDIR")

# shellcheck source=conf/chipseq.config
source "$CONFIG"
if (( ! DRY_RUN )); then
    command -v sbatch > /dev/null || die "sbatch not found (use --dry-run to test without SLURM)"
fi

# In a dry run, missing files are only reported as warnings
ERRORS=0
err() { echo "ERROR: $*" >&2; ERRORS=$(( ERRORS + 1 )); }
file_err() { if (( DRY_RUN )); then warn "$*"; else err "$*"; fi; }

# ---- references ---------------------------------------------------------------
if [[ ! -f $BOWTIE2_INDEX.1.bt2 && ! -f $BOWTIE2_INDEX.1.bt2l ]]; then
    file_err "bowtie2 index not found: $BOWTIE2_INDEX(.1.bt2|.1.bt2l)"
fi
[[ -f $CHROM_SIZES ]] || file_err "CHROM_SIZES not found: $CHROM_SIZES"
if [[ -n $BLACKLIST && ! -f $BLACKLIST ]]; then
    file_err "BLACKLIST not found: $BLACKLIST (set BLACKLIST=\"\" to skip blacklist filtering)"
fi

# ---- parse samplesheet ------------------------------------------------------------
# Normalise to one record per row, fields separated by \x1f (not whitespace, so
# empty fields survive `read`), in a fixed column order.
PARSED=$(mktemp)
trap 'rm -f "$PARSED"' EXIT
awk -F',' '
    { sub(/\r$/, "") }
    NR == 1 {
        for (i = 1; i <= NF; i++) { h = tolower($i); gsub(/^[ \t"]+|[ \t"]+$/, "", h); col[h] = i }
        n = split("sample fastq_1 replicate", req, " ")
        for (i = 1; i <= n; i++) if (!(req[i] in col)) { print "missing required column: " req[i] > "/dev/stderr"; bad = 1 }
        if (bad) exit 1
        nf = split("sample fastq_1 fastq_2 replicate antibody control control_replicate peak_type", fields, " ")
        next
    }
    /^[ \t,]*$/ || /^#/ { next }
    {
        for (i = 1; i <= NF; i++) gsub(/^[ \t"]+|[ \t"]+$/, "", $i)
        out = NR
        for (i = 1; i <= nf; i++) out = out "\x1f" ((fields[i] in col) ? $(col[fields[i]]) : "")
        print out
    }' "$SAMPLESHEET" > "$PARSED" || die "could not parse $SAMPLESHEET"

LIBS=()           # library ids in samplesheet order
RUNS=()           # run ids in samplesheet order
declare -A LIB_SAMPLE LIB_REP LIB_PAIRED LIB_AB LIB_CTRL_SAMPLE LIB_CTRL_REP LIB_PTYPE LIB_RUNS LIB_NRUNS
declare -A RUN_LIB RUN_FQ1 RUN_FQ2
declare -A SEEN_FASTQ

while IFS=$'\x1f' read -r line sample fq1 fq2 rep ab ctrl ctrl_rep ptype; do
    where="samplesheet line $line"
    if [[ ! $sample =~ ^[A-Za-z0-9._-]+$ ]]; then
        err "$where: sample '$sample' must be non-empty and use only letters, numbers, '.', '_' or '-'"
        continue
    fi
    [[ $rep =~ ^[1-9][0-9]*$ ]] || { err "$where: replicate '$rep' must be a positive integer"; continue; }
    [[ -n $fq1 ]] || { err "$where: fastq_1 is empty"; continue; }
    if [[ -z $ab && ( -n $ctrl || -n $ctrl_rep ) ]]; then
        err "$where: control/control_replicate given for '$sample' but antibody is empty"
    fi
    ptype=${ptype,,}
    if [[ -n $ab ]]; then
        ptype=${ptype:-$PEAK_TYPE_DEFAULT}
        [[ $ptype == narrow || $ptype == broad ]] || err "$where: peak_type '$ptype' must be narrow or broad"
    else
        ptype=""
    fi

    # FASTQ checks; relative paths are relative to the current directory
    paired=0
    fq1=$(realpath -m "$fq1")
    [[ -n $fq2 ]] && { fq2=$(realpath -m "$fq2"); paired=1; }
    for fq in "$fq1" ${fq2:+"$fq2"}; do
        [[ $fq =~ \.f(ast)?q\.gz$ ]] || err "$where: $fq must end in .fastq.gz or .fq.gz"
        [[ -f $fq ]] || file_err "$where: file not found $fq"
        [[ -n ${SEEN_FASTQ[$fq]:-} ]] && err "$where: $fq is also used on line ${SEEN_FASTQ[$fq]}"
        SEEN_FASTQ[$fq]=$line
    done

    lib="${sample}_REP${rep}"
    if [[ -z ${LIB_SAMPLE[$lib]:-} ]]; then
        LIBS+=("$lib")
        LIB_SAMPLE[$lib]=$sample; LIB_REP[$lib]=$rep; LIB_PAIRED[$lib]=$paired
        LIB_AB[$lib]=$ab; LIB_CTRL_SAMPLE[$lib]=$ctrl; LIB_CTRL_REP[$lib]=$ctrl_rep; LIB_PTYPE[$lib]=$ptype
        LIB_NRUNS[$lib]=0; LIB_RUNS[$lib]=""
    else
        [[ ${LIB_PAIRED[$lib]} == "$paired" ]] || err "$where: $lib mixes single-end and paired-end runs"
        [[ ${LIB_AB[$lib]} == "$ab" ]] || err "$where: $lib has different antibody values across rows"
        [[ ${LIB_CTRL_SAMPLE[$lib]} == "$ctrl" && ${LIB_CTRL_REP[$lib]} == "$ctrl_rep" ]] \
            || err "$where: $lib has different control values across rows"
        [[ ${LIB_PTYPE[$lib]} == "$ptype" ]] || err "$where: $lib has different peak_type values across rows"
    fi
    LIB_NRUNS[$lib]=$(( LIB_NRUNS[$lib] + 1 ))
    run="${lib}_T${LIB_NRUNS[$lib]}"
    RUNS+=("$run")
    RUN_LIB[$run]=$lib; RUN_FQ1[$run]=$fq1; RUN_FQ2[$run]=$fq2
    LIB_RUNS[$lib]+="${LIB_RUNS[$lib]:+ }$run"
done < "$PARSED"

(( ${#LIBS[@]} > 0 )) || die "no samples found in $SAMPLESHEET"

# ---- resolve design: controls, samples, antibodies ---------------------------------
declare -A LIB_CTRL SAMPLE_LIBS SAMPLE_PTYPE SAMPLE_AB AB_LIBS AB_PTYPE
SAMPLES=() ANTIBODIES=() IP_LIBS=()
for lib in "${LIBS[@]}"; do
    s=${LIB_SAMPLE[$lib]} ab=${LIB_AB[$lib]}
    # a sample is either a control or an IP for one antibody
    if [[ -n ${SAMPLE_AB[$s]+x} && ${SAMPLE_AB[$s]} != "$ab" ]]; then
        err "sample $s has inconsistent antibody values across replicates"
    fi
    SAMPLE_AB[$s]=$ab
    [[ -z $ab ]] && continue

    IP_LIBS+=("$lib")
    ctrl=${LIB_CTRL_SAMPLE[$lib]} crep=${LIB_CTRL_REP[$lib]}
    if [[ -z $ctrl ]]; then
        warn "$lib has no control; peaks will be called without an input"
        LIB_CTRL[$lib]="none"
    else
        if [[ -z $crep ]]; then
            # default: same replicate number if it exists, otherwise replicate 1
            crep=${LIB_REP[$lib]}
            [[ -n ${LIB_SAMPLE[${ctrl}_REP${crep}]:-} ]] || crep=1
        fi
        cl="${ctrl}_REP${crep}"
        if [[ -z ${LIB_SAMPLE[$cl]:-} ]]; then
            err "$lib: control $ctrl replicate $crep ($cl) is not in the samplesheet"
        elif [[ -n ${LIB_AB[$cl]} ]]; then
            err "$lib: control $cl is itself an IP (antibody '${LIB_AB[$cl]}')"
        fi
        LIB_CTRL[$lib]=$cl
    fi

    if [[ -z ${SAMPLE_LIBS[$s]:-} ]]; then
        SAMPLES+=("$s"); SAMPLE_PTYPE[$s]=${LIB_PTYPE[$lib]}
    elif [[ ${SAMPLE_PTYPE[$s]} != "${LIB_PTYPE[$lib]}" ]]; then
        err "sample $s has replicates with different peak types"
    fi
    SAMPLE_LIBS[$s]+="${SAMPLE_LIBS[$s]:+ }$lib"

    if [[ -z ${AB_LIBS[$ab]:-} ]]; then
        ANTIBODIES+=("$ab"); AB_PTYPE[$ab]=${LIB_PTYPE[$lib]}
    elif [[ ${AB_PTYPE[$ab]} != "${LIB_PTYPE[$lib]}" ]]; then
        err "antibody $ab is used with both narrow and broad peak types"
    fi
    AB_LIBS[$ab]+="${AB_LIBS[$ab]:+ }$lib"
done
(( ${#IP_LIBS[@]} > 0 )) || err "no IP libraries (rows with an antibody) in the samplesheet"

(( ERRORS == 0 )) || die "$ERRORS problem(s) found in the samplesheet / config, nothing submitted"

echo "Samplesheet OK: ${#RUNS[@]} runs, ${#LIBS[@]} libraries (${#IP_LIBS[@]} IP), ${#SAMPLES[@]} IP samples, ${#ANTIBODIES[@]} antibodies"

# ---- job submission ---------------------------------------------------------------
RUN_TAG=$(date '+%Y%m%d_%H%M%S')
DONE_DIR=$OUTDIR/pipeline_info/done
RUN_CONFIG=$CONFIG
MANIFEST=/dev/null
if (( ! DRY_RUN )); then
    mkdir -p "$OUTDIR/pipeline_info/design" "$DONE_DIR"
    # jobs read a frozen copy of the config and samplesheet
    RUN_CONFIG=$OUTDIR/pipeline_info/chipseq_${RUN_TAG}.config
    cp "$CONFIG" "$RUN_CONFIG"
    cp "$SAMPLESHEET" "$OUTDIR/pipeline_info/samplesheet_${RUN_TAG}.csv"
    MANIFEST=$OUTDIR/pipeline_info/jobs_${RUN_TAG}.tsv
    printf 'job_id\tstep\tid\tdependencies\tlog\n' > "$MANIFEST"
fi

N_SUBMITTED=0 N_SKIPPED=0 DRY_COUNTER=0
ALL_JOBS=()
JOBID=""

# submit <step> <id> <RESOURCE_PREFIX> <afterok|afterany> "<dep job ids>" <script> [args...]
# Sets JOBID ("" when the step is skipped by --resume).
submit() {
    local step=$1 id=$2 res=$3 dep_type=$4 deps=$5 script=$6
    shift 6
    local cpus_var="${res}_CPUS" mem_var="${res}_MEM" time_var="${res}_TIME"
    JOBID=""
    if (( RESUME )) && [[ $step != multiqc && -f $DONE_DIR/$step.$id.done ]]; then
        echo "  [skip] $step $id (already done)"
        N_SKIPPED=$(( N_SKIPPED + 1 ))
        return
    fi
    local logdir=$OUTDIR/logs/$step
    local opts=(
        --parsable
        --job-name="${JOB_PREFIX}.${step}.${id}"
        --partition="$SLURM_PARTITION"
        --nodes=1 --ntasks=1
        --cpus-per-task="${!cpus_var}"
        --mem="${!mem_var}"
        --time="${!time_var}"
        --output="$logdir/${step}.${id}.%j.log"
        --export="ALL,CHIPSEQ_PIPELINE_DIR=$PIPELINE_DIR,CHIPSEQ_CONFIG=$RUN_CONFIG,CHIPSEQ_OUTDIR=$OUTDIR"
    )
    [[ -n $SLURM_ACCOUNT ]] && opts+=(--account="$SLURM_ACCOUNT")
    [[ -n $SLURM_QOS ]] && opts+=(--qos="$SLURM_QOS")
    deps=$(echo "$deps" | xargs)
    if [[ -n $deps ]]; then
        opts+=(--dependency="${dep_type}:${deps// /:}")
        # cancel instead of pending forever when an upstream job fails
        [[ $dep_type == afterok ]] && opts+=(--kill-on-invalid-dep=yes)
    fi
    opts+=("${SLURM_EXTRA_ARGS[@]}")

    local cmd=(sbatch "${opts[@]}" "$PIPELINE_DIR/bin/$script" "$@")
    local ndeps shown
    read -ra ndeps <<< "$deps"
    shown=${deps:-none}
    (( ${#ndeps[@]} > 6 )) && shown="${#ndeps[@]} jobs"
    if (( DRY_RUN )); then
        DRY_COUNTER=$(( DRY_COUNTER + 1 ))
        JOBID="job${DRY_COUNTER}"
        printf '  [%s] %-9s %-32s after: %s\n' "$JOBID" "$step" "$id" "$shown"
        if [[ -n ${VERBOSE_DRY_RUN:-} ]]; then printf '        '; printf '%q ' "${cmd[@]}"; echo; fi
    else
        mkdir -p "$logdir"
        JOBID=$("${cmd[@]}") || die "sbatch failed for $step $id"
        JOBID=${JOBID%%;*}
        printf '  [%s] %-9s %-32s after: %s\n' "$JOBID" "$step" "$id" "$shown"
        printf '%s\t%s\t%s\t%s\t%s\n' "$JOBID" "$step" "$id" "${deps:-none}" \
            "$logdir/${step}.${id}.${JOBID}.log" >> "$MANIFEST"
    fi
    N_SUBMITTED=$(( N_SUBMITTED + 1 ))
    ALL_JOBS+=("$JOBID")
}

declare -A ALIGN_JOB FILTER_JOB BAMQC_JOB PEAKS_JOB

echo "Submitting runs (trim -> align)"
for run in "${RUNS[@]}"; do
    lib=${RUN_LIB[$run]}
    submit trim "$run" TRIM afterok "" 01_trim.sh "$run" "${RUN_FQ1[$run]}" ${RUN_FQ2[$run]:+"${RUN_FQ2[$run]}"}
    submit align "$run" ALIGN afterok "$JOBID" 02_align.sh "$run" "$lib" "${LIB_PAIRED[$lib]}"
    ALIGN_JOB[$run]=$JOBID
done

echo "Submitting libraries (filter -> bamqc)"
for lib in "${LIBS[@]}"; do
    deps=""
    # shellcheck disable=SC2086
    for run in ${LIB_RUNS[$lib]}; do deps+=" ${ALIGN_JOB[$run]}"; done
    # shellcheck disable=SC2086
    submit filter "$lib" FILTER afterok "$deps" 03_filter.sh "$lib" "${LIB_PAIRED[$lib]}" ${LIB_RUNS[$lib]}
    FILTER_JOB[$lib]=$JOBID
    submit bamqc "$lib" BAMQC afterok "$JOBID" 04_bam_qc.sh "$lib" "${LIB_PAIRED[$lib]}"
    BAMQC_JOB[$lib]=$JOBID
done

echo "Submitting peak calling"
for lib in "${IP_LIBS[@]}"; do
    ctrl=${LIB_CTRL[$lib]}
    deps="${BAMQC_JOB[$lib]}"
    [[ $ctrl != none ]] && deps+=" ${FILTER_JOB[$ctrl]}"
    submit peaks "$lib" PEAKS afterok "$deps" 05_call_peaks.sh "$lib" "$ctrl" "${LIB_PAIRED[$lib]}" "${LIB_PTYPE[$lib]}"
    PEAKS_JOB[$lib]=$JOBID
done

if is_true "$RUN_REPRODUCIBILITY"; then
    echo "Submitting replicate reproducibility"
    for s in "${SAMPLES[@]}"; do
        design=$OUTDIR/pipeline_info/design/$s.tsv
        deps=""
        if (( ! DRY_RUN )); then : > "$design"; fi
        # shellcheck disable=SC2086
        for lib in ${SAMPLE_LIBS[$s]}; do
            deps+=" ${PEAKS_JOB[$lib]}"
            (( DRY_RUN )) || printf '%s\t%s\t%s\n' "$lib" "${LIB_CTRL[$lib]}" "${LIB_PAIRED[$lib]}" >> "$design"
        done
        submit repro "$s" REPRO afterok "$deps" 06_reproducibility.sh "$s" "${SAMPLE_PTYPE[$s]}" "$design"
    done
fi

if is_true "$RUN_CONSENSUS"; then
    echo "Submitting consensus peaks"
    for ab in "${ANTIBODIES[@]}"; do
        deps="" all_pe=1
        # shellcheck disable=SC2086
        for lib in ${AB_LIBS[$ab]}; do
            deps+=" ${PEAKS_JOB[$lib]}"
            [[ ${LIB_PAIRED[$lib]} == 1 ]] || all_pe=0
        done
        # shellcheck disable=SC2086
        submit consensus "$ab" CONSENSUS afterok "$deps" 07_consensus.sh "$ab" "${AB_PTYPE[$ab]}" "$all_pe" ${AB_LIBS[$ab]}
    done
fi

echo "Submitting MultiQC"
submit multiqc all MULTIQC afterany "${ALL_JOBS[*]}" 08_multiqc.sh

echo
if (( DRY_RUN )); then
    echo "Dry run: $N_SUBMITTED jobs would be submitted ($N_SKIPPED skipped). Nothing was submitted."
    echo "Set VERBOSE_DRY_RUN=1 to print the full sbatch commands."
else
    echo "Submitted $N_SUBMITTED jobs ($N_SKIPPED skipped as already done)."
    echo "Job list : $MANIFEST"
    echo "Logs     : $OUTDIR/logs/<step>/"
    echo "Status   : $PIPELINE_DIR/pipeline_status.sh $OUTDIR"
fi
