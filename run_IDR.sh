#!/bin/bash 

#SBATCH --partition=bluemoon
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=21G
#SBATCH --time=30:00:00

#=============
# Simple script to run IDR for all pairs of repliaces in directory 
# Based on ENCODE
#==================

PREFIX=$1
REP1_PEAK_FILE=$2
REP2_PEAK_FILE=$3
POOLED_PEAK_FILE=$4

REP1_VS_REP2="${PREFIX}_R1vR2"

IDR_THRESH=0.05
IDR_OUTPUT="${REP1_VS_REP2}.IDR0.05.value.txt"

# =============================
# Perform IDR analysis.
# Generate a plot and IDR output with additional columns including IDR scores.
# =============================
idr --samples ${REP1_PEAK_FILE} ${REP2_PEAK_FILE} --peak-list ${POOLED_PEAK_FILE} --input-file-type narrowPeak --output-file ${IDR_OUTPUT} --rank signal.value --soft-idr-threshold ${IDR_THRESH} --plot --use-best-multisummit-IDR

# =============================
# Get peaks passing IDR threshold of 5%
# =============================
IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}')

awk 'BEGIN{OFS="\t"} $12>='"${IDR_THRESH_TRANSFORMED}"' {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' ${IDR_OUTPUT} | sort | uniq | sort -k7n,7n | gzip -nc > ${REP1_VS_REP2}.IDR0.05.narrowPeak.gz

NPEAKS_IDR=$(zcat ${REP1_VS_REP2}.IDR0.05.narrowPeak.gz | wc -l)
