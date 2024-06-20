#!/bin/bash

#SBATCH --partition=bluemoon
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=21G
#SBATCH --time=30:00:00

#================================
# Simple script to call broadPeaks with MACS2
# Generates broadPeaks and bdg
#================================

INPUT_BAM=$1
TREAT_BAM=$2
OUT_DIR=$3
QVAL="5e-3"
BROADCUTOFF="1e-1"

PREFIX="$(basename -- ${TREAT_BAM/".bam"/"_broad"})"

# call broad peaks

macs2 callpeak -t $TREAT_BAM -c $INPUT_BAM -g "hs" --outdir $OUT_DIR -n $PREFIX -q $QVAL --broad --broad-cutoff $BROADCUTOFF --bdg

