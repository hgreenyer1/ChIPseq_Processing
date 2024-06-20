#!/bin/bash

#SBATCH --partition=bluemoon
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=21G
#SBATCH --time=30:00:00

#================================
# Simple script to call loose peaks with MACS2 for IDR and rep comparison 
# Generates loose narrowPeaks and bdg
#================================

INPUT_BAM=$1
TREAT_BAM=$2
OUT_DIR=$3
QVAL="2e-2"

PREFIX="$(basename -- ${TREAT_BAM/".bam"/"_loose"})"

# call loose narrow peaks

macs2 callpeak -t $TREAT_BAM -c $INPUT_BAM -g "hs" --outdir $OUT_DIR -n $PREFIX -q $QVAL --bdg

