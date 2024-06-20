#!/bin/bash 

IN_DIR=$1

for PEAKS in $IN_DIR/*.narrowPeak
do
        PREFIX="$(basename -- ${PEAKS/"pooled_R1_001_sorted.filt.nodup_peaks.narrowPeak"/""})"
        OUT_BED="$PREFIX_pooledPeaks.bed"

        sbatch /users/h/g/hgreenye/scripts/chip_scripts/narrowPeaktoBED.sh $PEAKS $OUT_BED
done


