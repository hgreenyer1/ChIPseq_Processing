#!/bin/bash

IN_DIR=$1
OUT_DIR=$2

REP1_SUFF="_rep1.Aligned.sortedByCoord.out_loose_peaks.narrowPeak"
REP2_SUFF="_rep2.Aligned.sortedByCoord.out_loose_peaks.narrowPeak"
POOL_SUFF="_pooled.Aligned.sortedByCoord.out_loose_peaks.narrowPeak"
	
for REP1_PEAKS in $IN_DIR/*$REP1_SUFF
do
	PREFIX="$(basename -- ${REP1_PEAKS/"$REP1_SUFF"/""})"
	REP2_PEAKS="$(basename -- ${REP1_PEAKS/"$REP1_SUFF"/"$REP2_SUFF"})"
	POOL_PEAKS="$(basename -- ${REP1_PEAKS/"$REP1_SUFF"/"$POOL_SUFF"})"
	
	sbatch /users/h/g/hgreenye/scripts/chip_scripts/run_IDR.sh $PREFIX $REP1_PEAKS $REP2_PEAKS $POOL_PEAKS
done 	
