#!/bin/bash 

#SBATCH --partition=bluemoon
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=21G
#SBATCH --time=30:00:00

REP1_PEAKS=$1
REP2_PEAKS=$2
POOLED_PEAKS=$3

cat $REP1_PEAKS $REP2_PEAKS > $POOLED_PEAKS

sort-bed --max-mem 19G $POOLED_PEAKS > $POOLED_PEAKS.tmp

rm $POOLED_PEAKS

mv $POOLED_PEAKS.tmp $POOLED_PEAKS
