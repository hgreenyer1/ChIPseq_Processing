#!/bin/bash

IN_DIR=$1
OUT_DIR=$2

IN_BAM=$3
SUFF=$4


for TREAT_BAM in $IN_DIR/*$SUFF
do
        sbatch /users/h/g/hgreenye/scripts/chip_scripts/run_macs2_broad.sh $IN_BAM $TREAT_BAM $OUT_DIR
done

