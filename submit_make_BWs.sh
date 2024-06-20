#!/bin/bash 

#=========
#generate FE bigwigs for all peaks in a given directory
#==========

IN_DIR=$1
OUT_DIR=$2
SUFF=$3


for PEAKS in $IN_DIR/*$SUFF
do
	PREFIX=$(basename $PEAKS $SUFF)
        sbatch /users/h/g/hgreenye/scripts/chip_scripts/run_make_BWs.sh $PREFIX $IN_DIR $OUT_DIR
done
