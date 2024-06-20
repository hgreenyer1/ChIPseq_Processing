#!/bin/bash 

#filter all ChIPseq bam files in a given directory 

#suf_in should be .bam or _sorted.bam depending on if merging was needed 

suf_in=$3
suf_out="_filtered.bam"

in_dir=$1
out_dir=$2

for F in $in_dir/*$suf_in
do
        F_out="$out_dir/$(basename -- ${F/".bam"/$suf_out})"
        sbatch /users/h/g/hgreenye/scripts/chip_scripts/filter_bam.sh $F $F_out
done



