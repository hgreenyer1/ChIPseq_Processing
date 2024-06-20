#!/bin/bash 

# simple script to run MACS2 with no control input 

in_bam=$1
out_dir=$2

prefix="$(basename -- ${in_bam/"rm_multi.bam"/""})"

macs2 callpeak -t $in_bam --outdir $out_dir -n $prefix
