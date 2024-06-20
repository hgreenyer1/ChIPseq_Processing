#!/bin/bash 

#script to filter ChIPseq alignments 

#remove low quality mappings 

in_bam=$1
out_bam=$2

samtools view -b -q 10 $in_bam > $out_bam
