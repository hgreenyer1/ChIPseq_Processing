#!/bin/bash 

#call peaks for all files in directory 

# need to include code for the type and then set the command based on the input (broad narrow or loose will call a different script) 
# hopefully eventually I can turn this into a pipeline with dependencies 
# right now set to run for entire directory of bam files 

in_dir=$1
out_dir=$2

for f in $in_dir/*filtered.bam
do 
	sbatch /users/h/g/hgreenye/scripts/chip_scripts/simple_macs2.sh $f $out_dir
done 

