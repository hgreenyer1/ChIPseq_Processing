#!/bin/bash

#=============================
# Simple script to generate FE bigwigs from macs2 bdgs 
#==============================

#SBATCH --partition=bluemoon
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=21G
#SBATCH --time=30:00:00

PREFIX=$1
IN_DIR=$2
OUT_DIR=$3

CHR_SIZES="/users/h/g/hgreenye/INDEXES/HG38/BOWTIE_INDEX/GRCh38_noalt_as/hg38NameLength.txt"

TREATMENT=$IN_DIR/$PREFIX"_macs2_treat_pileup.bdg"
CONTROL=$IN_DIR/$PREFIX"_macs2_control_lambda.bdg"

METHOD=FE

echo $TREATMENT
echo $CONTROL
	
CMP_BDG=$OUT_DIR/$(basename $TREATMENT)
CMP_BDG="${CMP_BDG/_treat_pileup.bdg/}"_"$METHOD".bdg

CMP_BDG_CLIP="${CMP_BDG}.clip"
CMP_BDG_TMP="${CMP_BDG}.tmp"
CMP_BW=${CMP_BDG/.bdg/.bw}

echo "$CMP_BDG"
echo "$CMP_BW"

#get FE bed
macs2 bdgcmp -t $TREATMENT -c $CONTROL -m $METHOD -o $CMP_BDG

#echo "finished bdgcmp"

#clip bed
bedClip $CMP_BDG $CHR_SIZES $CMP_BDG_CLIP

sort-bed --max-mem 19G $CMP_BDG_CLIP > $CMP_BDG_TMP

rm $CMP_BDG
rm $CMP_BDG_CLIP

mv $CMP_BDG_TMP $CMP_BDG

#echo "finished bed-sort"

bedGraphToBigWig $CMP_BDG $CHR_SIZES $CMP_BW

echo "DONE!"


