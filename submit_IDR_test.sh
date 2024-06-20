IN_DIR=$1

REP1_SUFF="R1_CKDL230020107-1A_pooled_R1_001_sorted.filt.nodup_summits.bed"
REP2_SUFF="R2_CKDL230020107-1A_pooled_R1_001_sorted.filt.nodup_summits.bed"
POOL_SUFF="pooled_summits.bed"

for REP1_PEAKS in $IN_DIR/*$REP1_SUFF
do
        PREFIX="$(basename -- ${REP1_PEAKS/"$REP1_SUFF"/""})"
        REP2_PEAKS="$(basename -- ${REP1_PEAKS/"$REP1_SUFF"/"$REP2_SUFF"})"
        POOL_PEAKS="$(basename -- ${REP1_PEAKS/"$REP1_SUFF"/"$POOL_SUFF"})"

        sbatch /users/h/g/hgreenye/scripts/chip_scripts/run_IDR_test.sh $PREFIX $REP1_PEAKS $REP2_PEAKS $POOL_PEAKS
done
