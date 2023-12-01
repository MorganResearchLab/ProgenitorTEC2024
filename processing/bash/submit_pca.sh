#! /usr/bin/bash

## Make the combined UMAP after batch correction
module load r-4.1.0-gcc-9.3.0-wvnko7v

BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SRC=$(echo $BASE_DIR"/src")

OUTLOG=$(eval 'echo "$BASE_DIR"/logs/Perinatal_PCA.out')
ERROR=$(eval 'echo "$BASE_DIR"/logs/Perinatal_PCA.err')

SCE=$(echo $BASE_DIR"/SCE.dir/Perinatal_SCE-norm.RDS")
OUTPREF=$(echo $BASE_DIR"/SCE.dir/Perinatal")

JOB="bsub -q standard -M 120000 -R "rusage[mem=60000]" -J Pertinatal_PCA -n 1 -T 1500 -o $OUTLOG -e $ERROR  Rscript $SRC/find_hvgs.R --SCE $SCE --dimensions 50 --output $OUTPREF --assay logcounts --center --scale"

echo $JOB
eval $JOB
