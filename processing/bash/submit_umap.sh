#! /usr/bin/bash

## Make the combined UMAP after batch correction
module load r-4.1.0-gcc-9.3.0-wvnko7v

BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SRC=$(echo $BASE_DIR"/src")

OUTLOG=$(eval 'echo "$BASE_DIR"/logs/Perinatal_UMAP.out')
ERROR=$(eval 'echo "$BASE_DIR"/logs/Perinatal_UMAP.err')

#SCE="package"
SCE=$(echo $BASE_DIR"/SCE.dir/Perinatal-PCA.RDS")
OUTPREF=$(echo $BASE_DIR"/SCE.dir/Perinatal")

JOB="bsub -q standard -M 42000 -R "rusage[mem=32000]" -J Perinatal_UMAP -n 1 -T 1500 -o $OUTLOG -e $ERROR Rscript $SRC/make_SCE_umap.R --SCE $SCE --dimensions 50 --reduced PCA --output $OUTPREF"

echo $JOB
eval $JOB
