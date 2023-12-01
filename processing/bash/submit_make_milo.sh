#! /usr/bin/bash

module load r-4.1.0-gcc-9.3.0-wvnko7v

## Make the combined UMAP after batch correction
BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SRC=$(echo $BASE_DIR"/src")

OUTLOG=$(eval 'echo "$BASE_DIR"/logs/Perinatal_milo.out')
ERROR=$(eval 'echo "$BASE_DIR"/logs/Perinatal_milo.err')

# The input SCE object will come from the MouseThymusAgeing package
SCE=$(echo $BASE_DIR"/SCE.dir/Perinatal_SCE-FDL.RDS")
OUTPREF=$(echo $BASE_DIR"/milo.dir/Perinatal_Milo.RDS")

JOB="bsub -q standard -M 60000 -R "rusage[mem=42000]" -J Perinatal_milo -n 1 -T 1500 -o $OUTLOG -e $ERROR Rscript $SRC/make_milo.R --SCE $SCE --sample HTO --reduced 'PCA' --knn 50 --props 0.05 --dimensions 50 --overlap 10 --output $OUTPREF"

echo $JOB
eval $JOB
