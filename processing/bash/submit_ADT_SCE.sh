#! /usr/bin/bash

## Construct an ADT SCE object
module load r-4.1.0-gcc-9.3.0-wvnko7v

BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SRC=$(echo $BASE_DIR"/src")

MATRIX=$(echo $BASE_DIR"/ALL_counts-ADT.tsv")
CELLS=$(echo $BASE_DIR"/ALL_counts_colnames-ADT.tsv")
GENES=$(echo $BASE_DIR"/ALL_counts_rownames-ADT.tsv")

ERROR=$(echo $BASE_DIR"/logs/ADT_SCE.err")
LOG=$(echo $BASE_DIR"/logs/ADT_SCE.out")

OUTPUT=$(echo $BASE_DIR"/meta.dir/Perinatal_ADT_SCE.RDS")

JOB="bsub -q standard -T 12000 -M 24000 -R "rusage[mem=21000]" -J ADT_SCE -e $ERROR -o $LOG Rscript $SRC/ADT_SCE.R --matrix $MATRIX --cells $CELLS --features $GENES --output $OUTPUT"

echo $JOB
eval $JOB


