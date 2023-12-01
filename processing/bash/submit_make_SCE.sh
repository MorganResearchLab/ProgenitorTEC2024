#! /usr/bin/bash

module load r-4.1.0-gcc-9.3.0-wvnko7v

BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SRC=$(echo $BASE_DIR"/src")

## Submit jobs for normalising Pertinatal data
ERROR=$(echo $BASE_DIR"/logs/Perinatal_make_SCE.err")
OUTLOG=$(echo $BASE_DIR"/logs/Perinatal_make_SCE.out")

# define input and output files
COUNTS=$(echo $BASE_DIR"/ALL_counts-HTO.tsv")
ID=$(echo $BASE_DIR"/ALL_counts_colnames-HTO.tsv")
GENES=$(echo $BASE_DIR"/ALL_counts_rownames-HTO.tsv")
HTO=$(echo $BASE_DIR"/HTO_class.tsv")

OUTPUT=$(echo $BASE_DIR"/SCE.dir/Perinatal_SCE.RDS")

HTO_TASK="R CMD $SRC/make_SCE.R --counts $COUNTS --ids $ID --features $GENES --hto $HTO --output $OUTPUT"
HTO_JOB="bsub -R "rusage[mem=80000]" -M 170000 -T 5000 -q standard -J Perinatal_makeSCE -e $ERROR -o $OUTLOG $HTO_TASK"

echo $HTO_JOB
eval $HTO_JOB

