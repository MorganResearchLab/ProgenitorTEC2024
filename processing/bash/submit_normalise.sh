#! /usr/bin/bash

module load r-4.1.0-gcc-9.3.0-wvnko7v

BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SRC=$(echo $BASE_DIR"/src")

## Submit jobs for normalising ZsG data
ERROR=$(echo $BASE_DIR"/logs/Perinatal_norm_SCE.err")
OUTLOG=$(echo $BASE_DIR"/logs/Perinatal_norm_SCE.out")

# define input and output files
SCE=$(echo $BASE_DIR"/SCE.dir/Perinatal_SCE.RDS")
OUTPUT=$(echo $BASE_DIR"/SCE.dir/Perinatal_SCE-norm.RDS")

HTO_TASK="R CMD $SRC/norm_counts.R --SCE $SCE --breaks --output $OUTPUT"
HTO_JOB="bsub -R "rusage[mem=80000]" -M 170000 -T 5000 -q standard -J Perinatal_norm -e $ERROR -o $OUTLOG $HTO_TASK"

echo $HTO_JOB
eval $HTO_JOB

