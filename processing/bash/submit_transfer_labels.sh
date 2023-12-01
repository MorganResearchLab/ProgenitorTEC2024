#! /usr/bin/bash

module load r-4.1.0-gcc-9.3.0-wvnko7v

BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SRC=$(echo $BASE_DIR"/src")
#QUERY="package"
QUERY=$(echo $BASE_DIR"/SCE.dir/Perinatal_SCE-UMAP.RDS")
REFERENCE=$(echo $BASE_DIR"/Integration/AgeThymus_SCE.RDS")

ERROR=$(echo $BASE_DIR"/logs/Transfer_label.err")
LOG=$(echo $BASE_DIR"/logs/Transfer_label.out")

OUTPUT=$(echo $BASE_DIR"/Integration/Perinatal_SMARTcombined")
COL="Cluster"
HVG=$(echo $BASE_DIR"/Integration/Thymus_HVG.tsv")

JOB="bsub -q standard -T 12000 -M 38000 -R "rusage[mem=36000]" -J Transfer_label -e $ERROR -o $LOG Rscript $SRC/combine_data.R --reference $REFERENCE --query $QUERY --column $COL --dimensions 50 --knn 30 --l2norm --hvg $HVG --output $OUTPUT"

echo $JOB
eval $JOB
