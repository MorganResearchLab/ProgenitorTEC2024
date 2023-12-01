#! /usr/bin/bash

module load r-4.1.0-gcc-9.3.0-wvnko7v

BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SHORT_SCE="/nfs/research/marioni/mdmorgan/Thymus_ShortKGF/SCE.dir"
LONG_SCE="/nfs/research/marioni/mdmorgan/Thymus_KGF/SCE.dir"
ZSG_SCE="/nfs/research/marioni/mdmorgan/Thymus_ZsG/SCE.dir"

SRC=$(echo $BASE_DIR"/src")
SCE_LIST=$(echo $BASE_DIR"/SCE.dir/Perinatal_MultiSCE.RDS,"$SHORT_SCE"/ShortKGF_MultiSCE.RDS,"$LONG_SCE"/KGF_MultiSCE.RDS,"$ZSG_SCE"/ZsG_MultiSCE.RDS")

ERROR=$(echo $BASE_DIR"/logs/Combine_data.err")
LOG=$(echo $BASE_DIR"/logs/Combine_data.out")

OUTPUT=$(echo $BASE_DIR"/Integration/All_Thymus_Combined")

JOB="bsub -q bigmem -T 12000 -M 118000 -R "rusage[mem=106000]" -J Combine_data -e $ERROR -o $LOG Rscript $SRC/integrate_data.R --SCEList $SCE_LIST --dimensions 50 --knn 30 --l2norm --output $OUTPUT"

echo $JOB
eval $JOB
