#! /usr/bin/bash

module load r-4.1.0-gcc-9.3.0-wvnko7v

BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_All"
SRC=$(echo $BASE_DIR"/src")

MILO="/nfs/research/marioni/mdmorgan/Thymus_Newborn/milo.dir/Perinatal_Milo.RDS,/nfs/research/marioni/mdmorgan/Thymus_ZsG/milo.dir/ZsG_Milo.RDS"

ERR=$(echo $BASE_DIR"/logs/NhoodCorr.err")
OUTLOG=$(echo $BASE_DIR"/logs/NhoodCorr.out")

OUTPREF=$(echo $BASE_DIR"/nhood.dir/Perinatal_Ageing")

JOB="bsub -q standard -M 56000 -R "rusage[mem=52000]" -T 1200 -e $ERR -o $OUTLOG -J Newborn_Age_NhoodCor Rscript $SRC/comparse_nhoods.R --milo $MILO --assay 'logcounts' --output $OUTPREF"

echo $JOB
eval $JOB
