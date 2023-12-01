#! /usr/bin/bash

## Run mylo on a subset of the data
module load r-4.1.0-gcc-9.3.0-wvnko7v

BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SRC=$(echo $BASE_DIR"/src")

OUTLOG=$(eval 'echo "$BASE_DIR"/logs/Pertinatal_GOI_genes.out')
ERROR=$(eval 'echo "$BASE_DIR"/logs/Pertinatal_GOI_genes.err')

SCE=$(echo $BASE_DIR"/SCE.dir/Perinatal_MultiSCE.RDS")
OUTPREF=$(echo $BASE_DIR"/meta.dir/Perinatal")

JOB="bsub -q standard -M 160000 -R "rusage[mem=50000]" -J Perinatal_GOI_genes -n 1 -T 1500 -o $OUTLOG -e $ERROR Rscript $SRC/select_genes.R --SCE $SCE --assay logcounts --output $OUTPREF"

echo $JOB
eval $JOB
