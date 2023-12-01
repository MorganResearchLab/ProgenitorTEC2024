#! /usr/bin/bash

module load r-4.1.1-gcc-9.3.0-jkdw35f 

## loop over sample and run emptyDrops + cell filtering on each sample separately
BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SRCDIR=$BASE_DIR"/src"
LOGDIR=$BASE_DIR"/logs"
SAMPDIR=$BASE_DIR

for fle in `ls $SAMPDIR | grep "^Newborn*" | grep "HTO" | grep -v "tsv"`;
do
    echo $fle
    LOGERROR=$LOGDIR/$fle-emptyDrops.err
    LOGOUT=$LOGDIR/$fle-emptyDrops.out

    JOB="bsub -R "rusage[mem=54000]" -M 74000 -T 15000 -q standard  -J $fle -e $LOGERROR -o $LOGOUT Rscript $SRCDIR/ZsG_emptyDrops_HTO.R $fle"
    echo $JOB
    eval $JOB
done
