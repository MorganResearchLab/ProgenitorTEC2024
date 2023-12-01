#! /usr/bin/bash

## loop over sample and run emptyDrops + cell filtering on each sample separately
SRCDIR="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/src"
LOGDIR="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/logs"
SAMPDIR="/nfs/research1/marioni/mdmorgan/Thymus_Newborn"

for fle in `ls $SAMPDIR | grep "^Newborn*" | grep "ADT" | grep -v "tsv"`;
do
    echo $fle
    LOGERROR=$LOGDIR/$fle-emptyDrops.err
    LOGOUT=$LOGDIR/$fle-emptyDrops.out

    JOB="bsub -R "rusage[mem=32000]" -M 64000 -T 15000 -q research-rh74 -J $fle -e $LOGERROR -o $LOGOUT R CMD $SRCDIR/ZsG_emptyDrops_ADT.R $fle"
    echo $JOB
    eval $JOB
done
