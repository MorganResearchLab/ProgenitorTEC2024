#! /usr/bin/bash

## Submit jobs for normalising ZsG data
SRCDIR="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/src"

HTO_TASK="R CMD $SRCDIR/ZsG_10X_normalisation_HTO.R"
HTO_JOB="bsub -R "rusage[mem=96000]" -M 180000 -T 50000 -q research-rh74 -J HTO_norm -e logs/HTO_norm.err -o logs/HTO_norm.out $HTO_TASK"

#ADT_TASK="R CMD $SRCDIR/ZsG_10X_normalisation_ADT.R"
#ADT_JOB="bsub -R "rusage[mem=150000]" -M 250000 -T 50000 -q research-rh74 -J ADT_norm -e logs/ADT_norm.err -o logs/ADT_norm.out $ADT_TASK"

echo $HTO_JOB
#echo $ADT_JOB
eval $HTO_JOB
#$ADT_JOB
