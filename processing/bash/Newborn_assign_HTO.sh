#! /usr/bin/bash
## Assign HTOs and doublet status to each cell.

SRCDIR="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/src"
LOGFILE="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/logs/Newborn_HTO_assign"

JOB="bsub -R "rusage[mem=40000]" -T 1500 -q research-rh74 -M 80000 -e $LOGFILE.err -o $LOGFILE.out Rscript $SRCDIR/ZsG_10X_BCassign.R"

echo $JOB
eval $JOB
