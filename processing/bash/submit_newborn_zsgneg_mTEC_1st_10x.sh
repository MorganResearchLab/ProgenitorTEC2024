#! /usr/bin/bash

## submit a 10X cell ranger job for the ageing ZsGreen-lineage tracing experiment 1st sample
LOG="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/logs/"
SRCDIR="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/src"

HTO_LIBFILE="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/lib_files/Newborn_ZsGneg_1stRun_mTEC_HTO.tsv"  
ADT_LIBFILE="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/lib_files/Newborn_ZsGneg_1stRun_mTEC_ADT.tsv"

ADT_REFLIB="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/refs/Newborn_ZsGneg_1stRun_mTEC_ADT_refs.tsv"
HTO_REFLIB="/nfs/research1/marioni/mdmorgan/Thymus_Newborn/refs/Newborn_ZsGneg_1stRun_mTEC_HTO_refs.tsv"

# I need to run the feature barcoding for the CITE-seq antibodies and hashtags separately

ADT_ID="Newborn_ZsGneg_1stRun_mTEC_ADT"
HTO_ID="Newborn_ZsGneg_1stRun_mTEC_HTO"

HTO_JOB="bsub -R "rusage[mem=24000]" -M 36000 -W 168:00  -T 50000 -q research-rh74 -J $HTO_ID -e $LOG.$HTO_ID.err -o $LOG.$HTO_ID.out bash $SRCDIR/submit_10X_cellranger.sh $HTO_LIBFILE $HTO_REFLIB $HTO_ID"
echo $HTO_JOB
eval $HTO_JOB

ADT_JOB="bsub -R "rusage[mem=24000]" -M 36000 -W 168:00 -T 50000 -q research-rh74 -J $ADT_ID -e $LOG.$ADT_ID.err -o $LOG.$ADT_ID.out bash $SRCDIR/submit_10X_cellranger.sh $ADT_LIBFILE $ADT_REFLIB $ADT_ID"

echo $ADT_JOB
eval $ADT_JOB
