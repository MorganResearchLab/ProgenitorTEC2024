#! /usr/bin/bash

## Submit a 10X Genomics cellranger job using feature barcoding, i.e. CITE-seq
## 1) input libraries CSV
## 2) feature barcode reference CSV
## 3) Run ID

SRCDIR="/nfs/research1/marioni/software/cellranger-3.1.0"
REFDATA="/nfs/research1/marioni/genomes/refdata-cellranger-mm10-3.0.0"

# path to LSF cluster mode template file
TEMPLATEPATH="/nfs/research1/marioni/software/cellranger-3.1.0/martian-cs/v3.2.3/jobmanagers/lsf.template"

$SRCDIR/cellranger count --id=$3  --transcriptome=$REFDATA --libraries=$1 --feature-ref=$2 --localcores=180 --localmem=24 --chemistry=SC3Pv3 --nosecondary  #--jobmode=$TEMPLATEPATH

