#! /usr/bin/bash

module load r-4.1.1-gcc-9.3.0-jkdw35f

## Make the combined UMAP after batch correction
BASE_DIR="/nfs/research/marioni/mdmorgan/Thymus_Newborn"
SRC=$(echo $BASE_DIR"/src")

OUTLOG=$(eval 'echo "$BASE_DIR"/logs/Perinatal_noise.out')
ERROR=$(eval 'echo "$BASE_DIR"/logs/Perinatal_noise.err')

# The input SCE object will come from the MouseThymusAgeing package
SCE=$(echo $BASE_DIR"/milo.dir/iTEC_Milo.RDS")
OUTPREF=$(echo $BASE_DIR"/milo.dir/Perinatal_Milo-noise.txt")

META=$(echo $BASE_DIR"/milo.dir/iTEC_meta.txt")
NHOOD=$(echo $BASE_DIR"/milo.dir/iTEC_Nhood_map.tsv")

GENES="ENSMUSG00000054263,ENSMUSG00000034634,ENSMUSG00000025432,ENSMUSG00000000682,ENSMUSG00000072423,ENSMUSG00000075602,ENSMUSG00000061353,ENSMUSG00000000731,ENSMUSG00000002057,ENSMUSG00000023235,ENSMUSG00000027314,ENSMUSG00000015396,ENSMUSG00000006179,ENSMUSG00000052759,ENSMUSG00000039691,ENSMUSG00000020096"

JOB="bsub -q standard -M 60000 -R "rusage[mem=42000]" -J Perinatal_noise -n 1 -T 1500 -o $OUTLOG -e $ERROR Rscript $SRC/noise_analysis.R --SCE $SCE --meta $META --nhoods $NHOOD --genes $GENES --output $OUTPREF"

echo $JOB
eval $JOB
