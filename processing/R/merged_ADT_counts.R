#! /usr/bin/env Rscript

#load libraries

library(SingleCellExperiment) 
library(DropletUtils)  
library(scran)
library(scater)
library(Matrix)
library(BiocParallel)
library(biomaRt)

args <- commandArgs(trailingOnly=TRUE) # this should be the regex pattern

base_path <- "/nfs/research/marioni/mdmorgan/Thymus_Newborn/"

# find files
sample_list <- list.files(base_path, pattern="^Newborn_ZsG")
sample_list <- sample_list[grepl(sample_list, pattern="ADT$")]
out_locs <- paste0(base_path, sample_list, "_counts.tsv")
bc_locs <- paste0(base_path, sample_list, "_cellids.tsv")
genes_locs <- paste0(base_path, sample_list, "_features.tsv")

print(out_locs)
print(bc_locs)
print(genes_locs)

intest_matrices <- bplapply(out_locs, readMM)
intest_bcs <- bplapply(bc_locs, function(X) read.table(X, header=TRUE, stringsAsFactors=FALSE)[, 1])
intest_features <- bplapply(genes_locs, function(G) read.table(G, sep="\t", header=FALSE, stringsAsFactors=FALSE)[, 1])

intest_barcode <- list()
for(i in seq_along(sample_list)){
      samp <- sample_list[i]
      intest_barcode[[i]] <- paste0(samp, "_", intest_bcs[[i]])
}

intest_list <- list()
for(i in seq_along(1:length(intest_bcs))){
      mat <- intest_matrices[[i]]
      print(dim(mat))
      print(intest_features[[i]])
      rownames(mat) <- intest_features[[i]][2:5]
      colnames(mat) <- intest_bcs[[i]]
      intest_list[[i]] <- mat
}

# smoosh all matrices together
intest_counts <- do.call(cbind, intest_list)
intest_ids <- do.call(c, intest_bcs)

print(dim(intest_counts))
print(length(intest_ids))

writeMM(intest_counts,
	file="/nfs/research/marioni/mdmorgan/Thymus_Newborn/ALL_counts-ADT.tsv")

write.table(intest_ids,
	file="/nfs/research/marioni/mdmorgan/Thymus_Newborn/ALL_counts_colnames-ADT.tsv",
	sep="\t", row.names=FALSE, quote=FALSE)

write.table(rownames(intest_counts),
	file="/nfs/research/marioni/mdmorgan/Thymus_Newborn/ALL_counts_rownames-ADT.tsv",
	sep="\t", row.names=FALSE, quote=FALSE)