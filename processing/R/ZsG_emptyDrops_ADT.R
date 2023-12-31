#! /usr/bin/env Rscript

# This script takes in all of the 10X experiments calls cells and normalizes them all together
# Hopefully this will help to mitigate against any major batch effects between the different experiments.
# Due to the size of the data I need to run the cell calling and filtering on each sample separately.

#load libraries
library(SingleCellExperiment)
library(DropletUtils)
library(scran)
library(scater)
library(Matrix)
library(BiocParallel)
library(biomaRt)
source("/nfs/research1/marioni/mdmorgan/Thymus_Newborn/src/single_cell_functions.R")

# should only be 1 command line argument - sample name
args <- commandArgs(trailingOnly=TRUE)

message(paste0("Running cell calling on sample: ", args[1]))
# rather than relying on the CellRanger to call cells I'll use emptyDrops and remove poor quality cells at a later point.
ncores <- 4
mcparam <- MulticoreParam(workers=ncores)
register(mcparam)

base_path <- "/nfs/research1/marioni/mdmorgan/Thymus_Newborn/"
sample_list <- list.files(base_path, pattern="^Newborn")
sample_list <- sample_list[grepl(sample_list, pattern=args[1])]
sample_list <- sample_list[!grepl(sample_list, pattern="tsv")]

# check that only a single sample found
message(paste0("Found ", length(sample_list), " samples"))

# locations of hdf5 molecule information across all samples
mol_locs <- paste0(base_path, sample_list, "/outs/molecule_info.h5")
out_locs <- paste0(base_path, sample_list, "/outs/raw_feature_bc_matrix/matrix.mtx.gz")
bc_locs <- paste0(base_path, sample_list, "/outs/raw_feature_bc_matrix/barcodes.tsv.gz")
genes_locs <- paste0(base_path, sample_list, "/outs/raw_feature_bc_matrix/features.tsv.gz")

message("Reading in unfiltered matrices")
intest_matrices <- bplapply(out_locs, readMM)
intest_bcs <- bplapply(bc_locs, function(X) read.table(X, header=FALSE, stringsAsFactors=FALSE)[, 1])
# extract the genes and features

message("Reading in features")
intest_features<- bplapply(genes_locs, function(G) read.table(G, sep="\t", header=FALSE, stringsAsFactors=FALSE))
intest_genes <- bplapply(intest_features, function(G) G[grepl(G[, 3], pattern="Gene Expression"), 1])
intest_hto <- bplapply(intest_features, function(G) G[grepl(G[, 3], pattern="Antibody Capture"), 1])

message("Adding sample ID to cell barcodes")
# add sample ID to barcodes
intest_barcode <- list()
for(i in seq_along(sample_list)){
  samp <- sample_list[i]
  intest_barcode[[i]] <- paste0(samp, "_", intest_bcs[[i]])
}

message("Running EmptyDrops to estimate cells")
# call cells with emptyDrops - these matrices include the HTO as well - will this be a problem?
lower.umi <- 50
intest_calls <- lapply(intest_matrices, emptyDrops, niters=20000, ignore=4999, BPPARAM=mcparam, lower=lower.umi, retain=Inf)

message("Identifying cells at 1% FDR")
# identify cells using a FDR 1%
sig_cells <- lapply(intest_calls, function(X) X$FDR <= 0.01 & !is.na(X$FDR))

# subset the called cells - the numbers are very similar to what CellRanger gives
intest_cells <- lapply(1:length(intest_barcode), function(i) intest_matrices[[i]][, sig_cells[[i]]])
cell_barcodes <- lapply(1:length(intest_barcode), function(x) intest_barcode[[x]][sig_cells[[x]]])

# apply gene names to matrix rows and subset the common genes
#common.genes <- intersect(intest_genes[[1]], intest_genes[[2]])
common.genes <- intest_genes[[1]]
hto_tags <- unique(unlist(intest_hto))

message("Subsetting matrices to called cells")
intest_list <- list()
for(i in seq_along(1:length(cell_barcodes))){
  mat <- intest_cells[[i]]
  rownames(mat) <- intest_features[[i]][, 1]
  intest_list[[i]] <- mat[common.genes, ]
}

message("Writing out HTO feature matrices")
# need to figure out a way to call a cell as one sample or the other based on the relative counts after normalisation
# I'll just write out the HTO/ADT tags for each sample separately
intest_feat.mat <- list()
for(i in seq_along(1:length(cell_barcodes))){
  mat <- intest_cells[[i]]
  rownames(mat) <- intest_features[[i]][, 1]
  colnames(mat) <- cell_barcodes[[i]]
  
  # write out the HTO matrix - grab the sample name from first barcode
  sname <- unlist(lapply(strsplit(cell_barcodes[[i]][1], fixed=TRUE, split="_"), FUN=function(C) paste(C[1:5], collapse="_")))
  fname <- paste0(base_path, sname)
  writeMM(mat[rownames(mat) %in% hto_tags, ], 
          file=paste0(fname, "_counts.tsv"))
  write.table(rownames(mat[rownames(mat) %in% hto_tags, ]),
              file=paste0(fname, "_features.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  write.table(colnames(mat),
              file=paste0(fname, "_cellids.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE)
  intest_feat.mat[[i]] <- mat[rownames(mat) %in% hto_tags, ]
}

message("Combining together counts matrices")
# smoosh all matrices together
intest_counts <- do.call(cbind, intest_list)
intest_ids <- do.call(c, cell_barcodes)

# remove supefluous matrices
rm(list=c("intest_matrices", "intest_cells"))
gc()

######################
## Single-cell QC
######################
## Remove very low complexity cell libraries, i.e. ones that express < 1000 genes
lib.sizes <- Matrix::colSums(intest_counts)
nonzero.genes <- Matrix::colSums(intest_counts > 0)

# remove cells with < 1000 genes
#sum(nonzero.genes > 1000)

message("Selecting cells with > 1000 UMIs")
intest_counts <- intest_counts[, nonzero.genes > 1000]
intest_ids <- intest_ids[nonzero.genes > 1000]

## remove cells with excessively high mitochondrial content
biomaRt.connection <- useMart("ensembl", "mmusculus_gene_ensembl")
gene.df <- getBM(attributes = c("ensembl_gene_id", "chromosome_name"), filters="ensembl_gene_id", 
                 values=rownames(intest_counts), mart=biomaRt.connection)
rownames(gene.df) <- gene.df$ensembl_gene_id

mt.counts <- intest_counts[which(common.genes %in% gene.df$ensembl_gene_id[gene.df$chromosome_name == "MT"]), ]
mt.fraction <- Matrix::colSums(mt.counts)/Matrix::colSums(intest_counts)

message("Filtering cells with high MT content")
# fit a median-centred, MAD variance model normal to dervive p-values for outlier proportions
# this isn't the best calibrated model, so might be a little liberal
mt.p <- pnorm(mt.fraction, mean=median(mt.fraction), sd=mad(mt.fraction)*2, lower.tail=FALSE)
mt.lim <- min(mt.fraction[which(p.adjust(mt.p, method="fdr") < 0.05)])

# cells with a high mitochondrial fraction > `mt.lim` are removed as outliers.
message(paste0("Removing ", sum(mt.fraction >= mt.lim), " cells for high MT content"))
intest_counts <- intest_counts[, mt.fraction < mt.lim]
intest_ids <- intest_ids[mt.fraction < mt.lim]

message("Writing unnormalised counts matrix to file")
# Now we have fairly decent cells I'll proceed to the normalisation using deconvolution-estimated size factors
writeMM(intest_counts,
        file=paste0("/nfs/research1/marioni/mdmorgan/Thymus_Newborn/", args[1], "_counts-ADT.tsv"))
write.table(intest_ids,
            file=paste0("/nfs/research1/marioni/mdmorgan/Thymus_Newborn/", args[1], "_counts_colnames-ADT.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)
write.table(rownames(intest_counts),
            file=paste0("/nfs/research1/marioni/mdmorgan/Thymus_Newborn/", args[1], "_counts_rownames-ADT.tsv"),
            sep="\t", row.names=FALSE, quote=FALSE)

