#! /usr/bin/env Rscript

## Given an input dataset of 2 batches - downsample the batches to the same median sequencing depth

library(DropletUtils)
library(SingleCellExperiment)
library(scran)
library(scater)
library(irlba)
library(optparse)
library(Matrix)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the combined SCE object")

parser <- add_option(parser, c("-d", "--dimensions"), type="numeric",
                     help="The number of dimensions to use for batch correction")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files denoised data and combined SCE object")

parser <- add_option(parser, c("-c", "--center"), action="store_true",
                     help="Flag to center data before PCA")

parser <- add_option(parser, c("-l", "--scale"), action="store_true",
                     help="Flag to scale data before PCA")

parser <- add_option(parser, c("-a", "--assay"), type="character",
                     help="Gene expression measurements to use for PCA")

opt <- parse_args(parser)

message(paste0("Reading in ", opt$SCE))
big.sce <- readRDS(opt$SCE)

message("Computing highly variable genes")
hvg.stats <- modelGeneVar(big.sce)
hvg.stats$FDR[is.na(hvg.stats$FDR)] <- 1
hvg.stats$GeneID <- rownames(hvg.stats)
write.table(hvg.stats, paste0(opt$output, "_hvgs.txt"), sep="\t", quote=FALSE, row.names=FALSE)  

common.hvg <- hvg.stats$FDR < 0.1

## define based on HVGs
message(paste0("Performing joint PCA with ", sum(common.hvg), " HVGs"))
all.pca <- prcomp_irlba(t(assay(big.sce[common.hvg ,], opt$assay)), n=opt$dimensions + 1, .scale=opt$scale, center=opt$center)

out.pca <- as.data.frame(all.pca$x)
out.pca$CellID <- colnames(big.sce)
write.table(out.pca, file=paste0(opt$output, "-PCA.tsv"), sep="\t", quote=FALSE)

reducedDim(big.sce, "PCA") <- all.pca$x

ofile.sce <- paste0(opt$output, "-PCA.RDS")

message(paste0("Saving SCE to: ", ofile.sce))
saveRDS(big.sce, file=ofile.sce)

message("All done")
