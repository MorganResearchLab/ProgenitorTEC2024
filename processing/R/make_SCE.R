#! /usr/bin/env Rscript

## Given an input dataset of multiple batches - normalize across batches

library(DropletUtils)
library(SingleCellExperiment)
library(scran)
library(scater)
library(optparse)
library(Matrix)
library(biomaRt)


parser <- OptionParser()
parser <- add_option(parser, c("-c", "--counts"), type="character",
       	             help="The path to a mtx format counts matrix")

parser <- add_option(parser, c("-i", "--ids"), type="character",
                     help="The path to a file containing cell IDs")

parser <- add_option(parser, c("-f", "--features"), type="character",
                     help="The path to a file containing feature names")

parser <- add_option(parser, c("-t", "--hto"), type="character",
                     help="The path to a file containing HTO data per cell")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output combined SCE object")

opt <- parse_args(parser)

message(paste0("Reading in counts matrix: ", opt$counts))
count.mat <- readMM(opt$counts)

message(paste0("Reading in cell IDs: ", opt$ids))
bcs <- read.table(opt$ids, header=TRUE, stringsAsFactors=FALSE, sep="\t")
colnames(count.mat) <- bcs[, 1]

message(paste0("Reading in feature names: ", opt$features))
features <- read.table(opt$features, header=TRUE, stringsAsFactors=FALSE, sep="\t")
rownames(count.mat) <- features[, 1]

message("Making SCE object")
sce <- SingleCellExperiment(list(counts=count.mat))

message(paste0("Reading in HTO data: ", opt$hto))
hto <- read.table(opt$hto, header=TRUE, sep="\t", stringsAsFactors=FALSE)
rownames(hto) <- hto$ID
hto <- hto[colnames(sce), ]

colData(sce)$HTO <- hto$HTO
colData(sce)$Doublet.Class <- hto$Class

biomaRt.connection <- useMart("ensembl", "mmusculus_gene_ensembl", host="uswest.ensembl.org")

gene.df <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 mart=biomaRt.connection)
rownames(gene.df) <- gene.df$ensembl_gene_id

# also find doublets based on T cell genes
tcell.genes <- c("Cd3g", "Cd3e", "Lck", "Zap70", "Fyn", "Trgv2")
tcell.ensembl <- gene.df[gene.df$external_gene_name %in% tcell.genes,]$ensembl_gene_id
any.tcells <- colSums2(counts(sce)[rownames(sce) %in% tcell.ensembl, ] > 0) > 0

message("Dropping ", sum(any.tcells), " suspected thymocyte doublets")
sce <- sce[, !any.tcells]

n.doublets <- sum(!grepl(colData(sce)$Doublet.Class, pattern="Singlet"))
message("Droping ", n.doublets, " hash-tag doublets")
sce <- sce[, grepl(colData(sce)$Doublet.Class, pattern="Singlet")]

message(paste0("Writing SCE to file: ", opt$output))
saveRDS(sce, file=opt$output)

message("All done")
