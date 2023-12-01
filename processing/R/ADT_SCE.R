#! /usr/bin/env Rscript

## Take in a ADT counts matrix, feature names and cell IDs
# cpm normalise and create an SCE object

library(optparse)
library(SingleCellExperiment)
library(Matrix)

parser <- OptionParser()
parser <- add_option(parser, c("-m", "--matrix"), type="character",
       	             help="A counts matrix")

parser <- add_option(parser, c("-c", "--cells"), type="character",
                     help="Files of cell IDs")

parser <- add_option(parser, c("-f", "--features"), type="character",
                     help="A file of feature names")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file name for SCE object")

opt <- parse_args(parser)

message("Reading in ", opt$matrix)
in.matrix <- readMM(opt$matrix)

if(!class(in.matrix) %in% c("dgCMatrix")){
    in.matrix <- as(in.matrix, "dgCMatrix")
}

print(colnames(in.matrix))
print(rownames(in.matrix))

cell.ids <- gsub(read.table(opt$cells, sep="\t", header=TRUE, stringsAsFactors=FALSE)[, 1],
	         pattern="ADT",  replacement="HTO")
features <- read.table(opt$features, sep="\t", header=TRUE, stringsAsFactors=FALSE)[, 1]

print(dim(in.matrix))
print(length(features))
print(length(cell.ids))

print(head(features))
print(tail(cell.ids))

rownames(in.matrix) <- features
colnames(in.matrix) <- cell.ids

message("Making SCE")
adt.sce <- SingleCellExperiment(assay=list(counts=in.matrix))

message("Computing CPMs")
intest_cpm <- as.data.frame(apply(in.matrix, 2, FUN=function(X) log10((X+1)/((sum(X)+1)/(1e6+1)))))
rownames(intest_cpm) <- rownames(adt.sce)
colnames(intest_cpm) <- colnames(adt.sce)

assay(adt.sce, "cpm") <- intest_cpm

saveRDS(adt.sce, file=opt$output)
message("All done")