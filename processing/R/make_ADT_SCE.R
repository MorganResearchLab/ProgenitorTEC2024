#! /usr/bin/env Rscript

## make an SCE object of multi-modal data from 2 separate SCE objects
library(SingleCellExperiment)
library(scuttle)
library(Matrix)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCElist"), type="character",
       	             help="A comma-separated list of input SCE objects")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file name for multi-experiment object")

parser <- add_option(parser, c("-i", "--index"), type="numeric",
                     help="Index of SCE object that should be the main experiment data")

parser <- add_option(parser, c("-n", "--names"), type="character",
                     help="A comma-separated list of N-1 names to add to the additional experiments")

opt <- parse_args(parser)

sce.file.list <- unlist(strsplit(opt$SCElist, split=",", fixed=TRUE))
message("Found ", length(sce.file.list), " SCE objects to combine")

message("Reading in SCE objects")
sce.list <- lapply(sce.file.list,
	           FUN=function(X) readRDS(X))

# fix ADT vs. HTO naming issue
print(lapply(sce.list, assayNames))


# check they have the same colnames and ordering
comm.cells <- Reduce(x=sce.list,
	             f=function(x, y) intersect(colnames(x), colnames(y)))

if(length(comm.cells) == 0){
  print(head(colnames(sce.list[[1]])))
  print(head(colnames(sce.list[[2]])))
  stop("No cells in common - check cell barcodes")
}

message("Identified ", length(comm.cells), " common cell barcodes")
main.sce <- sce.list[[opt$index]]
main.sce <- main.sce[, colnames(main.sce) %in% comm.cells]
#rownames(main.sce) <- rowData(main.sce)$ensembl_gene_id
#rowData(main.sce) <- NULL

alt.names <- unlist(strsplit(opt$names, split=",", fixed=TRUE))
for(x in seq_along(alt.names)){
  x.name <- alt.names[x]
  message("Adding ", x.name, " to SCE object")
  altExp(main.sce, x.name) <- sce.list[[x+1]][, colnames(main.sce)]
}

saveRDS(main.sce, file=opt$output)

message("All done")