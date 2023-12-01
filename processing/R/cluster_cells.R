#! /usr/bin/env Rscript

## apply graph clustering to kNN-graph of single-cells

library(igraph)
library(SingleCellExperiment)
library(scran)
library(scater)
library(scuttle)
library(optparse)
library(Matrix)

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the combined SCE object")

parser <- add_option(parser, c("-d", "--dimensions"), type="numeric",
                     help="The number of dimensions to use for graph building")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files denoised data and combined SCE object")

parser <- add_option(parser, c("-k", "--knn"), type="numeric",
       	             help="k value to use for graph building")

parser <- add_option(parser, c("-r", "--reduced"), type="character",
       	  	     help="Name of reducedDim slot in SCE object to use")

parser <- add_option(parser, c("-c", "--cluster"), type="character",
                     default="walktrap", help="igraph clustering algorithm to use")

opt <- parse_args(parser)

message(paste0("Reading in ", opt$SCE))
big.sce <- readRDS(opt$SCE)

message("Building ", opt$knn, "-NN graph using ", opt$reduced)
nn.graph <- buildKNNGraph(reducedDim(big.sce, opt$reduced)[, c(1:opt$dimensions)], transposed=TRUE, k=opt$knn)

message("Clustering cells using ", opt$cluster)
if(opt$cluster %in% c("walktrap")){
  nn.cluster <- cluster_walktrap(nn.graph, membership=TRUE)
} else if(opt$cluster %in% c("louvain")){
  nn.cluster <- cluster_louvain(nn.graph)
} else {
  stop(opt$cluster, " algorithm not recognised")
}

nn.clusts <- membership(nn.cluster)
message("Found ", length(unique(nn.clusts)), " clusters")
clust.df <- data.frame("CellID"=colnames(big.sce), "Cluster"=as.character(nn.clusts))

clust.ofile <- paste0(opt$output, "_clusters.tsv")
message("Saving clustering to ", clust.ofile)
write.table(clust.df, file=clust.ofile, sep="\t", quote=FALSE, row.names=FALSE)

message("All done")

