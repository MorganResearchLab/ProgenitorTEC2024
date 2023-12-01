#! /usr/bin/env Rscript

## Compute the mean adjusted noise (CV2) over single cells in each nhood
## Use these to perform differential variability testing between conditions

### Build a Milo object ready for DA testing
library(SingleCellExperiment)
library(miloR)
library(Matrix)
library(MatrixGenerics)
library(optparse)


# within an nhood group, compute the mean-adjusted gene expression noise across nhoods - how does one account for the sharing of cells here?
meanAdjustedNoise <- function(gene.matrix, nhood.matrix, nhood.groups, design.df){
    # compute the noise (CV^2) in each group
    # regress the mean expression across groups on the CV^2
    # group.assign is a named vector
    unique.nhoods <- colnames(nhood.matrix) # these are the nhoods
    unique.groups <- unique(nhood.groups$NhoodGroup) # these the nhood groups, i.e. clusters
    
    gene.noise <- list()
    unique.genes <- rownames(gene.matrix)
    for(q in seq_along(unique.groups)){
        q.group <- unique.groups[q]
	message("Extracting nhoods from nhood group ", q.group)
        q.nhoods <- nhood.matrix[, nhood.groups[nhood.groups$NhoodGroup %in% q.group, ]$Nhood, drop=FALSE]
        
        q.mean.list <- list()
        q.cv2.list <- list()
	message("Nhood group, ", q.group, " contains ", ncol(q.nhoods), " nhoods")
        
        for(l in seq_len(ncol(q.nhoods))){
            l.cells <- rownames(q.nhoods)[rowSums(q.nhoods[, l, drop=FALSE]) > 0]
            B <- model.matrix( ~ 0 + HTO, design.df[design.df$CellID %in% l.cells, ])
            rownames(B) <- design.df$CellID[design.df$CellID %in% l.cells]
            QB <- gene.matrix[, l.cells, drop=FALSE] %*% B[l.cells, , drop=FALSE]
            q.agg <- apply(QB, 1, FUN=function(QX) QX/colSums(B)) # mean matrix
            
            q.var <- t(sapply(seq_len(ncol(B)), FUN=function(BX) {
                BX.cells <- rownames(B)[sum(B[, BX]) > 0]
                rowVars(gene.matrix[, BX.cells, drop=FALSE])
            }))
            
            colnames(q.var) <- rownames(gene.matrix)
            rownames(q.var) <- paste(colnames(q.nhoods)[l], colnames(B), sep=".")
            
            colnames(q.agg) <- rownames(gene.matrix)
            rownames(q.agg) <- paste(colnames(q.nhoods)[l], colnames(B), sep=".")
            
            q.cv2 <- q.var/q.agg
            q.mean.list[[l]] <- t(q.agg) # gene X sample.nhood
            q.cv2.list[[l]] <- t(q.cv2)
        }

	message("Concatenating mean and CV2 matrices over nhoods")
        # concatenate over nhoodsXSample
        q.cv2.mat <- do.call(cbind, q.cv2.list)
        q.mean.mat <- do.call(cbind, q.mean.list)

	message("Removing non-expressed and 0-variance genes")
        is.exprs <- apply(!is.na(q.mean.mat), 1, all) # ||  rowMeans(q.mean.mat) > 0
        is.var <- apply(!is.na(q.cv2.mat), 1, all) # || rowMeans(q.cv2.mat) > 0

	# drop the infinites, or set to zero?
        for(i in seq_len(ncol(q.cv2.mat))){            
            q.cv2.mat[is.infinite(q.cv2.mat[, i]), i] <- 0
        }
        
        keep.genes <- unique.genes[is.exprs & is.var]
	message("Keeping ", length(keep.genes), " genes")
        
	message("Fitting a local polynomial regression for each gene")
	# fit a local polynomial regression for each gene - gives a vector of adjusted noise values for each nhood
        q.etaz <- sapply(seq_along(keep.genes), FUN=function(GX){
            gx.gene <- keep.genes[GX]
            gx.df <- data.frame("CV2"=q.cv2.mat[gx.gene, ], "Mean"=q.mean.mat[gx.gene, ])
            
            gx.fit <- loess(CV2 ~ Mean, data=gx.df, span=0.75)
            gx.eta <- residuals(gx.fit)
            gx.z <- (gx.eta - mean(gx.eta))/var(gx.eta)
            return(gx.z)
        }) # gene X sample.nhood
        
        rownames(q.etaz) <- paste(q.group, rownames(q.etaz), sep=".")
        colnames(q.etaz) <- keep.genes
	q.etaz <- as.data.frame(t(q.etaz))
	q.etaz$Gene <- keep.genes
        gene.noise[[q.group]] <- q.etaz

	message("Mean-adjusted noise estimated for ", ncol(q.etaz), " genes and ", nrow(q.etaz), " nhood X sample combinations")
        # break
        sink(file="/dev/null")
        rm(list=c("q.etaz", "q.mean.mat", "q.cv2.mat", "q.mean.list", "q.cv2.list"))
        gc()
        sink(file=NULL)
    }
    
    noise.df <- Reduce(x=gene.noise, f=function(x,y) merge(x, y, by="Gene", all=TRUE))
    sink(file="/dev/null")
    rm(list=c("gene.noise"))
    gc()
    sink(file=NULL)
    
    return(noise.df)
}

parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the combined SCE object")

parser <- add_option(parser, c("-m", "--meta"), type="character",
       	             help="Meta-data file, tab separated")

parser <- add_option(parser, c("-n", "--nhoods"), type="character",
                     help="File containing mapping of nhoods to nhood groups")

parser <- add_option(parser, c("-g", "--genes"), type="character", default=NULL,
       	  	     help="Comma-separated list of genes to subset noise values")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Output file for noise values")

opt <- parse_args(parser)

message("Reading SCE object ", opt$SCE)
itec.milo <- readRDS(opt$SCE)

message("Reading meta data file ", opt$meta)
meta.proj.merge <- read.table(opt$meta, sep="\t", header=TRUE, stringsAsFactors=FALSE)

message("Reading in nhood group mapping file ", opt$nhoods)
nhood.group.df <- read.table(opt$nhoods, sep="\t", header=TRUE, stringsAsFactors=FALSE)

message("Estimating noise over all genes, nhoods and nhood groups")
noise.data <- meanAdjustedNoise(logcounts(itec.milo), nhoods(itec.milo), nhood.group.df, meta.proj.merge)

if(!is.null(opt$genes)){
    genes <- unlist(strsplit(opt$genes, split=",", fixed=TRUE))
    genes <- intersect(genes, noise.data$Gene)
    message("Found ", length(genes), " genes to subset noise data")

    print(sum(noise.data$Gene %in% genes))
    print(intersect(noise.data$Gene, genes))

    sub.file <- gsub(opt$output, pattern="\\.txt", replacement="SubsetGenes.txt")
    message("Writing subset noise data to ", sub.file)
    write.table(noise.data[noise.data$Gene %in% genes, ],
                file=sub.file, sep="\t", row.names=FALSE, quote=FALSE)
}


# this is a MASSIVE file - I'll also write a subset of genes
message("Writing noise data to ", opt$output)
write.table(noise.data, file=opt$output, sep="\t", row.names=FALSE, quote=FALSE)

message("All done")