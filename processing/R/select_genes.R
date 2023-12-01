#! /usr/bin/env Rscript

library(SingleCellExperiment)
library(miloR)
library(ggplot2)
library(optparse)
library(biomaRt)
parser <- OptionParser()
parser <- add_option(parser, c("-s", "--SCE"), type="character",
       	             help="The path to the combined SCE object")

parser <- add_option(parser, c("-o", "--output"), type="character",
                     help="Prefix for output files denoised data and combined SCE object")

parser <- add_option(parser, c("-a", "--assay"), type="character",
                     help="Expression to use in the SCE assays slot")

parser <- add_option(parser, c("-r", "--remove"), action="store_true", default=FALSE,
                     help="Switch to remove all non-Singlet cell barcodes")
                  

opt <- parse_args(parser)

message(paste0("Reading in ", opt$SCE))
all.sce <- readRDS(opt$SCE)
print(head(rownames(all.sce)))
if(is.null(rownames(all.sce))){
    all.sce <- all.sce[!is.na(rowData(all.sce)$ensembl_gene_id), ]
    rownames(all.sce) <- rowData(all.sce)$ensembl_gene_id
}

#print(head(rowData(all.sce)))
print(head(rownames(all.sce)))
rowData(all.sce) <- NULL

if(opt$remove){
  all.sce <- all.sce[, all.sce$Doublet.Class %in% c("Singlet")]
  message("Keeping ", ncol(all.sce), " singlets")
}

message("Extracting marker genes of interest")
goi.genes <- c("Cd3e", "Cd3g", "Cd14", "Epcam", "Prss16", "Ly6a", "Aire", "Cd52",
	       "Psmb11", "Zap70", "Cdk1", "Top2a", "Ube2c", "Mki67", "Tyms", "Foxn1", "Ccl25", "Cxcl12",
	       "Dll4", "Cd83", "Prss16", "Gpr25", "Tspan10", "Tbata", "Ajap1", "Apba2", "Dtx4", "Pds5b",
	       "Stc2", "Zar1", "Snap91", "Srxn1", "Ctnna2", "Lhx3", "Cecr6", "Ankmy1", "Tmie", "Abca15",
	       "Gmfg", "Arhgef10l", "Ptprn2", "Lmx1a", "Endou", "Mfsd2", "Akap2", "Tasp1", "Psmb10", "Psmb9",
	       "Psmb4", "Psma4", "Mcam", "Myc", "Mycn", "Icam2", "Mdk", "Cd205", "Itga6", "Plet1", "Gas1",
	       "Tbx21", "Fezf2", "Dsg3", "Lifr", "H2-Ab1", "H2-Ob", "H2-Dma", "H2-Dmb1", "H2-K1", 
	       "H2-Ea", "H2-Eb1", "H2-Eb2", "H2-Aa", "H2-Oa", "Krt5", "Pdpn", "Ccl21a", "Sod3", "Dpt",
	       "Trpm5", "Avil", "Krt80", "Spink5", "Pba2", "Ccna2", "Cxcl12", "Syngr1", "Gper1", "Cd177",
	       "Car8", "Tnfrsf11a", "Tnfsf11", "Cd40", "Cd40lg", "Fgf7", "Fst", "Fgfr1", "Fgfr2", "Bmpr1a",
	       "Bmp4a", "Inhba", "Fgfr3", "Bmp2")

# need to convert the gene symbols back to ensembl IDs
biomaRt.connection <- useMart("ensembl", "mmusculus_gene_ensembl", host="uswest.ensembl.org")

gene.df <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                 mart=biomaRt.connection)
rownames(gene.df) <- gene.df$ensembl_gene_id
#gene.df <- gene.df[!duplicated(gene.df$external_gene_name), ]

goi.ensembl.id <- gene.df$ensembl_gene_id[gene.df$external_gene_name %in% unique(goi.genes)]
if(any(rownames(all.sce) %in% c("ENSGZsGreen"))){
  goi.gene.id <- intersect(rownames(all.sce), c(goi.ensembl.id, "ENSGZsGreen"))
} else{
  goi.gene.id <- intersect(rownames(all.sce), goi.ensembl.id)
}

print(length(goi.gene.id))
print(opt$assay)
print(dim(all.sce[goi.gene.id, ]))

print(assayNames(all.sce))

goi.gene.exprs <-  as.data.frame(as.matrix(t(assay(all.sce[goi.gene.id, ], opt$assay))))

if(any(rownames(all.sce) %in% "ENSGZsGreen")){
  message("ENSGZsGreen")
  colnames(goi.gene.exprs) <- c(gene.df[setdiff(colnames(goi.gene.exprs), "ENSGZsGreen"), ]$external_gene_name, "ENSGZsGreen")
} else{
  colnames(goi.gene.exprs) <- gene.df[colnames(goi.gene.exprs), ]$external_gene_name
}
print(head(goi.gene.exprs))
goi.gene.exprs$CellID <- colnames(all.sce)
rownames(goi.gene.exprs) <- goi.gene.exprs$CellID

write.table(goi.gene.exprs,
            file=paste0(opt$output, "_GeneGOI.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

message("Calculating mean expression in HTOs")
## calculate average expression of genes in each HTO

tab.form <- as.formula("~ 0 + HTO")
clust.model.matrix <- model.matrix(tab.form, data=colData(all.sce))
rownames(clust.model.matrix) <- colnames(all.sce)

meta.df <- as.data.frame(colData(all.sce))
meta.df$CellID <- colnames(all.sce)
goi.merge <- merge(goi.gene.exprs, meta.df, by='CellID')
rownames(goi.merge) <- goi.merge$CellID

common.cells <- intersect(rownames(clust.model.matrix), rownames(goi.merge))
gene.cols <- setdiff(colnames(goi.merge), colnames(meta.df))

goi.by.cluster <- t(clust.model.matrix) %*% as.matrix(goi.merge[, gene.cols])
goi.by.cluster.bin <- t(clust.model.matrix) %*% as.matrix(goi.merge[, gene.cols] > 0)
ave.by.label <- as.data.frame(apply(goi.by.cluster, 2, function(X) X/colSums(clust.model.matrix)))
print(head(ave.by.label))

rownames(ave.by.label) <- gsub(rownames(ave.by.label), pattern="HTO", replacement="")

ave.by.label$GeneID <- rownames(ave.by.label)
write.table(ave.by.label,
            file=paste0(opt$output, "-HTO_gene-Ave.tsv"),
	    quote=FALSE, row.names=FALSE, sep="\t")

### extract proteins too
# switch experiment
cite.sce <- altExp(all.sce, "CITE")
protein.goi.exprs <- as.data.frame(as.matrix(t(assay(cite.sce, "cpm"))))
protein.goi.exprs$CellID <- colnames(cite.sce)

write.table(protein.goi.exprs, file=paste0(opt$output, "_ProteinGOI.tsv"),
	    sep="\t", quote=FALSE, row.names=FALSE)