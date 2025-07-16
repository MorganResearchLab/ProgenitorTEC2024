rm(list = ls())

library(Signac)
library(Seurat)
library(S4Vectors)
library(ggplot2)
library(GenomicRanges)
library(GenomeInfoDb)
library(harmony)
library(SeuratWrappers)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnsDb.Mmusculus.v79)
library(Biostrings)
library(rtracklayer)
library(SeuratData)
library(SeuratObject)
library(dplyr)
library(FNN)
library(patchwork)
library(tidyr)
library(future)
plan(sequential)
options(future.globals.maxSize = 100 * 1024^3)
set.seed(1234)

data <- readRDS("/Volumes/T7/Xenium_1/20241016__175903__3FD076_HOLLANDER/QCed/output-XETG00160__0037486__YH2__20241016__180103_seurat_object_qced.rds")
metadata <- readRDS("/Volumes/T7/Ly6d_Project/Final_obj_annotation_to_use_for_export.rds")
data <- AddMetaData(data, meta = metadata)
addmetadata <- read.csv("/Volumes/T7/Ly6d_Project/refined_nh_ids.csv")
rownames(addmetadata) <- addmetadata$barcodes
addmetadata <- addmetadata %>% dplyr::select(-barcodes)
data <- AddMetaData(data, meta = addmetadata)

data$x = as.numeric(data@images$fov@boundaries$centroids@coords[, 1])
data$y = as.numeric(data@images$fov@boundaries$centroids@coords[, 2])

df_full <- data.frame(
  row.names = colnames(data),
  cellbarcode = colnames(data),
  celltype = data$final_celltypes_xenium,
  x = as.numeric(data@images$fov@boundaries$centroids@coords[, 1]),
  y = as.numeric(data@images$fov@boundaries$centroids@coords[, 2])
)

cortex_celltypes <- c(
  "Xenium_DP (Q1)", "Xenium_DP (Q2)", "Xenium_DP (Sig.)", "Xenium_DP (P)",
  "Xenium_DN", "Xenium_cTEC", "Xenium_tc_stromal_mix(cor)",
  "Xenium_CapFb", "Xenium_Adipo"
)

medulla_celltypes <- c(
  "Xenium_mTEC", "Xenium_TlTEC", "Xenium_MedFb", "Xenium_Mature CD4",
  "Xenium_tc_stromal_mix(med)", "Xenium_SP4_resting / Treg-like"
)

# library(DropletUtils)
# DropletUtils::write10xCounts(path = "/Volumes/T7/Ly6d_Project/raw_data_full_exp/", x = as(data@assays$RNA@counts, "dgCMatrix"))
# metainfo <- data@meta.data
# rownames(metainfo) <- colnames(data)
# valid_celltypes <- c(cortex_celltypes, medulla_celltypes)
# metainfo$final_celltypes_xenium <- as.character(metainfo$final_celltypes_xenium)
# metainfo$final_celltypes_xenium[!(metainfo$final_celltypes_xenium %in% valid_celltypes)] <- NA
# metainfo <- metainfo[, c("nCount_RNA", "nFeature_RNA", "final_celltypes_xenium", "x", "y")]
# write.csv(metainfo, file = "/Volumes/T7/Ly6d_Project/raw_data_full_exp/meta.csv", row.names = TRUE)

data$final_celltypes_xenium <- as.character(data$final_celltypes_xenium)
valid_celltypes <- c(cortex_celltypes, medulla_celltypes)
Idents(data) <- "final_celltypes_xenium"
data_rep <- subset(data, idents=valid_celltypes)
Idents(data_rep) <- "final_celltypes_xenium"

grep("Prss", rownames(data), value = TRUE)

features_sel <- c("Ptprc", "Cd3e", "Cd8a", "Cd4", "Cd5", "Cd69", "Il2ra", "Foxp3", "Col6a1", "Col14a1", "Adipoq", "Epcam", "Cd80", "Cd40", "Aire", "Trpm5", "Psmb11")

dot_plot = DotPlot(data_rep, features = features_sel, dot.scale = 4)
dot_plot <- dot_plot +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dot_plot

ggsave(filename = "/Volumes/T7/Ly6d_Project/Feature_Plot_Xen_Celltypes.pdf", plot = dot_plot, dpi = 400, width = 9, height = 6)

df_full$compartment <- ifelse(
  df_full$celltype %in% cortex_celltypes, "cortex",
  ifelse(df_full$celltype %in% medulla_celltypes, "medulla", NA)
)

confident <- df_full %>% filter(celltype %in% c(cortex_celltypes, medulla_celltypes))
ambiguous <- df_full %>% filter(!celltype %in% c(cortex_celltypes, medulla_celltypes))

conf_coords <- as.matrix(confident[, c("x", "y")])
ambig_coords <- as.matrix(ambiguous[, c("x", "y")])

nn_11 <- get.knnx(data = conf_coords, query = ambig_coords, k = 11)
nn_labels_11 <- matrix(confident$compartment[nn_11$nn.index], ncol = 11)

majority_vote_11 <- apply(nn_labels_11, 1, function(labels) {
  tab <- table(labels)
  winner <- names(tab)[which.max(tab)]
  if (tab[winner] >= 6) return(winner)
  return(NA)
})

ambig_idx <- match(ambiguous$cellbarcode, df_full$cellbarcode)
df_full$compartment[ambig_idx] <- majority_vote_11

coords_all <- as.matrix(df_full[, c("x", "y")])
nn_501 <- get.knn(coords_all, k = 501)
nn_idx_501 <- nn_501$nn.index

smoothed_labels <- sapply(1:nrow(df_full), function(i) {
  neighbor_labels <- df_full$compartment[nn_idx_501[i, ]]
  neighbor_labels <- neighbor_labels[!is.na(neighbor_labels)]
  if (length(neighbor_labels) == 0) return(sample(c("cortex", "medulla"), 1))
  tab <- table(neighbor_labels)
  names(tab)[which.max(tab)]
})

df_full$compartment_final <- smoothed_labels

df_full$compartment_final <- as.character(df_full$compartment_final)

df_plot <- df_full %>% slice_sample(n = 50000)

p1 <- ggplot(df_plot, aes(x = x, y = y)) +
  geom_point(aes(color = compartment_final), size = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = c("cortex" = "lightblue", "medulla" = "plum"),
    na.value = "grey80",
    guide = guide_legend(override.aes = list(size = 4))  # increase legend marker size
  ) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "Compartment Definition", color = "Compartment")

celltypes_unique <- unique(df_plot$celltype)
celltype_colors <- setNames(
  alpha(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(celltypes_unique)), 1),  # 1 = no transparency
  celltypes_unique
)

p2 <- ggplot(df_plot, aes(x = x, y = y)) +
  geom_point(aes(color = celltype), size = 0.5, alpha = 0.8) +
  scale_color_manual(values = celltype_colors, na.value = "grey80") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Cell Type Annotation")

p1 <- p1 #/ p2
p1

ggsave(filename = "/Volumes/T7/Ly6d_Project/compartments_wo_cmj_region.pdf", plot = p1, dpi = 400, width = 12, height = 7)

df_full$compartment_cmj <- df_full$compartment_final

medulla_cells <- df_full %>% filter(compartment_final == "medulla")
cortex_cells  <- df_full %>% filter(compartment_final == "cortex")

nn_dist <- get.knnx(
  data = medulla_cells[, c("x", "y")],
  query = cortex_cells[, c("x", "y")],
  k = 1
)$nn.dist

# 50um into the cortex
cmj_indices <- which(nn_dist <= 50)

df_full$compartment_cmj[rownames(df_full) %in% rownames(cortex_cells)[cmj_indices]] <- "cmj"

# subsample for representation
df_plot <- df_full %>% slice_sample(n = 50000)


# df_full$compartment_cmj <- df_full$compartment_final  # make editable copy
# 
# medulla_cells <- df_full %>% filter(compartment_final == "medulla")
# non_medulla <- df_full %>% filter(compartment_final != "medulla")
# 
# nn1 <- get.knnx(
#   data = medulla_cells[, c("x", "y")],
#   query = non_medulla[, c("x", "y")],
#   k = 1
# )
# 
# within_50_1 <- which(nn1$nn.dist <= 50)
# medulla_edge_barcodes <- rownames(non_medulla)[within_50_1]
# 
# medulla_edge_cells <- df_full[medulla_edge_barcodes, ]
# cortex_cells <- df_full %>% filter(compartment_final == "cortex")
# 
# nn2 <- get.knnx(
#   data = medulla_edge_cells[, c("x", "y")],
#   query = cortex_cells[, c("x", "y")],
#   k = 1
# )
# 
# within_50_2 <- which(nn2$nn.dist <= 50)
# cortex_edge_barcodes <- rownames(cortex_cells)[within_50_2]
# 
# cmj_all_barcodes <- unique(c(medulla_edge_barcodes, cortex_edge_barcodes))
# df_full$compartment_cmj[rownames(df_full) %in% cmj_all_barcodes] <- "cmj"
# 
# df_plot <- df_full %>% slice_sample(n = 50000)


p1 <- ggplot(df_plot, aes(x = x, y = y)) +
  geom_point(aes(color = compartment_cmj), size = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = c("cortex" = "lightblue", "medulla" = "plum", "cmj" = "moccasin"),
    na.value = "grey80",
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "Compartment Assignment with CMJ", color = "Compartment")

celltypes_unique <- unique(df_plot$celltype)
celltype_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(celltypes_unique)),
  celltypes_unique
)

p2 <- ggplot(df_plot, aes(x = x, y = y)) +
  geom_point(aes(color = celltype), size = 0.5, alpha = 0.9) +
  scale_color_manual(values = celltype_colors, na.value = "grey80") +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Cell Type Annotation")

p1 <- p1 #/ p2
p1

ggsave(filename = "/Volumes/T7/Ly6d_Project/compartments_w_cmj_region.pdf", plot = p1, dpi = 400, width = 12, height = 7)

df_clean <- df_full %>% filter(!is.na(compartment_cmj))

compartment_counts <- df_clean %>%
  group_by(compartment_cmj) %>%
  summarise(count = n()) %>%
  ungroup()

total_cells <- sum(compartment_counts$count)
compartment_counts <- compartment_counts %>%
  mutate(fraction = count / total_cells)

compartment_fractions <- setNames(compartment_counts$fraction, compartment_counts$compartment_cmj)

compartment_fractions <- compartment_fractions[c("cortex", "cmj", "medulla")]

print(compartment_fractions)

addmetadata <- df_full$compartment_cmj

data <- AddMetaData(data, metadata = addmetadata, col.name = "compartments")

DefaultAssay(data) <- "RNAPreSCT"

grep("^Krt", rownames(data), value = TRUE)

# Only cells having transcripts for both Ccl21a and Lifr
expr <- FetchData(data, vars = c("Ccl21a", "Lifr"))

keep_cells <- rowSums(expr > 0) == ncol(expr)

data <- subset(data, cells = rownames(expr)[keep_cells])

data[["RNAPreSCT"]] <- data[["RNA"]]

data <- SCTransform(data, assay="RNAPreSCT", new.assay.name = "RNA", conserve.memory=TRUE, verbose=TRUE)
data <- FindVariableFeatures(data, nfeatures = 3000)
data <- RunPCA(data, assay = "RNA")

degs_up <- read.csv("/Volumes/T7/Ly6d_Project/Ly6dPositive_TEPC_Allmarkers.tsv", sep = " ")
degs_up <- degs_up$x
degs_down <- read.csv("/Volumes/T7/Ly6d_Project/Ly6dNegative_TEPC_Allmarkers.tsv", sep = " ")
degs_down <- degs_down$x

#degs <- as.data.frame(read.csv("/Volumes/T7/Ly6d_Project/iTEC_Ly6d_DEGs.csv"))

# top_up_genes <- degs[degs$FDR < 0.05 & degs$logFC > 0.5, ]
# top_up_genes <- top_up_genes[order(-top_up_genes$logFC), ]
# top_up_gene_names <- as.character(head(top_up_genes$external_gene_name, 200))
top_up_gene_names <- intersect(rownames(data), degs_up)

write.csv(top_up_gene_names, "/Volumes/T7/Ly6d_Project/Ly6d_Pos_Intersect.csv")

# top_down_genes <- degs[degs$FDR < 0.05 & degs$logFC < -0.5, ]
# top_down_genes <- top_down_genes[order(top_down_genes$logFC), ]
# top_down_gene_names <- as.character(head(top_down_genes$external_gene_name, 200))
top_down_gene_names <- intersect(rownames(data), degs_down)

write.csv(top_down_gene_names, "/Volumes/T7/Ly6d_Project/Ly6d_Neg_Intersect.csv")

data <- AddModuleScore(data, features = list(top_down_gene_names), name = "Ly6d_neg_It_TEC", nbin = 18)

data <- AddModuleScore(data, features = list(top_up_gene_names), name = "Ly6d_pos_It_TEC", nbin = 18)

scores_df <- data@meta.data[, c("Ly6d_neg_It_TEC1", "Ly6d_pos_It_TEC1")]

high_score <- data$Ly6d_pos_It_TEC1  
low_score <- data$Ly6d_neg_It_TEC1

############# Pick Cells ################

score_diff <- high_score - low_score

abs_max <- max(abs(score_diff))

scaled_score <- score_diff / abs_max

data$Ly6d_bias_score <- scaled_score

# # Number based assignment
# ranked_cells_tp3 <- rank(-high_score, ties.method = "min")
# ranked_cells_tp1 <- rank(-low_score, ties.method = "min")
# top_300_tp1 <- which(ranked_cells_tp1 <= 2000)
# top_300_tp3 <- which(ranked_cells_tp3 <= 2000)
# data$Cell_Type_Score <- NA
# data$Cell_Type_Score[top_300_tp1] <- "It_TEC_Ly6d_neg"
# data$Cell_Type_Score[top_300_tp3] <- "It_TEC_Ly6d_pos"

# Score based assignment
data$Cell_Type_Score <- NA  # or rename as needed
data$Cell_Type_Score[data$Ly6d_bias_score >= 0.2] <- "It_TEC_Ly6d_pos"
data$Cell_Type_Score[data$Ly6d_bias_score <= -0.2] <- "It_TEC_Ly6d_neg"

library(DropletUtils)
DropletUtils::write10xCounts(path = "/Volumes/T7/Ly6d_Project/raw_data_Ccl21_and_Lifr_pos_ItTEC/", x = as(data@assays$RNA@counts, "dgCMatrix"))
metainfo <- data@meta.data
rownames(metainfo) <- colnames(data)
valid_celltypes <- c(cortex_celltypes, medulla_celltypes)
metainfo$final_celltypes_xenium <- as.character(metainfo$final_celltypes_xenium)
metainfo$final_celltypes_xenium[!(metainfo$final_celltypes_xenium %in% valid_celltypes)] <- NA
metainfo <- metainfo[, c("nCount_RNA", "nFeature_RNA", "Ly6d_bias_score", "Cell_Type_Score", "x", "y")]
write.csv(metainfo, file = "/Volumes/T7/Ly6d_Project/raw_data_Ccl21_and_Lifr_pos_ItTEC/meta.csv", row.names = TRUE)

##########################################

p1 <- ImageDimPlot(data, group.by = "Cell_Type_Score", size=1.5, na.value = "gray30", flip_xy = FALSE, cols = c("blue", "red"), border.color = "black")
p1 = p1 + theme_void() +  # Remove axes/grid
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
    )
p1

ggsave(filename = "/Volumes/T7/Ly6d_Project/It_TEC_Ly6d_pos_neg_prediction_in_space.pdf" , plot = p1, dpi = 400, width = 15, height = 8)

# Subset Seurat object to cells with non-NA Cell_Type_Score
cells_to_plot <- colnames(data)[!is.na(data$Cell_Type_Score)]
data_subset <- subset(data, cells = cells_to_plot)

# Now plot only those
p1 <- ImageDimPlot(
  data_subset,
  group.by = "Cell_Type_Score",
  size = 1.5,
  cols = c("blue", "red"),
  border.color = "black",
  flip_xy = FALSE
)

# Beautify the plot
p1 <- p1 + theme_void() +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

p1

table(data_subset$Cell_Type_Score)

# Save
ggsave(
  filename = "/Volumes/T7/Ly6d_Project/It_TEC_Ly6d_pos_neg_prediction_in_space_no_nas.pdf",
  plot = p1,
  dpi = 400,
  width = 15,
  height = 8
)

print(table(data$Cell_Type_Score))

p1 <- ImageFeaturePlot(
  data,
  features = "Ly6d_bias_score",
  size = 1.5,
  cols = c("blue", "red"),
  min.cutoff = -1,
  max.cutoff = 1,
  border.color = "black",
  
) + 
  theme_void() +  # Remove axes/grid
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )
p1

#ggsave(filename = "/Volumes/T7/Ly6d_Project/It_TEC_Ly6d_pos_neg_continous_prediction_in_space.pdf" , plot = p1, dpi = 400, width = 8, height = 15)

data$Ly6d_bias_score <- as.numeric(as.character(data$Ly6d_bias_score))

p1 <- ImageFeaturePlot(
  data,
  features = "Ly6d_bias_score",
  size = 2,
  cols = c("blue", "gray30", "red"),
  border.color = "black"
) + 
  theme_void() +  # Remove axes/grid
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )
p1

p1 <- ImageFeaturePlot(data, features = c(top_down_gene_names[18:27]))
p1

p1 <- ImageFeaturePlot(data, features = c(top_up_gene_names[1:18]))
p1

p1 <- ImageFeaturePlot(data, features = c("Scn1a", "Ackr4", "Spon1", "Trem2", "Itgb5", "Runx2"), cols = c("darkgray", "firebrick"), size = 0.7)+ 
  theme_void() +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )
p1

features <- c("Scn1a", "Ackr4", "Spon1", "Trem2", "Itgb5", "Runx2")
plots <- ImageFeaturePlot(data, features = features, cols = c("darkgray", "firebrick"), size = 0.7)

plots <- lapply(plots, function(p) {
  p + 
    theme_void() +
    theme(
      plot.background = element_rect(fill = "black", color = NA),
      panel.background = element_rect(fill = "black", color = NA),
      legend.background = element_rect(fill = "black", color = NA),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white")
    )
})

p1 <- wrap_plots(plots)
p1

ggsave(filename = "/Volumes/T7/Ly6d_Project/It_TEC_Ly6d_pos_neg_6_features_in_space.pdf" , plot = p1, dpi = 400, width = 12, height = 17)

table(data$Cell_Type_Score)

df <- data@meta.data[, c("Cell_Type_Score", "compartments")]
df <- df %>% filter(!is.na(Cell_Type_Score), !is.na(compartments))
df$compartment_cmj <- factor(df$compartments, levels = c("cortex", "cmj", "medulla"))

cell_counts <- df %>%
  group_by(compartment_cmj, Cell_Type_Score) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Cell_Type_Score) %>%
  mutate(fraction_within_compartment = count / sum(count)) %>%
  ungroup()

cell_density <- cell_counts %>%
  mutate(compartment_fraction = compartment_fractions[as.character(compartment_cmj)],
         density = count * compartment_fractions)

density_fraction <- cell_counts %>%
  group_by(Cell_Type_Score) %>%
  mutate(prop_density = fraction_within_compartment / compartment_fractions)

p1 <- ggplot(cell_counts, aes(x = compartment_cmj, y = count, fill = Cell_Type_Score)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(
    "It_TEC_Ly6d_neg" = "blue", 
    "It_TEC_Ly6d_pos" = "red" 
  )) +
  labs(
    x = "Region",
    y = "Cell Count",
    fill = "Cell Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(filename = "/Volumes/T7/Ly6d_Project/barchart_cell_type_counts_in_region.pdf", plot = p1, dpi = 400, width = 10, height = 7)


p1 <- ggplot(cell_counts, aes(x = compartment_cmj, y = fraction_within_compartment, fill = Cell_Type_Score)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(
    "It_TEC_Ly6d_neg" = "blue", 
    "It_TEC_Ly6d_pos" = "red" 
  )) +
  labs(
    x = "Region",
    y = "Cell Category Fraction",
    fill = "Cell Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(filename = "/Volumes/T7/Ly6d_Project/barchart_cell_type_fraction_in_region.pdf", plot = p1, dpi = 400, width = 10, height = 7)


# p1 <- ggplot(density_fraction, aes(x = compartment_cmj, y = prop_density, fill = Cell_Type_Score)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   scale_fill_manual(values = c(
#     "It_TEC_Ly6d_neg" = "blue", 
#     "It_TEC_Ly6d_pos" = "red" 
#   )) +
#   labs(
#     x = "Region",
#     y = "Cell Category Fraction per region area",
#     fill = "Cell Category"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p1
# 
# ggsave(filename = "/Volumes/T7/Ly6d_Project/barchart_region_fractions_in_comp_per_comp_area.pdf", plot = p1, dpi = 400, width = 12, height = 7)

p1 <- ggplot(cell_density, aes(x = compartment_cmj, y = density, fill = Cell_Type_Score)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c(
    "It_TEC_Ly6d_neg" = "blue",
    "It_TEC_Ly6d_pos" = "red"
  )) +
  labs(
    x = "Region",
    y = "Cell Density in Region",
    fill = "Cell Category"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(filename = "/Volumes/T7/Ly6d_Project/barchart_cell_type_counts_multiplied_with_region_fraction.pdf", plot = p1, dpi = 400, width = 10, height = 7)

# p1 <- ggplot(density_fraction, aes(x = compartment_cmj, y = prop_density, fill = Cell_Type_Score)) +
#   geom_bar(stat = "identity", position = "stack") +
#   scale_fill_manual(values = c(
#     "It_TEC_Ly6d_neg" = "blue",
#     "It_TEC_Ly6d_pos" = "red"
#   )) +
#   labs(
#     title = "Density-Based Ratio",
#     x = "Region",
#     y = "Density Fraction of It_TECs",
#     fill = "It TEC Category"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# p1
# 
# ggsave("/Volumes/T7/Ly6d_Project/stacked_ratio_density.pdf", p1, dpi = 400, width = 7, height = 8)
