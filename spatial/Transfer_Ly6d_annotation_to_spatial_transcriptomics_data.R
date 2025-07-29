rm(list = ls())

library(Signac)
library(Seurat)
library(S4Vectors)
library(ggplot2)
library(GenomicRanges)
library(GenomeInfoDb)
library(harmony)
library(Matrix)
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

# Generate rds object from raw files ("data_full_exp") -> run once

# input_folder <- "E:/Ly6d_Project/codes/data/"
# output_folder <- "E:/Ly6d_Project/codes/testfolder/"

input_folder <- "/Volumes/T7/Ly6d_Project/codes/data/"
output_folder <- "/Volumes/T7/Ly6d_Project/codes/testfolder/"

########################
expression_matrix <- readMM(paste0(input_folder, "matrix.mtx"))
genes <- read.delim(paste0(input_folder, "genes.tsv"), header = FALSE)
genes <- genes$V1
barcodes <- read.delim(paste0(input_folder, "barcodes.tsv"), header = FALSE)
barcodes <- barcodes$V1
meta <- read.csv(paste0(input_folder, "meta.csv"))
rownames(expression_matrix) <- genes
colnames(expression_matrix) <- barcodes
expression_matrix <- CreateAssayObject(expression_matrix)
data <- CreateSeuratObject(counts = expression_matrix, assay = "RNAPreSCT")
data$nCount_RNA <- data$nCount_RNAPreSCT
data$nFeature_RNA <- data$nFeature_RNAPreSCT
data$nCount_RNAPreSCT <- NULL
data$nFeature_RNAPreSCT <- NULL
data <- AddMetaData(data, metadata = meta)
meta_compartment <- read.csv(paste0(input_folder, "refined_nh_ids.csv"), row.names = 1)
data <- AddMetaData(data, metadata = meta_compartment)
data$barcodes <- data$X
data@meta.data$X <- NULL
centroids_df <- data.frame(
  x = data$x,
  y = data$y,
  cell = colnames(data)
)
centroids <- CreateCentroids(centroids_df)
fov_img <- CreateFOV(
  coords = list(centroids = centroids),
  type = "centroids",
  assay = "RNAPreSCT",       
  key = "RNAPreSCT_",       
  name = "fov"
)
data@images[["fov"]] <- fov_img
rm(expression_matrix, barcodes, genes, meta)
saveRDS(data, paste0(input_folder, "rds_object.rds"))
#########################

data <- readRDS(paste0(input_folder, "rds_object.rds"))

df_full <- data.frame(
  row.names = colnames(data),
  cellbarcode = colnames(data),
  celltype = data$final_celltypes_xenium,
  x = data$x,
  y = data$y
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

data_rep <- data

Idents(data_rep) <- "final_celltypes_xenium"
data_rep <- subset(data_rep, subset = !is.na(final_celltypes_xenium))
Idents(data_rep) <- "final_celltypes_xenium"

data_rep <- SCTransform(data_rep, assay="RNAPreSCT", new.assay.name = "RNA", conserve.memory=TRUE, verbose=TRUE)

grep("Prss", rownames(data_rep), value = TRUE)

features_sel <- c("Ptprc", "Cd3e", "Cd8a", "Cd4", "Cd5", "Cd69", "Il2ra", "Foxp3", "Col6a1", "Col14a1", "Adipoq", "Epcam", "Cd80", "Cd40", "Aire", "Trpm5", "Psmb11")

dot_plot = DotPlot(data_rep, features = features_sel, dot.scale = 4)
dot_plot <- dot_plot +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dot_plot

ggsave(filename = paste0(output_folder, "Feature_Plot_Xen_Celltypes.pdf"), plot = dot_plot, dpi = 400, width = 9, height = 6)

rm(data_rep)

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
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "Compartment Definition", color = "Compartment")

celltypes_unique <- unique(df_plot$celltype)
celltype_colors <- setNames(
  alpha(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(celltypes_unique)), 1),
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

ggsave(filename = paste0(output_folder, "compartments_wo_cmj_region.pdf"), plot = p1, dpi = 400, width = 12, height = 7)

df_full$compartment_cmj <- df_full$compartment_final

medulla_cells <- df_full %>% filter(compartment_final == "medulla")
cortex_cells  <- df_full %>% filter(compartment_final == "cortex")

nn_dist <- get.knnx(
  data = cortex_cells[, c("x", "y")],
  query = medulla_cells[, c("x", "y")],
  k = 1
)$nn.dist

# 50um into the medulla

cmj_indices <- which(nn_dist <= 50)

df_full$compartment_cmj[rownames(df_full) %in% rownames(medulla_cells)[cmj_indices]] <- "cmj"

# subsample for representation
df_plot <- df_full %>% slice_sample(n = 50000)

p1 <- ggplot(df_plot, aes(x = x, y = y)) +
  geom_point(aes(color = compartment_cmj), size = 0.5, alpha = 0.8) +
  scale_color_manual(
    values = c("cortex" = "lightblue", "cmj" = "moccasin", "medulla" = "plum"),
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

ggsave(filename = paste0(output_folder, "compartments_w_cmj_region.pdf"), plot = p1, dpi = 400, width = 12, height = 7)

df_clean <- df_full %>% filter(!is.na(compartment_cmj))

compartment_counts <- df_clean %>%
  group_by(compartment_cmj) %>%
  summarise(count = n()) %>%
  ungroup()

total_cells <- sum(compartment_counts$count)
compartment_counts <- compartment_counts %>%
  mutate(fraction = count / total_cells)

compartment_fractions <- setNames(compartment_counts$fraction, compartment_counts$compartment_cmj)

compartments <- c("cortex", "cmj", "medulla")

compartment_fractions <- compartment_fractions[compartments]

print(compartment_fractions)

### Estimate area fractions ###

area_estimates <- data.frame(compartment = compartments, n_cells = NA, avg_distance = NA, area = NA)

# Assume hexagonal packing
hex_factor <- sqrt(3) / 2

for (i in seq_along(compartments)) {
  comp <- compartments[i]
  
  df_sub <- df_full %>% filter(compartment_cmj == comp)
  coords <- as.matrix(df_sub[, c("x", "y")])
  
  nn <- get.knn(coords, k = 2)
  d_avg <- mean(nn$nn.dist[, 2 - 1])  # nearest non-self

  n_cells <- nrow(coords)
  area_est <- n_cells * hex_factor * d_avg^2
  
  area_estimates[i, ] <- list(comp, n_cells, d_avg, area_est)
}

area_estimates <- area_estimates %>%
  mutate(area_fraction = area / sum(area, na.rm = TRUE))

area_fractions <- area_estimates$area_fraction
names(area_fractions) <- area_estimates$compartment
print(area_fractions)

###############################

addmetadata <- df_full$compartment_cmj

data <- AddMetaData(data, metadata = addmetadata, col.name = "compartments")

dim(data)

DefaultAssay(data) <- "RNAPreSCT"

# Only cells having transcripts for both Ccl21a and Lifr
expr <- FetchData(data, vars = c("Ccl21a", "Lifr"))

keep_cells <- rowSums(expr > 0) == ncol(expr)

data <- subset(data, cells = rownames(expr)[keep_cells])

dim(data)

DefaultAssay(data) <- "RNAPreSCT"

grep("Cd8", rownames(data), value = TRUE)

genes_to_filter <- c("Ptprc", "Cd8a", "Cd8b1", "Cd4", "Cd3d", "Cd3e", "Cd3g", "Rag1", "Cd5", "Cd69", "Cd68", "Sirpa", "Lck", "Tcf7", "Gzma", "Ighm", "Ighd", "Itgax")

expr <- GetAssayData(data, assay = "RNAPreSCT", slot = "counts")[genes_to_filter, ]

cells_to_remove <- colnames(expr)[Matrix::colSums(expr > 0) >= 4]

data <- subset(data, cells = setdiff(colnames(data), cells_to_remove))

dim(data)

data <- SCTransform(data, assay="RNAPreSCT", new.assay.name = "RNA", conserve.memory=TRUE, verbose=TRUE)
data <- FindVariableFeatures(data, nfeatures = 3000)
data <- RunPCA(data, assay = "RNA")

degs_up <- read.csv(paste0(input_folder, "Ly6dPositive_TEPC_Allmarkers.tsv"), sep = " ")
degs_up <- degs_up$x
degs_down <- read.csv(paste0(input_folder, "Ly6dNegative_TEPC_Allmarkers.tsv"), sep = " ")
degs_down <- degs_down$x

top_up_gene_names <- intersect(rownames(data), degs_up)

write.csv(top_up_gene_names, paste0(output_folder, "Ly6d_Pos_Intersect.csv"))

top_down_gene_names <- intersect(rownames(data), degs_down)

write.csv(top_down_gene_names, paste0(output_folder, "Ly6d_Neg_Intersect.csv"))

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

# Score based assignment
data$Cell_Type_Score <- NA

meta_df <- data@meta.data

ggplot(meta_df, aes(x = Ly6d_bias_score)) +
  geom_histogram(fill = "skyblue", color = "black", bins = 30) +
  labs(x = "Ly6d Bias Score", y = "Cell Count", title = "Distribution of Ly6d Bias Score") +
  theme_minimal()

data$Cell_Type_Score[data$Ly6d_bias_score >= 0.1] <- "It_TEC_Ly6d_pos"
data$Cell_Type_Score[data$Ly6d_bias_score <= -0.1] <- "It_TEC_Ly6d_neg"

table(data$Cell_Type_Score)

##########################################

p1 <- ImageDimPlot(data, group.by = "Cell_Type_Score", size=1.5, na.value = "gray30", flip_xy = FALSE, cols = c("blue", "red"), border.color = "black")
p1 = p1 + theme_void() +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
    )
p1

ggsave(filename = paste0(output_folder, "It_TEC_Ly6d_pos_neg_prediction_in_space.pdf") , plot = p1, dpi = 400, width = 15, height = 8)

cells_to_plot <- colnames(data)[!is.na(data$Cell_Type_Score)]
data_subset <- subset(data, cells = cells_to_plot)

p1 <- ImageDimPlot(
  data_subset,
  group.by = "Cell_Type_Score",
  size = 1.5,
  cols = c("blue", "red"),
  border.color = "black",
  flip_xy = FALSE
)

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

ggsave(filename = paste0(output_folder, "It_TEC_Ly6d_pos_neg_prediction_in_space_no_nas_dark.pdf"), plot = p1, dpi = 400, width = 15, height = 8)

p1 <- ImageDimPlot(
  data_subset,
  group.by = "Cell_Type_Score",
  size = 1.5,
  cols = c("blue", "red"),
  border.color = "black",
  flip_xy = FALSE
)

p1 <- p1 + theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )

p1

table(data_subset$Cell_Type_Score)

ggsave(filename = paste0(output_folder, "It_TEC_Ly6d_pos_neg_prediction_in_space_no_nas_bright.pdf"), plot = p1, dpi = 400, width = 15, height = 8)

p1 <- ImageDimPlot(
  data_subset,
  group.by = "Cell_Type_Score",
  size = 1.5,
  cols = c("blue", "red"),
  border.color = "black",
  flip_xy = TRUE
)

p1 <- p1 + theme_void() +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

p2 <- ImageDimPlot(
  data_subset,
  group.by = "compartments",
  size = 1.5,
  cols = c("cortex" = "lightblue", "cmj" = "moccasin", "medulla" = "plum"),
  border.color = "black",
  flip_xy = TRUE
)

p2 <- p2 + theme_void() +
  theme(
    plot.background = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    legend.background = element_rect(fill = "black", color = NA),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white")
  )

p1 <- p1 + p2
p1

ggsave(filename = paste0(output_folder, "It_TEC_Ly6d_pos_neg_prediction_and_compartment_in_space_no_nas_dark.pdf"), plot = p1, dpi = 400, width = 24, height = 13)


p1 <- ImageDimPlot(
  data_subset,
  group.by = "Cell_Type_Score",
  size = 1.5,
  cols = c("blue", "red"),
  border.color = "black",
  flip_xy = TRUE
)

p1 <- p1 + theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )

p2 <- ImageDimPlot(
  data_subset,
  group.by = "compartments",
  size = 1.5,
  cols = c("cortex" = "lightblue", "cmj" = "moccasin", "medulla" = "plum"),
  border.color = "black",
  flip_xy = TRUE
)

p2 <- p2 + theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )

p1 <- p1 + p2
p1

ggsave(filename = paste0(output_folder, "It_TEC_Ly6d_pos_neg_prediction_and_compartment_in_space_no_nas_bright.pdf"), plot = p1, dpi = 400, width = 24, height = 13)

print(table(data$Cell_Type_Score))

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

ggsave(filename = paste0(output_folder, "It_TEC_Ly6d_pos_neg_6_features_in_space.pdf") , plot = p1, dpi = 400, width = 12, height = 17)

table(data$Cell_Type_Score)

df <- data@meta.data[, c("Cell_Type_Score", "compartments")]
df <- df %>% filter(!is.na(Cell_Type_Score), !is.na(compartments))
df$compartment_cmj <- factor(df$compartments, levels = c("cortex", "cmj", "medulla"))

global_fractions <- df %>%
  group_by(Cell_Type_Score) %>%
  summarise(global_count = n(), .groups = "drop") %>%
  mutate(global_fraction = global_count / sum(global_count))

cell_counts <- df %>%
  group_by(compartment_cmj, Cell_Type_Score) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(global_fractions, by = "Cell_Type_Score") %>%
  mutate(
    fraction_within_compartment = count / global_count
  )

cell_density <- cell_counts %>%
  mutate(
    area_fraction = area_fractions[as.character(compartment_cmj)],
    compartment_fraction = compartment_fractions[as.character(compartment_cmj)],
    density_counts = fraction_within_compartment * compartment_fraction,
    density_area = fraction_within_compartment * area_fraction
  )

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

ggsave(filename = paste0(output_folder, "barchart_cell_type_counts_in_region.pdf"), plot = p1, dpi = 400, width = 10, height = 7)

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

ggsave(filename = paste0(output_folder, "barchart_cell_type_fraction_in_region.pdf"), plot = p1, dpi = 400, width = 10, height = 7)

p1 <- ggplot(cell_density, aes(x = compartment_cmj, y = density_counts, fill = Cell_Type_Score)) +
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

ggsave(filename = paste0(output_folder, "barchart_normalized_cell_type_counts_multiplied_with_cell_count_fraction_in_region.pdf"), plot = p1, dpi = 400, width = 10, height = 7)


p1 <- ggplot(cell_density, aes(x = compartment_cmj, y = density_area, fill = Cell_Type_Score)) +
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

ggsave(filename = paste0(output_folder, "barchart_normalized_cell_type_counts_multiplied_with_area_fraction_in_region.pdf"), plot = p1, dpi = 400, width = 10, height = 7)

writeLines(capture.output(sessionInfo()), paste0(output_folder, "SessionInfo.txt"))
