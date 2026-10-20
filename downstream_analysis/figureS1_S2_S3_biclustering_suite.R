# ==============================================================================
# MANUSCRIPT FIGURES S1, S2, S3: Biclustering Heatmap Suite
# 
# This script generates the full co-expression heatmap suite:
#   - Figure S1: Biclustering (100 genes, phase-structured)
#   - Figure S2: Biclustering (2,000 genes, fully unsupervised)
#   - Figure S3: Biclustering (2,000 genes, phase-structured)
#
# Output:
#   - Figure_S1_Biclustering_100genes_phase.pdf / .png
#   - Figure_S2_Biclustering_2000genes_unsupervised.pdf / .png
#   - Figure_S3_Biclustering_2000genes_phase.pdf / .png
#   - Clusters_100genes_unsupervised.csv (Table S3 source)
#   - Clusters_2000genes_unsupervised.csv (Table S4 source)
#   - Clusters_2000genes_phase.csv
# ==============================================================================

if (file.exists("utils.R")) {
  source("utils.R")
} else if (file.exists(file.path("downstream_analysis", "utils.R"))) {
  source(file.path("downstream_analysis", "utils.R"))
}

library(DESeq2)
library(sva)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(matrixStats)

log_header("Figures S1, S2, S3: Biclustering Heatmap Suite")

all_species <- c("maculatum", "laevis", "officinalis")
species_full_names <- c(
  "maculatum"   = "A. maculatum",
  "laevis"      = "P. laevis",
  "officinalis" = "A. officinalis"
)

# --- 1. LOAD & MERGE RAW COUNTS ACROSS SPECIES ---
get_og_data <- function(species_id) {
  cat("Loading counts for:", species_id, "...\n")
  dds_path <- get_data_path(species_id, "dge_analysis", "dds_group_design.rds")
  dds <- readRDS(dds_path)
  meta <- as.data.frame(colData(dds)) %>%
    dplyr::select(sample, tissue, phase, species) %>%
    dplyr::mutate(unique_sample_id = paste0(species, "_", sample))
  
  counts <- as.data.frame(counts(dds))
  colnames(counts) <- meta$unique_sample_id
  counts$gene_id <- rownames(counts)
  
  ann <- load_annotations(species_id) %>%
    dplyr::select(gene_id, ODB_OG) %>%
    dplyr::filter(!is.na(ODB_OG) & ODB_OG != "" & ODB_OG != "-")
  
  og_counts <- counts %>%
    dplyr::inner_join(ann, by = "gene_id") %>%
    dplyr::select(-gene_id) %>%
    dplyr::group_by(ODB_OG) %>%
    dplyr::summarise(across(everything(), sum), .groups = "drop")
  return(list(counts = og_counts, metadata = meta))
}

cat("Gathering shared orthologous counts...\n")
all_data <- lapply(all_species, get_og_data)
counts_list <- lapply(all_data, function(x) x$counts)
combined_matrix_raw <- Reduce(function(x, y) dplyr::inner_join(x, y, by = "ODB_OG"), counts_list)

meta_list <- lapply(all_data, function(x) x$metadata)
combined_meta <- dplyr::bind_rows(meta_list) %>%
  dplyr::filter(unique_sample_id %in% colnames(combined_matrix_raw))
combined_meta <- combined_meta[match(colnames(combined_matrix_raw)[-1], combined_meta$unique_sample_id), ]
combined_meta$phase <- factor(combined_meta$phase, levels = c("intermoult", "premoult", "intramoult", "postmoult"))

raw_counts_mat <- as.matrix(combined_matrix_raw[, -1])
rownames(raw_counts_mat) <- combined_matrix_raw$ODB_OG

# --- 2. COMBAT-SEQ BATCH CORRECTION & VST ---
cat("Applying ComBat-seq...\n")
corrected_counts <- ComBat_seq(raw_counts_mat, batch = combined_meta$species, group = combined_meta$phase)

cat("Performing VST...\n")
dds_corrected <- DESeqDataSetFromMatrix(countData = corrected_counts, colData = combined_meta, design = ~ phase + species)
vsd_corrected <- vst(dds_corrected, blind = FALSE)
debatched_vst <- assay(vsd_corrected)

# Load master annotation mapping
ann_master <- load_annotations("maculatum") %>%
  dplyr::select(ODB_OG, Preferred_name) %>%
  dplyr::distinct(ODB_OG, .keep_all = TRUE)

out_dir <- get_output_path("Supplementary_Figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Shared Column Annotation
ha_col_base <- HeatmapAnnotation(
  Phase = combined_meta$phase,
  Species = species_full_names[combined_meta$species],
  `Body Part` = combined_meta$tissue,
  col = list(
    Phase = moult_colors,
    Species = c("A. maculatum" = "#E41A1C", "P. laevis" = "#377EB8", "A. officinalis" = "#4DAF4A"),
    `Body Part` = c("f.legs" = "#F781BF", "h.legs" = "#FF7F00", "head" = "#A65628", "thorax" = "#999999")
  ),
  annotation_name_side = "left"
)

# ==============================================================================
# FIGURE S1: 100 Genes, Phase-Structured Biclustering
# ==============================================================================
cat("\nGenerating Figure S1: Biclustering (100 genes, phase-structured)...\n")
topVar100 <- head(order(rowVars(debatched_vst), decreasing = TRUE), 100)
mat100 <- debatched_vst[topVar100, ]
mat100_scaled <- t(scale(t(mat100)))

set.seed(42)
ht_s1 <- Heatmap(
  mat100_scaled,
  name = "Z-score",
  column_title = "Figure S1: Conserved Co-Expression Modules (100 Genes, Phase-Structured)",
  top_annotation = ha_col_base,
  row_km = 8,
  column_split = combined_meta$phase,
  row_gap = unit(1.5, "mm"),
  cluster_row_slices = TRUE,
  cluster_column_slices = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  heatmap_legend_param = list(title = "Rel. Exp\n(Z-score)")
)

pdf(file.path(out_dir, "Figure_S1_Biclustering_100genes_phase.pdf"), width = 12, height = 14)
ht_s1_drawn <- draw(ht_s1, merge_legend = TRUE)
dev.off()
png(file.path(out_dir, "Figure_S1_Biclustering_100genes_phase.png"), width = 12, height = 14, units = "in", res = 300)
draw(ht_s1, merge_legend = TRUE)
dev.off()

# Export 100-gene unsupervised cluster table (Table S3 source)
set.seed(42)
ht_s1_unsupervised <- Heatmap(
  mat100_scaled, row_km = 8, column_km = 4, cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE
)
ht_s1_un_drawn <- draw(ht_s1_unsupervised)
row_orders_100_un <- row_order(ht_s1_un_drawn)
cluster_100_list <- list()
for (i in seq_along(row_orders_100_un)) {
  df <- data.frame(ODB_OG = rownames(mat100_scaled)[row_orders_100_un[[i]]], Cluster = paste0("Cluster_", i)) %>%
    left_join(ann_master, by = "ODB_OG")
  cluster_100_list[[i]] <- df
}
write_csv(bind_rows(cluster_100_list), file.path(out_dir, "Clusters_100genes_unsupervised.csv"))

# ==============================================================================
# FIGURE S2: 2,000 Genes, Fully Unsupervised Biclustering
# ==============================================================================
cat("\nGenerating Figure S2: Biclustering (2,000 genes, unsupervised)...\n")
topVar2000 <- head(order(rowVars(debatched_vst), decreasing = TRUE), 2000)
mat2000 <- debatched_vst[topVar2000, ]
mat2000_scaled <- t(scale(t(mat2000)))

set.seed(42)
ht_s2 <- Heatmap(
  mat2000_scaled,
  name = "Z-score",
  column_title = "Figure S2: Universal Isopod Moult Biclustering (2,000 Genes, Unsupervised)",
  top_annotation = ha_col_base,
  row_km = 8,
  column_km = 4,
  row_gap = unit(1.5, "mm"),
  column_gap = unit(1.5, "mm"),
  cluster_row_slices = TRUE,
  cluster_column_slices = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  heatmap_legend_param = list(title = "Rel. Exp\n(Z-score)")
)

pdf(file.path(out_dir, "Figure_S2_Biclustering_2000genes_unsupervised.pdf"), width = 12, height = 14)
ht_s2_drawn <- draw(ht_s2, merge_legend = TRUE)
dev.off()
png(file.path(out_dir, "Figure_S2_Biclustering_2000genes_unsupervised.png"), width = 12, height = 14, units = "in", res = 300)
draw(ht_s2, merge_legend = TRUE)
dev.off()

# Export 2,000-gene cluster table (Table S4 source)
row_orders_2000_un <- row_order(ht_s2_drawn)
cluster_2000_list <- list()
for (i in seq_along(row_orders_2000_un)) {
  df <- data.frame(ODB_OG = rownames(mat2000_scaled)[row_orders_2000_un[[i]]], Cluster = paste0("Cluster_", i)) %>%
    left_join(ann_master, by = "ODB_OG")
  cluster_2000_list[[i]] <- df
}
write_csv(bind_rows(cluster_2000_list), file.path(out_dir, "Clusters_2000genes_unsupervised.csv"))

# ==============================================================================
# FIGURE S3: 2,000 Genes, Phase-Structured Biclustering
# ==============================================================================
cat("\nGenerating Figure S3: Biclustering (2,000 genes, phase-structured)...\n")
set.seed(42)
ht_s3 <- Heatmap(
  mat2000_scaled,
  name = "Z-score",
  column_title = "Figure S3: Universal Isopod Moult Biclustering (2,000 Genes, Phase-Structured)",
  top_annotation = ha_col_base,
  row_km = 8,
  column_split = combined_meta$phase,
  row_gap = unit(1.5, "mm"),
  cluster_row_slices = TRUE,
  cluster_column_slices = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  heatmap_legend_param = list(title = "Rel. Exp\n(Z-score)")
)

pdf(file.path(out_dir, "Figure_S3_Biclustering_2000genes_phase.pdf"), width = 12, height = 14)
ht_s3_drawn <- draw(ht_s3, merge_legend = TRUE)
dev.off()
png(file.path(out_dir, "Figure_S3_Biclustering_2000genes_phase.png"), width = 12, height = 14, units = "in", res = 300)
draw(ht_s3, merge_legend = TRUE)
dev.off()

row_orders_2000_st <- row_order(ht_s3_drawn)
cluster_2000_st_list <- list()
for (i in seq_along(row_orders_2000_st)) {
  df <- data.frame(ODB_OG = rownames(mat2000_scaled)[row_orders_2000_st[[i]]], Cluster = paste0("Phase_Module_", i)) %>%
    left_join(ann_master, by = "ODB_OG")
  cluster_2000_st_list[[i]] <- df
}
write_csv(bind_rows(cluster_2000_st_list), file.path(out_dir, "Clusters_2000genes_phase.csv"))

cat("\nSUCCESS: Figures S1, S2, and S3 generated successfully in:", out_dir, "\n")
