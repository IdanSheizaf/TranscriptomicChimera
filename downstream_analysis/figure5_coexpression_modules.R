# ==============================================================================
# MANUSCRIPT FIGURE 5: Conserved Co-Expression Modules (Unsupervised, 100 Genes)
# 
# This script identifies and visualizes the conserved co-expression modules 
# across all species, tissues, and moult phases using unsupervised ComplexHeatmap 
# on ComBat-seq batch-corrected VST expression of top 100 shared OGs.
#
# Output:
#   - Figure_5_Coexpression_Modules.pdf / .png (Unsupervised primary figure)
#   - Clusters_100genes_unsupervised.csv (Table S3 source)
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

log_header("Figure 5: Conserved Co-Expression Modules (Unsupervised, 100 Genes)")

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

# Select top 100 variable OGs for the primary manuscript figure (Figure 5)
topVarOGs <- head(order(matrixStats::rowVars(debatched_vst), decreasing = TRUE), 100)
mat <- debatched_vst[topVarOGs, ]
mat_scaled <- t(scale(t(mat)))

# --- 3. COMPLEXHEATMAP VISUALIZATION ---
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

out_dir <- get_output_path("Biclustering")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Unsupervised Biclustering Heatmap (Primary Manuscript Figure 5)
cat("Rendering Figure 5: Unsupervised Biclustering (100 Genes)...\n")
set.seed(42)
ht_unsupervised <- Heatmap(
  mat_scaled, 
  name = "Z-score",
  column_title = "Figure 5: Conserved Co-Expression Modules (Unsupervised Biclustering, 100 Genes)",
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
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(title = "Rel. Exp\n(Z-score)")
)

fig5_pdf <- file.path(out_dir, "Figure_5_Coexpression_Modules.pdf")
fig5_png <- file.path(out_dir, "Figure_5_Coexpression_Modules.png")

pdf(fig5_pdf, width = 13, height = 15)
ht_un_drawn <- draw(ht_unsupervised, merge_legend = TRUE)
dev.off()

png(fig5_png, width = 13, height = 15, units = "in", res = 300)
draw(ht_unsupervised, merge_legend = TRUE)
dev.off()

# Export Module Gene Lists & Table S3
ann_master <- load_annotations("maculatum") %>%
  dplyr::select(ODB_OG, Preferred_name) %>%
  dplyr::distinct(ODB_OG, .keep_all = TRUE)

row_orders_un <- row_order(ht_un_drawn)
cluster_100_list <- list()
for (i in seq_along(row_orders_un)) {
  module_df <- data.frame(ODB_OG = rownames(mat_scaled)[row_orders_un[[i]]]) %>%
    left_join(ann_master, by = "ODB_OG") %>%
    mutate(Cluster = paste0("Cluster_", i))
  cluster_100_list[[i]] <- module_df
  write_csv(module_df, file.path(out_dir, paste0("Unsupervised_Module_", i, "_Genes.csv")))
}

table_s3_df <- bind_rows(cluster_100_list)
table_s3_file <- file.path(out_dir, "Clusters_100genes_unsupervised.csv")
write_csv(table_s3_df, table_s3_file)

csv_out <- get_output_path()
if (dir.exists(csv_out)) {
  write_csv(table_s3_df, file.path(csv_out, "Clusters_100genes_unsupervised.csv"))
}
cat("  - Exported Table S3:", table_s3_file, "(", nrow(table_s3_df), " OGs)\n")

cat("\nSUCCESS: Figure 5 generated and saved to:\n")
cat("  -", fig5_pdf, "\n")
cat("  -", fig5_png, "\n")

