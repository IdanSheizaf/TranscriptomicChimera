# === UNIVERSAL BICLUSTERING VIA COMPLEXHEATMAP ===
# This script performs "biclustering-like" partitioning on batch-corrected 
# (ComBat-seq) VST counts to identify conserved co-expression modules.

# === SETUP ===
project_path <- "G:/My Drive/PhD/Projects/Comparative Expressional Changes During the Moult Cycle in Land Isopods/local_scripts"
setwd(project_path)
source("utils.R")
library(DESeq2)
library(sva)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

log_header("Biclustered Heatmap Analysis")

# --- 1. PREPARE DE-BATCHED DATA ---
all_species <- c("maculatum", "laevis", "officinalis")
species_full_names <- c(
  "maculatum"   = "A. maculatum",
  "laevis"      = "P. laevis",
  "officinalis" = "A. officinalis"
)

get_og_data <- function(species_id) {
  cat("Processing:", species_id, "...\n")
  dds_path <- get_data_path(species_id, "dge_analysis", "dds_group_design.rds")
  dds <- readRDS(dds_path)
  meta <- as.data.frame(colData(dds)) %>%
    dplyr::select(sample, tissue, phase, species) %>%
    dplyr::mutate(unique_sample_id = paste0(species, "_", sample))
  
  # Use RAW counts for ComBat-seq
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

cat("Gathering and merging raw counts across species...\n")
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

cat("\nApplying ComBat-seq batch correction...\n")
corrected_counts <- ComBat_seq(raw_counts_mat, batch=combined_meta$species, group=combined_meta$phase)

cat("Performing VST on corrected counts...\n")
dds_corrected <- DESeqDataSetFromMatrix(countData = corrected_counts, colData = combined_meta, design = ~ phase + species)
vsd_corrected <- vst(dds_corrected, blind=FALSE)
debatched_vst <- assay(vsd_corrected)

# Filter for Top 2000 high-variance OGs for a clean yet comprehensive heatmap
topVarOGs <- head(order(rowVars(debatched_vst), decreasing = TRUE), 500)
mat <- debatched_vst[topVarOGs, ]

# Z-score scaling
mat_scaled <- t(scale(t(mat)))

# --- 2. GENERATE HEATMAPS ---

# Shared Metadata Annotation
ha_col_base <- HeatmapAnnotation(
  Phase = combined_meta$phase,
  Species = species_full_names[combined_meta$species],
  Tissue = combined_meta$tissue,
  col = list(
    Phase = c("intermoult"="#440154", "premoult"="#31688E", "intramoult"="#35B779", "postmoult"="#FDE725"),
    Species = c("A. maculatum"="#E41A1C", "P. laevis"="#377EB8", "A. officinalis"="#4DAF4A"),
    Tissue = c("f.legs"="#F781BF", "h.legs"="#FF7F00", "head"="#A65628", "thorax"="#999999")
  ),
  annotation_name_side = "left"
)

out_dir <- "output_csv/Biclustering"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ann_master <- load_general_annotations() %>%
  dplyr::select(ODB_OG, Preferred_name) %>%
  dplyr::distinct(ODB_OG, .keep_all = TRUE)


# --- VERSION A: FULLY UNSUPERVISED PARTITIONED ---
cat("\nGenerating Version A: Unsupervised Partitioned Heatmap...\n")
set.seed(42)
ht_unsupervised <- Heatmap(mat_scaled, 
                          name = "Z-score",
                          column_title = "Universal Isopod Moult Biclustering (Fully Unsupervised)",
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
                          column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                          heatmap_legend_param = list(title = "Rel. Exp"))

pdf(file.path(out_dir, "Biclustering_VersionA_Unsupervised.pdf"), width = 12, height = 14)
ht_un_drawn <- draw(ht_unsupervised, merge_legend = TRUE)
dev.off()

# Export Unsupervised Modules
row_orders_un <- row_order(ht_un_drawn)
for (i in seq_along(row_orders_un)) {
  module_df <- data.frame(ODB_OG = rownames(mat_scaled)[row_orders_un[[i]]]) %>%
    left_join(ann_master, by = "ODB_OG") %>%
    mutate(Module = paste0("Unsupervised_Module_", i))
  write_csv(module_df, file.path(out_dir, paste0("Unsupervised_Module_", i, "_Genes.csv")))
}


# --- VERSION B: PHASE-STRUCTURED PARTITIONED ---
cat("\nGenerating Version B: Phase-Structured Heatmap...\n")
set.seed(42)
# Here we split columns by the biological factor 'phase'
ht_structured <- Heatmap(mat_scaled, 
                        name = "Z-score",
                        column_title = "Universal Isopod Moult Biclustering (Phase-Structured)",
                        top_annotation = ha_col_base,
                        row_km = 8, 
                        column_split = combined_meta$phase,
                        row_gap = unit(1.5, "mm"),
                        cluster_row_slices = TRUE, 
                        cluster_column_slices = FALSE, # Keep biological order
                        show_row_names = FALSE,
                        show_column_names = FALSE,
                        cluster_rows = TRUE,
                        cluster_columns = TRUE,
                        show_row_dend = TRUE,
                        show_column_dend = TRUE,
                        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                        heatmap_legend_param = list(title = "Rel. Exp"))

pdf(file.path(out_dir, "Biclustering_VersionB_PhaseStructured.pdf"), width = 12, height = 14)
ht_st_drawn <- draw(ht_structured, merge_legend = TRUE)
dev.off()

# Export Structured Modules
row_orders_st <- row_order(ht_st_drawn)
for (i in seq_along(row_orders_st)) {
  module_df <- data.frame(ODB_OG = rownames(mat_scaled)[row_orders_st[[i]]]) %>%
    left_join(ann_master, by = "ODB_OG") %>%
    mutate(Module = paste0("Structured_Module_", i))
  write_csv(module_df, file.path(out_dir, paste0("Structured_Module_", i, "_Genes.csv")))
}

cat("\nAnalysis complete! Both versions saved in:", out_dir, "\n")

