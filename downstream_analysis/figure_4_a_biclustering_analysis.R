# === MASS BICLUSTERING ANALYSIS ===
# This script automates 32 variations of biclustering heatmaps
# (8 gene counts x 4 ordering strategies) as requested by the user.
# Tissues included: All (including h.legs)

# === SETUP ===
project_path <- "G:/My Drive/PhD/Projects/Comparative Expressional Changes During the Moult Cycle in Land Isopods/local_scripts/bicluster_23.6"
setwd(project_path)
source("utils.R")
library(DESeq2)
library(sva)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(matrixStats)

log_header("Mass Biclustering Analysis (32 Variations - All Tissues)")

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
    dplyr::mutate(unique_sample_id = paste0(species_id, "_", sample))
  
  # Use RAW counts for ComBat-seq
  # We use the entire counts matrix because we are including ALL samples
  counts_raw <- counts(dds)
  counts <- as.data.frame(counts_raw)
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
all_data_raw <- lapply(all_species, get_og_data)

counts_list <- lapply(all_data_raw, function(x) x$counts)
combined_matrix_raw <- Reduce(function(x, y) dplyr::inner_join(x, y, by = "ODB_OG"), counts_list)

meta_list <- lapply(all_data_raw, function(x) x$metadata)
combined_meta <- dplyr::bind_rows(meta_list) %>%
  dplyr::filter(unique_sample_id %in% colnames(combined_matrix_raw))

# Reorder meta to match matrix columns
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

# --- 2. DEFINE VARIATIONS ---
gene_counts <- c(50, 100, 200, 500, 1000, 2000, 5000, "all")
strategies <- c("unsupervised", "species", "tissue", "phase")

out_base_dir <- "output_csv/Mass_Biclustering_All_Tissues"
if (!dir.exists(out_base_dir)) dir.create(out_base_dir, recursive = TRUE)

ann_master <- load_general_annotations() %>%
  dplyr::select(ODB_OG, Preferred_name) %>%
  dplyr::distinct(ODB_OG, .keep_all = TRUE)

# Shared Heatmap Annotation
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

# --- 3. LOOP AND GENERATE ---
for (n in gene_counts) {
  cat("\n--- Processing Gene Count:", n, "---\n")
  
  # Select Genes
  if (n == "all") {
    current_mat <- debatched_vst
  } else {
    n_val <- as.numeric(n)
    if (n_val > nrow(debatched_vst)) n_val <- nrow(debatched_vst)
    topVarOGs <- head(order(rowVars(debatched_vst), decreasing = TRUE), n_val)
    current_mat <- debatched_vst[topVarOGs, ]
  }
  
  # Z-score scaling
  mat_scaled <- t(scale(t(current_mat)))
  
  for (strat in strategies) {
    cat("  Strategy:", strat, "...")
    
    file_name <- paste0("Biclustering_", n, "genes_", strat, ".pdf")
    title_text <- paste0("Biclustering (", n, " genes, ", strat, ")")
    
    # Configure Heatmap Parameters based on Strategy
    if (strat == "unsupervised") {
      ht <- Heatmap(mat_scaled, 
                    name = "Z-score",
                    column_title = title_text,
                    top_annotation = ha_col_base,
                    row_km = min(8, nrow(mat_scaled)), 
                    column_km = 5
                    ,
                    row_gap = unit(1, "mm"),
                    column_gap = unit(1, "mm"),
                    cluster_row_slices = TRUE, 
                    cluster_column_slices = TRUE,
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    cluster_rows = TRUE,
                    cluster_columns = TRUE)
    } else {
      # For species, tissue, phase: split columns by the factor
      col_split_factor <- combined_meta[[strat]]
      ht <- Heatmap(mat_scaled, 
                    name = "Z-score",
                    column_title = title_text,
                    top_annotation = ha_col_base,
                    row_km = min(8, nrow(mat_scaled)), 
                    column_split = col_split_factor,
                    cluster_row_slices = TRUE, 
                    cluster_column_slices = FALSE, # Keep biological order
                    show_row_names = FALSE,
                    show_column_names = FALSE,
                    cluster_rows = TRUE,
                    cluster_columns = TRUE)
    }
    
    # Save PDF
    pdf_path <- file.path(out_base_dir, file_name)
    pdf(pdf_path, width = 12, height = 14)
    set.seed(42)
    ht_drawn <- draw(ht, merge_legend = TRUE)
    dev.off()
    
    # --- EXTRACT AND SAVE CLUSTERS ---
    row_orders <- row_order(ht_drawn)
    cluster_list <- list()
    for (i in seq_along(row_orders)) {
      cluster_list[[i]] <- data.frame(
        ODB_OG = rownames(mat_scaled)[row_orders[[i]]],
        Cluster = i
      )
    }
    
    all_clusters_df <- dplyr::bind_rows(cluster_list) %>%
      dplyr::left_join(ann_master, by = "ODB_OG")
    
    csv_name <- paste0("Clusters_", n, "genes_", strat, ".csv")
    write_csv(all_clusters_df, file.path(out_base_dir, csv_name))
    
    cat(" DONE.\n")
  }
}

cat("\nMass Biclustering Analysis Complete! All 32 heatmaps saved in:", out_base_dir, "\n")
