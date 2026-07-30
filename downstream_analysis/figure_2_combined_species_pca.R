# === COMBINED SPECIES PCA ANALYSIS ===
# This script merges VST counts from all species using Orthologous Groups (OGs)
# to create a single, unified PCA plot with full statistical validation.

# === SETUP WORKING DIRECTORY ===
project_path <- "G:/My Drive/PhD/Projects/Comparative Expressional Changes During the Moult Cycle in Land Isopods/local_scripts"
setwd(project_path)

# === LIBRARIES ===
source("utils.R")
library(DESeq2)
library(ggrepel)
library(patchwork)
library(grid)
library(vegan) # For PERMANOVA
library(sva)   # For ComBat-seq

# === SETUP ===
log_header("Combined Species PCA (ComBat-seq)")

all_species <- c("maculatum", "laevis", "officinalis")
gene_thresholds <- c(500, 1000, 2000, "all") 

# Mapping species to full binomial names
species_full_names <- c(
  "maculatum"   = "Armadillidium maculatum",
  "laevis"      = "Porcellio laevis",
  "officinalis" = "Armadillo officinalis"
)

# === FUNCTIONS ===

#' Load Raw Counts and Metadata, and Translate to OGs
get_og_data <- function(species_id) {
  cat("Processing:", species_id, "...\n")
  
  dds_path <- get_data_path(species_id, "dge_analysis", "dds_group_design.rds")
  if (!file.exists(dds_path)) stop("DDS not found for: ", species_id)
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
  
  # Sum raw counts by OG (integers required)
  og_counts <- counts %>%
    dplyr::inner_join(ann, by = "gene_id") %>%
    dplyr::select(-gene_id) %>%
    dplyr::group_by(ODB_OG) %>%
    dplyr::summarise(across(everything(), sum), .groups = "drop")
  
  return(list(counts = og_counts, metadata = meta))
}

# === 1. GATHER & MERGE DATA ===

all_data <- list()
for (sp in all_species) {
  all_data[[sp]] <- get_og_data(sp)
}

cat("\nMerging species counts by shared Orthologous Groups...\n")
counts_list <- lapply(all_data, function(x) x$counts)
combined_matrix_raw <- Reduce(function(x, y) dplyr::inner_join(x, y, by = "ODB_OG"), counts_list)

cat("Combining metadata...\n")
meta_list <- lapply(all_data, function(x) x$metadata)
combined_meta <- dplyr::bind_rows(meta_list) %>%
  dplyr::filter(unique_sample_id %in% colnames(combined_matrix_raw))

# Synchronize order
combined_meta <- combined_meta[match(colnames(combined_matrix_raw)[-1], combined_meta$unique_sample_id), ]
combined_meta$phase <- factor(combined_meta$phase, levels = c("intermoult", "premoult", "intramoult", "postmoult"))
combined_meta$species_full <- species_full_names[combined_meta$species]

# Extract raw matrix for ComBat-seq
raw_counts_mat <- as.matrix(combined_matrix_raw[, -1])
rownames(raw_counts_mat) <- combined_matrix_raw$ODB_OG

cat("Shared OGs for ComBat-seq:", nrow(raw_counts_mat), "\n")

# === 2. APPLY COMBAT-SEQ BATCH CORRECTION ===
cat("\nApplying ComBat-seq to remove Species batch effects...\n")
# batch: Species
# group: Phase (to preserve the biological signal)
corrected_counts <- ComBat_seq(raw_counts_mat, batch=combined_meta$species, group=combined_meta$phase)

# === 3. VST TRANSFORMATION ON CORRECTED COUNTS ===
cat("Performing VST on corrected counts...\n")
# Create a dummy DESeq2 object for VST
# Note: ComBat-seq outputs counts, so we use DESeqDataSetFromMatrix
dds_corrected <- DESeqDataSetFromMatrix(countData = corrected_counts,
                                        colData = combined_meta,
                                        design = ~ phase + species)
vsd_corrected <- vst(dds_corrected, blind=FALSE)
pca_input_full <- assay(vsd_corrected)

# === 4. GLOBAL STATISTICAL VALIDATION ===
cat("\n--- Running Global PERMANOVA on Corrected Data ---\n")
global_perm <- adonis2(t(pca_input_full) ~ phase + species, data = combined_meta, permutations = 999)
print(global_perm)

perm_out_dir <- file.path(get_project_root(), "PCA", "Combined")
dir.create(perm_out_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(print(global_perm)), file.path(perm_out_dir, "Combined_PERMANOVA_Results.txt"))

# === 2. LOOP THROUGH GENE THRESHOLDS ===

for (n_genes in gene_thresholds) {
  cat("\n--- Analyzing", n_genes, "Genes ---\n")
  
  if (n_genes == "all") {
    pca_data <- pca_input_full
    actual_n <- nrow(pca_input_full)
  } else {
    n_genes_num <- as.numeric(n_genes)
    actual_n <- min(n_genes_num, nrow(pca_input_full))
    topVarOGs <- head(order(rowVars(pca_input_full), decreasing = TRUE), actual_n)
    pca_data <- pca_input_full[topVarOGs, ]
  }
  
  pca_full <- prcomp(t(pca_data))
  pca_var  <- round(100 * summary(pca_full)$importance[2,], 1)
  
  out_dir <- file.path(get_project_root(), "PCA", "Combined", paste0("top_", n_genes))
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
  }
  
  # Visualization helper
  plot_any_pca <- function(pca_res, pca_var_pct, x_pc = 1, y_pc = 2, color_var = "phase", title_prefix = "Combined PCA") {
    plot_df <- data.frame(
      PCx = pca_res$x[, x_pc],
      PCy = pca_res$x[, y_pc],
      species = combined_meta$species_full,
      phase = combined_meta$phase,
      tissue = combined_meta$tissue
    )
    
    sub_text <- paste0("Included OGs: ", actual_n)
    if (exists("global_perm")) {
       sub_text <- paste0(sub_text, " | Global PERMANOVA (Phase): R² = ", round(global_perm$R2[1], 2), 
                         ", p = ", format.pval(global_perm$`Pr(>F)`[1], digits = 3))
    }

    p <- ggplot(plot_df, aes(x = PCx, y = PCy, color = !!sym(color_var), shape = species)) +
      geom_point(size = 4, alpha = 0.7) +
      theme_bw() +
      labs(
        title = paste0(title_prefix, ": PC", x_pc, " vs PC", y_pc),
        subtitle = sub_text,
        x = paste0("PC", x_pc, ": ", pca_var_pct[x_pc], "% variance"),
        y = paste0("PC", y_pc, ": ", pca_var_pct[y_pc], "% variance"),
        shape = "Species"
      ) +
      theme(plot.title = element_text(face = "bold"),
            legend.text = element_text(face = "italic"))
    
    if (color_var == "phase") {
      p <- p + scale_color_manual(values = c("intermoult"="#440154", "premoult"="#31688E", "intramoult"="#35B779", "postmoult"="#FDE725"))
    }
    return(p)
  }

  ggsave(file.path(out_dir, "Combined_PC1v2_phase.png"),  plot_any_pca(pca_full, pca_var, 1, 2, "phase"),  width = 10, height = 7, create.dir = TRUE)
  ggsave(file.path(out_dir, "Combined_PC1v3_phase.png"),  plot_any_pca(pca_full, pca_var, 1, 3, "phase"),  width = 10, height = 7, create.dir = TRUE)
  ggsave(file.path(out_dir, "Combined_PC2v3_phase.png"),  plot_any_pca(pca_full, pca_var, 2, 3, "phase"),  width = 10, height = 7, create.dir = TRUE)

  # Added: Plots by Tissue
  ggsave(file.path(out_dir, "Combined_PC1v2_tissue.png"), plot_any_pca(pca_full, pca_var, 1, 2, "tissue"), width = 10, height = 7, create.dir = TRUE)
  ggsave(file.path(out_dir, "Combined_PC1v3_tissue.png"), plot_any_pca(pca_full, pca_var, 1, 3, "tissue"), width = 10, height = 7, create.dir = TRUE)
  ggsave(file.path(out_dir, "Combined_PC2v3_tissue.png"), plot_any_pca(pca_full, pca_var, 2, 3, "tissue"), width = 10, height = 7, create.dir = TRUE)

  # PC Boxplots with Stats
  pc_boxplots <- list()
  for (i in 1:min(7, ncol(pca_full$x))) {
    plot_df_pc <- data.frame(score = pca_full$x[, i], phase = combined_meta$phase, species = combined_meta$species_full)
    kw_test <- kruskal.test(score ~ phase, data = plot_df_pc)
    
    pc_boxplots[[i]] <- ggplot(plot_df_pc, aes(x = phase, y = score, fill = phase)) +
      geom_violin(alpha = 0.3, color = NA) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.2) + 
      geom_jitter(aes(shape = species), width = 0.1, size = 1.2, alpha = 0.5) +
      theme_bw() + scale_fill_manual(values = c("intermoult"="#440154", "premoult"="#31688E", "intramoult"="#35B779", "postmoult"="#FDE725")) +
      labs(title = paste0("PC", i, " (", pca_var[i], "%)"), 
           subtitle = paste0("Kruskal-Wallis p = ", format.pval(kw_test$p.value, digits = 3)),
           x = NULL, y = NULL) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 7), plot.subtitle = element_text(size = 8))
  }
  ggsave(file.path(out_dir, "Combined_PC1_7_Boxplots.png"), wrap_plots(pc_boxplots, ncol = 4), width = 14, height = 8, create.dir = TRUE)

  # --- DRIVER EXTRACTION (New & Robust) ---
  # We save drivers for PC1-5 to identify different waves of the toolkit.
  cat("  Saving PC1-5 Drivers (ComBat Corrected)...\n")
  
  master_ann_loop <- load_general_annotations() %>%
    dplyr::select(ODB_OG, Preferred_name) %>%
    dplyr::distinct(ODB_OG, .keep_all = TRUE)
  
  # Extract top 100 drivers for each of the first 5 PCs
  driver_list <- list()
  for (i in 1:min(5, ncol(pca_full$rotation))) {
    loadings <- pca_full$rotation[, i]
    top_indices <- order(abs(loadings), decreasing = TRUE)[1:100]
    
    driver_list[[i]] <- data.frame(
      ODB_OG = rownames(pca_full$rotation)[top_indices],
      Loading = loadings[top_indices],
      PC = paste0("PC", i)
    ) %>%
      dplyr::left_join(master_ann_loop, by = "ODB_OG") %>%
      dplyr::select(ODB_OG, Preferred_name, Loading, PC)
  }
  
  all_drivers <- dplyr::bind_rows(driver_list)
  write.csv(all_drivers, file.path(out_dir, "Global_Drivers_ComBat_Corrected.csv"), row.names = FALSE)
}


cat("\nAnalysis complete!\n")
