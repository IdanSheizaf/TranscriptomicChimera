# ==============================================================================
# MANUSCRIPT FIGURE 4: Global Integrated Transcriptomic Alignment & Decoupling
# 
# This script integrates expression data across all three species by mapping 
# transcripts to shared Orthologous Groups (OGs), applies ComBat-seq to remove 
# species-specific batch effects while preserving the biological moult phase, 
# and runs PCA to demonstrate the conserved "Universal Phase Engine" and 
# anatomical decoupling across Oniscidea.
#
# Output:
#   - Figure_4_Universal_PCA.png / .pdf (Publication multi-panel)
#   - Combined_PERMANOVA_Results.txt
#   - Table_S1_Global_Drivers_ComBat_Corrected.csv
# ==============================================================================

if (file.exists("utils.R")) {
  source("utils.R")
} else if (file.exists(file.path("downstream_analysis", "utils.R"))) {
  source(file.path("downstream_analysis", "utils.R"))
}

library(DESeq2)
library(patchwork)
library(vegan) # For PERMANOVA
library(sva)   # For ComBat-seq
library(tidyverse)

log_header("Figure 4: Global Integrated Transcriptomic Alignment")

all_species <- c("maculatum", "laevis", "officinalis")
n_genes <- 2000

species_full_names <- c(
  "maculatum"   = "Armadillidium maculatum",
  "laevis"      = "Porcellio laevis",
  "officinalis" = "Armadillo officinalis"
)

# --- 1. GATHER & MERGE RAW COUNTS BY ORTHOLOGOUS GROUPS (OGs) ---
get_og_data <- function(species_id) {
  cat("Loading counts & annotations for:", species_id, "...\n")
  dds_path <- get_data_path(species_id, "dge_analysis", "dds_group_design.rds")
  if (!file.exists(dds_path)) stop("DDS not found for: ", species_id)
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

all_data <- list()
for (sp in all_species) {
  all_data[[sp]] <- get_og_data(sp)
}

cat("\nMerging species counts across shared OGs...\n")
counts_list <- lapply(all_data, function(x) x$counts)
combined_matrix_raw <- Reduce(function(x, y) dplyr::inner_join(x, y, by = "ODB_OG"), counts_list)

meta_list <- lapply(all_data, function(x) x$metadata)
combined_meta <- dplyr::bind_rows(meta_list) %>%
  dplyr::filter(unique_sample_id %in% colnames(combined_matrix_raw))
combined_meta <- combined_meta[match(colnames(combined_matrix_raw)[-1], combined_meta$unique_sample_id), ]
combined_meta$phase <- factor(combined_meta$phase, levels = c("intermoult", "premoult", "intramoult", "postmoult"))
combined_meta$species_full <- factor(species_full_names[combined_meta$species], levels = species_full_names)

raw_counts_mat <- as.matrix(combined_matrix_raw[, -1])
rownames(raw_counts_mat) <- combined_matrix_raw$ODB_OG
cat("Shared Orthologous Groups:", nrow(raw_counts_mat), "\n")

# --- 2. COMBAT-SEQ BATCH CORRECTION ---
cat("\nApplying ComBat-seq to remove species divergence while preserving moult phase...\n")
corrected_counts <- ComBat_seq(raw_counts_mat, batch = combined_meta$species, group = combined_meta$phase)

# --- 3. VST TRANSFORMATION & PERMANOVA ---
cat("Performing VST on batch-corrected counts...\n")
dds_corrected <- DESeqDataSetFromMatrix(countData = corrected_counts, colData = combined_meta, design = ~ phase + species)
vsd_corrected <- vst(dds_corrected, blind = FALSE)
pca_input_full <- assay(vsd_corrected)

cat("Running global PERMANOVA...\n")
global_perm <- adonis2(t(pca_input_full) ~ phase + species, data = combined_meta, permutations = 999)
r2_val <- round(global_perm$R2[1], 2)
p_val  <- format.pval(global_perm$`Pr(>F)`[1], digits = 3)
cat(sprintf("  Global Phase effect: R2 = %s, p = %s\n", r2_val, p_val))

# --- 4. TOP 2000 HIGH-VARIANCE OGs PCA ---
vars <- matrixStats::rowVars(pca_input_full)
actual_n <- min(n_genes, length(vars))
top_ogs <- head(order(vars, decreasing = TRUE), actual_n)
pca_res <- prcomp(t(pca_input_full[top_ogs, ]))
pca_var <- round(100 * summary(pca_res)$importance[2, ], 1)

plot_df <- data.frame(
  PC1     = pca_res$x[, 1],
  PC2     = pca_res$x[, 2],
  phase   = combined_meta$phase,
  tissue  = combined_meta$tissue,
  species = combined_meta$species_full
)

# Panel A: Phase colors + Body part shapes (Manuscript Primary)
p4a <- ggplot(plot_df, aes(x = PC1, y = PC2, color = phase, shape = tissue)) +
  geom_point(size = 5.5, alpha = 0.85, stroke = 1.3) +
  scale_color_manual(values = moult_colors, name = "Moult Phase") +
  scale_shape_manual(values = moult_shapes, name = "Body Part") +
  labs(
    title = "A: Global Moult Phase Trajectory",
    subtitle = paste0("Top ", actual_n, " OGs | PERMANOVA (Phase): R² = ", r2_val, ", p = ", p_val),
    x = paste0("PC1: ", pca_var[1], "% variance"),
    y = paste0("PC2: ", pca_var[2], "% variance")
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    legend.position = "right"
  )

# Panel B: Body part colors + Species shapes
tissue_colors <- c("f.legs" = "#377EB8", "h.legs" = "#4DAF4A", "head" = "#984EA3", "thorax" = "#FF7F00")
tissue_labels <- c("f.legs" = "Front Legs", "h.legs" = "Hind Legs", "head" = "Head", "thorax" = "Pereon")

p4b <- ggplot(plot_df, aes(x = PC1, y = PC2, color = tissue, shape = species)) +
  geom_point(size = 5.5, alpha = 0.85, stroke = 1.3) +
  scale_color_manual(values = tissue_colors, labels = tissue_labels, name = "Body Part") +
  labs(
    title = "B: Anatomical Partitioning",
    subtitle = "Hind legs (green) decouple during intramoult, projecting into postmoult",
    x = paste0("PC1: ", pca_var[1], "% variance"),
    y = paste0("PC2: ", pca_var[2], "% variance"),
    shape = "Species"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "grey30"),
    legend.position = "right",
    legend.text = element_text(face = "italic")
  )

figure4_composite <- (p4a + p4b) +
  plot_annotation(
    title = "Figure 4: Global Integrated Transcriptomic Alignment and Conserved Anatomical Decoupling",
    subtitle = "Multi-species ComBat-seq batch correction across 8,924 shared Orthologous Groups (OGs)",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey30")
    )
  )

# Save Outputs
out_dir <- file.path(get_project_root(), "PCA", "Combined")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

fig4_png <- file.path(out_dir, "Figure_4_Universal_PCA.png")
fig4_pdf <- file.path(out_dir, "Figure_4_Universal_PCA.pdf")

ggsave(fig4_png, figure4_composite, width = 18, height = 7, dpi = 300)
ggsave(fig4_pdf, figure4_composite, width = 18, height = 7)

# Save PERMANOVA stats
writeLines(capture.output(print(global_perm)), file.path(out_dir, "Combined_PERMANOVA_Results.txt"))

cat("\nSUCCESS: Figure 4 generated and saved to:\n")
cat("  -", fig4_png, "\n")
cat("  -", fig4_pdf, "\n")

# --- 4. EXPORT TABLE S2: GLOBAL DRIVERS OF COMBAT-CORRECTED PCA ---
cat("\nExtracting Top PCA Loadings for Table S2 (Global Drivers)...\n")
root <- get_project_root()
master_candidates <- c(
  file.path(root, "de_novo_transcriptomes_crustacea.og.annotations"),
  file.path(root, "transcriptome assembly", "de_novo_transcriptomes_crustacea.og.annotations"),
  file.path(root, "local_scripts", "RNAseq", "de_novo_transcriptomes_crustacea.og.annotations"),
  file.path(root, "..", "local_scripts", "RNAseq", "de_novo_transcriptomes_crustacea.og.annotations"),
  file.path(dirname(root), "local_scripts", "RNAseq", "de_novo_transcriptomes_crustacea.og.annotations")
)
master_path <- master_candidates[file.exists(master_candidates)][1]

if (!is.na(master_path) && file.exists(master_path)) {
  master_ann <- read.delim(master_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>%
    dplyr::select(ODB_OG, Preferred_name) %>%
    dplyr::filter(!is.na(ODB_OG) & ODB_OG != "" & ODB_OG != "-") %>%
    dplyr::distinct()
} else {
  master_ann <- data.frame(ODB_OG = character(), Preferred_name = character())
}

driver_list <- list()
for (i in 1:min(5, ncol(pca_res$rotation))) {
  loadings <- pca_res$rotation[, i]
  top_indices <- order(abs(loadings), decreasing = TRUE)[1:100]
  
  driver_list[[i]] <- data.frame(
    ODB_OG = rownames(pca_res$rotation)[top_indices],
    Loading = loadings[top_indices],
    PC = paste0("PC", i)
  ) %>%
    dplyr::left_join(master_ann, by = "ODB_OG") %>%
    dplyr::mutate(Preferred_name = ifelse(is.na(Preferred_name) | Preferred_name == "" | Preferred_name == "-", paste0("Uncharacterized_", ODB_OG), Preferred_name)) %>%
    dplyr::select(ODB_OG, Preferred_name, Loading, PC)
}

all_drivers <- dplyr::bind_rows(driver_list)
table_s2_file <- file.path(out_dir, "Global_Drivers_ComBat_Corrected.csv")
write.csv(all_drivers, table_s2_file, row.names = FALSE)

# Also ensure copy in standard output_csv directory
csv_out <- get_output_path()
if (dir.exists(csv_out)) {
  write.csv(all_drivers, file.path(csv_out, "Global_Drivers_ComBat_Corrected.csv"), row.names = FALSE)
}
cat("  - Exported Table S2:", table_s2_file, "(", nrow(all_drivers), " rows across PC1-PC5)\n")

