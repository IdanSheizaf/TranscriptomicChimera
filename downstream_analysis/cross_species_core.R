# === CROSS-SPECIES CORE ANALYSIS ===
# This script identifies genes that are significantly regulated across
# all three isopod species for multiple phase transitions.

# === SETUP ===
if (file.exists("utils.R")) {
  source("utils.R")
} else if (file.exists(file.path("downstream_analysis", "utils.R"))) {
  source(file.path("downstream_analysis", "utils.R"))
}
library(tidyverse)
library(ggVennDiagram)
library(patchwork)
library(DESeq2)
library(matrixStats)

log_header("Comprehensive Cross-Species Core")

all_species <- c("maculatum", "laevis", "officinalis")

# Abbreviated names for the lists
sp_names <- c("maculatum" = "A. maculatum", "laevis" = "P. laevis", "officinalis" = "A. officinalis")

target_comparisons <- c(
  "type1d_premoult_vs_intermoult_no_hlegs",
  "type1d_intramoult_vs_intermoult_no_hlegs",
  "type1d_postmoult_vs_intermoult_no_hlegs"
)

# === FUNCTIONS ===

get_sig_genes <- function(species, comparison) {
  cat("  [", species, "] ", sep="")
  # Using Type 1D results (no h.legs) as per project refinement
  dge_path <- get_data_path(species, "dge_analysis/dge_type1d_phase_vs_intermoult_no_hlegs/full", paste0(comparison, ".csv"))
  if (!file.exists(dge_path)) return(NULL)
  
  dge <- read_csv(dge_path, show_col_types = FALSE)
  annotations <- load_annotations(species)
  
  sig_set <- dge %>%
    inner_join(annotations, by = c("transcript_id" = "gene_id")) %>%
    filter(padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
    filter(!is.na(ODB_OG) & ODB_OG != "" & ODB_OG != "-") %>%
    pull(ODB_OG) %>% unique()
  
  cat("Found", length(sig_set), "sig OGs.\n")
  return(sig_set)
}

# === EXECUTION ===

# Load Master Mapping for ODB_OG to Preferred_name
cat("\nLoading master OG mapping...\n")
root <- get_project_root()
master_candidates <- c(
  file.path(root, "de_novo_transcriptomes_crustacea.og.annotations"),
  file.path(root, "transcriptome assembly", "de_novo_transcriptomes_crustacea.og.annotations"),
  file.path(root, "local_scripts", "RNAseq", "de_novo_transcriptomes_crustacea.og.annotations")
)
master_ann_path <- master_candidates[file.exists(master_candidates)][1]
master_mapping <- read.delim(master_ann_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>%
  select(ODB_OG, Preferred_name) %>%
  filter(!is.na(ODB_OG) & ODB_OG != "" & ODB_OG != "-") %>%
  distinct(ODB_OG, .keep_all = TRUE)

all_venn_plots <- list()
master_core_list <- list()

for (comp in target_comparisons) {
  cat("\nAnalyzing Comparison:", comp, "\n")
  
  og_lists <- list()
  for (sp in all_species) {
    sig <- get_sig_genes(sp, comp)
    if (!is.null(sig)) {
      og_lists[[sp_names[sp]]] <- sig
    } else {
      cat("    Warning: Missing data for", sp, "\n")
    }
  }
  
  if (length(og_lists) < 3) {
    cat("  SKIPPING: Less than 3 species have data for this comparison.\n")
    next
  }
  
  phase_name <- str_extract(comp, "(?<=type1[bd]_)[a-z]+(?=_vs)")
  conserved_core_ogs <- Reduce(intersect, og_lists)
  
  cat("  Found", length(conserved_core_ogs), "OGs in the 3-way conserved core.\n")
  
  # Map OGs back to names for the CSV output
  core_df <- data.frame(ODB_OG = conserved_core_ogs) %>%
    left_join(master_mapping, by = "ODB_OG") %>%
    mutate(Phase = phase_name) %>%
    select(ODB_OG, Preferred_name, Phase)
  
  master_core_list[[phase_name]] <- core_df
  
  # --- CREATE VENN PLOT ---
  # We use label = "both" to show Count and Percentage
  p <- ggVennDiagram(og_lists, label = "both", label_alpha = 0, label_size = 5, edge_size = 0.8) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1) + 
    labs(title = str_to_title(phase_name)) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 28, margin = margin(b=10, t=10)),
      plot.margin = margin(t=20, b=20, l=50, r=50) 
    ) +
    coord_fixed(clip = "off")
  
  # FORCING ITALICS via layer modification (This part worked)
  for(i in 1:length(p$layers)){
    if (class(p$layers[[i]]$geom)[1] == "GeomText" && 
        (!is.null(p$layers[[i]]$mapping$label) && rlang::as_label(p$layers[[i]]$mapping$label) != "label") || 
        is.null(p$layers[[i]]$mapping$label)) {
      p$layers[[i]]$aes_params$fontface <- "italic"
      p$layers[[i]]$aes_params$size <- 7 
    }
  }
  
  all_venn_plots[[comp]] <- p
}

# === FINAL OUTPUTS ===

out_dir_path <- normalizePath(file.path(getwd(), "output_csv"), mustWork = FALSE)
if (!dir.exists(out_dir_path)) dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)

if (length(master_core_list) > 0) {
  write_csv(bind_rows(master_core_list), file.path(out_dir_path, "conserved_core_master_all_phases.csv"))
}

if (length(all_venn_plots) > 0) {
  # Dynamically build the combined plot to avoid "subscript out of bounds"
  combined_figure <- wrap_plots(all_venn_plots) +
    plot_layout(ncol = min(3, length(all_venn_plots))) +
    plot_annotation(
      title = "Conserved Moult Toolkit Across Three Isopod Species",
      subtitle = paste("Thresholds: padj <", PADJ_CUTOFF, ", |log2FC| >", LFC_CUTOFF),
      theme = theme(
        plot.title = element_text(size = 36, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 24, hjust = 0.5, margin = margin(b=10))
      )
    )
  
  output_img <- file.path(out_dir_path, "combined_venn_moult_core.png")
  ggsave(output_img, combined_figure, width = 10 * length(all_venn_plots), height = 9, bg = "white", dpi = 300)
  cat("\nSUCCESS: Saved Venn diagram with italics and high-contrast colors.\n")
}

# === 3B. HEATMAP OF CONSERVED CORE ===
cat("\nGenerating Conserved Core Heatmap (Z-scores)...\n")
library(ComplexHeatmap)
library(circlize)

# Gather unique OGs from the master core list
core_ogs_all <- bind_rows(master_core_list) %>% pull(ODB_OG) %>% unique()

# Helper to get averaged Z-scores for a species
get_species_zscores <- function(species_id) {
  cat("  Processing:", species_id, "...\n")
  
  # Load DDS for normalization
  dds_path <- get_data_path(species_id, "dge_analysis", "dds_group_design.rds")
  if (!file.exists(dds_path)) {
    cat("    Warning: DDS not found for", species_id, "\n")
    return(NULL)
  }
  dds <- readRDS(dds_path)
  
  # Normalize metadata and EXCLUDE h.legs
  meta <- as.data.frame(colData(dds)) %>%
    mutate(sample_id = rownames(.)) %>%
    filter(tissue != "h.legs") %>%
    dplyr::select(sample_id, phase)
  
  if (nrow(meta) == 0) {
    cat("    Warning: No samples left after excluding h.legs for", species_id, "\n")
    return(NULL)
  }

  # VST Transform
  vsd <- vst(dds, blind=FALSE)
  counts <- as.data.frame(assay(vsd))
  # Subset counts to keep only the non-h.legs samples
  counts <- counts[, meta$sample_id, drop=FALSE]
  counts$gene_id <- rownames(counts)
  
  # Load mapping for this species
  ann <- load_annotations(species_id) %>%
    dplyr::select(gene_id, ODB_OG) %>%
    filter(ODB_OG %in% core_ogs_all)
  
  # Filter counts for core OGs and aggregate by ODB_OG
  sp_core_counts <- counts %>%
    inner_join(ann, by = "gene_id") %>%
    dplyr::select(-gene_id) %>%
    group_by(ODB_OG) %>%
    summarise(across(everything(), mean), .groups = "drop")
  
  # Pivot to long, join metadata, and average by phase (ignoring tissue for clarity)
  sp_long <- sp_core_counts %>%
    pivot_longer(cols = -ODB_OG, names_to = "sample_id", values_to = "vst") %>%
    left_join(meta, by = "sample_id") %>%
    group_by(ODB_OG, phase) %>%
    summarise(vst = mean(vst), .groups = "drop") %>%
    # Scale WITHIN SPECIES to highlight conserved waves
    group_by(ODB_OG) %>%
    mutate(vst_scaled = (vst - mean(vst)) / sd(vst)) %>%
    ungroup() %>%
    mutate(species = species_id)
  
  return(sp_long)
}

# Combine all species data
all_sp_data <- bind_rows(lapply(all_species, get_species_zscores))

# --- DATA RESHAPING FOR DIFFERENT VIEWS ---

# 1. Species-Split (Grouped by Species)
heatmap_df_split <- all_sp_data %>%
  mutate(phase = factor(phase, levels = c("intermoult", "premoult", "intramoult", "postmoult"))) %>%
  arrange(species, phase) %>%
  mutate(col_name = paste0(species, "_", phase)) %>%
  select(ODB_OG, col_name, vst_scaled) %>%
  pivot_wider(names_from = col_name, values_from = vst_scaled) %>%
  column_to_rownames("ODB_OG")

# 2. Synchrony Grid (Interleaved by Phase)
heatmap_df_sync <- all_sp_data %>%
  mutate(phase = factor(phase, levels = c("intermoult", "premoult", "intramoult", "postmoult"))) %>%
  arrange(phase, species) %>%
  mutate(col_name = paste0(species, "_", phase)) %>%
  select(ODB_OG, col_name, vst_scaled) %>%
  pivot_wider(names_from = col_name, values_from = vst_scaled) %>%
  column_to_rownames("ODB_OG")

# 3. Consensus Wave (Averaged across Species)
heatmap_df_consensus <- all_sp_data %>%
  group_by(ODB_OG, phase) %>%
  summarise(vst_consensus = mean(vst_scaled, na.rm = TRUE), .groups = "drop") %>%
  mutate(phase = factor(phase, levels = c("intermoult", "premoult", "intramoult", "postmoult"))) %>%
  arrange(phase) %>%
  pivot_wider(names_from = phase, values_from = vst_consensus) %>%
  column_to_rownames("ODB_OG")

# --- PLOTTING ---

# Global Row Labels Mapping
master_row_labels <- data.frame(ODB_OG = unique(c(rownames(heatmap_df_split), rownames(heatmap_df_consensus)))) %>%
  left_join(master_mapping, by = "ODB_OG") %>%
  mutate(label = ifelse(is.na(Preferred_name) | Preferred_name == "", ODB_OG, Preferred_name)) %>%
  select(ODB_OG, label)

get_label <- function(ogs) {
  master_row_labels$label[match(ogs, master_row_labels$ODB_OG)]
}

plot_and_save_heatmap <- function(mat, file_name, split_by = NULL, title = "") {
  # Filter for top 100 most variable OGs within this specific matrix
  top_vars <- head(order(matrixStats::rowVars(as.matrix(mat), na.rm = TRUE), decreasing = TRUE), 100)
  mat_sub <- mat[top_vars, ]
  labels_sub <- get_label(rownames(mat_sub))
  
  # Column Annotations logic
  if (any(grepl("_", colnames(mat_sub)))) {
    col_info <- data.frame(col_name = colnames(mat_sub)) %>%
      separate(col_name, into = c("species", "phase"), sep = "_", remove = FALSE)
    
    col_ann <- columnAnnotation(
      Species = sp_names[col_info$species],
      Phase = factor(col_info$phase, levels = c("intermoult", "premoult", "intramoult", "postmoult")),
      col = list(
        Species = setNames(RColorBrewer::brewer.pal(3, "Set1"), sp_names),
        Phase = c("intermoult"="#440154", "premoult"="#31688E", "intramoult"="#35B779", "postmoult"="#FDE725")
      ),
      show_legend = TRUE
    )
    col_split <- if(!is.null(split_by)) col_info[[split_by]] else NULL
  } else {
    col_ann <- columnAnnotation(
      Phase = factor(colnames(mat_sub), levels = c("intermoult", "premoult", "intramoult", "postmoult")),
      col = list(Phase = c("intermoult"="#440154", "premoult"="#31688E", "intramoult"="#35B779", "postmoult"="#FDE725"))
    )
    col_split <- NULL
  }

  p <- Heatmap(
    as.matrix(mat_sub),
    name = "Z-score",
    column_title = title,
    top_annotation = col_ann,
    row_labels = labels_sub,
    cluster_columns = FALSE, 
    column_split = col_split,
    row_names_gp = gpar(fontsize = 7),
    show_column_names = FALSE,
    heatmap_legend_param = list(title = "Z-score")
  )
  
  pdf(file.path(out_dir_path, file_name), width = 10, height = 12)
  draw(p)
  dev.off()
  cat("  Saved:", file_name, "\n")
}

# Save
plot_and_save_heatmap(heatmap_df_split, "conserved_core_heatmap_species_split.pdf", split_by = "species", title = "Species-Specific Moult Waves")
plot_and_save_heatmap(heatmap_df_sync, "conserved_core_heatmap_synchrony.pdf", title = "Universal Synchrony Grid")
plot_and_save_heatmap(heatmap_df_consensus, "conserved_core_heatmap_consensus.pdf", title = "Consensus Isopod Moult Wave")

# === 4. STATISTICAL VALIDATION ===
cat("\n--- PERFORMING STATISTICAL VALIDATION ---\n")
stats_log <- c("--- CROSS-SPECIES STATISTICAL VALIDATION ---", 
               paste("Date:", Sys.time()), "")

# A. Hypergeometric Test for Overlap Significance
universe_size <- nrow(master_mapping)
stats_log <- c(stats_log, "A. Overlap Significance (Hypergeometric Test)",
               paste("Universe Size (Shared OGs):", universe_size))

for (comp in target_comparisons) {
  phase_name <- str_extract(comp, "(?<=type1[bd]_)[a-z]+(?=_vs)")
  og_lists <- list()
  for (sp in all_species) {
    sig <- get_sig_genes(sp, comp)
    if (!is.null(sig)) og_lists[[sp]] <- sig
  }

  if (length(og_lists) == 3) {
    core_3way <- Reduce(intersect, og_lists)
    intersect_12 <- length(intersect(og_lists[[1]], og_lists[[2]]))
    p_val <- phyper(length(core_3way) - 1, length(og_lists[[3]]), universe_size - length(og_lists[[3]]), intersect_12, lower.tail = FALSE)
    
    res_line <- sprintf("  [%s] 3-way Overlap Significance: p = %.2e", str_to_title(phase_name), p_val)
    cat(res_line, "\n")
    stats_log <- c(stats_log, res_line)
  }
}

# B. Inter-Species Correlation (Wave Synchrony)
cat("\nCalculating Inter-Species Wave Correlation (R)...\n")
stats_log <- c(stats_log, "", "B. Inter-Species Wave Correlation (Pearson R)")

cor_data <- all_sp_data %>%
  select(ODB_OG, species, phase, vst_scaled) %>%
  pivot_wider(names_from = species, values_from = vst_scaled)

cor_ml <- cor(cor_data$maculatum, cor_data$laevis, use = "complete.obs")
cor_mo <- cor(cor_data$maculatum, cor_data$officinalis, use = "complete.obs")
cor_lo <- cor(cor_data$laevis, cor_data$officinalis, use = "complete.obs")
mean_cor <- mean(c(cor_ml, cor_mo, cor_lo))

res_lines <- c(
  sprintf("  A. maculatum  vs  P. laevis:      R = %.3f", cor_ml),
  sprintf("  A. maculatum  vs  A. officinalis: R = %.3f", cor_mo),
  sprintf("  P. laevis     vs  A. officinalis: R = %.3f", cor_lo),
  sprintf("  MEAN CROSS-SPECIES CORRELATION:   R = %.3f", mean_cor)
)

for(line in res_lines) cat(line, "\n")
stats_log <- c(stats_log, res_lines)

# Save to file
writeLines(stats_log, file.path(out_dir_path, "fig_3 statistical_validation_results.txt"))
cat("\nSUCCESS: Saved statistical results to output_csv/statistical_validation_results.txt\n")

cat("\nAnalysis complete! All metrics and heatmaps finalized.\n")

