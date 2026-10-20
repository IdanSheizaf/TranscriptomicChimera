# ==============================================================================
# MANUSCRIPT FIGURE 3: Surge in Appendicular Transcriptomic Asymmetry
# 
# This script calculates the degree of appendicular transcriptomic asymmetry 
# (Front legs vs Back legs) across the four moult cycle phases, performing 
# Wilcoxon rank sum testing to quantify the statistically significant ~34-fold 
# surge during the Intramoult interval (The Transcriptomic Chimera).
#
# Output:
#   - Figure_3_Asymmetry_Surge.png / .pdf
#   - intramoult_Asymmetry_Statistics.txt
# ==============================================================================

if (file.exists("utils.R")) {
  source("utils.R")
} else if (file.exists(file.path("downstream_analysis", "utils.R"))) {
  source(file.path("downstream_analysis", "utils.R"))
}

library(tidyverse)
library(patchwork)

log_header("Figure 3: Appendicular Asymmetry Surge")

all_species <- c("maculatum", "laevis", "officinalis")
phases <- c("intermoult", "premoult", "intramoult", "postmoult")
out_dir <- get_output_path()
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 1. DATA COLLECTION ---
cat("Collecting Front vs Hind Leg DEG counts across species and phases...\n")
deg_counts_list <- list()

for (sp in all_species) {
  for (ph in phases) {
    file_name <- paste0("type4_", ph, "_f.legs_vs_h.legs.csv")
    dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", file_name)
    
    if (!file.exists(dge_path)) {
      # Check reverse naming
      file_name <- paste0("type4_", ph, "_h.legs_vs_f.legs.csv")
      dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", file_name)
    }
    
    if (file.exists(dge_path)) {
      dge <- read_csv(dge_path, show_col_types = FALSE)
      sig_count <- dge %>%
        dplyr::filter(padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
        nrow()
      deg_counts_list[[paste0(sp, "_", ph)]] <- data.frame(species = sp, phase = ph, deg_count = sig_count)
      cat("  [", sp, "] Phase:", ph, "->", sig_count, "DEGs\n")
    }
  }
}

deg_summary <- bind_rows(deg_counts_list)
deg_summary$phase <- factor(deg_summary$phase, levels = phases)
sp_names <- c("maculatum" = "Armadillidium maculatum", "laevis" = "Porcellio laevis", "officinalis" = "Armadillo officinalis")
deg_summary$species_name <- factor(sp_names[deg_summary$species], levels = sp_names)

# --- 2. STATISTICAL VALIDATION (The Intramoult Surge) ---
cat("\nCalculating Wilcoxon significance of the Intramoult Surge...\n")
intramoult_vals <- deg_summary %>% dplyr::filter(phase == "intramoult") %>% pull(deg_count)
baseline_vals   <- deg_summary %>% dplyr::filter(phase == "intermoult") %>% pull(deg_count)
other_vals      <- deg_summary %>% dplyr::filter(phase != "intramoult") %>% pull(deg_count)

w_res         <- wilcox.test(intramoult_vals, other_vals)
mean_intra    <- mean(intramoult_vals)
mean_baseline <- mean(baseline_vals)
fold_increase <- mean_intra / max(1, mean_baseline)

cat(sprintf("  Intramoult Mean DEGs: %.1f\n", mean_intra))
cat(sprintf("  Intermoult Mean DEGs: %.1f\n", mean_baseline))
cat(sprintf("  Fold Increase: %.1f-fold (Wilcoxon p = %.6f)\n\n", fold_increase, w_res$p.value))

stat_file <- file.path(out_dir, "intramoult_Asymmetry_Statistics.txt")
writeLines(
  c(
    "=== APPENDICULAR ASYMMETRY SURGE STATISTICS ===",
    paste("Date:", Sys.time()),
    "",
    sprintf("Wilcoxon Rank Sum Test (Intramoult vs Other Phases): p-value = %.6f", w_res$p.value),
    sprintf("Mean DEGs Intramoult: %.1f", mean_intra),
    sprintf("Mean DEGs Intermoult: %.1f", mean_baseline),
    sprintf("Fold Increase: %.2f-fold", fold_increase)
  ),
  stat_file
)

# --- 3. PUBLICATION VISUALIZATION ---
sig_label <- ifelse(w_res$p.value < 0.001, "p < 0.001", 
                    ifelse(w_res$p.value < 0.01, "p < 0.01", 
                           ifelse(w_res$p.value < 0.05, "p < 0.05", "ns")))

# Custom palette for species
species_colors <- c(
  "Armadillidium maculatum" = "#E41A1C", 
  "Porcellio laevis"        = "#377EB8", 
  "Armadillo officinalis"   = "#4DAF4A"
)

p_asym <- ggplot(deg_summary, aes(x = phase, y = deg_count, fill = species_name)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7, linewidth = 0.8) +
  scale_fill_manual(values = species_colors, name = "Species") +
  labs(
    title = "Figure 3: Rhythmic Surge in Transcriptomic Asymmetry During the Intramoult Phase",
    subtitle = paste0("Appendicular divergence (Front vs Hind Legs): ~", round(fold_increase), "-fold surge during intramoult (Wilcoxon, ", sig_label, ")"),
    x = "Moult Phase",
    y = "Significant DEGs (|log2FC| > 1, padj < 0.05)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.text = element_text(face = "italic")
  )

fig3_png <- file.path(out_dir, "Figure_3_Asymmetry_Surge.png")
fig3_pdf <- file.path(out_dir, "Figure_3_Asymmetry_Surge.pdf")

ggsave(fig3_png, p_asym, width = 11, height = 7, dpi = 300)
ggsave(fig3_pdf, p_asym, width = 11, height = 7)

cat("SUCCESS: Figure 3 generated and saved to:\n")
cat("  -", fig3_png, "\n")
cat("  -", fig3_pdf, "\n")

# --- 4. EXPORT TABLE S1: CONSERVED INTRAMOULT POLARITY LIST ---
cat("\nAnalyzing Conserved Functional Polarity for Intramoult (Table S1)...\n")
intramoult_genes <- list()
for (sp in all_species) {
  file_name <- "type4_intramoult_f.legs_vs_h.legs.csv"
  dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", file_name)
  
  if (!file.exists(dge_path)) {
    file_name <- "type4_intramoult_h.legs_vs_f.legs.csv"
    dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", file_name)
  }
  
  if (file.exists(dge_path)) {
    ann <- load_annotations(sp)
    intramoult_genes[[sp]] <- read_csv(dge_path, show_col_types = FALSE) %>%
      dplyr::inner_join(ann, by = c("transcript_id" = "gene_id")) %>%
      dplyr::filter(padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
      dplyr::select(Preferred_name, ODB_OG, log2FoldChange) %>%
      dplyr::distinct()
  }
}

if (length(intramoult_genes) == 3) {
  og_list_for_int <- lapply(intramoult_genes, function(df) df$ODB_OG)
  common_3way_ogs <- Reduce(intersect, og_list_for_int)
  
  df_polarity <- bind_rows(intramoult_genes) %>%
    dplyr::filter(ODB_OG %in% common_3way_ogs) %>%
    dplyr::group_by(ODB_OG) %>%
    dplyr::summarise(
      Preferred_name = dplyr::first(na.omit(Preferred_name[Preferred_name != "" & Preferred_name != "-"])),
      avg_lfc = mean(log2FoldChange),
      consistent_direction = all(log2FoldChange > 0) | all(log2FoldChange < 0),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Preferred_name = ifelse(is.na(Preferred_name), paste0("Uncharacterized_", ODB_OG), Preferred_name)) %>%
    dplyr::mutate(Bias = ifelse(avg_lfc > 0, "Anterior (Front)", "Posterior (Back)")) %>%
    dplyr::arrange(desc(abs(avg_lfc)))
  
  table_s1_file <- file.path(out_dir, "Conserved_Asymmetry_Polarity_List.csv")
  write_csv(df_polarity, table_s1_file)
  cat("  - Exported Table S1:", table_s1_file, "(", nrow(df_polarity), " OGs)\n")
}

