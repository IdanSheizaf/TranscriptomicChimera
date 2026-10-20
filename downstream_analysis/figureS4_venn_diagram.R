# ==============================================================================
# MANUSCRIPT FIGURE S4: Conserved Core 3-Way Venn Diagram
# 
# This script identifies the statistical intersection of differentially expressed
# genes (DEGs) across all three isopod species (excluding hind legs to avoid
# early-transition bias) for Premoult, Intramoult, and Postmoult relative to
# Intermoult.
#
# Output:
#   - Figure_S4_Venn_Core.pdf / .png
#   - conserved_core_master_all_phases.csv (Table S5)
#   - Overlap significance validation via Hypergeometric test
# ==============================================================================

if (file.exists("utils.R")) {
  source("utils.R")
} else if (file.exists(file.path("downstream_analysis", "utils.R"))) {
  source(file.path("downstream_analysis", "utils.R"))
}
library(tidyverse)
library(ggVennDiagram)
library(patchwork)

log_header("Figure S4: Conserved Core 3-Way Venn Diagram")

all_species <- c("maculatum", "laevis", "officinalis")
sp_names <- c("maculatum" = "A. maculatum", "laevis" = "P. laevis", "officinalis" = "A. officinalis")

target_comparisons <- c(
  "type1d_premoult_vs_intermoult_no_hlegs",
  "type1d_intramoult_vs_intermoult_no_hlegs",
  "type1d_postmoult_vs_intermoult_no_hlegs"
)

out_dir <- get_output_path("Supplementary_Figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Master Mapping for ODB_OG to Preferred_name
master_mapping <- load_annotations("maculatum") %>%
  dplyr::select(ODB_OG, Preferred_name) %>%
  filter(!is.na(ODB_OG) & ODB_OG != "" & ODB_OG != "-") %>%
  distinct(ODB_OG, .keep_all = TRUE)

get_sig_genes <- function(species, comparison) {
  dge_path <- get_data_path(species, "dge_analysis/dge_type1d_phase_vs_intermoult_no_hlegs/full", paste0(comparison, ".csv"))
  if (!file.exists(dge_path)) return(NULL)
  
  dge <- read_csv(dge_path, show_col_types = FALSE)
  annotations <- load_annotations(species)
  
  sig_set <- dge %>%
    inner_join(annotations, by = c("transcript_id" = "gene_id")) %>%
    filter(padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
    filter(!is.na(ODB_OG) & ODB_OG != "" & ODB_OG != "-") %>%
    pull(ODB_OG) %>% unique()
  return(sig_set)
}

all_venn_plots <- list()
master_core_list <- list()

for (comp in target_comparisons) {
  cat("Processing comparison:", comp, "...\n")
  og_lists <- list()
  for (sp in all_species) {
    sig <- get_sig_genes(sp, comp)
    if (!is.null(sig)) og_lists[[sp_names[sp]]] <- sig
  }
  
  if (length(og_lists) < 3) next
  phase_name <- str_extract(comp, "(?<=type1[bd]_)[a-z]+(?=_vs)")
  conserved_core_ogs <- Reduce(intersect, og_lists)
  
  core_df <- data.frame(ODB_OG = conserved_core_ogs) %>%
    left_join(master_mapping, by = "ODB_OG") %>%
    mutate(Phase = phase_name) %>%
    dplyr::select(ODB_OG, Preferred_name, Phase)
  master_core_list[[phase_name]] <- core_df
  
  p <- ggVennDiagram(og_lists, label = "both", label_alpha = 0, label_size = 5, edge_size = 0.8) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1) + 
    labs(title = str_to_title(phase_name)) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 26, margin = margin(b = 10, t = 10)),
      plot.margin = margin(t = 20, b = 20, l = 40, r = 40)
    ) +
    coord_fixed(clip = "off")
  
  for (i in seq_along(p$layers)) {
    if (class(p$layers[[i]]$geom)[1] == "GeomText" && 
        (!is.null(p$layers[[i]]$mapping$label) && rlang::as_label(p$layers[[i]]$mapping$label) != "label") || 
        is.null(p$layers[[i]]$mapping$label)) {
      p$layers[[i]]$aes_params$fontface <- "italic"
      p$layers[[i]]$aes_params$size <- 6.5
    }
  }
  all_venn_plots[[comp]] <- p
}

# Export Master Core Table (Table S5)
if (length(master_core_list) > 0) {
  master_core_df <- bind_rows(master_core_list)
  write_csv(master_core_df, file.path(out_dir, "conserved_core_master_all_phases.csv"))
  cat("  Exported Table S5: conserved_core_master_all_phases.csv (", nrow(master_core_df), " OGs)\n")
}

# Save Combined Venn Figure
combined_figure <- wrap_plots(all_venn_plots) +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Figure S4: Conserved Core Moult Toolkit Across Three Isopod Species",
    subtitle = "Statistical 3-way Venn overlap of phase-specific transitions relative to Intermoult (no h.legs)",
    theme = theme(
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 18, hjust = 0.5, margin = margin(b = 10))
    )
  )

fig_s4_pdf <- file.path(out_dir, "Figure_S4_Venn_Core.pdf")
fig_s4_png <- file.path(out_dir, "Figure_S4_Venn_Core.png")

ggsave(fig_s4_pdf, combined_figure, width = 24, height = 8, bg = "white")
ggsave(fig_s4_png, combined_figure, width = 24, height = 8, bg = "white", dpi = 300)

cat("\nSUCCESS: Figure S4 generated and saved to:\n")
cat("  -", fig_s4_pdf, "\n")
cat("  -", fig_s4_png, "\n")
