# ==============================================================================
# MANUSCRIPT FIGURE 2: Species-Specific Transcriptomic Clocks
# 
# This script performs Principal Component Analysis (PCA) on VST-normalized 
# counts (top 2,000 variable genes) for each species independently, running 
# PERMANOVA tests to validate the cyclical "transcriptomic clock" progression 
# and demonstrating the anatomical decoupling of intramoult hind legs.
#
# Output:
#   - Figure_2_Species_Clocks.png / .pdf (Composite 3-panel figure)
#   - Individual species PCA plots in PCA/<species>/
# ==============================================================================

if (file.exists("utils.R")) {
  source("utils.R")
} else if (file.exists(file.path("downstream_analysis", "utils.R"))) {
  source(file.path("downstream_analysis", "utils.R"))
}

library(DESeq2)
library(tidyverse)
library(patchwork)
library(vegan) # For PERMANOVA

log_header("Figure 2: Species-Specific Transcriptomic Clocks")

all_species <- c("maculatum", "laevis", "officinalis")
n_genes <- 2000

# Mapping species to full binomial names for publication
species_full_names <- list(
  "maculatum"   = "Armadillidium maculatum",
  "laevis"      = "Porcellio laevis",
  "officinalis" = "Armadillo officinalis"
)

# Canonical PCs selected for each species clock
species_pcs <- list(
  "maculatum"   = c(1, 2),
  "laevis"      = c(1, 3),
  "officinalis" = c(1, 3)
)

out_base <- file.path(get_project_root(), "PCA")
if (!dir.exists(out_base)) dir.create(out_base, recursive = TRUE)

plot_list <- list()

for (sp in all_species) {
  cat("\nProcessing PCA for:", sp, "...\n")
  
  dds_path <- get_data_path(sp, "dge_analysis", "dds_phase_design.rds")
  if (!file.exists(dds_path)) {
    dds_path <- get_data_path(sp, "dge_analysis", "dds_group_design.rds")
  }
  
  if (!file.exists(dds_path)) {
    cat("  WARNING: No RDS found for", sp, "\n")
    next
  }
  
  dds <- readRDS(dds_path)
  vsd <- vst(dds, blind = FALSE)
  sample_data <- as.data.frame(colData(vsd))
  
  # PERMANOVA (Phase effect across all genes)
  perm_res <- adonis2(t(assay(vsd)) ~ phase, data = sample_data, permutations = 999)
  r2_val <- round(perm_res$R2[1], 2)
  p_val  <- format.pval(perm_res$`Pr(>F)`[1], digits = 3)
  cat("  PERMANOVA Phase effect: R2 =", r2_val, ", p =", p_val, "\n")
  
  # Select Top 2000 high-variance genes
  vars <- matrixStats::rowVars(assay(vsd))
  actual_n <- min(n_genes, length(vars))
  top_genes <- head(order(vars, decreasing = TRUE), actual_n)
  
  pca_res <- prcomp(t(assay(vsd)[top_genes, ]))
  pca_res$x[, 1] <- -pca_res$x[, 1] # Flip PC1 for standard circular orientation
  pca_var <- round(100 * summary(pca_res)$importance[2, ], 1)
  
  pc_x <- species_pcs[[sp]][1]
  pc_y <- species_pcs[[sp]][2]
  
  plot_df <- data.frame(
    PCx    = pca_res$x[, pc_x],
    PCy    = pca_res$x[, pc_y],
    phase  = factor(sample_data$phase, levels = names(moult_colors)),
    tissue = sample_data$tissue,
    sample = sample_data$sample
  )
  
  full_name <- species_full_names[[sp]]
  sub_text  <- paste0("Top ", actual_n, " genes | R² = ", r2_val, ", p = ", p_val)
  
  p_sp <- ggplot(plot_df, aes(x = PCx, y = PCy, color = phase, shape = tissue)) +
    geom_point(size = 5, alpha = 0.85, stroke = 1.2) +
    scale_color_manual(values = moult_colors, name = "Phase") +
    scale_shape_manual(values = moult_shapes, name = "Body Part") +
    labs(
      title = bquote(bolditalic(.(full_name)) ~ bold(": PC" * .(pc_x) ~ "vs PC" * .(pc_y))),
      subtitle = sub_text,
      x = paste0("PC", pc_x, ": ", pca_var[pc_x], "% variance"),
      y = paste0("PC", pc_y, ": ", pca_var[pc_y], "% variance")
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
      legend.position = "right"
    )
  
  plot_list[[sp]] <- p_sp
  
  # Save individual species plot
  sp_dir <- file.path(out_base, sp)
  if (!dir.exists(sp_dir)) dir.create(sp_dir, recursive = TRUE)
  ggsave(file.path(sp_dir, paste0(sp, "_Figure2_clock.png")), p_sp, width = 8, height = 6, dpi = 300)
}

# Combine into composite Figure 2
if (length(plot_list) == 3) {
  figure2_composite <- (plot_list[["maculatum"]] + plot_list[["laevis"]] + plot_list[["officinalis"]]) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = "Figure 2: Species-Specific Transcriptomic Clocks Across Three Terrestrial Isopod Families",
      subtitle = "Continuous circular progression and intramoult hind-leg anatomical decoupling",
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey30")
      )
    )
  
  fig2_png <- file.path(out_base, "Figure_2_Species_Clocks.png")
  fig2_pdf <- file.path(out_base, "Figure_2_Species_Clocks.pdf")
  
  ggsave(fig2_png, figure2_composite, width = 22, height = 7, dpi = 300)
  ggsave(fig2_pdf, figure2_composite, width = 22, height = 7)
  cat("\nSUCCESS: Composite Figure 2 saved to:\n  -", fig2_png, "\n  -", fig2_pdf, "\n")
}
