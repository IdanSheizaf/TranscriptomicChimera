# === INDIVIDUAL SPECIES PCA ANALYSIS ===
# This script generates PC1/2, PC1/3, and PC2/3 plots for each species
# across multiple gene selection thresholds (500, 1000, 2000, 5000 genes).

# === SETUP WORKING DIRECTORY ===
project_path <- "G:/My Drive/PhD/Projects/Comparative Expressional Changes During the Moult Cycle in Land Isopods/local_scripts"
setwd(project_path)

# === LIBRARIES ===
source("utils.R")
library(DESeq2)
library(ggrepel)
library(vegan) # For PERMANOVA

# === SETUP ===
log_header("Individual Species PCA")

all_species <- c("maculatum", "laevis", "officinalis")
gene_thresholds <- c(500, 1000, 2000, 5000) 

# Standard Phase Colors (The "Clock" Colors)
moult_colors <- c(
  "intermoult" = "#440154", 
  "premoult"   = "#31688E", 
  "intramoult" = "#35B779", 
  "postmoult"  = "#FDE725"
)

# Standard Tissue Shapes
# thorax is set to shape 5 (diamond) as per user request
moult_shapes <- c(
  "f.legs" = 16, # Closed circle
  "h.legs" = 17, # Closed triangle
  "head"   = 15, # Closed diamond (filled)
  "thorax" = 18   # Open diamond
)

# Mapping species to full binomial names for publication
species_full_names <- list(
  "maculatum"   = "Armadillidium maculatum",
  "laevis"      = "Porcellio laevis",
  "officinalis" = "Armadillo officinalis"
)

# === PLOTTING FUNCTION ===
plot_pca_species <- function(pca_full, pca_var, sample_data, species_name, n_genes, pc_x = 1, pc_y = 2, permanova_res = NULL) {

  full_name <- species_full_names[[species_name]]

  plot_data <- data.frame(
    PCx = pca_full$x[, pc_x],
    PCy = pca_full$x[, pc_y],
    phase = factor(sample_data$phase, levels = names(moult_colors)),
    tissue = sample_data$tissue,
    sample = sample_data$sample
  )

  # Create subtitle with stats if provided
  sub_text <- paste0("Top ", n_genes, " genes")
  if (!is.null(permanova_res)) {
    sub_text <- paste0(sub_text, " | PERMANOVA (Phase): R² = ", round(permanova_res$R2[1], 2), 
                       ", p = ", format.pval(permanova_res$`Pr(>F)`[1], digits = 3))
  }

  p <- ggplot(plot_data, aes(x=PCx, y=PCy, color=phase, shape=tissue, label=sample)) +
    geom_point(size=5, alpha=0.8) +
    scale_color_manual(values = moult_colors) +
    scale_shape_manual(values = moult_shapes) +
    labs(
      title = bquote(italic(.(full_name)) ~ ": PC" * .(pc_x) ~ "vs PC" * .(pc_y)),
      subtitle = sub_text,
      x = paste0("PC", pc_x, ": ", pca_var[pc_x], "% variance"),
      y = paste0("PC", pc_y, ": ", pca_var[pc_y], "% variance"),
      color = "Phase", shape = "Tissue"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size=16, face="bold", hjust=0.5),
      plot.subtitle = element_text(size=10, hjust=0.5, color="grey30"),
      legend.position = "right"
    )

  return(p)
}


# === PROCESSING LOOP ===

for (sp in all_species) {
  cat("\nProcessing PCA and Stats for:", sp, "...\n")
  
  # 1. Load data
  dds_path <- get_data_path(sp, "dge_analysis", "dds_phase_design.rds")
  if (!file.exists(dds_path)) {
    cat("  Note: Phase design RDS not found, falling back to group design...\n")
    dds_path <- get_data_path(sp, "dge_analysis", "dds_group_design.rds")
  }
  
  if (!file.exists(dds_path)) {
    cat("  SKIPPING: No RDS found for", sp, "\n")
    next
  }
  
  dds <- readRDS(dds_path)
  
  # 2. VST Transformation
  cat("  Performing VST...\n")
  vsd <- vst(dds, blind=FALSE)
  sample_data <- as.data.frame(colData(vsd))
  
  # --- STATISTICAL VALIDATION (PERMANOVA) ---
  cat("  Running PERMANOVA (Phase effect)... ")
  # Use all genes for the most robust statistical test of the global state
  perm_res <- adonis2(t(assay(vsd)) ~ phase, data = sample_data, permutations = 999)
  cat("p =", perm_res$`Pr(>F)`[1], "\n")
  
  # Save Stats to text file
  stat_file <- file.path(get_project_root(), "PCA", sp, "permanova_results.txt")
  if (!dir.exists(dirname(stat_file))) dir.create(dirname(stat_file), recursive = TRUE)
  writeLines(capture.output(print(perm_res)), stat_file)

  # Pre-calculate variances once per species for plotting
  vars <- rowVars(assay(vsd))
  
  # 3. Loop through different gene counts for plotting
  for (n_genes in gene_thresholds) {
    
    # Ensure we don't exceed actual gene count
    actual_n <- min(n_genes, length(vars))
    cat("  Analyzing Top", actual_n, "genes... ")
    
    topVarGenes <- head(order(vars, decreasing = TRUE), actual_n)
    pca_data <- assay(vsd)[topVarGenes, ]
    
    pca_full <- prcomp(t(pca_data))
    pca_full$x[, 1] <- -pca_full$x[, 1] # Flip PC1
    # pca_full$x[, 2] <- -pca_full$x[, 2] # Flip PC2
    # pca_full$x[, 3] <- -pca_full$x[, 3] # Flip PC3
    pca_var  <- round(100 * summary(pca_full)$importance[2,], 1)
    cat("PC1 Var:", pca_var[1], "%\n")
    
    # 4. Create Specific Output Directory
    out_dir_path <- normalizePath(file.path(get_project_root(), "PCA", sp, paste0("top_", actual_n)), mustWork = FALSE)
    if (!dir.exists(out_dir_path)) dir.create(out_dir_path, recursive = TRUE, showWarnings = FALSE)
    
    # 5. Save Plots
    ggsave(file.path(out_dir_path, paste0(sp, "_PCA_PC1v2_top_reversed", actual_n, ".png")), 
           plot_pca_species(pca_full, pca_var, sample_data, sp, actual_n, 1, 2, perm_res), width = 8, height = 6)
    
    ggsave(file.path(out_dir_path, paste0(sp, "_PCA_PC1v3_top_reversed", actual_n, ".png")), 
           plot_pca_species(pca_full, pca_var, sample_data, sp, actual_n, 1, 3, perm_res), width = 8, height = 6)
    
    ggsave(file.path(out_dir_path, paste0(sp, "_PCA_PC2v3_top_reversed", actual_n, ".png")), 
           plot_pca_species(pca_full, pca_var, sample_data, sp, actual_n, 2, 3, perm_res), width = 8, height = 6)
  }
}

cat("\nAll species PCA plots generated with PERMANOVA stats!\n")
