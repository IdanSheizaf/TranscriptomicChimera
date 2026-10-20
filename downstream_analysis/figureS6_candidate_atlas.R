# ==============================================================================
# MANUSCRIPT FIGURE S6: Conserved Candidate Gene Atlas
# 
# This script generates standardized expression trajectories for the top 20
# conserved moult toolkit candidate orthologs across all three terrestrial isopod
# species (Armadillidium maculatum, Porcellio laevis, Armadillidium officinalis).
#
# Output:
#   - Figure_S6_Candidate_Gene_Atlas.pdf
#   - Figure_S6_Candidate_Gene_Atlas.png
#   - Universal_Moult_Toolkit_Candidates.csv (Table S10 source)
# ==============================================================================

if (file.exists("utils.R")) {
  source("utils.R")
} else if (file.exists(file.path("downstream_analysis", "utils.R"))) {
  source(file.path("downstream_analysis", "utils.R"))
}
library(tidyverse)
library(patchwork)
library(DESeq2)

log_header("Figure S6: Candidate Gene Atlas (Supplementary)")

out_dir <- get_output_path("Supplementary_Figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- 1. LOAD CANDIDATE LIST (Table S10 Source) ---
candidate_candidates <- c(
  file.path(get_output_path(), "Universal_Moult_Toolkit_Candidates.csv"),
  file.path(get_project_root(), "local_scripts", "output_csv", "Universal_Moult_Toolkit_Candidates.csv"),
  file.path(dirname(get_project_root()), "local_scripts", "output_csv", "Universal_Moult_Toolkit_Candidates.csv"),
  file.path(get_project_root(), "figures_final", "Supplementary", "Table S10 Universal_Moult_Toolkit_Candidates.csv"),
  file.path(dirname(get_project_root()), "figures_final", "Supplementary", "Table S10 Universal_Moult_Toolkit_Candidates.csv")
)
candidate_path <- candidate_candidates[file.exists(candidate_candidates)][1]

if (is.na(candidate_path) || !file.exists(candidate_path)) {
  stop("Universal_Moult_Toolkit_Candidates.csv not found! Please run universal_moult_toolkit_synthesis.R first.")
}


# Define Publication Order (Top 20 candidates)
atlas_order <- c(
  "Cuticle protein",
  "cuticle protein AM1274-like",
  "cuticle protein CP14.6-like",
  "cuticle protein CP1499-like",
  "cuticle protein AM1159-like",
  "cuticle protein AMP1A-like",
  "cuticle protein 7-like",
  "collagen alpha-1(XI) chain-like",
  "ice-structuring glycoprotein-like",
  "U-scoloptoxin(01)-Cw1a-like",
  "U-scoloptoxin(01)-Cw1a",
  "C-type lectin BfL-1-like",
  "venom phosphodiesterase 2-like",
  "mucin-2-like",
  "5'-nucleotidase",
  "Glutaredoxin domain-containing cysteine-rich protein",
  "oleosin-B6-like",
  "formin-like protein 20",
  "pneumococcal serine-rich repeat protein-like isoform X1",
  "methyl farnesoate epoxidase-like"
)

# --- 2. LOAD COUNTS FOR ALL SPECIES ---
cat("Loading normalized counts for all species...\n")

get_species_counts <- function(species_id) {
  cat("  Processing:", species_id, "...\n")
  
  dds_path <- get_data_path(species_id, "dge_analysis", "dds_group_design.rds")
  dds <- readRDS(dds_path)
  
  meta <- as.data.frame(colData(dds)) %>%
    mutate(sample_id = rownames(.)) %>%
    dplyr::select(sample_id, tissue, phase) %>%
    mutate(species = case_when(
      species_id == "maculatum" ~ "A. maculatum",
      species_id == "laevis" ~ "P. laevis",
      species_id == "officinalis" ~ "A. officinalis",
      TRUE ~ species_id
    ))
  
  vsd <- vst(dds, blind = FALSE)
  counts <- as.data.frame(assay(vsd))
  sample_names <- colnames(counts)
  counts$gene_id <- rownames(counts)
  
  ann <- load_annotations(species_id) %>%
    dplyr::select(gene_id, Preferred_name) %>%
    filter(!is.na(Preferred_name) & Preferred_name != "" & Preferred_name != "-")
  
  counts_mapped <- counts %>%
    inner_join(ann, by = "gene_id") %>%
    dplyr::select(-gene_id) %>%
    pivot_longer(cols = all_of(sample_names), names_to = "sample_id", values_to = "vst_count") %>%
    left_join(meta, by = "sample_id")
  
  return(counts_mapped)
}

all_counts <- bind_rows(lapply(c("maculatum", "laevis", "officinalis"), get_species_counts))

# --- 3. PLOTTING ENGINE ---
plot_candidate_gene <- function(gene_name) {
  gene_data <- all_counts %>%
    dplyr::filter(Preferred_name == gene_name)
  
  if (nrow(gene_data) == 0) return(NULL)
  
  # Format Phase Order
  gene_data <- gene_data %>%
    mutate(phase = factor(tolower(phase), 
                          levels = c("intermoult", "premoult", "intramoult", "postmoult"),
                          labels = c("Intermoult", "Premoult", "Intramoult", "Postmoult")))
  
  # Scale Z-score by species
  gene_data <- gene_data %>%
    group_by(species) %>%
    mutate(vst_scaled = (vst_count - mean(vst_count)) / sd(vst_count)) %>%
    ungroup()
  
  gene_summary <- gene_data %>%
    group_by(phase, species) %>%
    summarise(mean_scaled = mean(vst_scaled), .groups = "drop")
  
  sp_palette <- c("A. maculatum" = "#E41A1C", "P. laevis" = "#377EB8", "A. officinalis" = "#4DAF4A")
  
  p <- ggplot() +
    geom_line(data = gene_summary, aes(x = phase, y = mean_scaled, color = species, group = species), 
              linewidth = 1.0, alpha = 0.8) +
    geom_point(data = gene_summary, aes(x = phase, y = mean_scaled, color = species), size = 2.0) +
    theme_bw() +
    scale_color_manual(values = sp_palette) +
    labs(title = gene_name, x = NULL, y = "Rel. Exp (Z-score)", color = "Species") +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9),
      axis.title.y = element_text(size = 9),
      legend.position = "none",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# --- 4. ASSEMBLE COMBINED ATLAS ---
cat("Generating 20-Gene Candidate Atlas Grid...\n")
atlas_plots <- list()
for (gene in atlas_order) {
  p <- plot_candidate_gene(gene)
  if (!is.null(p)) atlas_plots[[gene]] <- p
}

subtitle_text <- "Cross-species relative expression (Z-score) across the biphasic moult cycle (averaged across tissues)"

combined_atlas <- wrap_plots(atlas_plots, ncol = 4, guides = "collect") +
  plot_annotation(
    title = "Figure S6: The Universal Isopod Moult Toolkit Candidate Atlas",
    subtitle = subtitle_text,
    theme = theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  ) & theme(legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 13, face = "bold"))

fig_s6_pdf <- file.path(out_dir, "Figure_S6_Candidate_Gene_Atlas.pdf")
fig_s6_png <- file.path(out_dir, "Figure_S6_Candidate_Gene_Atlas.png")

ggsave(fig_s6_pdf, combined_atlas, width = 16, height = 20)
ggsave(fig_s6_png, combined_atlas, width = 16, height = 20, dpi = 300)

cat("\nSUCCESS: Figure S6 generated and saved to:\n")
cat("  -", fig_s6_pdf, "\n")
cat("  -", fig_s6_png, "\n")

# --- 4. EXPORT TABLE S10: CANDIDATE TOOLKIT TABLE ---
candidates_df <- read_csv(candidate_path, show_col_types = FALSE)
table_s10_out <- file.path(out_dir, "Universal_Moult_Toolkit_Candidates.csv")
write_csv(candidates_df, table_s10_out)

csv_out <- get_output_path()
if (dir.exists(csv_out)) {
  write_csv(candidates_df, file.path(csv_out, "Universal_Moult_Toolkit_Candidates.csv"))
}
cat("  - Exported Table S10:", table_s10_out, "(", nrow(candidates_df), " candidates)\n")

