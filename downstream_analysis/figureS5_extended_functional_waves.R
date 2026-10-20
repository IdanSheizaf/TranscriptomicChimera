# ==============================================================================
# MANUSCRIPT FIGURE S5: Extended Functional Asymmetry & Cephalic Regulation
# 
# This script generates the supplementary functional enrichment figures:
#   1. Hind Legs vs Front Legs: Molecular Function (MF) (Figure S5A)
#   2. Head vs Front Legs: Biological Process (BP) (Figure S5B)
#   3. Head vs Front Legs: Molecular Function (MF) (Figure S5C)
#
# Output:
#   - Figure_S5A_HLegs_vs_FLegs_MF.pdf
#   - Figure_S5B_Head_vs_FLegs_BP.pdf
#   - Figure_S5C_Head_vs_FLegs_MF.pdf
#   - Tables S7, S8, S9 consensus CSV tables
# ==============================================================================

if (file.exists("utils.R")) {
  source("utils.R")
} else if (file.exists(file.path("downstream_analysis", "utils.R"))) {
  source(file.path("downstream_analysis", "utils.R"))
}
library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(GO.db)

log_header("Figure S5: Extended Functional Waves (Supplementary)")

all_species <- c("maculatum", "laevis", "officinalis")
out_dir <- get_output_path("Supplementary_Figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("Pre-loading species annotations...\n")
annotations_list <- setNames(lapply(all_species, load_annotations), all_species)

# Shared Enrichment Function
get_enrichment_from_csv <- function(species, dge_path, go_type = "mf", ann = NULL) {
  if (!file.exists(dge_path)) return(NULL)
  dge <- read_csv(dge_path, show_col_types = FALSE) %>% 
    filter(padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF)
  if (nrow(dge) < 5) return(NULL)
  
  if (is.null(ann)) ann <- load_annotations(species)
  target_col <- if(go_type == "mf") "GOs_mf" else "GOs_bp"
  
  ann_processed <- ann %>%
    mutate(all_gos_raw = replace_na(as.character(!!sym(target_col)), "")) %>%
    mutate(GO_list = str_extract_all(all_gos_raw, "GO:[0-9]{7}")) %>%
    unnest(GO_list) %>%
    filter(!is.na(GO_list) & GO_list != "") %>%
    distinct(gene_id, GO_list, .keep_all = TRUE)
  
  term2gene <- ann_processed %>%
    dplyr::select(GO_term = GO_list, Preferred_name) %>%
    filter(!is.na(Preferred_name) & Preferred_name != "" & Preferred_name != "-") %>%
    distinct()
  
  query_genes <- dge %>%
    inner_join(ann, by = c("transcript_id" = "gene_id")) %>%
    pull(Preferred_name) %>% unique()
  
  res <- enricher(gene = query_genes, TERM2GENE = term2gene, pvalueCutoff = 0.05)
  if (is.null(res) || nrow(res@result) == 0) return(NULL)
  
  res_df <- as.data.frame(res@result)
  go_mapping <- AnnotationDbi::select(GO.db, keys = res_df$ID, columns = c("TERM"), keytype = "GOID")
  
  res_df %>%
    dplyr::select(-Description) %>%
    left_join(go_mapping %>% dplyr::rename(ID = GOID, Description = TERM), by = "ID") %>%
    mutate(species = species)
}

phases <- c("intermoult", "premoult", "intramoult", "postmoult")

render_go_dotplot <- function(consensus_df, title_text, subtitle_text, out_pdf_path) {
  consensus_filtered <- consensus_df %>%
    filter(n_species >= 2) %>%
    mutate(phase = factor(phase, levels = phases))
  
  if (nrow(consensus_filtered) == 0) {
    cat("Warning: No consensus terms found for:", title_text, "\n")
    return(NULL)
  }
  
  plot_data <- consensus_filtered %>%
    mutate(Conservation = ifelse(n_species == 3, "Universal (3-way)", "Conserved (2-way)")) %>%
    mutate(Domain = case_when(
      grepl("neural|synapse|axon|sensory|brain|neuron|receptor|signaling", Description, ignore.case=T) ~ "Neural &\nSignaling",
      grepl("metabolic|oxid|energy|redox|lipid|carbohydrate|glycoly", Description, ignore.case=T) ~ "Metabolism",
      grepl("protease|peptidase|hydroly|proteoly|serine|endopeptidase", Description, ignore.case=T) ~ "Proteolysis",
      grepl("chitin|cuticle|structural|matrix|peritrophic", Description, ignore.case=T) ~ "Cuticle &\nStructure",
      TRUE ~ "Other\nFunctional"
    ))
  
  top_plot_data <- plot_data %>% 
    group_by(phase, Domain) %>% 
    slice_min(avg_padj, n = 5) %>% 
    ungroup()
  
  p <- ggplot(top_plot_data, aes(x = phase, y = reorder(format_go_term(Description), phase, FUN = function(x) min(as.numeric(x))), 
                                 size = max_count, color = -log10(avg_padj))) +
    geom_point(aes(shape = Conservation), stroke = 3.0) +
    scale_size_continuous(range = c(9, 30), name = "Max Gene Count") +
    scale_color_viridis_c(option = "plasma", name = "-log10(avg p.adjust)") +
    scale_shape_manual(values = c("Universal (3-way)" = 17, "Conserved (2-way)" = 16)) +
    scale_x_discrete(drop = FALSE) + 
    facet_grid(Domain ~ ., scales = "free_y", space = "free_y") +
    theme_bw(base_size = 40) +
    theme(
      axis.text.y = element_text(size = 48, face = "bold", lineheight = 0.78, color = "black"),
      axis.text.x = element_text(size = 52, face = "bold", color = "black", angle = 45, hjust = 1, vjust = 1),
      axis.title.x = element_text(size = 54, face = "bold", color = "black", margin = margin(t = 24)),
      axis.ticks.y = element_line(color = "grey30", linewidth = 1.5),
      axis.ticks.length.y = unit(0.6, "cm"),
      strip.text.y = element_text(size = 46, face = "bold", color = "black", angle = 0, lineheight = 0.85),
      plot.title = element_text(size = 56, face = "bold"),
      plot.subtitle = element_text(size = 42),
      legend.text = element_text(size = 36),
      legend.title = element_text(size = 40, face = "bold"),
      legend.key.size = unit(4.5, "lines"),
      panel.spacing = unit(2.5, "lines"),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 1.0),
      plot.margin = margin(30, 30, 30, 30)
    ) +
    guides(shape = guide_legend(override.aes = list(size = 20))) +
    labs(title = title_text, subtitle = subtitle_text, x = "Moult Phase", y = NULL, shape = "Conservation")
  
  ggsave(out_pdf_path, p, width = 38, height = 70, limitsize = FALSE, bg = "white")
  cat("  Saved plot to:", out_pdf_path, "\n")
  return(p)
}

# --- 1. HEAD VS FRONT LEGS (BP & MF) ---
for (gt in c("bp", "mf")) {
  cat("\nProcessing Head vs Front Legs (", toupper(gt), ")...\n")
  results_list <- list()
  for (ph in phases) {
    for (sp in all_species) {
      dge_file <- paste0("type4_", ph, "_f.legs_vs_head.csv")
      dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", dge_file)
      if (file.exists(dge_path)) {
        dge_head_up <- read_csv(dge_path, show_col_types = FALSE) %>% filter(log2FoldChange < 0)
        temp_file <- tempfile(fileext = ".csv")
        write_csv(dge_head_up, temp_file)
        res <- get_enrichment_from_csv(sp, temp_file, gt, ann = annotations_list[[sp]])
        if (!is.null(res)) results_list[[paste0(sp, "_", ph)]] <- res %>% mutate(phase = ph, go_type = gt)
      }
    }
  }
  if (length(results_list) > 0) {
    consensus <- bind_rows(results_list) %>%
      group_by(ID, Description, phase) %>%
      summarise(
        n_species = n(), min_padj = min(p.adjust), avg_padj = mean(p.adjust), 
        max_count = max(as.numeric(str_extract(Count, "^[0-9]+"))), .groups = "drop"
      )
    csv_file <- file.path(out_dir, paste0("Head_vs_FLegs_Consensus_", toupper(gt), ".csv"))
    write_csv(consensus, csv_file)
    
    sub_letter <- if (gt == "bp") "B" else "C"
    pdf_file <- file.path(out_dir, paste0("Figure_S5", sub_letter, "_Head_vs_FLegs_", toupper(gt), ".pdf"))
    render_go_dotplot(
      consensus, 
      paste0("Figure S5", sub_letter, ": Head vs Front Legs (", toupper(gt), ")"),
      "Conserved functional signatures enriched in the Head across the moult cycle",
      pdf_file
    )
  }
}

# --- 2. HIND LEGS VS FRONT LEGS (MF ONLY) ---
cat("\nProcessing Hind Legs vs Front Legs (MF)...\n")
results_list_hl_mf <- list()
for (ph in phases) {
  for (sp in all_species) {
    dge_file <- paste0("type4_", ph, "_f.legs_vs_h.legs.csv")
    dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", dge_file)
    if (file.exists(dge_path)) {
      dge_hlegs_up <- read_csv(dge_path, show_col_types = FALSE) %>% filter(log2FoldChange < 0)
      temp_file <- tempfile(fileext = ".csv")
      write_csv(dge_hlegs_up, temp_file)
      res <- get_enrichment_from_csv(sp, temp_file, "mf", ann = annotations_list[[sp]])
      if (!is.null(res)) results_list_hl_mf[[paste0(sp, "_", ph)]] <- res %>% mutate(phase = ph, go_type = "mf")
    }
  }
}
if (length(results_list_hl_mf) > 0) {
  consensus_hl_mf <- bind_rows(results_list_hl_mf) %>%
    group_by(ID, Description, phase) %>%
    summarise(
      n_species = n(), min_padj = min(p.adjust), avg_padj = mean(p.adjust), 
      max_count = max(as.numeric(str_extract(Count, "^[0-9]+"))), .groups = "drop"
    )
  csv_file <- file.path(out_dir, "Asymmetry_Consensus_All_Phases_MF.csv")
  write_csv(consensus_hl_mf, csv_file)
  
  pdf_file <- file.path(out_dir, "Figure_S5A_HLegs_vs_FLegs_MF.pdf")
  render_go_dotplot(
    consensus_hl_mf,
    "Figure S5A: Locomotory Asymmetry Waves - Hind Legs vs Front Legs (MF)",
    "Conserved molecular functions enriched in Hind Legs across the moult cycle",
    pdf_file
  )
}

cat("\nSUCCESS: Figure S5 suite generated successfully in:", out_dir, "\n")
