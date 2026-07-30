# === HIND-LEG VS FRONT-LEG TEMPORAL ASYMMETRY (MF & BP) ===
# This script identifies functional asymmetry by contrasting 
# Hind Legs against Front Legs across ALL moult phases (Type 4 DGE).

# === SETUP ===
project_path <- "G:/My Drive/PhD/Projects/Comparative Expressional Changes During the Moult Cycle in Land Isopods/local_scripts"
setwd(project_path)
source("utils.R")
library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(GO.db)

log_header("H.legs vs F.legs All Phases")

all_species <- c("maculatum", "laevis", "officinalis")
out_dir <- "output_csv/Asymmetry_Waves_All_Phases"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Pre-load annotations for all species ONCE to prevent redundant disk I/O in loops
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

for (gt in c("mf", "bp")) {
  cat("\n--- PROCESSING GO TYPE:", toupper(gt), "---\n")
  results_list <- list()
  
  for (ph in phases) {
    for (sp in all_species) {
      cat("  Processing:", sp, "@", ph, "...\n")
      dge_file <- paste0("type4_", ph, "_f.legs_vs_h.legs.csv")
      dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", dge_file)
      
      if (file.exists(dge_path)) {
        dge_h_up <- read_csv(dge_path, show_col_types = FALSE) %>% filter(log2FoldChange < 0)
        if (nrow(dge_h_up) >= 5) {
          temp_file <- tempfile(fileext = ".csv")
          write_csv(dge_h_up, temp_file)
          res <- get_enrichment_from_csv(sp, temp_file, gt, ann = annotations_list[[sp]])
          if (!is.null(res)) results_list[[paste0(sp, "_", ph)]] <- res %>% mutate(phase = ph, go_type = gt)
        }
      }
    }
  }
  
  if (length(results_list) == 0) {
    p_empty <- ggplot() + 
      scale_x_discrete(limits = phases, drop = FALSE) +
      theme_bw(base_size = 32) + 
      theme(
        axis.text.x = element_text(size = 38, face = "bold", color = "black", angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(size = 42, face = "bold")
      ) +
      labs(title = paste("Hind-Leg Specific Asymmetry Waves (", toupper(gt), ") - NO CONSERVED TERMS", sep=""), x = "Moult Phase")
    ggsave(file.path(out_dir, paste0("Asymmetry_Waves_All_Phases_", toupper(gt), ".png")), p_empty, width = 26, height = 12, bg = "white")
    next
  }
  
  # Consensus Synthesis - Strictly 2+ species
  consensus <- bind_rows(results_list) %>%
    group_by(ID, Description, phase) %>%
    summarise(
      n_species = n(), 
      min_padj = min(p.adjust), 
      avg_padj = mean(p.adjust), 
      max_count = max(as.numeric(str_extract(Count, "^[0-9]+"))), 
      .groups = "drop"
    ) %>%
    filter(n_species >= 2) %>%
    mutate(phase = factor(phase, levels = phases))

  if (nrow(consensus) == 0) {
     p_empty <- ggplot() + 
      scale_x_discrete(limits = phases, drop = FALSE) +
      theme_bw(base_size = 32) + 
      theme(
        axis.text.x = element_text(size = 38, face = "bold", color = "black", angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(size = 42, face = "bold")
      ) +
      labs(title = paste("Hind-Leg Specific Asymmetry Waves (", toupper(gt), ") - NO CONSERVED TERMS", sep=""), x = "Moult Phase")
    ggsave(file.path(out_dir, paste0("Asymmetry_Waves_All_Phases_", toupper(gt), ".png")), p_empty, width = 26, height = 12, bg = "white")
    next
  }

  # Domain Categorization
  plot_data <- consensus %>%
    mutate(Conservation = ifelse(n_species == 3, "Universal (3-way)", "Conserved (2-way)")) %>%
    mutate(Domain = case_when(
      grepl("neural|synapse|axon|sensory|brain|neuron|receptor|signaling", Description, ignore.case=T) ~ "Neural &\nSignaling",
      grepl("metabolic|oxid|energy|redox", Description, ignore.case=T) ~ "Metabolism",
      grepl("protease|peptidase|hydroly|proteoly", Description, ignore.case=T) ~ "Proteolysis",
      grepl("chitin|cuticle", Description, ignore.case=T) ~ "Cuticle",
      grepl("\\bion\\b|\\bions\\b|transport|calcif|calcium", Description, ignore.case=T) ~ "Mineral &\nIon Transport",
      TRUE ~ "Other\nFunctional"
    ))

  top_plot_data <- plot_data %>% group_by(phase, Domain) %>% slice_min(avg_padj, n = 5) %>% ungroup()

  # Visualization with massive font sizes, giant shape points, line wrapping, and horizontal domain headers
  p <- ggplot(top_plot_data, aes(x = phase, y = reorder(stringr::str_wrap(Description, width = 45), phase, FUN = function(x) min(as.numeric(x))), 
                               size = max_count, color = -log10(avg_padj))) +
    geom_point(aes(shape = Conservation), stroke = 2.5) +
    scale_size_continuous(range = c(8, 30), name = "Max Gene Count") +
    scale_color_viridis_c(option = "plasma", name = "-log10(avg p.adj)") +
    scale_shape_manual(values = c("Universal (3-way)" = 17, "Conserved (2-way)" = 16)) +
    scale_x_discrete(drop = FALSE) + 
    facet_grid(Domain ~ ., scales = "free_y", space = "free_y") +
    theme_bw(base_size = 32) +
    theme(
      axis.text.y = element_text(size = 34, face = "bold", lineheight = 0.85, color = "black"),
      axis.text.x = element_text(size = 38, face = "bold", color = "black", angle = 45, hjust = 1, vjust = 1),
      strip.text.y = element_text(size = 28, face = "bold", color = "black", angle = 0, lineheight = 0.85),
      plot.title = element_text(size = 42, face = "bold"),
      plot.subtitle = element_text(size = 34),
      legend.text = element_text(size = 30),
      legend.title = element_text(size = 32, face = "bold"),
      legend.key.size = unit(3, "lines"),
      panel.spacing = unit(2.5, "lines")
    ) +
    guides(shape = guide_legend(override.aes = list(size = 18))) +
    labs(title = paste("Hind-Leg Specific Asymmetry Waves (", toupper(gt), ")", sep=""), 
         subtitle = "Terms enriched in H.legs relative to F.legs across all phases", 
         x = "Moult Phase", y = NULL, shape = "Conservation")

  ggsave(file.path(out_dir, paste0("Asymmetry_Waves_All_Phases_", toupper(gt), ".png")), p, width = 26, height = 48, bg = "white")
  write_csv(consensus, file.path(out_dir, paste0("Asymmetry_Consensus_All_Phases_", toupper(gt), ".csv")))
}
