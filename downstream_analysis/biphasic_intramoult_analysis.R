# === BIPHASIC INTRAMOULT ANALYSIS ===
# This script analyzes anatomical asymmetry across all moult phases
# to highlight the Intramoult "Chimera" state.
# Now expanded to multiple contrasts (Front vs Back legs, Head vs Front legs).

# === SETUP ===
if (file.exists("utils.R")) {
  source("utils.R")
} else if (file.exists(file.path("downstream_analysis", "utils.R"))) {
  source(file.path("downstream_analysis", "utils.R"))
}
output_dir <- file.path("output_csv")
library(tidyverse)
library(patchwork)
library(DESeq2)
library(matrixStats)

log_header("Biphasic Asymmetry Analysis")

all_species <- c("maculatum", "laevis", "officinalis")
phases <- c("intermoult", "premoult", "intramoult", "postmoult")

# # Define the contrasts to analyze
# Each list item: c(Contrast_ID, Plot_Title, File_Pattern)
contrasts <- list(
  list(id = "f.legs_vs_h.legs", title = "Front vs Back Legs", file = "f.legs_vs_h.legs")
)

# Function to run the analysis for a given contrast
run_asymmetry_analysis <- function(contrast_info) {
  cid <- contrast_info$id
  ctitle <- contrast_info$title
  cfile <- contrast_info$file
  
  cat("\n", rep("=", 40), "\n", sep = "")
  cat(" ANALYZING CONTRAST:", ctitle, "\n")
  cat(rep("=", 40), "\n\n", sep = "")
  
  # --- 1. DATA COLLECTION ---
  cat("Collecting DEG counts...\n")
  deg_counts_list <- list()
  for (sp in all_species) {
    for (ph in phases) {
      cat("  [", sp, "] Phase:", ph, "... ")
      file_name <- paste0("type4_", ph, "_", cfile, ".csv")
      dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", file_name)
      
      # Handle potential reverse naming (head_vs_f.legs vs f.legs_vs_head)
      if (!file.exists(dge_path)) {
        parts <- str_split(cfile, "_vs_")[[1]]
        rev_file <- paste0(parts[2], "_vs_", parts[1])
        file_name <- paste0("type4_", ph, "_", rev_file, ".csv")
        dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", file_name)
      }
      
      if (file.exists(dge_path)) {
        dge <- read_csv(dge_path, show_col_types = FALSE)
        sig_count <- dge %>%
          dplyr::filter(padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
          nrow()
        deg_counts_list[[paste0(sp, "_", ph)]] <- data.frame(species = sp, phase = ph, deg_count = sig_count)
        cat("Found", sig_count, "DEGs\n")
      } else {
        cat("FAILED (File not found)\n")
      }
    }
  }
  
  if (length(deg_counts_list) == 0) return(NULL)
  
  deg_summary <- bind_rows(deg_counts_list)
  deg_summary$phase <- factor(deg_summary$phase, levels = phases)
  sp_names <- c("maculatum" = "A. maculatum", "laevis" = "P. laevis", "officinalis" = "A. officinalis")
  deg_summary$species_name <- sp_names[deg_summary$species]
  deg_summary$species_name <- factor(deg_summary$species_name, levels = sp_names)
  
  # --- 2. STATISTICAL VALIDATION (The Intramoult Surge) ---
  cat("\nCalculating Significance of the Intramoult Surge...\n")
  intramoult_vals <- deg_summary %>% dplyr::filter(phase == "intramoult") %>% pull(deg_count)
  baseline_vals <- deg_summary %>% dplyr::filter(phase == "intermoult") %>% pull(deg_count)
  other_vals <- deg_summary %>% dplyr::filter(phase != "intramoult") %>% pull(deg_count)
  
  w_res <- wilcox.test(intramoult_vals, other_vals)
  mean_intra <- mean(intramoult_vals)
  mean_baseline <- mean(baseline_vals)
  fold_increase <- mean_intra / max(1, mean_baseline)
  
  stat_res_lines <- c(
    paste0("--- ", toupper(cid), " ASYMMETRY SURGE STATISTICS ---"),
    paste("Date:", Sys.time()),
    "",
    sprintf("Wilcoxon Rank Sum Test (Intramoult vs Other Phases):"),
    sprintf("  p-value = %.6f", w_res$p.value),
    "",
    sprintf("Biological Significance (Relative to Intermoult):"),
    sprintf("  Mean DEGs Intramoult: %.1f", mean_intra),
    sprintf("  Mean DEGs Intermoult: %.1f", mean_baseline),
    sprintf("  Fold Increase: %.2f-fold", fold_increase)
  )
  
  writeLines(stat_res_lines, file.path(output_dir, paste0("intramoult_Assymetry_", cid, "_statistics.txt")))
  
  # --- 3. VISUALIZATION ---
  cat("Generating Summary Plot...\n")
  # Determine significance label
  sig_label <- ifelse(w_res$p.value < 0.01, "** p < 0.01", 
                      ifelse(w_res$p.value < 0.05, "* p < 0.05", "ns"))
  
  p_summary <- ggplot(deg_summary, aes(x = phase, y = deg_count, fill = species_name)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    theme_bw() +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title = paste("Asymmetry:", ctitle),
      subtitle = paste0("Intramoult surge: ", round(fold_increase), "-fold increase (", sig_label, ")"),
      x = "Moult Phase",
      y = "Significant DEGs (|LFC|>1, padj<0.05)",
      fill = "Species"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      legend.text = element_text(face = "italic"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  if (w_res$p.value < 0.05) {
    max_y <- max(deg_summary$deg_count) * 1.1
    p_summary <- p_summary + 
      # annotate("segment", x = 1.6, xend = 4.4, y = max_y * 1.1, yend = max_y * 1.1, size = 1) +
      annotate("text", x = 3, y = max_y * 0.95, label = sig_label, fontface = "bold", size = 5)
  }
  
  ggsave(file.path(output_dir, "intramoult_Assymetry.png"), p_summary, width = 8, height = 6, bg = "white")
  
  # --- 4. FUNCTIONAL POLARITY ANALYSIS (Conserved Intramoult) ---
  cat("Analyzing Conserved Functional Polarity...\n")
  intramoult_genes <- list()
  for (sp in all_species) {
    file_name <- paste0("type4_intramoult_", cfile, ".csv")
    dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", file_name)
    
    if (!file.exists(dge_path)) {
      parts <- str_split(cfile, "_vs_")[[1]]
      rev_file <- paste0(parts[2], "_vs_", parts[1])
      file_name <- paste0("type4_intramoult_", rev_file, ".csv")
      dge_path <- get_data_path(sp, "dge_analysis/dge_type4_within_phase_tissues/full", file_name)
    }
    
    if (file.exists(dge_path)) {
      ann <- load_annotations(sp)
      intramoult_genes[[sp]] <- read_csv(dge_path, show_col_types = FALSE) %>%
        dplyr::inner_join(ann, by = c("transcript_id" = "gene_id")) %>%
        dplyr::filter(padj < PADJ_CUTOFF & abs(log2FoldChange) > LFC_CUTOFF) %>%
        # REMOVED strict Preferred_name filter to include uncharacterized OGs
        dplyr::select(Preferred_name, ODB_OG, log2FoldChange) %>%
        dplyr::distinct()
    }
  }
  
  if (length(intramoult_genes) == 3) {
    og_list_for_int <- lapply(intramoult_genes, function(df) df$ODB_OG)
    common_3way_ogs <- Reduce(intersect, og_list_for_int)
    
    parts <- str_split(ctitle, " vs ")[[1]]
    
    df_polarity <- bind_rows(intramoult_genes) %>%
      dplyr::filter(ODB_OG %in% common_3way_ogs) %>%
      dplyr::group_by(ODB_OG) %>%
      dplyr::summarise(
        # Robust naming: use first available Preferred_name, fallback to OG ID
        Preferred_name = dplyr::first(na.omit(Preferred_name[Preferred_name != "" & Preferred_name != "-"])),
        avg_lfc = mean(log2FoldChange),
        consistent_direction = all(log2FoldChange > 0) | all(log2FoldChange < 0),
        .groups = "drop"
      ) %>%
      dplyr::mutate(Preferred_name = ifelse(is.na(Preferred_name), paste0("Uncharacterized_", ODB_OG), Preferred_name)) %>%
      dplyr::mutate(Bias = ifelse(avg_lfc > 0, parts[1], parts[2])) %>%
      dplyr::arrange(desc(abs(avg_lfc)))
    
    write_csv(df_polarity, file.path(output_dir, paste0("intramoult_Assymetry_", cid, "_Polarity_List_Full.csv")))
  }
  
  return(p_summary)
}

# Run all contrasts
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
all_plots <- list()
for (cnt in contrasts) {
  p <- run_asymmetry_analysis(cnt)
  if (!is.null(p)) all_plots[[cnt$id]] <- p
}

cat("\nGlobal Analysis complete! All asymmetry outputs saved to 'output_csv/'.\n")
