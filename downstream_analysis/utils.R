#' Project Utilities for Isopod RNA-seq Analysis
#' 
#' This script contains shared functions for path management, data loading, 
#' and ID conversion across the different analysis scripts.

library(tidyverse)

# Source global configuration
# This ensures PROJECT_ROOT, ACTIVE_SPECIES, etc. are available
if (file.exists("config.R")) {
  source("config.R")
} else {
  stop("config.R not found! Please create it based on the template.")
}

# --- PATH MANAGEMENT ---

#' Get the Project Root Directory
get_project_root <- function() {
  return(PROJECT_ROOT)
}

#' Build a Path to a Species Data Directory
#' Build a Path to a Species Data Directory
get_data_path <- function(species = ACTIVE_SPECIES, sub_dir = NULL, filename = NULL) {
  root <- get_project_root()
  
  # Candidate base directories
  candidates <- c(
    file.path(root, species),
    file.path(root, "cluster DGE", species),
    file.path(root, "local_scripts", "RNAseq", species),
    file.path(root, "..", "cluster DGE", species)
  )
  
  base_path <- candidates[dir.exists(candidates)][1]
  if (is.na(base_path)) {
    base_path <- file.path(root, species)
  }
  
  path <- base_path
  if (!is.null(sub_dir)) {
    path <- file.path(path, sub_dir)
  }
  if (!is.null(filename)) {
    path <- file.path(path, filename)
  }
  
  return(path)
}

# --- DATA LOADING ---

#' Load Comprehensive Hybrid Annotations
#' 
#' This function merges the species-specific KEGG/GO annotations with the 
#' master Crustacea Orthologous Group (OG) annotations. 
#' It prioritizes the descriptive names from the Crustacea file while 
#' keeping the rich pathway data from the species file.
#' 
#' @param species Character. Species name.
#' @return A data frame with merged annotations.
load_annotations <- function(species = ACTIVE_SPECIES) {
  
  # 1. Load Master Crustacea OG file
  root <- get_project_root()
  master_candidates <- c(
    file.path(root, "de_novo_transcriptomes_crustacea.og.annotations"),
    file.path(root, "transcriptome assembly", "de_novo_transcriptomes_crustacea.og.annotations"),
    file.path(root, "local_scripts", "RNAseq", "de_novo_transcriptomes_crustacea.og.annotations")
  )
  master_path <- master_candidates[file.exists(master_candidates)][1]
  
  if (is.na(master_path) || !file.exists(master_path)) {
    stop("Master Crustacea annotation file not found in project root or subdirectories.")
  }
  
  cat("Loading Master Crustacea annotations...\n")
  master_ann <- read.delim(master_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  # 2. Load Species-Specific file
  # (Contains: gene_id, KEGG_ko, KEGG_Pathway, GO_names, etc.)
  spec_path <- get_data_path(species, filename = "combined_annotations_with_GO_from_KO_try.tsv")
  if (!file.exists(spec_path)) {
    warning("Species-specific file not found at: ", spec_path, ". Using Master file only.")
    return(master_ann %>% filter(str_starts(gene_id, species)))
  }
  
  cat("Loading Species-specific (", species, ") annotations...\n")
  spec_ann <- read.delim(spec_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  # 3. Hybrid Merge (Left Join on gene_id)
  # We use the Master file as the base, filtered for the current species
  cat("Merging annotation sources...\n")
  
  hybrid_ann <- master_ann %>%
    filter(str_starts(gene_id, species)) %>%
    left_join(
      spec_ann %>% dplyr::select(gene_id, KEGG_ko, KEGG_Pathway, KEGG_Module, BRITE, PFAMs, GO_names, KEGG_Pathway_names), 
      by = "gene_id", 
      suffix = c("_master", "_spec")
    )
  
  # Cleanup: Ensure we have a single set of columns and reasonable names
  # If a gene is missing from the master but present in the spec, we might lose it.
  # Let's check for "orphan" genes in the spec file.
  orphan_count <- nrow(spec_ann) - sum(spec_ann$gene_id %in% master_ann$gene_id)
  if (orphan_count > 0) {
    cat("  Note:", orphan_count, "genes in species file were not found in Master Crustacea file.\n")
  }
  
  return(hybrid_ann)
}

#' Load General Crustacea Annotations (Legacy)
#' Deprecated: Use load_annotations() instead for the hybrid view.
load_general_annotations <- function() {
  master_path <- file.path(get_project_root(), "local_scripts", "RNAseq", "de_novo_transcriptomes_crustacea.og.annotations")
  read.delim(master_path, sep = "\t", stringsAsFactors = FALSE)
}

#' Load Averaged Normalized Counts (VST)
#' 
#' This function first attempts to load the global multi-species heatmap file. 
#' If it fails to find the relevant genes for the current species, it 
#' generates them directly from the species-specific RDS objects.
#' 
#' @param species Character. Species name.
#' @return A data frame with gene_id in the first column and averaged sample groups in subsequent columns.
load_normalized_counts <- function(species = ACTIVE_SPECIES) {
  
  # 1. Try Global File first (Legacy/Optimization)
  path <- file.path(get_project_root(), "heatmaps", "vst_normalized_counts_averaged_ordered.csv")
  
  if (file.exists(path)) {
    cat("Checking global counts file for", species, "...\n")
    # Read just the first column to check IDs efficiently
    first_col <- read_csv(path, col_select = 1, n_max = 500, show_col_types = FALSE)
    
    if (any(str_starts(pull(first_col, 1), species))) {
      cat("  Found species in global file. Loading...\n")
      counts <- read_csv(path, show_col_types = FALSE) %>%
        rename(gene_id = 1)
      return(counts)
    }
  }
  
  # 2. Fallback: Generate Averaged Counts from DDS
  cat("  Species not found in global file. Generating from", species, "RDS...\n")
  
  dds_path <- get_data_path(species, "dge_analysis", "dds_group_design.rds")
  if (!file.exists(dds_path)) {
    stop("Could not find counts file or DDS object for species: ", species)
  }
  
  library(DESeq2)
  dds <- readRDS(dds_path)
  
  cat("  Performing VST transformation...\n")
  vsd <- vst(dds, blind = FALSE)
  counts_mat <- assay(vsd)
  
  # Get group information for averaging
  meta <- as.data.frame(colData(dds))
  groups <- meta$group
  
  # Average by group
  cat("  Averaging replicates by sample group...\n")
  avg_counts <- do.call(cbind, lapply(unique(groups), function(g) {
    rowMeans(counts_mat[, groups == g, drop = FALSE])
  }))
  colnames(avg_counts) <- unique(groups)
  
  # Convert to tidy data frame
  avg_counts_df <- as.data.frame(avg_counts) %>%
    rownames_to_column("gene_id")
  
  return(avg_counts_df)
}

# --- LOGGING ---

#' Print a Standard Header for Scripts
log_header <- function(script_name) {
  cat(rep("=", 60), "\n", sep = "")
  cat(" SCRIPT:  ", script_name, "\n")
  cat(" SPECIES: ", ACTIVE_SPECIES, "\n")
  cat(" DATE:    ", as.character(Sys.time()), "\n")
  cat(rep("=", 60), "\n\n", sep = "")
}
