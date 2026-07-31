# Custom R packages path (uncomment if custom R packages directory is needed)
# .libPaths("/path/to/R_packages")

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  project_dir <- args[1]
} else {
  project_dir <- ".."
}

cat("=== Complete Multi-Factor DGE Analysis (DESeq2 Only) ===\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("Project directory:", project_dir, "\n\n")

# === CONFIGURATION ===
quant_dir   <- file.path(project_dir, "quantification")
sample_file <- file.path(project_dir, "samples.tsv")

print(project_dir)

# =============================================================================
# === GLOBAL: FIXED PHASE COMPARISON DIRECTION (A vs B) ===
# =============================================================================

desired_phase_contrasts <- tribble(
  ~A,           ~B,
  "premoult",   "intermoult",
  "intramoult", "premoult",
  "postmoult",  "intramoult",
  "intermoult", "postmoult",
  "premoult",   "postmoult",
  "intramoult", "intermoult"
)

intermoult_ref_contrasts <- tribble(
  ~A,           ~B,
  "premoult",   "intermoult",
  "intramoult", "intermoult",
  "postmoult",  "intermoult"
)

phase_pair_map <- list()
for (k in seq_len(nrow(desired_phase_contrasts))) {
  A <- desired_phase_contrasts$A[k]
  B <- desired_phase_contrasts$B[k]
  phase_pair_map[[paste(A, B, sep="|")]] <- c(A=A, B=B)
  phase_pair_map[[paste(B, A, sep="|")]] <- c(A=A, B=B)
}

orient_phase_pair <- function(p1, p2) {
  key <- paste(p1, p2, sep="|")
  if (!key %in% names(phase_pair_map)) {
    stop("Missing desired phase pair mapping for: ", key,
         "\nAdd it to desired_phase_contrasts.")
  }
  phase_pair_map[[key]]
}

# =============================================================================
# === 1. LOAD AND PREPARE METADATA ===
# =============================================================================

cat("1. Loading samples...\n")
samples <- read_tsv(sample_file, show_col_types = FALSE)

cat("\n=== DIAGNOSTIC: Sample Data ===\n")
cat("Total rows:", nrow(samples), "\n")
cat("Column names:", paste(colnames(samples), collapse=", "), "\n")
cat("\nFirst 5 rows:\n")
print(head(samples, 5))

cat("\nChecking for NA values:\n")
cat("  sample NA:", sum(is.na(samples$sample)), "\n")
cat("  tissue NA:", sum(is.na(samples$tissue)), "\n")
cat("  phase NA:",  sum(is.na(samples$phase)),  "\n")

samples$tissue <- trimws(samples$tissue)
samples$phase  <- trimws(samples$phase)
samples$sample <- trimws(samples$sample)

samples$phase <- factor(samples$phase,
                        levels = c("premoult", "intermoult", "intramoult", "postmoult"))
samples$tissue  <- factor(samples$tissue)
samples$species <- factor(samples$species)

samples$group <- interaction(samples$tissue, samples$phase, sep="_", drop=TRUE)

cat("\nAfter factor conversion:\n")
cat("  phase NA:", sum(is.na(samples$phase)), "\n")
cat("  tissue NA:", sum(is.na(samples$tissue)), "\n")
cat("  group NA:", sum(is.na(samples$group)), "\n")

if (any(is.na(samples$group))) {
  cat("\n=== SAMPLES WITH NA IN GROUP ===\n")
  print(samples[is.na(samples$group), c("sample", "tissue", "phase")])
  stop("Fix NA values in samples.tsv")
}

cat("\nUnique phases found:", paste(unique(as.character(samples$phase)), collapse=", "), "\n")
cat("Unique tissues found:", paste(unique(as.character(samples$tissue)), collapse=", "), "\n")
cat("================================\n\n")

# Support both direct sample folder and legacy _merged_trimmed folder
sample_paths_direct  <- file.path(quant_dir, samples$sample, "quant.sf")
sample_paths_trimmed <- file.path(quant_dir, paste0(samples$sample, "_merged_trimmed"), "quant.sf")
samples$path <- ifelse(file.exists(sample_paths_direct), sample_paths_direct, sample_paths_trimmed)

files_exist <- file.exists(samples$path)
n_found <- sum(files_exist)
n_total <- nrow(samples)

cat("Files found:", n_found, "of", n_total, "\n")

if (n_found < n_total) {
  cat("\n=== ERROR: MISSING QUANT FILES ===\n")
  missing_samples <- samples[!files_exist, ]
  cat("Missing", nrow(missing_samples), "sample(s):\n\n")
  for (i in 1:nrow(missing_samples)) {
    cat("  Sample:", missing_samples$sample[i], "\n")
    cat("  Expected path:", missing_samples$path[i], "\n")
    cat("  Tissue:", missing_samples$tissue[i],
        " | Phase:", missing_samples$phase[i], "\n\n")
  }
  cat("=== ANALYSIS STOPPED ===\n")
  cat("Fix missing files and re-run.\n\n")
  stop("Missing quant.sf files detected. Cannot proceed.")
}

cat("All quant.sf files found! Proceeding...\n\n")
write_csv(samples, "samples_processed.csv")

# =============================================================================
# === 2. IMPORT SALMON DATA ===
# =============================================================================

cat("2. Importing Salmon quantifications...\n")
txi <- tximport(samples$path, type="salmon", txOut=TRUE)

# =============================================================================
# === 3. CREATE DESEQ2 OBJECTS ===
# =============================================================================

cat("3. Creating DESeq2 objects...\n")

dds <- DESeqDataSetFromTximport(txi, colData=samples, design=~group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat("  Main dds: Filtered to", nrow(dds), "genes with >=10 counts\n")

dds_phase <- DESeqDataSetFromTximport(txi, colData=samples, design=~phase)
keep_phase <- rowSums(counts(dds_phase)) >= 10
dds_phase <- dds_phase[keep_phase,]
cat("  Phase dds: Filtered to", nrow(dds_phase), "genes with >=10 counts\n")

cat("\nGroup levels:\n")
print(levels(dds$group))

cat("\n=== Checking replication ===\n")
group_counts <- table(dds$group)
single_rep_groups <- names(group_counts[group_counts == 1])

if (length(single_rep_groups) > 0) {
  cat("WARNING: Single-replicate groups detected:\n")
  print(single_rep_groups)
  cat("These groups will use shared dispersion estimates.\n")
  cat("Results involving these groups will be less statistically reliable.\n\n")
} else {
  cat("All groups have multiple replicates.\n\n")
}

# =============================================================================
# === 4. RUN DESEQ ===
# =============================================================================

cat("4. Running DESeq...\n")
dds       <- DESeq(dds)
dds_phase <- DESeq(dds_phase)

# Create a version excluding h.legs for the new categories (if h.legs tissue exists)
cat("  Creating dds_phase_no_hlegs (excluding h.legs if present)...\n")
if ("h.legs" %in% colData(dds_phase)$tissue) {
  dds_phase_no_hlegs <- dds_phase[, colData(dds_phase)$tissue != "h.legs"]
} else {
  dds_phase_no_hlegs <- dds_phase
}
# Filter out genes that might now have zero counts or very low counts in this subset
keep_no_hlegs <- rowSums(counts(dds_phase_no_hlegs)) >= 10
dds_phase_no_hlegs <- dds_phase_no_hlegs[keep_no_hlegs,]
dds_phase_no_hlegs <- DESeq(dds_phase_no_hlegs)

cat("  All DESeq2 objects completed.\n\n")

all_groups   <- levels(dds$group)
phase_levels <- levels(samples$phase)
tissues <- sort(unique(sapply(strsplit(as.character(all_groups), "_"), `[`, 1)))

# =============================================================================
# === TYPE 1: PHASE-VS-PHASE ACROSS ALL TISSUES (DESeq2 ~phase) ===
# =============================================================================

cat("\n=== TYPE 1: Phase-vs-phase across all tissues (DESeq2) ===\n")
dir.create("dge_type1_phase_vs_phase_alltissues",        showWarnings=FALSE)
dir.create("dge_type1_phase_vs_phase_alltissues/full",   showWarnings=FALSE)
dir.create("dge_type1_phase_vs_phase_alltissues/sig",    showWarnings=FALSE)

type1_summary <- data.frame()

for (k in seq_len(nrow(desired_phase_contrasts))) {
  test_phase <- desired_phase_contrasts$A[k]
  ref_phase  <- desired_phase_contrasts$B[k]

  cat("  ", test_phase, "vs", ref_phase, "\n")

  res <- results(dds_phase, contrast=c("phase", test_phase, ref_phase))

  res_full <- as.data.frame(res) %>%
    rownames_to_column("transcript_id")

  filename_base <- paste0("type1_", test_phase, "_vs_", ref_phase, "_alltissues")

  write_csv(res_full,
            file.path("dge_type1_phase_vs_phase_alltissues/full", paste0(filename_base, ".csv")))

  sig <- res_full %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(padj)

  write_csv(sig,
            file.path("dge_type1_phase_vs_phase_alltissues/sig", paste0(filename_base, ".csv")))

  type1_summary <- rbind(type1_summary,
                         data.frame(comparison    = paste0(test_phase, "_vs_", ref_phase, "_alltissues"),
                                    n_significant = nrow(sig)))
  cat("    Significant DEGs:", nrow(sig), "\n")
}

write_csv(type1_summary %>% arrange(desc(n_significant)),
          "type1_phase_alltissues_summary.csv")

# =============================================================================
# === TYPE 1B: PHASE-VS-INTERMOULT ACROSS ALL TISSUES (DESeq2 ~phase) ===
# =============================================================================

cat("\n=== TYPE 1B: Phase-vs-intermoult across all tissues (intermoult reference) ===\n")
dir.create("dge_type1b_phase_vs_intermoult_alltissues",       showWarnings=FALSE)
dir.create("dge_type1b_phase_vs_intermoult_alltissues/full",  showWarnings=FALSE)
dir.create("dge_type1b_phase_vs_intermoult_alltissues/sig",   showWarnings=FALSE)

type1b_summary <- data.frame()

for (k in seq_len(nrow(intermoult_ref_contrasts))) {
  test_phase <- intermoult_ref_contrasts$A[k]
  ref_phase  <- intermoult_ref_contrasts$B[k]  # always "intermoult"

  cat("  ", test_phase, "vs intermoult (reference)\n")

  res <- results(dds_phase, contrast=c("phase", test_phase, ref_phase))

  res_full <- as.data.frame(res) %>%
    rownames_to_column("transcript_id")

  filename_base <- paste0("type1b_", test_phase, "_vs_intermoult_alltissues")

  write_csv(res_full,
            file.path("dge_type1b_phase_vs_intermoult_alltissues/full", paste0(filename_base, ".csv")))

  sig <- res_full %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(padj)

  write_csv(sig,
            file.path("dge_type1b_phase_vs_intermoult_alltissues/sig", paste0(filename_base, ".csv")))

  type1b_summary <- rbind(type1b_summary,
                          data.frame(comparison    = paste0(test_phase, "_vs_intermoult_alltissues"),
                                     n_significant = nrow(sig)))
  cat("    Significant DEGs:", nrow(sig), "\n")
}

write_csv(type1b_summary %>% arrange(desc(n_significant)),
          "type1b_intermoult_ref_summary.csv")

# =============================================================================
# === TYPE 1C: PHASE-VS-PHASE EXCLUDING H.LEGS (DESeq2 ~phase) ===
# =============================================================================

cat("\n=== TYPE 1C: Phase-vs-phase excluding h.legs (DESeq2) ===\n")
dir.create("dge_type1c_phase_vs_phase_no_hlegs",        showWarnings=FALSE)
dir.create("dge_type1c_phase_vs_phase_no_hlegs/full",   showWarnings=FALSE)
dir.create("dge_type1c_phase_vs_phase_no_hlegs/sig",    showWarnings=FALSE)

type1c_summary <- data.frame()

for (k in seq_len(nrow(desired_phase_contrasts))) {
  test_phase <- desired_phase_contrasts$A[k]
  ref_phase  <- desired_phase_contrasts$B[k]

  cat("  ", test_phase, "vs", ref_phase, " (no h.legs)\n")

  res <- results(dds_phase_no_hlegs, contrast=c("phase", test_phase, ref_phase))

  res_full <- as.data.frame(res) %>%
    rownames_to_column("transcript_id")

  filename_base <- paste0("type1c_", test_phase, "_vs_", ref_phase, "_no_hlegs")

  write_csv(res_full,
            file.path("dge_type1c_phase_vs_phase_no_hlegs/full", paste0(filename_base, ".csv")))

  sig <- res_full %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(padj)

  write_csv(sig,
            file.path("dge_type1c_phase_vs_phase_no_hlegs/sig", paste0(filename_base, ".csv")))

  type1c_summary <- rbind(type1c_summary,
                         data.frame(comparison    = paste0(test_phase, "_vs_", ref_phase, "_no_hlegs"),
                                    n_significant = nrow(sig)))
  cat("    Significant DEGs:", nrow(sig), "\n")
}

write_csv(type1c_summary %>% arrange(desc(n_significant)),
          "type1c_phase_no_hlegs_summary.csv")

# =============================================================================
# === TYPE 1D: PHASE-VS-INTERMOULT EXCLUDING H.LEGS (DESeq2 ~phase) ===
# =============================================================================

cat("\n=== TYPE 1D: Phase-vs-intermoult excluding h.legs (intermoult reference) ===\n")
dir.create("dge_type1d_phase_vs_intermoult_no_hlegs",       showWarnings=FALSE)
dir.create("dge_type1d_phase_vs_intermoult_no_hlegs/full",  showWarnings=FALSE)
dir.create("dge_type1d_phase_vs_intermoult_no_hlegs/sig",   showWarnings=FALSE)

type1d_summary <- data.frame()

for (k in seq_len(nrow(intermoult_ref_contrasts))) {
  test_phase <- intermoult_ref_contrasts$A[k]
  ref_phase  <- intermoult_ref_contrasts$B[k]  # always "intermoult"

  cat("  ", test_phase, "vs intermoult (no h.legs)\n")

  res <- results(dds_phase_no_hlegs, contrast=c("phase", test_phase, ref_phase))

  res_full <- as.data.frame(res) %>%
    rownames_to_column("transcript_id")

  filename_base <- paste0("type1d_", test_phase, "_vs_intermoult_no_hlegs")

  write_csv(res_full,
            file.path("dge_type1d_phase_vs_intermoult_no_hlegs/full", paste0(filename_base, ".csv")))

  sig <- res_full %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
    arrange(padj)

  write_csv(sig,
            file.path("dge_type1d_phase_vs_intermoult_no_hlegs/sig", paste0(filename_base, ".csv")))

  type1d_summary <- rbind(type1d_summary,
                          data.frame(comparison    = paste0(test_phase, "_vs_intermoult_no_hlegs"),
                                     n_significant = nrow(sig)))
  cat("    Significant DEGs:", nrow(sig), "\n")
}

write_csv(type1d_summary %>% arrange(desc(n_significant)),
          "type1d_intermoult_ref_no_hlegs_summary.csv")

# =============================================================================
# === TYPE 2: WITHIN-TISSUE PHASE COMPARISONS (DESeq2 ~group) ===
# =============================================================================

cat("\n=== TYPE 2: Within-tissue phase comparisons (DESeq2) ===\n")
dir.create("dge_type2_within_tissue_phases",       showWarnings=FALSE)
dir.create("dge_type2_within_tissue_phases/full",  showWarnings=FALSE)
dir.create("dge_type2_within_tissue_phases/sig",   showWarnings=FALSE)

type2_contrasts <- list()

for (tis in tissues) {
  for (k in seq_len(nrow(desired_phase_contrasts))) {
    phA <- desired_phase_contrasts$A[k]
    phB <- desired_phase_contrasts$B[k]
    gA  <- paste0(tis, "_", phA)
    gB  <- paste0(tis, "_", phB)
    if (gA %in% all_groups && gB %in% all_groups) {
      name <- paste0("type2_", tis, "_", phA, "_vs_", phB)
      type2_contrasts[[name]] <- c(gA, gB)
    }
  }
}

cat("Total Type 2 contrasts:", length(type2_contrasts), "\n")

type2_summary <- data.frame()

for (i in seq_along(type2_contrasts)) {
  contrast_name <- names(type2_contrasts)[i]
  contrast_pair <- type2_contrasts[[contrast_name]]

  cat("  ", contrast_name, "\n")

  res <- tryCatch({
    results(dds, contrast=c("group", contrast_pair[1], contrast_pair[2]))
  }, error = function(e) {
    cat("    ERROR:", e$message, "\n"); NULL
  })

  if (!is.null(res)) {
    contrast_filename <- gsub(" ", "_", contrast_name)

    res_full <- as.data.frame(res) %>% rownames_to_column("transcript_id")
    write_csv(res_full,
              file.path("dge_type2_within_tissue_phases/full", paste0(contrast_filename, ".csv")))

    sig <- as.data.frame(res) %>%
      rownames_to_column("transcript_id") %>%
      filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
      arrange(padj)
    write_csv(sig,
              file.path("dge_type2_within_tissue_phases/sig", paste0(contrast_filename, ".csv")))

    type2_summary <- rbind(type2_summary,
                           data.frame(contrast=contrast_name, n_significant=nrow(sig)))
    cat("    Significant DEGs:", nrow(sig), "\n")
  }
}

write_csv(type2_summary %>% arrange(desc(n_significant)),
          "type2_contrasts_summary.csv")

# =============================================================================
# === TYPE 2B: WITHIN-TISSUE PHASE-VS-INTERMOULT (DESeq2 ~group) ===
# =============================================================================

cat("\n=== TYPE 2B: Within-tissue vs intermoult (intermoult reference) ===\n")
dir.create("dge_type2b_within_tissue_vs_intermoult",       showWarnings=FALSE)
dir.create("dge_type2b_within_tissue_vs_intermoult/full",  showWarnings=FALSE)
dir.create("dge_type2b_within_tissue_vs_intermoult/sig",   showWarnings=FALSE)

type2b_summary <- data.frame()

for (tis in tissues) {
  for (k in seq_len(nrow(intermoult_ref_contrasts))) {
    phA <- intermoult_ref_contrasts$A[k]
    phB <- intermoult_ref_contrasts$B[k]  # always "intermoult"
    gA  <- paste0(tis, "_", phA)
    gB  <- paste0(tis, "_", phB)

    if (gA %in% all_groups && gB %in% all_groups) {
      cat("  ", tis, "-", phA, "vs intermoult\n")

      res <- tryCatch({
        results(dds, contrast=c("group", gA, gB))
      }, error = function(e) {
        cat("    ERROR:", e$message, "\n"); NULL
      })

      if (!is.null(res)) {
        name <- paste0("type2b_", tis, "_", phA, "_vs_intermoult")
        contrast_filename <- gsub(" ", "_", name)

        res_full <- as.data.frame(res) %>% rownames_to_column("transcript_id")
        write_csv(res_full,
                  file.path("dge_type2b_within_tissue_vs_intermoult/full", paste0(contrast_filename, ".csv")))

        sig <- as.data.frame(res) %>%
          rownames_to_column("transcript_id") %>%
          filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
          arrange(padj)
        write_csv(sig,
                  file.path("dge_type2b_within_tissue_vs_intermoult/sig", paste0(contrast_filename, ".csv")))

        type2b_summary <- rbind(type2b_summary,
                                data.frame(contrast=name, n_significant=nrow(sig)))
        cat("    Significant DEGs:", nrow(sig), "\n")
      }
    }
  }
}

write_csv(type2b_summary %>% arrange(desc(n_significant)),
          "type2b_intermoult_ref_summary.csv")

# =============================================================================
# === TYPE 3: TISSUE-VS-ALL-OTHERS WITHIN PHASE (DESeq2 ~group) ===
# =============================================================================

cat("\n=== TYPE 3: Tissue-vs-all-others within each phase (DESeq2) ===\n")
dir.create("dge_type3_tissue_vs_allothers",       showWarnings=FALSE)
dir.create("dge_type3_tissue_vs_allothers/full",  showWarnings=FALSE)
dir.create("dge_type3_tissue_vs_allothers/sig",   showWarnings=FALSE)

type3_summary <- data.frame()

for (ph in phase_levels) {
  phase_groups  <- all_groups[grepl(paste0("_", ph, "$"), all_groups)]
  phase_tissues <- sort(unique(sapply(strsplit(as.character(phase_groups), "_"), `[`, 1)))

  for (test_tis in phase_tissues) {
    test_group   <- paste0(test_tis, "_", ph)
    if (!test_group %in% all_groups) next
    other_groups <- setdiff(phase_groups, test_group)
    if (length(other_groups) < 1) next

    cat("  ", test_tis, "@", ph, "vs all others\n")

    contrast_vec <- numeric(length(all_groups))
    names(contrast_vec) <- all_groups
    contrast_vec[test_group]   <- 1
    contrast_vec[other_groups] <- -1 / length(other_groups)

    res <- results(dds, contrast=contrast_vec)

    res_full      <- as.data.frame(res) %>% rownames_to_column("transcript_id")
    filename_base <- paste0("type3_", test_tis, "_vs_allothers_", ph)

    write_csv(res_full,
              file.path("dge_type3_tissue_vs_allothers/full", paste0(filename_base, ".csv")))

    sig <- res_full %>%
      filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
      arrange(padj)
    write_csv(sig,
              file.path("dge_type3_tissue_vs_allothers/sig", paste0(filename_base, ".csv")))

    type3_summary <- rbind(type3_summary,
                           data.frame(comparison    = paste0(test_tis, "_vs_allothers_", ph),
                                      n_significant = nrow(sig)))
    cat("    Significant DEGs:", nrow(sig), "\n")
  }
}

write_csv(type3_summary %>% arrange(desc(n_significant)),
          "type3_tissue_vs_allothers_summary.csv")

# =============================================================================
# === TYPE 4: WITHIN-PHASE TISSUE COMPARISONS (DESeq2 ~group) ===
# =============================================================================

cat("\n=== TYPE 4: Within-phase tissue comparisons (DESeq2) ===\n")
dir.create("dge_type4_within_phase_tissues",       showWarnings=FALSE)
dir.create("dge_type4_within_phase_tissues/full",  showWarnings=FALSE)
dir.create("dge_type4_within_phase_tissues/sig",   showWarnings=FALSE)

type4_contrasts <- list()

for (ph in phase_levels) {
  phase_groups  <- all_groups[grepl(paste0("_", ph, "$"), all_groups)]
  phase_tissues <- sort(unique(sapply(strsplit(as.character(phase_groups), "_"), `[`, 1)))

  if (length(phase_tissues) > 1) {
    for (i in 1:(length(phase_tissues)-1)) {
      for (j in (i+1):length(phase_tissues)) {
        tis1 <- phase_tissues[i]
        tis2 <- phase_tissues[j]
        g1   <- paste0(tis1, "_", ph)
        g2   <- paste0(tis2, "_", ph)
        if (g1 %in% all_groups && g2 %in% all_groups) {
          name <- paste0("type4_", gsub(" ", "", ph), "_", tis1, "_vs_", tis2)
          type4_contrasts[[name]] <- c(g1, g2)
        }
      }
    }
  }
}

cat("Total Type 4 contrasts:", length(type4_contrasts), "\n")

type4_summary <- data.frame()

for (i in seq_along(type4_contrasts)) {
  contrast_name <- names(type4_contrasts)[i]
  contrast_pair <- type4_contrasts[[contrast_name]]

  cat("  ", contrast_name, "\n")

  res <- tryCatch({
    results(dds, contrast=c("group", contrast_pair[1], contrast_pair[2]))
  }, error = function(e) {
    cat("    ERROR:", e$message, "\n"); NULL
  })

  if (!is.null(res)) {
    contrast_filename <- gsub(" ", "_", contrast_name)

    res_full <- as.data.frame(res) %>% rownames_to_column("transcript_id")
    write_csv(res_full,
              file.path("dge_type4_within_phase_tissues/full", paste0(contrast_filename, ".csv")))

    sig <- as.data.frame(res) %>%
      rownames_to_column("transcript_id") %>%
      filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
      arrange(padj)
    write_csv(sig,
              file.path("dge_type4_within_phase_tissues/sig", paste0(contrast_filename, ".csv")))

    type4_summary <- rbind(type4_summary,
                           data.frame(contrast=contrast_name, n_significant=nrow(sig)))
    cat("    Significant DEGs:", nrow(sig), "\n")
  }
}

write_csv(type4_summary %>% arrange(desc(n_significant)),
          "type4_within_phase_tissues_summary.csv")

# =============================================================================
# === TYPE 5: CROSS-TISSUE-CROSS-PHASE COMPARISONS (DESeq2 ~group) ===
# =============================================================================

cat("\n=== TYPE 5: Cross-tissue-cross-phase comparisons (DESeq2) ===\n")
dir.create("dge_type5_cross_tissue_cross_phase",       showWarnings=FALSE)
dir.create("dge_type5_cross_tissue_cross_phase/full",  showWarnings=FALSE)
dir.create("dge_type5_cross_tissue_cross_phase/sig",   showWarnings=FALSE)

type5_contrasts <- list()

for (i in 1:(length(all_groups)-1)) {
  for (j in (i+1):length(all_groups)) {
    group1 <- all_groups[i]
    group2 <- all_groups[j]

    parts1  <- strsplit(as.character(group1), "_")[[1]]
    parts2  <- strsplit(as.character(group2), "_")[[1]]
    tissue1 <- parts1[1];  phase1 <- paste(parts1[-1], collapse="_")
    tissue2 <- parts2[1];  phase2 <- paste(parts2[-1], collapse="_")

    if (tissue1 == tissue2) next
    if (phase1  == phase2)  next

    ab <- orient_phase_pair(phase1, phase2)

    if (phase1 == ab["A"]) {
      groupA <- group1; groupB <- group2
      tissueA <- tissue1; phaseA <- phase1
      tissueB <- tissue2; phaseB <- phase2
    } else {
      groupA <- group2; groupB <- group1
      tissueA <- tissue2; phaseA <- phase2
      tissueB <- tissue1; phaseB <- phase1
    }

    name <- paste0("type5_", tissueA, "_", gsub(" ", "", phaseA),
                   "_vs_", tissueB, "_", gsub(" ", "", phaseB))
    type5_contrasts[[name]] <- c(groupA, groupB)
  }
}

cat("Total Type 5 contrasts:", length(type5_contrasts), "\n")
cat("(Cross-tissue AND cross-phase only, excluding Type 2 and Type 4)\n\n")

type5_summary    <- data.frame()
comparison_count <- 0

for (i in seq_along(type5_contrasts)) {
  contrast_name <- names(type5_contrasts)[i]
  contrast_pair <- type5_contrasts[[contrast_name]]

  comparison_count <- comparison_count + 1
  if (comparison_count %% 10 == 0) {
    cat("  Processed", comparison_count, "of", length(type5_contrasts), "comparisons...\n")
  }

  res <- tryCatch({
    results(dds, contrast=c("group", contrast_pair[1], contrast_pair[2]))
  }, error = function(e) {
    cat("    ERROR:", contrast_name, "-", e$message, "\n"); NULL
  })

  if (!is.null(res)) {
    contrast_filename <- gsub(" ", "_", contrast_name)

    res_full <- as.data.frame(res) %>% rownames_to_column("transcript_id")
    write_csv(res_full,
              file.path("dge_type5_cross_tissue_cross_phase/full", paste0(contrast_filename, ".csv")))

    sig <- as.data.frame(res) %>%
      rownames_to_column("transcript_id") %>%
      filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
      arrange(padj)
    write_csv(sig,
              file.path("dge_type5_cross_tissue_cross_phase/sig", paste0(contrast_filename, ".csv")))

    type5_summary <- rbind(type5_summary,
                           data.frame(contrast=contrast_name, n_significant=nrow(sig)))
  }
}

cat("\nType 5 complete:", nrow(type5_summary), "comparisons\n")
write_csv(type5_summary %>% arrange(desc(n_significant)),
          "type5_cross_tissue_cross_phase_summary.csv")

# =============================================================================
# === TYPE 5B: CROSS-TISSUE VS INTERMOULT PHASE (DESeq2 ~group) ===
# =============================================================================

cat("\n=== TYPE 5B: Cross-tissue vs intermoult (intermoult reference phase) ===\n")
dir.create("dge_type5b_cross_tissue_vs_intermoult",       showWarnings=FALSE)
dir.create("dge_type5b_cross_tissue_vs_intermoult/full",  showWarnings=FALSE)
dir.create("dge_type5b_cross_tissue_vs_intermoult/sig",   showWarnings=FALSE)

type5b_summary     <- data.frame()
intermoult_groups  <- all_groups[grepl("_intermoult$", all_groups)]

for (test_group in setdiff(all_groups, intermoult_groups)) {
  test_parts  <- strsplit(test_group, "_")[[1]]
  test_tissue <- test_parts[1]
  test_phase  <- paste(test_parts[-1], collapse="_")

  for (int_group in intermoult_groups) {
    int_tissue <- strsplit(int_group, "_")[[1]][1]

    # Cross-tissue only; skip same tissue and skip if test is already intermoult
    if (test_tissue == int_tissue) next
    if (test_phase  == "intermoult") next

    cat("  ", test_tissue, "-", test_phase, "vs", int_tissue, "-intermoult\n")

    res <- tryCatch({
      results(dds, contrast=c("group", test_group, int_group))
    }, error = function(e) {
      cat("    ERROR:", e$message, "\n"); NULL
    })

    if (!is.null(res)) {
      name <- paste0("type5b_", test_tissue, "_", gsub(" ", "", test_phase),
                     "_vs_", int_tissue, "_intermoult")
      contrast_filename <- gsub(" ", "_", name)

      res_full <- as.data.frame(res) %>% rownames_to_column("transcript_id")
      write_csv(res_full,
                file.path("dge_type5b_cross_tissue_vs_intermoult/full", paste0(contrast_filename, ".csv")))

      sig <- as.data.frame(res) %>%
        rownames_to_column("transcript_id") %>%
        filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
        arrange(padj)
      write_csv(sig,
                file.path("dge_type5b_cross_tissue_vs_intermoult/sig", paste0(contrast_filename, ".csv")))

      type5b_summary <- rbind(type5b_summary,
                              data.frame(contrast=name, n_significant=nrow(sig)))
      cat("    Significant DEGs:", nrow(sig), "\n")
    }
  }
}

write_csv(type5b_summary %>% arrange(desc(n_significant)),
          "type5b_intermoult_ref_summary.csv")

# =============================================================================
# === SAVE DESEQ2 OBJECTS ===
# =============================================================================

cat("\n5. Saving DESeq2 objects...\n")
saveRDS(dds,       "dds_group_design.rds")
saveRDS(dds_phase, "dds_phase_design.rds")

# =============================================================================
# === FINAL SUMMARY ===
# =============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("ORIGINAL TYPES:\n")
cat("  Type 1  (phase-vs-phase all tissues):      ", nrow(type1_summary),  "\n")
cat("  Type 1C (phase-vs-phase excluding h.legs): ", nrow(type1c_summary), "\n")
cat("  Type 2  (within-tissue phases):            ", nrow(type2_summary),  "\n")
cat("  Type 3  (tissue-vs-others):                ", nrow(type3_summary),  "\n")
cat("  Type 4  (within-phase tissues):            ", nrow(type4_summary),  "\n")
cat("  Type 5  (cross-tissue-cross-phase):        ", nrow(type5_summary),  "\n")
cat("\nINTERMOULT-REFERENCE TYPES:\n")
cat("  Type 1B (phase vs intermoult all tissues): ", nrow(type1b_summary), "\n")
cat("  Type 1D (phase vs intermoult no h.legs):   ", nrow(type1d_summary), "\n")
cat("  Type 2B (within-tissue vs intermoult):     ", nrow(type2b_summary), "\n")
cat("  Type 5B (cross-tissue vs intermoult):      ", nrow(type5b_summary), "\n")
cat("\nOutput directories:\n")
cat("  - dge_type1_phase_vs_phase_alltissues/\n")
cat("  - dge_type1b_phase_vs_intermoult_alltissues/\n")
cat("  - dge_type1c_phase_vs_phase_no_hlegs/\n")
cat("  - dge_type1d_phase_vs_intermoult_no_hlegs/\n")
cat("  - dge_type2_within_tissue_phases/\n")
cat("  - dge_type2b_within_tissue_vs_intermoult/\n")
cat("  - dge_type3_tissue_vs_allothers/\n")
cat("  - dge_type4_within_phase_tissues/\n")
cat("  - dge_type5_cross_tissue_cross_phase/\n")
cat("  - dge_type5b_cross_tissue_vs_intermoult/\n")
