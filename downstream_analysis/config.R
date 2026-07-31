# Global Configuration for Isopod RNA-seq Analysis
# Change these settings here to update the entire pipeline.

# --- PROJECT SETTINGS ---
# Set the project root directory dynamically (from environment variable, current working directory, or relative parent)
default_root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (default_root == "" || !dir.exists(default_root)) {
  default_root <- getwd()
  if (basename(default_root) == "downstream_analysis") {
    PROJECT_ROOT <- dirname(default_root)
  } else {
    PROJECT_ROOT <- default_root
  }
} else {
  PROJECT_ROOT <- default_root
}

# --- ACTIVE ANALYSIS ---
# Current species under investigation: "maculatum", "laevis", or "officinalis"
ACTIVE_SPECIES <- "officinalis"

# Default DGE folder to look into
DEFAULT_DGE_FOLDER <- "dge_type6_f.legs_h.legs_within_phase"

# --- SIGNIFICANCE THRESHOLDS ---
PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1.0  # log2FoldChange threshold (1.0 = 2-fold change)

# --- VISUALIZATION DEFAULTS ---
# Standard colors for plots
COLOR_UP   <- "#E41A1C" # Red
COLOR_DOWN <- "#377EB8" # Blue
COLOR_NS   <- "grey80"  # Not Significant
