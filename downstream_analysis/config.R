# Global Configuration for Isopod RNA-seq Analysis
# Change these settings here to update the entire pipeline.

# --- PROJECT SETTINGS ---
# Set the absolute path to your project root
PROJECT_ROOT <- "G:/My Drive/PhD/Projects/Comparative Expressional Changes During the Moult Cycle in Land Isopods"

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
