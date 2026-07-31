### Everything until CONFIGURATION is specific for HUJI CLUSTER environment
#!/bin/bash
#SBATCH --job-name=ko_to_go
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --mail-user=user@email.com
#SBATCH --output=/path/to/output_%j.out
#SBATCH --error=/path/to/error_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END

# ============================================================================
# Enrich combined_annotations.tsv with GO terms from KEGG KO
# ============================================================================
# Usage: sbatch 03_KO_2_GO.sh /path/to/project/folder [script_dir]
# ============================================================================

# 1. Environment Setup
. /etc/profile.d/huji-lmod.sh
module load spack python 
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh
conda activate base

set -e

if [ $# -lt 1 ]; then
    echo "Error: Project folder required"
    echo "Usage: sbatch 03_KO_2_GO.sh /path/to/project/folder [script_dir]"
    exit 1
fi

# ============================================================================
# CONFIGURATION
# ============================================================================

SCRIPT_DIR="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
INPUT_FILE="$1/transcriptome/combined_annotations.tsv"
OUTPUT_FILE="${INPUT_FILE%.tsv}_with_GO_from_KO.tsv"
PROGRESS_FILE="${INPUT_FILE%.tsv}_ko2go_progress.txt"

# ============================================================================
# SLURM ENVIRONMENT
# ============================================================================

echo "========================================="
echo "SLURM Job Information"
echo "========================================="
echo "Job ID:        $SLURM_JOB_ID"
echo "Job Name:      $SLURM_JOB_NAME"
echo "Node:          $SLURMD_NODENAME"
echo "Start Time:    $(date)"
echo "========================================="
echo ""

# Find Python
PYTHON_CMD=$(which python3 2>/dev/null || which python 2>/dev/null || echo "")

if [ -z "$PYTHON_CMD" ]; then
    echo "Error: Python not found in PATH"
    exit 1
fi

echo "Using Python: $PYTHON_CMD"
$PYTHON_CMD --version
echo ""

# ============================================================================
# CHECK INPUT FILE
# ============================================================================

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File not found: $INPUT_FILE"
    exit 1
fi

if [ ! -f "${SCRIPT_DIR}/ko_to_go.py" ]; then
    echo "Error: ko_to_go.py not found in $SCRIPT_DIR"
    exit 1
fi

echo "Input file:  $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "Python script: ${SCRIPT_DIR}/ko_to_go.py"
echo ""

# ============================================================================
# INSTALL DEPENDENCIES
# ============================================================================

echo "Installing Python dependencies..."
$PYTHON_CMD -m pip install --user requests --quiet --no-warn-script-location 2>/dev/null || true
echo "Dependencies installed."
echo ""

# ============================================================================
# RUN ENRICHMENT
# ============================================================================

echo "Starting KO-to-GO enrichment..."
echo "This may take 1-4 hours depending on the number of unique KO terms."
echo ""

$PYTHON_CMD "${SCRIPT_DIR}/ko_to_go.py" "$INPUT_FILE" "$OUTPUT_FILE" "$PROGRESS_FILE"

EXITCODE=$?

# ============================================================================
# SUMMARY
# ============================================================================

if [ $EXITCODE -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "SUCCESS!"
    echo "========================================="
    echo ""
    echo "Enriched annotation file created:"
    echo "  $OUTPUT_FILE"
    echo ""
    echo "File statistics:"
    echo "  Original lines: $(wc -l < "$INPUT_FILE")"
    echo "  Enriched lines: $(wc -l < "$OUTPUT_FILE")"
    echo "  File size:      $(ls -lh "$OUTPUT_FILE" | awk '{print $5}')"
    echo ""
    echo "========================================="
    echo "Job completed at: $(date)"
    echo "========================================="
    
    rm -f "$PROGRESS_FILE"
else
    echo ""
    echo "========================================="
    echo "ERROR: Enrichment failed"
    echo "========================================="
    exit 1
fi
