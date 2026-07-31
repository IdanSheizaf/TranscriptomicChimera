### Everything until CONFIGURATION is specific for HUJI CLUSTER environment
#!/bin/bash
#SBATCH --job-name=dge_complete
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --mail-user=user@email.com
#SBATCH --output=/path/to/output_%j.out
#SBATCH --error=/path/to/error_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END

. /etc/profile.d/huji-lmod.sh
module load spack R4/4.4.1

set -e

if [ $# -lt 1 ]; then
    echo "Error: Project folder required"
    echo "Usage: sbatch 04_dge.sh /path/to/project/folder [script_dir]"
    exit 1
fi

# --- CONFIGURATION ---
PROJECT_DIR="$1"
SCRIPT_DIR="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# Create and navigate to analysis folder
mkdir -p "$PROJECT_DIR/dge_analysis"
cd "$PROJECT_DIR/dge_analysis"

echo "Starting DGE analysis..."
echo "Project directory: $PROJECT_DIR"
echo "Working directory: $(pwd)"

# Pass project directory to R script
Rscript "$SCRIPT_DIR/run_deseq2.R" "$PROJECT_DIR"

echo ""
echo "Analysis complete!"
