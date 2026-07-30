#!/bin/bash
#SBATCH --job-name=dge_complete
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --mail-user=idan.slurm@gmail.com
#SBATCH --output=/sci/labs/ariel.chipman/idansh/scripts/logs/deseq2_%j.out
#SBATCH --error=/sci/labs/ariel.chipman/idansh/scripts/errors/deseq2_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END

. /etc/profile.d/huji-lmod.sh
module load spack R4/4.4.1

PROJECT_DIR="$1"

# Create and navigate to analysis folder
mkdir -p "$PROJECT_DIR/dge_analysis"
cd "$PROJECT_DIR/dge_analysis"

echo "Starting DGE analysis..."
echo "Project directory: $PROJECT_DIR"
echo "Working directory: $(pwd)"

# Pass project directory to R script
Rscript /sci/labs/ariel.chipman/idansh/scripts/rnaseq/new_protocol/run_deseq2.R "$PROJECT_DIR"

echo ""
echo "Analysis complete!"
