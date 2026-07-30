#!/bin/bash
#SBATCH --job-name=01_salmon_index
#SBATCH --cpus-per-task=8
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-user=idan.slurm@gmail.com
#SBATCH --output=/sci/labs/ariel.chipman/idansh/scripts/logs/salmon_index_%j.out
#SBATCH --error=/sci/labs/ariel.chipman/idansh/scripts/errors/salmon_index_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END

# 2. Environment
. /etc/profile.d/huji-lmod.sh
module load spack salmon

# --- CONFIGURATION ---
INPUT_FASTA="$1/transcriptome/transcriptome_decontam.fasta"   # Your clustered assembly
INDEX_DIR="$1/transcriptome/salmon_index"            		# Output directory for index

# --- COMMAND ---
echo "Starting Salmon Indexing..."
echo "Input: $INPUT_FASTA"
echo "Output Dir: $INDEX_DIR"

salmon index \
    -t "$INPUT_FASTA" \
    -i "$INDEX_DIR" \
    -k 31 \
    --threads $SLURM_CPUS_PER_TASK

echo "Indexing complete."
