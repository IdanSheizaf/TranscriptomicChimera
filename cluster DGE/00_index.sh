### Everything until CONFIGURATION is specific for HUJI CLUSTER environment
#!/bin/bash
#SBATCH --job-name=01_salmon_index
#SBATCH --cpus-per-task=8
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-user=user@email.com
#SBATCH --output=/path/to/log_%j.out
#SBATCH --error=/path/to/error_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END

# --- Environment ---
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
