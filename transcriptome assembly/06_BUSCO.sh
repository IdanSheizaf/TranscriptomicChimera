#!/bin/bash
#SBATCH --job-name=BUSCO
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=idan.slurm@mail.huji.ac.il
#SBATCH --output=/sci/labs/ariel.chipman/idansh/scripts/logs/BUSCO_%j.out
#SBATCH --error=/sci/labs/ariel.chipman/idansh/scripts/errors/BUSCO_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=5G
#SBATCH --time=2:00:00

# Load conda and activate environment
. /etc/profile.d/huji-lmod.sh
module load spack miniconda3
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh
conda activate busco
# ==== USER SETTINGS ====
LINEAGE=crustacea_odb12
LINEAGE_PATH=/sci/labs/ariel.chipman/idansh/.config/busco/

# ==== Get folder from argument ====
FULL_FOLDER_PATH="$1"

if [[ -z "$FULL_FOLDER_PATH" ]]; then
    echo "No folder provided! Exiting."
    exit 1
fi

CDHIT_DIR="${FULL_FOLDER_PATH}/cd-hit"


# --- FIX: Use Array Assignment and Globbing ---
# 1. Enable nullglob for safety (if no files match, the array is empty instead of containing the literal pattern)
shopt -s nullglob
# 2. Populate array with all files ending in .fasta in the CD-HIT directory
INPUT_FILES=("$CDHIT_DIR"/*.fasta)
# 3. Assign the first file found to INPUT_FILE
INPUT_FILE="${INPUT_FILES[0]}"


# === INPUT CHECK ===
if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "CRITICAL ERROR: No .fasta file found in $CDHIT_DIR. Check CD-HIT output."
    exit 1
elif [[ ${#INPUT_FILES[@]} -gt 1 ]]; then
    echo "WARNING: Multiple .fasta files found in $CDHIT_DIR. Using the first one: $INPUT_FILE"
fi
# === INPUT CHECK END ===


mkdir -p "BUSCO"
	
if [[ -f "$INPUT_FILE" ]]; then
    echo "Input File found: $INPUT_FILE"
    
    # Change directory to the output folder before running BUSCO to prevent clutter
    cd "$FULL_FOLDER_PATH" || { echo "Failed to cd into $FULL_FOLDER_PATH"; exit 1; }

    busco -i "$INPUT_FILE" \
          -m transcriptome \
          -l "$LINEAGE" \
          -o "BUSCO" \
          --offline \
          --download_path "$LINEAGE_PATH" \
          -c 8

    # The 'cd - > /dev/null' command returns to the original directory.
    cd - > /dev/null
else
    # This should be caught by the array check above, but kept as a final safeguard.
    echo "Warning: Input file not found: $INPUT_FILE"
    exit 1
fi