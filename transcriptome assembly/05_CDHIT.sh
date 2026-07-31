### Everything until Environment is specific for HUJI CLUSTER environment
#!/bin/bash
#SBATCH --job-name=CDHIT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --output=/path/to/output_%j.out
#SBATCH --error=/path/to/error_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=user@email.com


# ==== Load environment (customize as needed) ====
. /etc/profile.d/huji-lmod.sh
module load spack miniconda3
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh
conda activate cd-hit_env

# ==== USER SETTINGS ====
CDHIT_IDENTITY=0.95
CDHIT_THREADS=8

# ==== Get folder from argument ====
INPUT_DIR="$1" # The argument passed by submit_all.sh is $1

if [[ -z "$INPUT_DIR" ]]; then
  echo "No folder provided! Exiting."
  exit 1
fi

# 1. Dynamically find the input CDS file using globbing
# This ensures we find the file regardless of the Trinity assembly's casing (Trinity.fasta or trinity.fasta).
shopt -s nullglob
CDS_FILES=("$INPUT_DIR"/TD2/*.cds)

if [[ ${#CDS_FILES[@]} -eq 0 ]]; then
  echo "CRITICAL ERROR: No final CDS file (*.transdecoder.cds) found in $INPUT_DIR. Check TD2 output."
  exit 1
fi

# Assuming TD2 refinement was run on one Trinity assembly, there should only be one CDS file.
INPUT_FILE="${CDS_FILES[0]}"


# 2. Set up paths
# Get the leaf directory name (e.g., 'trinity_output') for the output file prefix
FINAL_NAME=$(basename "$INPUT_DIR") 

OUTPUT_DIR="${INPUT_DIR}/cd-hit"
OUTPUT_PREFIX="${OUTPUT_DIR}/${FINAL_NAME}_CD-HIT.fasta"

# Make output directory
mkdir -p "$OUTPUT_DIR"

echo "Input File: $INPUT_FILE"
echo "Clustering Assembly: $FINAL_NAME"
echo "Output will be saved to: $OUTPUT_DIR"


# 3. Run CD-HIT
# -c 0.95: 95% sequence identity threshold
# -T 8: Use 8 threads
# -M 0: Use all available memory (recommended for large assemblies)
cd-hit-est -i "$INPUT_FILE" \
           -o "$OUTPUT_PREFIX" \
           -c "$CDHIT_IDENTITY" \
           -T "$CDHIT_THREADS" \
           -M 0 

echo "CD-HIT clustering complete. Final file: $OUTPUT_PREFIX"