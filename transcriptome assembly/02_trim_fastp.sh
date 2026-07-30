#!/bin/bash
#SBATCH --job-name=trimming_batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/sci/labs/ariel.chipman/idansh/scripts/logs/trim_%j.out
#SBATCH --error=/sci/labs/ariel.chipman/idansh/scripts/errors/trim_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=idan.slurm@gmail.com


INPUT_DIR="$1"

if [[ -z "$INPUT_DIR" ]]; then
    echo "Error: No input directory provided."
    exit 1
fi

. /etc/profile.d/huji-lmod.sh
module load spack miniconda3
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh
conda activate rnaseq

OUT_DIR="$INPUT_DIR/trimmed"
mkdir -p "$OUT_DIR"

echo "Processing directory: $INPUT_DIR"

# Enable nullglob
shopt -s nullglob
FILES=("$INPUT_DIR"/raw_reads/*.fastq.gz)

# Check if files exist
if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: No .fastq.gz files found in $INPUT_DIR"
    exit 1
fi

for FILE in "${FILES[@]}"; do
    # Define Sample Name (Remove .fastq.gz extension)
    SAMPLE_NAME=$(basename "$FILE" .fastq.gz)

    echo "  > Trimming Sample: $SAMPLE_NAME"

    # Single-End fastp command (removed -I, -O, and PE detection)
    fastp \
      -i "$FILE" \
      -o "$OUT_DIR/${SAMPLE_NAME}_trimmed.fastq.gz" \
      --thread 8 \
      --html "$OUT_DIR/${SAMPLE_NAME}_fastp.html" \
      --json "$OUT_DIR/${SAMPLE_NAME}_fastp.json"
done