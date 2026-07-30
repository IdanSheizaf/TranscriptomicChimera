#!/bin/bash
#SBATCH --job-name=trim_PE
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=/sci/labs/ariel.chipman/idansh/scripts/logs/trim_pe_%j.out
#SBATCH --error=/sci/labs/ariel.chipman/idansh/scripts/errors/trim_pe_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=idan.slurm@gmail.com

INPUT_DIR="$1/raw_reads"

if [[ -z "$INPUT_DIR" ]]; then
    echo "Error: No input directory provided."
    exit 1
fi

. /etc/profile.d/huji-lmod.sh
module load spack miniconda3
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh
conda activate rnaseq

OUT_DIR="$1/trimmed_PE"
mkdir -p "$OUT_DIR"

echo "Processing Paired-End files in: $INPUT_DIR"

# Enable nullglob to handle empty directories gracefully
shopt -s nullglob
R1_FILES=("$INPUT_DIR"/*_R1.fastq.gz)

if [ ${#R1_FILES[@]} -eq 0 ]; then
    echo "ERROR: No files ending in *_R1.fastq.gz found in $INPUT_DIR"
    exit 1
fi

for R1_FILE in "${R1_FILES[@]}"; do
    # 1. Identify the R2 partner file
    R2_FILE="${R1_FILE/_R1.fastq.gz/_R2.fastq.gz}"
    
    # 2. Get the clean Sample Name
    SAMPLE_NAME=$(basename "$R1_FILE" | sed 's/_R1.fastq.gz//')

    # 3. Check if R2 exists
    if [[ ! -f "$R2_FILE" ]]; then
        echo "WARNING: R2 file missing for $SAMPLE_NAME. Skipping..."
        continue
    fi

    echo "  > Trimming Pair: $SAMPLE_NAME"

    # 4. Run fastp in PAIRED mode
    # -i/-I : Input R1/R2
    # -o/-O : Output R1/R2
    # --detect_adapter_for_pe : Auto-detect PE adapters
    
    fastp \
      -i "$R1_FILE" \
      -I "$R2_FILE" \
      -o "$OUT_DIR/${SAMPLE_NAME}_R1_trimmed.fastq.gz" \
      -O "$OUT_DIR/${SAMPLE_NAME}_R2_trimmed.fastq.gz" \
      --detect_adapter_for_pe \
      --thread 8 \
      --html "$OUT_DIR/${SAMPLE_NAME}_fastp.html" \
      --json "$OUT_DIR/${SAMPLE_NAME}_fastp.json"
done