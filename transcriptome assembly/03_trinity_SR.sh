#!/bin/bash
#SBATCH --job-name=trinity
#SBATCH --output=/sci/labs/ariel.chipman/idansh/scripts/logs/trinity_%j.out
#SBATCH --error=/sci/labs/ariel.chipman/idansh/scripts/errors/trinity_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=idan.slurm@gmail.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=120:00:00

# --- 1. Environment Setup ---
. /etc/profile.d/huji-lmod.sh
module load spack
module load apptainer

# --- 2. USER SETTINGS ---
BASE_DIR="$(readlink -f "$1")"
DATA_DIR="$BASE_DIR/trimmed"
OUTPUT_DIR="$BASE_DIR/trinity_out"
SAMPLE_NAME="$(basename "$BASE_DIR")"

# --- CRITICAL: Apptainer Path ---
TRINITY_SIF="/sci/labs/ariel.chipman/idansh/shared_software/trinityrnaseq.v2.15.2.simg"

if [[ -z "$DATA_DIR" ]]; then
    echo "Error: No input directory provided."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Input Directory: $DATA_DIR"
echo "Output Directory: $OUTPUT_DIR"

# --- 3. Generate File Lists (Comma Separated) ---
# We look for all R1/R2 trimmed files and sort them to ensure pair synchronization.
FILES=$(ls -1 "$DATA_DIR"/*trimmed.fastq.gz 2>/dev/null | sort | paste -sd, -)


if [[ -z "$FILES" ]]; then
    echo "CRITICAL ERROR: No trimmed R1 files found in $DATA_DIR."
    exit 1
fi


# --- 4. Run Trinity via Apptainer ---
# We bind DATA_DIR (to read files) and OUTPUT_DIR (to write results)
apptainer exec \
    --bind "$BASE_DIR":"$BASE_DIR" \
    "$TRINITY_SIF" Trinity \
    --seqType fq \
    --single "$FILES" \
    --CPU 16 \
    --max_memory 200G \
    --output "$OUTPUT_DIR" \
    --full_cleanup

# --- 5. Rename TRINITY in FASTA headers ---
TRINITY_FASTA="$OUTPUT_DIR.Trinity.fasta"
TEMP_FASTA="$OUTPUT_DIR.Trinity.tmp.fasta"
GENE_MAP="$TRINITY_FASTA.gene_trans_map"
GENE_TEMP="$TRINITY_FASTA.gene_trans_map.tmp"

if [[ -f "$TRINITY_FASTA" ]]; then
    awk -v repl="$SAMPLE_NAME" '
        /^>/ {
            gsub(/TRINITY/, repl, $0)
            print
            next
        }
        { print }
    ' "$TRINITY_FASTA" > "$TEMP_FASTA"

    mv "$TEMP_FASTA" "$TRINITY_FASTA"

    echo "Headers in $TRINITY_FASTA updated to use \"$SAMPLE_NAME\" instead of TRINITY."
else
    echo "WARNING: Trinity FASTA not found at $TRINITY_FASTA, skipping header rename."
fi

if [[ -f "$GENE_MAP" ]]; then
    awk -v repl="$SAMPLE_NAME" '
        /^>/ {
            gsub(/TRINITY/, repl, $0)
            print
            next
        }
        { print }
    ' "$GENE_MAP" > "$GENE_TEMP"

    mv "$GENE_TEMP" "$GENE_MAP"

    echo "Headers in $GENE_MAP updated to use \"$SAMPLE_NAME\" instead of TRINITY."
else
    echo "WARNING: Trinity gene_map not found at $GENE_MAP, skipping header rename."
fi


rm -rf "$OUTPUT_DIR"
echo "TRINITY pipeline finished, and removed all temp files at " "$OUTPUT_DIR"