### Everything until USER SETTINGS is specific for HUJI CLUSTER environment
#!/bin/bash
#SBATCH --job-name=trinity
#SBATCH --mail-user=user@email.com
#SBATCH --output=/path/to/output_%j.out
#SBATCH --error=/path/to/error_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
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
DATA_DIR="$BASE_DIR/trimmed_PE"
OUTPUT_DIR="$BASE_DIR/trinity_out"
SAMPLE_NAME="$(basename "$BASE_DIR")"

# --- Apptainer Path ---
TRINITY_SIF="${2:-/path/to/trinityrnaseq.simg}"

if [[ -z "$DATA_DIR" ]]; then
    echo "Error: No input directory provided."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Input Directory: $DATA_DIR"
echo "Output Directory: $OUTPUT_DIR"

# --- 3. Generate File Lists (Comma Separated) ---
# We look for all R1/R2 trimmed files and sort them to ensure pair synchronization.
R1_FILES=$(ls -1 "$DATA_DIR"/*_R1*trimmed.fastq.gz 2>/dev/null | sort | paste -sd, -)
R2_FILES=$(ls -1 "$DATA_DIR"/*_R2*trimmed.fastq.gz 2>/dev/null | sort | paste -sd, -)

if [[ -z "$R1_FILES" ]]; then
    echo "CRITICAL ERROR: No trimmed R1 files found in $DATA_DIR."
    exit 1
fi

# Safety Check: Compare file counts (by counting commas + 1)
COUNT_R1=$(echo "$R1_FILES" | tr ',' '\n' | wc -l)
COUNT_R2=$(echo "$R2_FILES" | tr ',' '\n' | wc -l)

if [[ "$COUNT_R1" != "$COUNT_R2" ]]; then
    echo "CRITICAL ERROR: File counts mismatch! R1=$COUNT_R1, R2=$COUNT_R2. Cannot run assembly."
    exit 1
fi

echo "Detected $COUNT_R1 synchronized pairs. Starting Assembly..."

# --- 4. Run Trinity via Apptainer ---
# We bind DATA_DIR (to read files) and OUTPUT_DIR (to write results)
apptainer exec \
    --bind "$BASE_DIR":"$BASE_DIR" \
    "$TRINITY_SIF" Trinity \
    --seqType fq \
    --left "$R1_FILES" \
    --right "$R2_FILES" \
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