### Everything until Environment Setup is specific for HUJI CLUSTER environment
#!/bin/bash
#SBATCH --job-name=bowtie_BAM
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --mail-user=user@email.com
#SBATCH --output=/path/to/output_%j.out
#SBATCH --error=/path/to/error_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END

# --- 1. Environment Setup (Requires bowtie2, samtools) ---
. /etc/profile.d/huji-lmod.sh
module load spack miniconda3
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh
conda activate rnaseq # Using the general environment with Bowtie2

# --- 2. USER INPUTS ---
# $1: The path to the directory containing the final CD-HIT FASTA (Assembly folder)
INPUT_DIR="$1"
# $2: The path to the directory containing ALL split R1/R2 trimmed reads
TRIMMED_READS_DIR="$1"/trimmed_PE

if [[ -z "$TRIMMED_READS_DIR" ]]; then
    echo "Usage: sbatch script.sh <assembly_folder> <trimmed_reads_directory>"
    exit 1
fi

# --- 3. Generate Comma-Separated File Lists ---
echo "--- Finding R1 and R2 files in $TRIMMED_READS_DIR ---"
# Find all R1 files, sort them to ensure R2 lines up, and join them with commas.
LEFT_READS=$(ls -1 "${TRIMMED_READS_DIR}"/*R1*.fastq.gz 2>/dev/null | sort | paste -sd, -)
RIGHT_READS=$(ls -1 "${TRIMMED_READS_DIR}"/*R2*.fastq.gz 2>/dev/null | sort | paste -sd, -)

if [[ -z "$LEFT_READS" ]]; then
    echo "CRITICAL ERROR: No R1 trimmed files found in $TRIMMED_READS_DIR. Check file naming pattern (*R1*.fastq.gz)."
    exit 1
fi

# --- 4. Define Paths ---
# Find the final clustered reference FASTA (from CD-HIT)
shopt -s nullglob
REF_FASTA_FILE=("$INPUT_DIR"/cd-hit/*_CD-HIT.fasta)

if [[ ${#REF_FASTA_FILE[@]} -eq 0 ]]; then
    echo "CRITICAL ERROR: Could not find the clustered CD-HIT FASTA file in $INPUT_DIR/cd-hit/."
    exit 1
fi

FINAL_REF="${REF_FASTA_FILE[0]}"
INDEX_BASE="${INPUT_DIR}/blobtools_index/assembly_index"
OUTPUT_DIR="${INPUT_DIR}/blobtools_data"
BAM_OUT="${OUTPUT_DIR}/reads_mapped_to_cdhit.bam"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname "$INDEX_BASE")"

echo "Reference FASTA: $FINAL_REF"
echo "Mapping using the following R1 list: $LEFT_READS"

# --- 5. Bowtie2 Indexing ---
echo "--- Step A: Building Bowtie2 Index ---"
bowtie2-build --threads 16 "$FINAL_REF" "$INDEX_BASE"

# --- 6. Bowtie2 Alignment ---
echo "--- Step B: Aligning Reads (PE - Split Files) ---"
# Bowtie2 handles comma-separated lists directly via -1 and -2
bowtie2 -p 16 --sensitive-local \
    -x "$INDEX_BASE" \
    -1 "$LEFT_READS" \
    -2 "$RIGHT_READS" \
    | samtools view -bS - \
    | samtools sort -@ 16 -o "$BAM_OUT"

# --- 7. BAM Indexing ---
echo "--- Step C: Indexing BAM file ---"
samtools index "$BAM_OUT"

echo "Done. Final BAM file is ready for BlobTools at: $BAM_OUT"