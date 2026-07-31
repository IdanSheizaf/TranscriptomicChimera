### Everything until CONFIGURATION is specific for HUJI CLUSTER environment
#!/bin/bash
#SBATCH --job-name=01_salmon_quant
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --mail-user=user@email.com
#SBATCH --output=/path/to/output_%A_%a.out
#SBATCH --error=/path/to/error_%A_%a.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --array=1-49%10                # <-- CHANGE 49 to your number of samples

# Environment
. /etc/profile.d/huji-lmod.sh
module load spack salmon

# --- CONFIGURATION ---
INDEX_DIR="$1/transcriptome/salmon_index"
READS_DIR="$1/raw_reads/concat/trimmed"                        # Directory with your .fastq.gz files
OUTPUT_BASE="$1/quantification"

# Library Stats (CRITICAL FOR SINGLE-END)
FLD_MEAN=250    # Fragment length mean (adjust if you have Bioanalyzer data)
FLD_SD=25       # Fragment length standard deviation (adjust if you have Bioanalyzer data)

# --- GET SAMPLE FOR THIS ARRAY TASK ---
SAMPLE_FILE=$(ls ${READS_DIR}/*.fastq.gz | sed -n "${SLURM_ARRAY_TASK_ID}p")
SAMPLE_NAME=$(basename "$SAMPLE_FILE" .fastq.gz)

# --- COMMAND ---
echo "========================================"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing Sample: $SAMPLE_NAME"
echo "Input File: $SAMPLE_FILE"
echo "Output Dir: $OUTPUT_BASE/$SAMPLE_NAME"
echo "========================================"

mkdir -p "$OUTPUT_BASE/$SAMPLE_NAME"

salmon quant \
    -i "$INDEX_DIR" \
    -l A \
    -r "$SAMPLE_FILE" \
    --validateMappings \
    --fldMean $FLD_MEAN \
    --fldSD $FLD_SD \
    --gcBias \
    --seqBias \
    --numBootstraps 100 \
    -p $SLURM_CPUS_PER_TASK \
    -o "$OUTPUT_BASE/$SAMPLE_NAME"

echo "Quantification complete for $SAMPLE_NAME"
