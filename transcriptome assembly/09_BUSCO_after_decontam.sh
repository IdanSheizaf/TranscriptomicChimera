### Everything until USER SETTINGS is specific for HUJI CLUSTER environment
#!/bin/bash
#
#SBATCH --job-name=BUSCO_decontam
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=user@email.com
#SBATCH --output=/path/to/output_%j.out
#SBATCH --error=/path/to/error_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
#SBATCH --time=20:00:00

set -euo pipefail

echo "=== BUSCO on decontaminated transcriptome started ==="
date

############################################
# 0. ENVIRONMENT
############################################
# Load conda and activate environment
. /etc/profile.d/huji-lmod.sh
module load spack miniconda3
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh
conda activate busco

############################################
# 1. USER SETTINGS
############################################

# 1st argument: the *sample folder* you used as INPUTDIR for 08_mmseqs_decontam_ID_removal.sh
SAMPLE_DIR="$(readlink -f "${1:-}")"
BASE_NAME=$(basename "$SAMPLE_DIR") 

if [[ -z "${SAMPLE_DIR:-}" ]]; then
  echo "No sample folder provided! Exiting."
  echo "Usage: sbatch 09_BUSCO_after_decontam.sh /path/to/sample_folder [lineage_path]"
  exit 1
fi

# BUSCO lineage and offline DB path
LINEAGE="crustacea_odb12"
LINEAGE_PATH="${2:-/path/to/busco_downloads}"

# Decontaminated transcriptome from step 8
DECONTAM_DIR="${SAMPLE_DIR}/decontamination"
INPUT_FASTA="${SAMPLE_DIR}/${BASE_NAME}_decontam.fasta"

OUTDIR="${SAMPLE_DIR}/BUSCO_decontam"

echo "Sample dir:        ${SAMPLE_DIR}"
echo "Decontam dir:      ${DECONTAM_DIR}"
echo "Input FASTA:       ${INPUT_FASTA}"
echo "BUSCO output dir:  ${OUTDIR}"

############################################
# 2. INPUT CHECKS
############################################

if [[ ! -d "${DECONTAM_DIR}" ]]; then
  echo "CRITICAL ERROR: mmseqs_decontam directory not found:"
  echo "  ${DECONTAM_DIR}"
  echo "Run 08_mmseqs_decontam_ID_removal.sh first."
  exit 1
fi

if [[ ! -s "${INPUT_FASTA}" ]]; then
  echo "CRITICAL ERROR: decontaminated FASTA not found or empty:"
  echo "  ${INPUT_FASTA}"
  exit 1
fi

mkdir -p "${OUTDIR}"

############################################
# 3. RUN BUSCO
############################################

cd "${OUTDIR}" || { echo "Failed to cd into ${OUTDIR}"; exit 1; }


echo "Running BUSCO on decontaminated transcriptome..."
busco \
  -i "${INPUT_FASTA}" \
  -m transcriptome \
  -l "${LINEAGE}" \
  -o "BUSCO_decontam" \
  --offline \
  --download_path "${LINEAGE_PATH}" \
  -c 16

echo "BUSCO run finished."
cd - > /dev/null

date
echo "=== BUSCO on decontaminated transcriptome completed ==="
