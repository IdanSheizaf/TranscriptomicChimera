### Everything until USER SETTINGS is specific for HUJI CLUSTER environment
#!/bin/bash
#SBATCH --job-name=TD2
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=user@email.com
#SBATCH --output=/path/to/output_%j.out
#SBATCH --error=/path/to/error_%j.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=270G

# ==== Load environment ====
. /etc/profile.d/huji-lmod.sh
module load spack miniconda3
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh
conda activate TD2_env

# ==== USER SETTINGS ====
SAMPLE_DIR="$1"
MMSEQS_DB_DIR="${2:-/path/to/mmseq_db}"

# ==== Get folder from argument ====


INPUT_FASTA="${SAMPLE_DIR}/trinity_out.Trinity.fasta"
OUTPUT_DIR="${SAMPLE_DIR}/TD2"

if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "Missing input file in $SAMPLE_DIR. Skipping."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || { echo "Could not cd into $OUTPUT_DIR"; exit 1; }

echo "Running TD2.LongOrfs..."
TD2.LongOrfs -t "$INPUT_FASTA" > longorfs.log 2>&1

LONGEST_ORFS_PEP="$OUTPUT_DIR/trinity_out.Trinity/longest_orfs.pep"
if [[ ! -f "$LONGEST_ORFS_PEP" ]]; then
    echo "Missing longest_orfs.pep in $OUTPUT_DIR Skipping MMseqs2."
    exit 1
fi

echo "Running MMseqs2 homology searches..."
mmseqs easy-search "$LONGEST_ORFS_PEP" ${MMSEQS_DB_DIR}/swissprot swissprot_hits.m8 tmp_swissprot -s 7.0 --split-memory-limit 250G
mmseqs easy-search "$LONGEST_ORFS_PEP" ${MMSEQS_DB_DIR}/uniref90 uniref90_hits.m8 tmp_uniref90 -s 7.0 --split-memory-limit 250G
mmseqs easy-search "$LONGEST_ORFS_PEP" ${MMSEQS_DB_DIR}/pfam_full pfam_hits.m8 tmp_pfam -s 7.0 --split-memory-limit 250G
mmseqs easy-search "$LONGEST_ORFS_PEP" ${MMSEQS_DB_DIR}/crustome crustome_hits.m8 tmp_crustome -s 7.0 --split-memory-limit 250G

echo "Concatenating homology search results..."
cat swissprot_hits.m8 uniref90_hits.m8 pfam_hits.m8 crustome_hits.m8 > combined_hits.m8

rm "${OUTPUT_DIR}/crustome_hits.m8"
rm "${OUTPUT_DIR}/pfam_hits.m8"
rm "${OUTPUT_DIR}/uniref90_hits.m8"
rm "${OUTPUT_DIR}/swissprot_hits.m8"
rm -rf "${OUTPUT_DIR}/tmp_crustome"
rm -rf "${OUTPUT_DIR}/tmp_pfam"
rm -rf "${OUTPUT_DIR}/tmp_uniref90"
rm -rf "${OUTPUT_DIR}/tmp_swissprot"

echo "Running TD2.Predict with homology retention..."
TD2.Predict -t "$INPUT_FASTA" --retain-mmseqs-hits combined_hits.m8 > predict.log 2>&1


rm "${OUTPUT_DIR}/combined_hits.m8"


echo "Done with $1."
