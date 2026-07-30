###specific for SLURM managers, environment is user-specific and related to CLUSTER envirnment.
#!/bin/bash
#SBATCH --job-name=full_protein_annotations
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G  
#SBATCH --time=12:00:00
#SBATCH --output=/path/to/log.out
#SBATCH --error=/path/to/error.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=user@email.com

set -e
INPUT_DIR=$1
PROTEOME_FASTA=${INPUT_DIR}/mmseqsdecontamproteome.decontam.fasta
OUTPUT_DIR=${INPUT_DIR}/annotations

[[ ! -f $PROTEOME_FASTA ]] && { echo "ERROR: $PROTEOME_FASTA missing"; exit 1; }

# Environment
. /etc/profile.d/huji-lmod.sh
module load spack miniconda3 
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh

# Activate conda environment
conda activate bio_env

# Source annotation tool environment variables
source /sci/labs/ariel.chipman/idansh/env_setup.sh

mkdir -p $OUTPUT_DIR/{eggnog,diamond,signalp6,logs}

BASENAME=$(basename $PROTEOME_FASTA .fasta)
echo "=== 1. EggNOG-mapper ==="
emapper.py -i $PROTEOME_FASTA --output $BASENAME --outputdir $OUTPUT_DIR/eggnog \
  -m diamond --cpu $SLURM_CPUS_PER_TASK --datadir $EGGNOGDATADIR \
  --taxscope auto --go_evidence all --target_orthologs all [...]


echo "=== 2. SignalP6 ==="
signalp6 --fastafile $PROTEOME_FASTA --organism euk --mode fast \
  --outputdir $OUTPUT_DIR/signalp6

echo "=== 3. MERGE ALL ? combinedannotations.tsv ==="
# ===== YOUR 04_annotation_merger.sh LOGIC HERE =====
# Copy-paste the core merging/joining code from 04_annotation_merger.sh
# Example (adapt from your file):
cd $INPUT_DIR
# Assume 04 has function like: merge_annotations proteome_dir
# e.g.:
# join_eggnog_diamond_signalp.py $OUTPUT_DIR  # Or bash joins
# Produces: transcriptomecombinedannotations.tsv

echo "=== SUCCESS: transcriptomecombinedannotations.tsv ready for KO2GO ==="
ls -lh transcriptomecombinedannotations.tsv
