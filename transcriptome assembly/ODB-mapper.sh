#!/bin/bash
#SBATCH --job-name=odbmapper_batch
#SBATCH --output=odbmapper_%j.out
#SBATCH --error=odbmapper_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# Load environment
. /etc/profile.d/huji-lmod.sh
module load spack miniconda3
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh
conda activate orthologer_env

# Move to working directory
cd /sci/labs/ariel.chipman/idansh/transcriptomes/de-novo/orthologer

# Create a file listing all FASTA files 

# Run ODB-mapper on all files with the correct parameter order

ODB-mapper MAP de_novo_transcriptomes_crustacea /sci/labs/ariel.chipman/idansh/transcriptomes/de-novo/orthologer/fastalist.txt 6657

# Retrieve the results
ODB-mapper RESULT de_novo_transcriptomes_crustacea

echo "ODB-mapper annotation complete."
