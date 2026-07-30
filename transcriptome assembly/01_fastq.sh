#!/bin/bash
#SBATCH --job-name=fastqc_batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/sci/labs/ariel.chipman/idansh/scripts/logs/FASTQC_%j.out
#SBATCH --error=/sci/labs/ariel.chipman/idansh/scripts/errors/FASTQC_%j.err
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=idan.slurm@gmail.com

# 1. Input: Path to directory containing fastq files
INPUT_DIR="$1"

if [[ -z "$INPUT_DIR" ]]; then
    echo "Error: No input directory provided."
    exit 1
fi

# 2. Environment
. /etc/profile.d/huji-lmod.sh
module load spack miniconda3 fastqc
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh


OUT_DIR="$INPUT_DIR/fastqc"
mkdir -p $OUT_DIR

echo "Processing directory: $INPUT_DIR/raw_reads"

# 3. Run FastQC on all fastq.gz files in the dir
# -t 8 means it processes 8 files simultaneously
fastqc -t 8 -o "$OUT_DIR" "$INPUT_DIR"/raw_reads/*.fastq.gz