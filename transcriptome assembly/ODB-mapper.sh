### Everything until ENVIRONMENT is specific for HUJI CLUSTER environment
#!/bin/bash
#SBATCH --job-name=odbmapper_batch
#SBATCH --mail-user=user@email.com
#SBATCH --output=/path/to/output_%j.out
#SBATCH --error=/path/to/error_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

if [ $# -lt 1 ]; then
    echo "Error: Working folder required"
    echo "Usage: sbatch ODB-mapper.sh /path/to/working/folder [fastalist_file]"
    exit 1
fi

# Load environment
. /etc/profile.d/huji-lmod.sh
module load spack miniconda3
source /usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh
conda activate orthologer_env

WORK_DIR="$1"
FASTALIST="${2:-$WORK_DIR/fastalist.txt}"

# Move to working directory
cd "$WORK_DIR"

# Run ODB-mapper on all files with the correct parameter order
ODB-mapper MAP de_novo_transcriptomes_crustacea "$FASTALIST" 6657

# Retrieve the results
ODB-mapper RESULT de_novo_transcriptomes_crustacea

echo "ODB-mapper annotation complete."
