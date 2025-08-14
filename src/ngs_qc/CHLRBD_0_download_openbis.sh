#!/bin/bash

#SBATCH --account=es_reddy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=16384
#SBATCH --job-name="download_openbis"
#SBATCH --mem-per-cpu=16384
#SBATCH --output="download_openbis_out"
#SBATCH --error="download_openbis_error"
#SBATCH --open-mode=truncate

echo "Loading Modules"

module load stack/2024-06 gcc/12.2.0  bbmap/39.01
module load stack/2024-06 gcc/12.2.0 python/3.9.18
module load stack/2024-06 gcc/12.2.0 openjdk/21.0.3_9

echo "Loading Venv"

VENV_PATH="${VENV_PATH:-/cluster/project/reddy/rakuhn/plmfit_venv}"
source "$VENV_PATH/bin/activate"

SEQPIPE_DIR="${SEQPIPE_DIR:-/cluster/project/reddy/rakuhn/seq-pipeline}"
OUT_DIR="${1:-data}"
POOL_ID="${2:-/BSSE_Reddy/POOL-1008}"

echo "Starting Download"

python "$SEQPIPE_DIR/openbis.py" --path_to_output "$OUT_DIR" --pool_id "$POOL_ID"
