#!/bin/bash

#SBATCH --account=es_reddy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16384
#SBATCH --job-name="A3-CDR3_preprocessing"
#SBATCH --mem-per-cpu=16384
#SBATCH --output=/cluster/project/reddy/marluca/NGS_pipeline/logs/%x_%A_%a.out
#SBATCH --error=/cluster/project/reddy/marluca/NGS_pipeline/logs/%x_%A_%a.err
#SBATCH --open-mode=truncate
#SBATCH --array=0-3  # number of samples minus one

# Minimal sample identifiers
SAMPLES=("59053" "59054" "59055" "59056")

# Pick the current sample based on SLURM_ARRAY_TASK_ID
SAMPLE_PREFIX="GFB-${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

echo "Loading Modules"
module load stack/2024-06 gcc/12.2.0 bbmap/39.01
module load stack/2024-06 gcc/12.2.0 python/3.9.18
module load stack/2024-06 gcc/12.2.0 openjdk/21.0.3_9

echo "Loading Venv"
VENV_PATH="${VENV_PATH:-/cluster/project/reddy/marluca/NGS_pipeline/venv}"
source "$VENV_PATH/bin/activate"

SEQPIPE_DIR="/cluster/project/reddy/marluca/NGS_pipeline/src/ngs_qc"
CONFIG_DIR="/cluster/project/reddy/marluca/NGS_pipeline/config"

CONFIG_FILE="$CONFIG_DIR/A3_CDR3_library-config_${SAMPLES[$SLURM_ARRAY_TASK_ID]}.yaml"


echo "Submitting Job"
echo "Using config: $CONFIG_FILE"

python "$SEQPIPE_DIR/pre_process2.py" --yaml_config "$CONFIG_FILE" --sample "$SAMPLE_PREFIX" 