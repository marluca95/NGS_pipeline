#!/bin/bash
#SBATCH --account=es_reddy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --mem-per-cpu=16384
#SBATCH --job-name="AA_extract"
#SBATCH --output=/cluster/project/reddy/marluca/NGS_pipeline/logs/%x_%A_%a.out
#SBATCH --error=/cluster/project/reddy/marluca/NGS_pipeline/logs/%x_%A_%a.err
#SBATCH --open-mode=truncate
#SBATCH --array=0-3  # adjust to number of samples minus one

echo "Loading Modules"
module load stack/2024-06 gcc/12.2.0 python/3.9.18

echo "Activating venv"
VENV_PATH="/cluster/project/reddy/marluca/NGS_pipeline/venv"
source "$VENV_PATH/bin/activate"

# Directories
SEQPIPE_DIR="/cluster/project/reddy/marluca/NGS_pipeline/src"
CONFIG_DIR="/cluster/project/reddy/marluca/NGS_pipeline/config"

# --- SAMPLE HANDLING ---
# Minimal identifiers (same as before in preprocessing)
SAMPLES=("59053" "59054" "59055" "59056")
SAMPLE_ID="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

# Config file for this pool/sample
CONFIG_FILE="$CONFIG_DIR/A3_library_AAextraction-config.yaml"

echo "Submitting extraction job for $SAMPLE_ID"
echo "Using config: $CONFIG_FILE"

python "$SEQPIPE_DIR/extract_codons.py" --config "$CONFIG_FILE"

echo "Finished sample $SAMPLE_ID"