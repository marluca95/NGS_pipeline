#!/bin/bash
#SBATCH --job-name=umi_process
#SBATCH --output=/cluster/project/reddy/marluca/NGS_pipeline/logs/umi_process_%A_%a.out
#SBATCH --error=/cluster/project/reddy/marluca/NGS_pipeline/logs/umi_process_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-3  # adjust based on number of samples

module load stack/2024-06 gcc/12.2.0 python/3.9.18
source /cluster/project/reddy/marluca/NGS_pipeline/venv/bin/activate

# Directory with sample-specific YAML configs
CONFIG_DIR="/cluster/project/reddy/marluca/NGS_pipeline/config"

# Gather all configs into an array
CONFIG_FILES=(${CONFIG_DIR}/*_umi_config.yaml)
NUM_SAMPLES=${#CONFIG_FILES[@]}

if [ ${SLURM_ARRAY_TASK_ID} -ge ${NUM_SAMPLES} ]; then
    echo "Error: SLURM_ARRAY_TASK_ID (${SLURM_ARRAY_TASK_ID}) >= NUM_SAMPLES (${NUM_SAMPLES})"
    exit 1
fi

# Select YAML for this array job
CONFIG_FILE=${CONFIG_FILES[$SLURM_ARRAY_TASK_ID]}

echo "Running UMI processing for config: ${CONFIG_FILE}"

# Call Python script
python /cluster/project/reddy/marluca/NGS_pipeline/src/ngs_qc/UMI-consensus2.py \
    --config ${CONFIG_FILE} --top_k 1 3 --min_count_for_report 1 2

echo "Finished ${CONFIG_FILE}"