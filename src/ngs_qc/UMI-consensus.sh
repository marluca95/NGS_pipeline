#!/bin/bash
#SBATCH --job-name=umi_process
#SBATCH --output=/cluster/project/reddy/marluca/NGS_pipeline/jobs/logs/umi_process_%A_%a.out
#SBATCH --error=/cluster/project/reddy/marluca/NGS_pipeline/jobs/logs/umi_process_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-9  # Adjust based on number of samples

# Load required modules / activate venv
module load python/3.11
source /cluster/project/reddy/marluca/NGS_pipeline/venv/bin/activate

# Base output directory containing per-sample folders
BASE_OUTPUT_DIR="/cluster/project/reddy/marluca/NGS_pipeline/output"

# Directory containing YAML configs for each sample
CONFIG_DIR="/cluster/project/reddy/marluca/NGS_pipeline/config/umi"

# List all sample directories (one per sample)
SAMPLE_DIRS=($(ls -d ${BASE_OUTPUT_DIR}/GFB-*))
NUM_SAMPLES=${#SAMPLE_DIRS[@]}

# Safety check: ensure the array index does not exceed number of samples
if [ ${SLURM_ARRAY_TASK_ID} -ge $NUM_SAMPLES ]; then
    echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} exceeds number of samples (${NUM_SAMPLES}), exiting."
    exit 1
fi

# Pick the sample directory corresponding to this array task
SAMPLE_DIR=${SAMPLE_DIRS[$SLURM_ARRAY_TASK_ID]}
SAMPLE_NAME=$(basename ${SAMPLE_DIR})

# Combined FASTQ file inside this sample folder
SAMPLE_FILE=$(ls ${SAMPLE_DIR}/*_combined_final.fastq.gz)

# Corresponding YAML config
YAML_CONFIG="${CONFIG_DIR}/${SAMPLE_NAME}_umi_config.yaml"

# Create output folder for UMI processing
OUTPUT_DIR="${SAMPLE_DIR}/umi"
mkdir -p ${OUTPUT_DIR}

echo "Processing sample: ${SAMPLE_NAME}"
echo "Input FASTQ: ${SAMPLE_FILE}"
echo "Config: ${YAML_CONFIG}"
echo "Output directory: ${OUTPUT_DIR}"

# Run the UMI Python script
python /cluster/project/reddy/marluca/NGS_pipeline/src/ngs_qc/umi_process.py \
    --yaml_config ${YAML_CONFIG} \
    --input_fastq ${SAMPLE_FILE} \
    --output_dir ${OUTPUT_DIR} \
    --sample ${SAMPLE_NAME}

echo "Finished processing ${SAMPLE_NAME}"