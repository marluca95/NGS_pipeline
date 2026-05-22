#!/bin/bash
# Paired-Chain Pipeline — Step 01: Paired-end quality trimming (BBDuk)
#
# Submit from the paired_chain_pipeline root:
#   cd /cluster/project/reddy/marluca/NGS_pipeline/paired_chain_pipeline
#   sbatch run/01_run_paired_preprocessing.sbatch

#SBATCH --job-name=01_paired_preprocessing
#SBATCH --output=/cluster/project/reddy/marluca/NGS_pipeline/logs/P3620_LUCA-TCRDMF5-abDMS/01_preprocessed/slurm_%j.out
#SBATCH --error=/cluster/project/reddy/marluca/NGS_pipeline/logs/P3620_LUCA-TCRDMF5-abDMS/01_preprocessed/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=06:00:00

set -euo pipefail

PROJECT_ROOT="/cluster/project/reddy/marluca/NGS_pipeline"
PIPELINE_ROOT="${PROJECT_ROOT}/paired_chain_pipeline"
YAML="${PIPELINE_ROOT}/config/P3620_abDMS/script1_preprocessing_config.yaml"
VENV="${PROJECT_ROOT}/venv/bin/activate"

mkdir -p "${PROJECT_ROOT}/logs/P3620_LUCA-TCRDMF5-abDMS/01_preprocessed"

module load jdk/8u141-b15

cd "${PIPELINE_ROOT}/scripts"
source "${VENV}"

echo "============================================="
echo "Job:    ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo "Step:   01_paired_preprocessing"
echo "Config: ${YAML}"
echo "============================================="

python -u script1_paired_preprocessing.py \
    --yaml_config "${YAML}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
