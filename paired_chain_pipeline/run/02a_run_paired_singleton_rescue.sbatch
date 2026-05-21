#!/bin/bash
# Paired-Chain Pipeline — Step 02a: Paired singleton rescue
# Runs after step 02 has completed for all samples.
#
# Submit from the paired_chain_pipeline root:
#   cd /cluster/project/reddy/marluca/NGS_pipeline/paired_chain_pipeline
#   sbatch run/02a_run_paired_singleton_rescue.sbatch

#SBATCH --job-name=02a_paired_singleton_rescue
#SBATCH --output=/cluster/project/reddy/marluca/NGS_pipeline/logs/P3620_LUCA-TCRDMF5-abDMS/02a_singleton_rescue/slurm_%j.out
#SBATCH --error=/cluster/project/reddy/marluca/NGS_pipeline/logs/P3620_LUCA-TCRDMF5-abDMS/02a_singleton_rescue/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00

set -euo pipefail

PROJECT_ROOT="/cluster/project/reddy/marluca/NGS_pipeline"
PIPELINE_ROOT="${PROJECT_ROOT}/paired_chain_pipeline"
YAML="${PIPELINE_ROOT}/config/P3620_abDMS/script2a_paired_singleton_rescue_config.yaml"
VENV="${PROJECT_ROOT}/venv/bin/activate"

mkdir -p "${PROJECT_ROOT}/logs/P3620_LUCA-TCRDMF5-abDMS/02a_singleton_rescue"

cd "${PIPELINE_ROOT}/scripts"
source "${VENV}"

echo "============================================="
echo "Job:    ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo "Step:   02a_paired_singleton_rescue"
echo "Config: ${YAML}"
echo "============================================="

python -u script2a_paired_singleton_rescue.py \
    --yaml_config "${YAML}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
