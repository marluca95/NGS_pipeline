#!/bin/bash
# Paired-Chain Pipeline — Step 04: Paired variant labeling
# Runs after step 03 has completed for all samples.
#
# Submit from the paired_chain_pipeline root:
#   cd /cluster/project/reddy/marluca/NGS_pipeline/paired_chain_pipeline
#   sbatch run/04_run_paired_variant_labeling.sbatch

#SBATCH --job-name=04_paired_variant_labeling
#SBATCH --output=/cluster/project/reddy/marluca/NGS_pipeline/logs/P3620_LUCA-TCRDMF5-abDMS/04_variant_labeling/slurm_%j.out
#SBATCH --error=/cluster/project/reddy/marluca/NGS_pipeline/logs/P3620_LUCA-TCRDMF5-abDMS/04_variant_labeling/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=02:00:00

set -euo pipefail

PROJECT_ROOT="/cluster/project/reddy/marluca/NGS_pipeline"
PIPELINE_ROOT="${PROJECT_ROOT}/paired_chain_pipeline"
YAML="${PIPELINE_ROOT}/config/P3620_abDMS/script4_paired_variant_labeling_config.yaml"
VENV="${PROJECT_ROOT}/venv/bin/activate"

mkdir -p "${PROJECT_ROOT}/logs/P3620_LUCA-TCRDMF5-abDMS/04_variant_labeling"

cd "${PIPELINE_ROOT}/scripts"
source "${VENV}"

echo "============================================="
echo "Job:    ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo "Step:   04_paired_variant_labeling"
echo "Config: ${YAML}"
echo "============================================="

python -u script4_paired_variant_labeling.py \
    --yaml_config "${YAML}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
