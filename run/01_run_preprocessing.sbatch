#!/bin/bash
# =============================================================================
# Step 01 — Quality Trimming (script1_preprocessing.py)
# =============================================================================
# Submit:
#   cd /cluster/project/reddy/katja/NGS_pipeline
#   sbatch run/01_run_preprocessing.sbatch
# =============================================================================
#SBATCH --job-name=01_preprocessing
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/01_preprocessing/slurm_preprocessing_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/01_preprocessing/slurm_preprocessing_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=04:00:00

set -euo pipefail


PROJECT_ROOT="/cluster/project/reddy/katja/NGS_pipeline"
YAML="${PROJECT_ROOT}/config/TCRA3_comblib/script1_preprocessing_config.yaml"
VENV="${PROJECT_ROOT}/venv/bin/activate"

# ---------------------------------------------------------------------------
echo "============================================="
echo "Job:   $SLURM_JOB_ID"
echo "Node:  $(hostname)"
echo "Start: $(date)"
echo "Step:  01_preprocessing"
echo "Config: ${YAML}"
echo "============================================="

mkdir -p "${PROJECT_ROOT}/logs/01_preprocessing"
cd "${PROJECT_ROOT}"
source "${VENV}"

python -u scripts/script1_preprocessing.py \
    --yaml_config "${YAML}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
