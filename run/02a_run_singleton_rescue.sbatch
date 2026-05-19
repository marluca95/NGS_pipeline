#!/bin/bash
# =============================================================================
# Step 02a — Singleton Rescue (Script2a_singleton_rescue.py)
# =============================================================================
# Runs after step 02. Processes all samples in a single job.
#
# Submit:
#   cd /cluster/project/reddy/katja/NGS_pipeline
#   sbatch run/02a_run_singleton_rescue.sbatch
# =============================================================================
#SBATCH --job-name=02a_singleton_rescue
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/02a_singleton_rescue/slurm_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/02a_singleton_rescue/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00

set -euo pipefail

PROJECT_ROOT="/cluster/project/reddy/katja/NGS_pipeline"
YAML="${PROJECT_ROOT}/config/TCRA3_KH157/script2a_singleton_rescue.yaml"
VENV="${PROJECT_ROOT}/venv/bin/activate"

# ---------------------------------------------------------------------------
mkdir -p "${PROJECT_ROOT}/logs/02a_singleton_rescue"
cd "${PROJECT_ROOT}"
source "${VENV}"

echo "============================================="
echo "Job:   $SLURM_JOB_ID"
echo "Node:  $(hostname)"
echo "Start: $(date)"
echo "Step:  02a_singleton_rescue"
echo "Config: ${YAML}"
echo "============================================="

python -u scripts/script2a_singleton_rescue.py \
    --yaml_config "${YAML}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
