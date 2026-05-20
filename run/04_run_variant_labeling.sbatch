#!/bin/bash
# =============================================================================
# Step 04 — Variant Labeling and Enrichment (script4_variant_labeling.py)
# =============================================================================
# Runs after step 03. Processes all peptide groups in a single job.
#
# Submit:
#   cd /cluster/project/reddy/katja/NGS_pipeline
#   sbatch run/04_run_variant_labeling.sbatch
# =============================================================================
#SBATCH --job-name=04_variant_labeling
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/04_variant_labeling/slurm_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/04_variant_labeling/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00

set -euo pipefail

PROJECT_ROOT="/cluster/project/reddy/katja/NGS_pipeline"
YAML="${PROJECT_ROOT}/config/04_variant_labeling.yaml" # make sure to sort out config structure and point to correct YAML
VENV="${PROJECT_ROOT}/venv/bin/activate"

# ---------------------------------------------------------------------------
mkdir -p "${PROJECT_ROOT}/logs/04_variant_labeling"
cd "${PROJECT_ROOT}"
source "${VENV}"

echo "============================================="
echo "Job:   $SLURM_JOB_ID"
echo "Node:  $(hostname)"
echo "Start: $(date)"
echo "Step:  04_variant_labeling"
echo "Config: ${YAML}"
echo "============================================="

python -u scripts/script4_variant_labeling.py \
    --yaml_config "${YAML}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
