#!/bin/bash
# =============================================================================
# Step 04 — KH157 Variant Labeling with tcr_population
# =============================================================================
# Submit:
#   cd /cluster/project/reddy/katja/NGS_pipeline
#   sbatch run/04_run_variant_labeling_KH157_tcr_population.sbatch
# =============================================================================
#SBATCH --job-name=KH157_variant_labeling
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/04_variant_labeling_KH157/slurm_KH157_combined_tcr_population_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/04_variant_labeling_KH157/slurm_KH157_combined_tcr_population_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=1:00:00

set -euo pipefail

PROJECT_ROOT="/cluster/project/reddy/katja/NGS_pipeline"
YAML="${PROJECT_ROOT}/config/TCRA3_KH157/script4_variant_labeling_22_05_2026_combined_with_tcr_population.yaml"
VENV="${PROJECT_ROOT}/venv/bin/activate"

# ---------------------------------------------------------------------------
mkdir -p "${PROJECT_ROOT}/logs/04_variant_labeling_KH157"
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
