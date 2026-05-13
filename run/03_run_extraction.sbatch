#!/bin/bash
# =============================================================================
# Step 03 — Library Filtering and Extraction (script3_library_filtering_and_extraction.py)
# =============================================================================
# This is an array job: one task per unique sample_id.
#
# Submit:
#   cd /cluster/project/reddy/katja/NGS_pipeline
#   TSV="scripts/P3408_LUCA-TCRA3_samples.tsv"
#   N=$(awk -F'\t' 'NR>1 {print $1}' "${TSV}" | sort -u | wc -l)
#   sbatch --array=0-$((N-1)) run/03_run_extraction.sbatch
# =============================================================================
#SBATCH --job-name=03_extraction
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/03_extraction/slurm_%A_%a.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/03_extraction/slurm_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00

set -euo pipefail

PROJECT_ROOT="/cluster/project/reddy/katja/NGS_pipeline"
YAML="${PROJECT_ROOT}/config/TCRA3_KH157/script3_library_filtering_and_extraction.yaml"
TSV="${PROJECT_ROOT}/data/P3569_LUCA-TCRA3_KH157/P3569_LUCA-TCRA3_KH157_samples.tsv"
VENV="${PROJECT_ROOT}/venv/bin/activate"

# ---------------------------------------------------------------------------
# Resolve sample_id from array task index
SAMPLE_ID=$(awk -F'\t' 'NR>1 {print $1}' "${TSV}" | sort -u | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

if [[ -z "${SAMPLE_ID:-}" ]]; then
    echo "ERROR: No sample_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

mkdir -p "${PROJECT_ROOT}/run/logs/03_extraction"
cd "${PROJECT_ROOT}"
source "${VENV}"

echo "============================================="
echo "Job:    ${SLURM_JOB_ID}, Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo "Step:   03_extraction"
echo "Sample: ${SAMPLE_ID}"
echo "Config: ${YAML}"
echo "============================================="

python -u scripts/script3_library_filtering_and_extraction.py \
    --yaml_config "${YAML}" \
    --sample_id   "${SAMPLE_ID}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
