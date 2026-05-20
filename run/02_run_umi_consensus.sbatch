#!/bin/bash
# =============================================================================
# Step 02 — UMI Consensus Calling (script2_umi_consensus.py)
# =============================================================================
# This is an array job: one task per unique sample_id.
#
# Submit:
#   cd /cluster/project/reddy/katja/NGS_pipeline
#   YAML="config/TCRA3_KH157/script2_umi_config_KH157.yaml"
#   TSV=$(awk -F': ' '/^sample_sheet:/ {gsub(/"/, "", $2); print $2; exit}' "${YAML}")
#   N=$(awk -F'\t' 'NR>1 {print $1}' "${TSV}" | sort -u | wc -l)
#   sbatch --array=0-$((N-1)) run/02_run_umi_consensus.sbatch
# =============================================================================
#SBATCH --job-name=02_umi_consensus
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/02_umi_consensus/slurm_%A_%a.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/02_umi_consensus/slurm_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=48:00:00

set -euo pipefail

PROJECT_ROOT="/cluster/project/reddy/katja/NGS_pipeline"
YAML="${YAML:-${PROJECT_ROOT}/config/TCRA3_KH157/script2_umi_config_KH157.yaml}"
VENV="${PROJECT_ROOT}/venv/bin/activate"

if [[ ! -f "${YAML}" ]]; then
    echo "ERROR: YAML config not found: ${YAML}" >&2
    exit 1
fi

# Keep array indexing in sync with the script config by reading sample_sheet from YAML.
TSV=$(awk -F': ' '/^sample_sheet:/ {gsub(/"/, "", $2); print $2; exit}' "${YAML}")
if [[ -z "${TSV:-}" ]]; then
    echo "ERROR: Could not resolve sample_sheet from YAML: ${YAML}" >&2
    exit 1
fi
if [[ ! -f "${TSV}" ]]; then
    echo "ERROR: sample_sheet file not found: ${TSV}" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Resolve sample_id from array task index
SAMPLE_ID=$(awk -F'\t' 'NR>1 {print $1}' "${TSV}" | sort -u | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

if [[ -z "${SAMPLE_ID:-}" ]]; then
    echo "ERROR: No sample_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

mkdir -p "${PROJECT_ROOT}/logs/02_umi_consensus"
cd "${PROJECT_ROOT}"
source "${VENV}"

echo "============================================="
echo "Job:    ${SLURM_JOB_ID}, Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo "Step:   02_umi_consensus"
echo "Sample: ${SAMPLE_ID}"
echo "Config: ${YAML}"
echo "Sample sheet: ${TSV}"
echo "============================================="

python -u scripts/script2_umi_consensus.py \
    --yaml_config "${YAML}" \
    --sample_id   "${SAMPLE_ID}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
