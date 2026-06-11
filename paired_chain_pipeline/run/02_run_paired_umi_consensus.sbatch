#!/bin/bash
# Paired-Chain Pipeline — Step 02: Paired UMI consensus (ARRAY JOB)
# One task per unique sample_id.
#
# Submit from the paired_chain_pipeline root:
#   cd /cluster/project/reddy/marluca/NGS_pipeline/paired_chain_pipeline
#   YAML="config/P3620_abDMS/script2_paired_umi_config.yaml"
#   TSV=$(awk -F': ' '/^sample_sheet:/ {gsub(/"/, "", $2); print $2; exit}' "${YAML}")
#   N=$(awk -F'\t' 'NR>1 {print $1}' "${TSV}" | sort -u | wc -l)
#   sbatch --array=0-$((N-1)) run/02_run_paired_umi_consensus.sbatch

#SBATCH --job-name=02_paired_umi
#SBATCH --output=/cluster/project/reddy/marluca/NGS_pipeline/logs/P3620_LUCA-TCRDMF5-abDMS/02_umi_consensus/slurm_%A_%a.out
#SBATCH --error=/cluster/project/reddy/marluca/NGS_pipeline/logs/P3620_LUCA-TCRDMF5-abDMS/02_umi_consensus/slurm_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16G
#SBATCH --time=48:00:00

PROJECT_ROOT="/cluster/project/reddy/marluca/NGS_pipeline"
PIPELINE_ROOT="${PROJECT_ROOT}/paired_chain_pipeline"
YAML="${YAML:-${PIPELINE_ROOT}/config/P3620_abDMS/script2_paired_umi_config.yaml}"
VENV="${PROJECT_ROOT}/venv/bin/activate"

# Auto-submit as array job if invoked without --array
if [[ -z "${SLURM_ARRAY_TASK_ID+x}" ]]; then
    TSV=$(awk -F': ' '/^sample_sheet:/ {gsub(/"/, "", $2); print $2; exit}' "${YAML}")
    N=$(awk -F'\t' 'NR>1 {print $1}' "${TSV}" | sort -u | wc -l)
    mkdir -p "${PROJECT_ROOT}/logs/P3620_LUCA-TCRDMF5-abDMS/02_umi_consensus"
    echo "Auto-submitting as array job: 0-$((N-1)) tasks"
    exec sbatch --array=0-$((N-1)) "${PIPELINE_ROOT}/run/02_run_paired_umi_consensus.sbatch"
fi

set -euo pipefail

mkdir -p "${PROJECT_ROOT}/logs/P3620_LUCA-TCRDMF5-abDMS/02_umi_consensus"

if [[ ! -f "${YAML}" ]]; then
    echo "ERROR: YAML config not found: ${YAML}" >&2
    exit 1
fi

TSV=$(awk -F': ' '/^sample_sheet:/ {gsub(/"/, "", $2); print $2; exit}' "${YAML}")
if [[ -z "${TSV:-}" || ! -f "${TSV}" ]]; then
    echo "ERROR: Could not resolve sample_sheet from YAML: ${YAML}" >&2
    exit 1
fi

SAMPLE_ID=$(awk -F'\t' 'NR>1 {print $1}' "${TSV}" | sort -u | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
if [[ -z "${SAMPLE_ID:-}" ]]; then
    echo "ERROR: No sample_id for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi


cd "${PIPELINE_ROOT}/scripts"
source "${VENV}"

echo "============================================="
echo "Job:    ${SLURM_JOB_ID}, Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo "Step:   02_paired_umi_consensus"
echo "Sample: ${SAMPLE_ID}"
echo "Config: ${YAML}"
echo "============================================="

python -u script2_paired_umi_consensus.py \
    --yaml_config "${YAML}" \
    --sample_id   "${SAMPLE_ID}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
