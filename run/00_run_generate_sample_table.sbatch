#!/bin/bash
# =============================================================================
# Step 00 — Sample Sheet Generation (script0_generate_sample_table.py)
# =============================================================================
# Submit:
#   cd /cluster/project/reddy/katja/NGS_pipeline
#   sbatch run/00_run_generate_sample_table.sbatch
# =============================================================================
#SBATCH --job-name=00_sample_table
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/00_sample_table/slurm_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/00_sample_table/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00

set -euo pipefail

PROJECT_ROOT="/cluster/project/reddy/katja/NGS_pipeline"
ROOT_DIR="/cluster/project/reddy/katja/data/raw/P3408_LUCA-TCRA3_KH157"
OUTPUT="${PROJECT_ROOT}/data/P3408_LUCA-TCRA3_KH157/P3408_LUCA-TCRA3-KH157_samples.tsv"
LOGS_DIR="${PROJECT_ROOT}/logs/P3408_LUCA-TCRA3_KH157/00_sample_table"
VENV="${PROJECT_ROOT}/venv/bin/activate"

mkdir -p "${PROJECT_ROOT}/logs/00_sample_table"
cd "${PROJECT_ROOT}"
source "${VENV}"

echo "============================================="
echo "Job:   ${SLURM_JOB_ID}"
echo "Node:  $(hostname)"
echo "Start: $(date)"
echo "Step:  00_sample_table"
echo "Root:  ${ROOT_DIR}"
echo "Out:   ${OUTPUT}"
echo "============================================="

python -u scripts/script0_generate_sample_table.py \
    --root_dir "${ROOT_DIR}" \
    --output   "${OUTPUT}" \
    --logs_dir "${LOGS_DIR}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
