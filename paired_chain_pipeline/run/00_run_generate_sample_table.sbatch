#!/bin/bash
# Paired-Chain Pipeline — Step 00: Generate sample table
# Reuses the original script0_generate_sample_table.py (no changes needed).
#
# Submit from the paired_chain_pipeline root:
#   cd /cluster/project/reddy/marluca/NGS_pipeline/paired_chain_pipeline
#   sbatch run/00_run_generate_sample_table.sbatch

#SBATCH --job-name=00_paired_sample_table
#SBATCH --output=/cluster/project/reddy/marluca/NGS_pipeline/logs/P3620_LUCA-TCRDMF5-abDMS/00_sample_table/slurm_%j.out
#SBATCH --error=/cluster/project/reddy/marluca/NGS_pipeline/logs/P3620_LUCA-TCRDMF5-abDMS/00_sample_table/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00

set -euo pipefail

PROJECT_ROOT="/cluster/project/reddy/marluca/NGS_pipeline"
PIPELINE_ROOT="${PROJECT_ROOT}/paired_chain_pipeline"
ROOT_DIR="${PROJECT_ROOT}/data/raw/P3620_LUCA-TCRDMF5-abDMS"
OUTPUT="${PROJECT_ROOT}/data/P3620_LUCA-TCRDMF5-abDMS/P3620_LUCA-TCRDMF5-abDMS_samples.tsv"
LOGS_DIR="${PROJECT_ROOT}/logs/P3620_LUCA-TCRDMF5-abDMS/00_sample_table"
VENV="${PROJECT_ROOT}/venv/bin/activate"

mkdir -p "${LOGS_DIR}"
mkdir -p "$(dirname "${OUTPUT}")"

cd "${PIPELINE_ROOT}/scripts"
source "${VENV}"

echo "============================================="
echo "Job:   ${SLURM_JOB_ID}"
echo "Node:  $(hostname)"
echo "Start: $(date)"
echo "Step:  00_sample_table (Paired-Chain Pipeline)"
echo "Root:  ${ROOT_DIR}"
echo "Out:   ${OUTPUT}"
echo "============================================="

python -u script0_generate_sample_table.py \
    --root_dir "${ROOT_DIR}" \
    --output   "${OUTPUT}" \
    --logs_dir "${LOGS_DIR}"

echo "============================================="
echo "Done: $(date)"
echo "============================================="
