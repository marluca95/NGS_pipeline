#!/bin/bash
#SBATCH --job-name=umi_consensus
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/umi/slurm_%A_%a.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/umi/slurm_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=24G
#SBATCH --time=4:00:00

set -euo pipefail

cd /cluster/project/reddy/katja/NGS_pipeline
source venv/bin/activate

YAML="config/script2_umi_config_DMF5.yaml"
TSV="/cluster/project/reddy/katja/NGS_pipeline/scripts/P3481_LUCA-TCRDMF5_samples.tsv"
LOGDIR="/cluster/project/reddy/katja/NGS_pipeline/logs/umi"

mkdir -p "$LOGDIR"

# Get sample ID for this array index
SAMPLE_ID=$(
  awk -F'\t' 'NR>1 {print $1}' "$TSV" \
  | sort -u \
  | sed -n "$((SLURM_ARRAY_TASK_ID+1))p"
)

if [[ -z "${SAMPLE_ID:-}" ]]; then
  echo "No SAMPLE_ID found for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID" >&2
  exit 1
fi

# Redirect all output AFTER this point to per-sample logs
exec > "${LOGDIR}/umi_${SAMPLE_ID}.out" \
     2> "${LOGDIR}/umi_${SAMPLE_ID}.err"

echo "========================================"
echo "Starting sample: $SAMPLE_ID"
echo "Date: $(date)"
echo "JobID: ${SLURM_JOB_ID}"
echo "ArrayTask: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "========================================"

python -u scripts/script2_umi_consensus.py \
  --yaml_config "$YAML" \
  --sample_id "$SAMPLE_ID"

echo "Finished sample: $SAMPLE_ID at $(date)"
