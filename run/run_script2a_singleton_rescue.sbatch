#!/bin/bash
#SBATCH --job-name=singleton_rescue
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/slurm_singleton_rescue_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/slurm_singleton_rescue_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00

# Stop on error
set -euo pipefail

echo "========================================="
echo "Starting singleton_rescue job"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Start time: $(date)"
echo "========================================="

# Go to project directory
cd /cluster/project/reddy/katja/NGS_pipeline

# Activate virtual environment
source venv/bin/activate

# Run singleton rescue script
python scripts/Script2a_singleton_rescue.py \
    --yaml_config config/script2a_singleton_rescueDMF5.yaml

echo "========================================="
echo "Finished singleton_rescue job"
echo "End time: $(date)"
echo "========================================="
