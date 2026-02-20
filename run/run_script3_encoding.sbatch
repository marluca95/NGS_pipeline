#!/bin/bash
#SBATCH --job-name=extract_codons
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/slurm_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00

# Stop on error
set -euo pipefail

echo "Starting extract_codons job"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "Start time: $(date)"

# Go to project directory
cd /cluster/project/reddy/katja/NGS_pipeline

# Activate virtual environment
source venv/bin/activate

# Run script
python scripts/script3_extract_codons.py \
    --config config/script3_codon_extraction.yaml

echo "Finished extract_codons job"
echo "End time: $(date)"
