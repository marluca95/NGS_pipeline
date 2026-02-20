#!/bin/bash
#SBATCH --job-name=variant_label
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/variant_label_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/variant_label_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G

set -euo pipefail

cd /cluster/project/reddy/katja/NGS_pipeline
mkdir -p logs

source /cluster/project/reddy/marluca/NGS_pipeline/venv/bin/activate

python scripts/script4_variant_labeling.py \
  --yaml_config /cluster/project/reddy/katja/NGS_pipeline/config/script4_variant_labelingDMF5.yaml
