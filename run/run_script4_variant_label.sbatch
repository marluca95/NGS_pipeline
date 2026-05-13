#!/bin/bash
#SBATCH --job-name=variant_label
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/04_variant_label/variant_label_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/04_variant_label/variant_label_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G

set -euo pipefail

cd /cluster/project/reddy/katja/NGS_pipeline
mkdir -p logs

source /cluster/project/reddy/katja/NGS_pipeline/venv/bin/activate

python scripts/script4_variant_labeling.py \
  --yaml_config /cluster/project/reddy/katja/NGS_pipeline/config/TCRA3_KH157/script4_variant_labeling.yaml
