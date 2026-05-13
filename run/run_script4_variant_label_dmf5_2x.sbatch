#!/bin/bash
#SBATCH --job-name=variant_label_dmf5_2x
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/variant_label_dmf5_2x_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/variant_label_dmf5_2x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G

set -euo pipefail

cd /cluster/project/reddy/katja/NGS_pipeline
mkdir -p logs

source /cluster/project/reddy/marluca/NGS_pipeline/venv/bin/activate

python scripts/script4_variant_labeling_dmf5_2x.py \
  --yaml_config /cluster/project/reddy/katja/NGS_pipeline/config/DMF5_3xNNK/script4_variant_labelingDMF5_2x.yaml
