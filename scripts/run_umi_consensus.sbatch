#!/bin/bash
#SBATCH --job-name=umi_consensus_FULLRUN
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/umi_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/umi_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=24G
#SBATCH --time=12:00:00

cd /cluster/project/reddy/katja/NGS_pipeline
source venv/bin/activate

python -u scripts/script2_umi_consensus.py \
  --yaml_config config/script2_umi_config.yaml \
  --sample_id GFB-61819
