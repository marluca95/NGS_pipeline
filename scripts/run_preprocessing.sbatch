#!/bin/bash
#SBATCH --job-name=bbduk_preproc
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/slurm_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/slurm_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=04:00:00


cd /cluster/project/reddy/katja/NGS_pipeline
source venv/bin/activate

python scripts/script1_preprocessing.py \
  --yaml_config config/script1_preprocessing_config.yaml \
  --sample_sheet scripts/P3408_LUCA-TCRA3_samples.tsv \
  --sample_id GFB-61811

