#!/bin/bash
#SBATCH --job-name=cdr3_eda_a100
#SBATCH --output=/cluster/project/reddy/marluca/NGS_pipeline/logs/cdr3_eda_%A_%a.out
#SBATCH --error=/cluster/project/reddy/marluca/NGS_pipeline/logs/cdr3_eda_%A_%a.err
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --partition=gpup
#SABTCH --gres=gpu:a100_80gb:1
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --time=24:00:00


# Load modules
module load stack/2024-06 gcc/12.2.0 python/3.9.18
module load stack/2024-06 cuda/12.8.0
source /cluster/project/reddy/marluca/NGS_pipeline/venv/bin/activate

OUTPUT_DIR="/cluster/project/reddy/marluca/NGS_pipeline/data/results"
mkdir -p $OUTPUT_DIR

python /cluster/project/reddy/marluca/NGS_pipeline/src/ml/cdr3b_eda.py --data_csv /cluster/project/reddy/marluca/NGS_pipeline/data/TCR-A3_CDR3b-t4x_MAGE-A3_specificity.csv --output_dir $OUTPUT_DIR