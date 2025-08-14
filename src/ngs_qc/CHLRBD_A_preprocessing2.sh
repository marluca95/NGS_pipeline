#!/bin/bash

#SBATCH --account=es_reddy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=16384
#SBATCH --job-name="CHLRBD_preprocessing"
#SBATCH --mem-per-cpu=16384
#SBATCH --output="CHLRBD_preprocessing_out"
#SBATCH --error="CHLRBD_preprocessing_error"
#SBATCH --open-mode=truncate

echo "Loading Modules"

module load stack/2024-06 gcc/12.2.0  bbmap/39.01
module load stack/2024-06 gcc/12.2.0 python/3.9.18
module load stack/2024-06 gcc/12.2.0 openjdk/21.0.3_9

echo "Loading Venv"

VENV_PATH="${VENV_PATH:-/cluster/project/reddy/rakuhn/plmfit_venv}"
source "$VENV_PATH/bin/activate"

SEQPIPE_DIR="${SEQPIPE_DIR:-/cluster/project/reddy/rakuhn/seq-pipeline}"
CONFIG_FILE="${1:-CHLRBD_A_config.yaml}"
SAMPLE_PREFIX="${2:-RBM}"

echo "Copying FASTQ files"

#cp /cluster/project/reddy/rakuhn/MammalianExpression/data/data/openbis-gfb.ethz.ch/*/original/*/raw_data/*.fastq.gz ./data/fastq/

echo "Submitting Job"

python "$SEQPIPE_DIR/pre_process2.py" --yaml_config "$CONFIG_FILE" --sample "$SAMPLE_PREFIX"
