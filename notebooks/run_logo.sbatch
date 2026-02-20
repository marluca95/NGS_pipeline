#!/bin/bash
#SBATCH --job-name=aa_logos
#SBATCH --output=/cluster/project/reddy/katja/NGS_pipeline/logs/aa_logos_%j.out
#SBATCH --error=/cluster/project/reddy/katja/NGS_pipeline/logs/aa_logos_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=01:00:00

set -euo pipefail

PROJECT=/cluster/project/reddy/katja/NGS_pipeline
SCRIPT=/cluster/project/reddy/katja/NGS_pipeline/notebooks/logomaker_aa.py
INPUT=/cluster/project/reddy/katja/NGS_pipeline/data/P3408_LUCA-TCRA3/extracted
OUTDIR=/cluster/project/reddy/katja/NGS_pipeline/data/P3408_LUCA-TCRA3/extracted/logos_aa

cd "$PROJECT"
source venv/bin/activate

python "$SCRIPT" --input_dir "$INPUT" --out_dir "$OUTDIR"
