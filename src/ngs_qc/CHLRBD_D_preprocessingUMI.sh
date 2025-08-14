#!/bin/bash

#SBATCH --account=es_reddy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=16384
#SBATCH --job-name="CHL_RBD_preprocessingUMI"
#SBATCH --mem-per-cpu=16384
#SBATCH --output="logs/CHL_RBD_preprocessingUMI_out"
#SBATCH --error="logs/CHL_RBD_preprocessingUMI_error"
#SBATCH --open-mode=truncate

#TODO: test mismatch workflow 

# Loop over all matching files
BASE_DIR="${1:-/cluster/project/reddy/rakuhn/MammalianExpression/output}"
CONFIG_FILE="${2:-CHLRBD_A_config.yaml}"
VENV_PATH="${VENV_PATH:-/cluster/project/reddy/rakuhn/plmfit_venv}"
SCRIPT_DIR="${1:-/cluster/project/reddy/rakuhn/MammalianExpression/CHL_RBD/script_python}"
OUT_DIR="${1:-/cluster/project/reddy/rakuhn/MammalianExpression/CHL_RBD}"
max_primer_mismatch=1


find "$BASE_DIR" -type f -name "*_combined.fastq.gz" | while read -r filepath; do

    # Extract just the filename
    filename=$(basename "$filepath")

    # Remove suffix to get sample name
    sample_name="${filename%%_combined.fastq.gz}"

    echo "Submitting job for $sample_name..."

    sbatch <<EOF
#!/bin/bash
#SBATCH --account=es_reddy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=16384
#SBATCH --job-name=$sample_name
#SBATCH --output=logs/CHLRBD_preprocessingUMI_${sample_name}.out
#SBATCH --error=logs/CHLRBD_preprocessingUMI_${sample_name}.err


echo "Loading Modules"

module load stack/2024-06 gcc/12.2.0  bbmap/39.01
module load stack/2024-06 gcc/12.2.0 python/3.9.18
module load stack/2024-06 gcc/12.2.0 openjdk/21.0.3_9

echo "Loading Venv"

source "$VENV_PATH/bin/activate"

# Run your actual analysis script

#python "$SCRIPT_DIR/script_python/CHLRBD_D_preprocessing_extract_UMI.py" --yaml_config "$CONFIG_FILE" --samplename "${sample_name}" --fastq "${filepath}" --max_primer_mismatch "${max_primer_mismatch}"

#python "$SCRIPT_DIR/script_python/CHLRBD_D_preprocessing_filter_UMI.py" --yaml_config "$CONFIG_FILE" --samplename "${sample_name}" --csv "UMI_extracted_${sample_name}.csv"

mkdir -p "$OUT_DIR/CHL_RBD_preprocessingUMI"
mkdir -p "$OUT_DIR/CHL_RBD_preprocessingUMI/Translation_summary"
mv ./umi_translation_*.csv CHL_RBD_preprocessingUMI/Translation_summary
mkdir -p "$OUT_DIR/CHL_RBD_preprocessingUMI/Extracted_UMI"
mv ./UMI_extracted_*.csv CHL_RBD_preprocessingUMI/Extracted_UMI
mkdir -p "$OUT_DIR/CHL_RBD_preprocessingUMI/Consensus_summary"
mv ./umi_consensus_summary_*.csv CHL_RBD_preprocessingUMI/Consensus_summary
mkdir -p "$OUT_DIR/CHL_RBD_preprocessingUMI/plots"

python "$SCRIPT_DIR/script_python/CHLRBD_D_preprossessingUMI_plot.py" --samplename "${sample_name}" --csv "CHL_RBD_preprocessingUMI/Translation_summary/umi_translation_${sample_name}.csv" --outpath "$OUT_DIR/CHL_RBD_preprocessingUMI/plots"

echo "end of script"

EOF

done
