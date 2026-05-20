# NGS Pipeline

This folder contains a multi-step NGS processing pipeline for TCR library experiments.

At a high level, the pipeline does this:

1. Builds a sample sheet from raw FASTQ files.
2. Trims reads with BBDuk and combines lanes per sample.
3. Collapses reads by UMI into consensus reads and keeps singletons.
4. Optionally rescues singleton-derived variants with a minimum UMI threshold.
5. Extracts and validates the variable region against library rules.
6. Computes per-peptide variant enrichment and specificity labels.
7. Optionally exports a final TCR-peptide-label CSV for ML.

## Pipeline Flow

```
raw FASTQ tree
	-> script0_generate_sample_table.py
sample_sheet.tsv
	-> script1_preprocessing.py
01_preprocessed/*_combined_trimmed.fastq.gz
	-> script2_umi_consensus.py
02_umi_consensus/*_consensus.fastq.gz + *_singletons.fastq.gz
	-> script2a_singleton_rescue.py (optional)
02_umi_consensus/*_singletons_rescued.fastq.gz
	-> script3_library_filtering_and_extraction.py
03_extracted/*.filtered.PASS.aa.tsv
	-> script4_variant_labeling.py (or script4_variant_labeling_dmf5_2x.py)
04_variant_labeling/*.variant_labeling.csv
	-> script5_create_tcr_peptide_specificity_csv.py (optional)
final tcr_peptide_specificity CSV
```

## Folder Structure

```
NGS_pipeline/
	config/                  # YAML configs per project/run and per step
	data/                    # pipeline outputs grouped by project
	logs/                    # script and slurm logs
	notebooks/               # exploratory and QC notebooks
	run/                     # sbatch launchers for steps 00-04
	scripts/                 # main pipeline scripts (documented below)
	src/                     # helper analysis/QC code
	requirements.txt         # python dependencies
```

Notes:

- `config/pipeline_config.yaml` is a full reference/template for all step parameters.
- Step-specific YAML files (for example in `config/TCRA3_KH157/`) are used for actual runs.
- `run/*.sbatch` scripts show the intended cluster execution order.

## Typical Run Order

From `NGS_pipeline/` on the cluster:

1. `sbatch run/00_run_generate_sample_table.sbatch`
2. `sbatch run/01_run_preprocessing.sbatch`
3. `sbatch --array=... run/02_run_umi_consensus.sbatch`
4. `sbatch run/02a_run_singleton_rescue.sbatch` (optional but typically used)
5. `sbatch --array=... run/03_run_extraction.sbatch`
6. `sbatch run/04_run_variant_labeling.sbatch`

## Scripts: quick reference

All scripts live in `scripts/`. Below is a compact summary — purpose, essential inputs, and main outputs. Configuration is provided via per-step YAML files in `config/` (see `run/` for example sbatch launchers).

- `script0_generate_sample_table.py` — Build a sample sheet from raw FASTQ files.
	- Input: raw `.fastq.gz` directory tree.
	- Output: tabular sample sheet TSV (sample_id, sample_name, fastq, ...).

- `script1_preprocessing.py` — Trim reads (BBDuk) and combine lanes per sample.
	- Input: sample sheet + FASTQ files.
	- Output: per-sample trimmed FASTQs and a combined `{sample}_combined_trimmed.fastq.gz` plus per-fastq summary CSV.

- `script2_umi_consensus.py` — Call UMI consensus reads and keep singletons.
	- Input: `{sample}_combined_trimmed.fastq.gz` files.
	- Output: `{sample}_consensus.fastq.gz`, `{sample}_singletons.fastq.gz`, QC TSV and per-sample summary.

- `Script2a_singleton_rescue.py` — (Optional) filter/retain singleton-derived variants meeting a UMI-count threshold.
	- Input: singleton FASTQs from step 02.
	- Output: rescued singleton FASTQs and per-sample variant TSVs, plus a run summary.

- `script3_library_filtering_and_extraction.py` — Extract variable region, validate library rules, write PASS/FAIL and AA TSV.
	- Input: consensus and/or rescued-singleton FASTQs.
	- Output: `.filtered.PASS.fastq.gz`, optional `.filtered.FAIL.fastq.gz`, `.PASS.aa.tsv`, and a summary TSV per sample.

- `script4_variant_labeling.py` (and `script4_variant_labeling_dmf5_2x.py`) — Count variants per condition, compute enrichment, assign specificity labels (0/1/2), and write per-peptide CSVs plus a run summary.
	- Input: `*.filtered.PASS.aa.tsv` files from step 03.
	- Output: `{peptide}.variant_labeling.csv` and `variant_labeling_summary.csv`.

- `script5_create_tcr_peptide_specificity_csv.py` — Combine labeled CSVs with a peptide mapping to produce a final ML-ready CSV.
	- Input: labeled CSVs and a peptide mapping file.
	- Output: consolidated TCR–peptide–label CSV.

## Script Utility Modules

`scripts/utils/` contains shared helpers used by the pipeline (config, logging, sample-sheet helpers, and metrics writing).

## Important Notes

- Use the per-step YAML configs in `config/` (examples in `config/TCRA3_KH157/`) when running the `run/*.sbatch` launchers.
- `script5_create_tcr_peptide_specificity_csv.py` currently uses fixed paths inside the script; you may want to parameterize it for reuse.
