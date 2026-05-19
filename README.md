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
	-> Script2a_singleton_rescue.py (optional)
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

## Scripts: Inputs and Outputs

All scripts below are in `scripts/`.

### `script0_generate_sample_table.py`

Purpose:
- Traverses a raw FASTQ directory tree and builds a sample sheet TSV.

CLI:
- `--root_dir` (required): root containing `.fastq.gz` files.
- `--output` (required): output sample table TSV.
- `--logs_dir` (optional): log directory.

Input:
- FASTQ files matching naming pattern like
	`{sample_id}_{...}_{sample_name}_S{...}_L{lane}_R{read}_{...}.fastq.gz`

Output:
- Sample sheet TSV with columns:
	`sample_id, sample_name, run_depth, lane, read, fastq`
- Run log file in the configured logs directory.

### `script0_generate_sample_table_P3569.py`

Purpose:
- Variant of step 00 for `P3569_LUCA-TCRA3_KH157`.

CLI:
- `--root_dir` (optional, has default path)
- `--output` (optional, has default path)

Input:
- FASTQ files under the provided root directory.

Output:
- Sample sheet TSV with the same columns as step 00.

### `script1_preprocessing.py`

Purpose:
- Quality-trims reads with BBDuk and combines trimmed FASTQs per sample.

CLI:
- `--yaml_config` (required): preprocessing YAML.
- `--sample_sheet` (optional): overrides YAML `sample_sheet`.
- `--sample_id` (optional): run one sample only.

Key YAML inputs:
- `sample_sheet`: TSV from step 00.
- `output_dir`: step 01 output directory.
- `bbmap_dir`: BBMap installation containing `bbduk.sh`.
- optional trimming params (`qtrim`, `quality_threshold`, `min_length`, `bbduk_xmx_gb`).

Input:
- Per-row FASTQ paths from sample sheet (`fastq` column).

Output per sample (`{output_dir}/{sample_id}_{sample_name}/`):
- `*.trimmed.fastq.gz` (trimmed lane files)
- `{sample_label}_combined_trimmed.fastq.gz`
- `{sample_label}.bbduk_summary.csv` (per-input FASTQ metrics)

Run-level output:
- `{output_dir}/bbduk_run_summary_{run_label_or_timestamp}.csv`
- `{output_dir}/per_sample_metrics.tsv`
- Logs in configured `logs_dir`.

### `script2_umi_consensus.py`

Purpose:
- Splits reads by UMI around anchor sequence, builds posterior/quality-aware consensus for multi-read UMIs, keeps singleton UMIs.

CLI:
- `--yaml_config` (required)
- `--sample_id` (optional): process one sample.

Key YAML inputs:
- `sample_sheet`
- `input_dir` (expects step 01 sample folders with `*_combined_trimmed.fastq.gz`)
- `output_dir`
- `umi_length`, `anchor_sequence`, `max_mismatches`, `spacer_min`, `spacer_max`
- `min_reads_per_umi`, `posterior_min_log_delta`
- `combined_suffix` (default `_combined_trimmed.fastq.gz`)

Output per sample (`output_dir`):
- `{sample_label}_consensus.fastq.gz`
- `{sample_label}_singletons.fastq.gz`
- `{sample_label}_consensus_qc.tsv`
- `{sample_label}_umi_summary.txt`

Run-level output:
- `{output_dir}/per_sample_metrics.tsv`
- Logs in configured `logs_dir`.

### `Script2a_singleton_rescue.py`

Purpose:
- From singleton inserts, extracts a variable region after a WT anchor and keeps variants observed in at least `min_umis` distinct UMIs.

CLI:
- `--yaml_config` (required)
- `--sample_id` (optional)

Key YAML inputs:
- `samples_tsv`
- `singleton_dir` (usually step 02 output dir)
- `output_dir`
- `wt_anchor`, `var_len`, `min_umis`
- `singleton_suffix` (default `_singletons.fastq.gz`)

Input:
- `{singleton_dir}/{sample_label}{singleton_suffix}`

Output per sample:
- `{output_dir}/{sample_label}_singletons_rescued.fastq.gz`
- `{output_dir}/{sample_label}_singleton_variants_{min_umis}umis.tsv`

Run-level output:
- `{output_dir}/singleton_rescue_summary.tsv`
- `{output_dir}/per_sample_metrics.tsv`
- Logs in configured `logs_dir`.

### `script3_library_filtering_and_extraction.py`

Purpose:
- Extracts the variable region after anchor, reverse-complements into design orientation, validates against library rules, writes PASS/FAIL FASTQ, and AA translation table.

CLI:
- `--yaml_config` (required)
- `--sample_id` (required): one sample per run/task.

Key YAML inputs:
- `samples_tsv`
- `input_dir`
- `input_suffixes_consensus`
- `input_suffixes_singletons`
- `output_dir`
- `anchor_sequence`, `anchor_max_mismatches`, `region_length_nt`
- `drop_if_contains_N`, `write_fail_fastq`, `write_aa_tsv`
- `library_mode` = `combinatorial` or `3xNNK`
- mode-specific settings:
	- combinatorial: `degenerate_codons` (11 codons)
	- 3xNNK: `WT_nt_33`, `max_mutated_codons`

Input:
- For each sample label, one or both:
	- `{input_dir}/{sample_label}{input_suffixes_consensus}`
	- `{input_dir}/{sample_label}{input_suffixes_singletons}`

Output per sample:
- `{output_dir}/{sample_label}.filtered.PASS.fastq.gz`
- `{output_dir}/{sample_label}.filtered.FAIL.fastq.gz` (if enabled)
- `{output_dir}/{sample_label}.filtered.PASS.aa.tsv`
- `{output_dir}/{sample_label}.filtered.summary.tsv`

Run-level output:
- `{output_dir}/per_sample_metrics.tsv`
- Logs in configured `logs_dir`.

### `script4_variant_labeling.py`

Purpose:
- Builds per-peptide enrichment tables and assigns specificity labels (`0/1/2`) using library, negative, and positive selections (`pos3x`, optional `pos1x`).

CLI:
- `--yaml_config` (required)

Key YAML inputs:
- `input_dir`, `output_dir`
- `input_glob` (default `*.filtered.PASS.aa.tsv`)
- `aa_column` (default `aa_seq`)
- file parsing suffixes: `lib_suffix`, `neg_suffix`, `pos3x_suffix`, `pos1x_suffix`, optional `file_token_prefix`
- thresholds: `pseudocount`, `status_up`, `status_down`
- `required_conditions`

Input:
- Step 03 AA TSV files grouped by peptide and condition.

Output:
- Per peptide: `{peptide_key}.variant_labeling.csv`
- Run summary: `variant_labeling_summary.csv`
- Logs in configured `logs_dir`.

### `script4_variant_labeling_dmf5_2x.py`

Purpose:
- Alternative variant-labeling workflow for datasets with one positive condition (`pos`) and one negative condition (`neg`) plus library.

CLI:
- `--yaml_config` (or legacy `--yaml_file`)

Input:
- Same type as step 04, but suffix schema is `lib/neg/pos` (2x design).

Output:
- Per peptide: `{peptide_key}.variant_labeling.csv`
- Run summary: `variant_labeling_summary.csv`
- Logs in configured `logs_dir`.

### `script5_create_tcr_peptide_specificity_csv.py`

Purpose:
- Combines variant-labeling CSV files and peptide-name mapping into a final ML-ready CSV.

CLI:
- No CLI arguments (uses hardcoded paths in the script).

Input:
- All `*.csv` from `DATA_DIR` (variant-labeling outputs).
- `MAPPING_PATH` containing peptide name to peptide sequence mapping.

Output:
- `OUTPUT_PATH` CSV with columns:
	- `tcr` (from `aa_seq`)
	- `peptide` (mapped peptide sequence)
	- `label` (from specificity; only 0/1 kept)

## Script Utility Modules

`scripts/utils/` contains shared helpers used by pipeline scripts:

- `config_utils.py`: YAML loading and validation.
- `logging_utils.py`: consistent run/sample log setup.
- `metrics_utils.py`: appending standardized per-sample metrics TSVs.
- `sample_utils.py`: sample-sheet loading and column checks.

## Important Notes

- `Script2a_singleton_rescue.py` has an uppercase `S` in filename; run scripts reference this exact name.
- `script5_create_tcr_peptide_specificity_csv.py` currently uses hardcoded absolute paths. If you want it reusable across projects, parameterize these via CLI arguments or YAML.
- Several datasets/config directories exist (`TCRA3_KH157`, `TCRA3_3xNNK`, `DMF5_3xNNK`, etc.). Use matching YAML files for your project when launching run steps.
