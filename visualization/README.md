# visualization/

This folder contains standalone plotting scripts for the NGS pipeline.
Each script reads pre-computed output files from the data/ directory and
saves figures here (or to a sub-directory you specify).

## Intended structure

```
visualization/
  qc/
    plot_preprocessing_qc.py   — reads 01_preprocessed/ BBDuk summary CSVs
    plot_umi_qc.py             — reads 02_umi_consensus/ UMI stats / QC TSVs
    plot_extraction_qc.py      — reads 03_extracted/ filtered.summary.tsv files
  analysis/
    plot_enrichment.py         — reads 04_variant_labeling/ variant_labeling.csv files
    plot_specificity.py        — specificity distributions across peptides
```

## Usage pattern

Import plotting functions from here in Jupyter notebooks to keep notebook cells
clean and avoid logic duplication between notebooks and scripts:

```python
# In a notebook:
import sys
sys.path.insert(0, "../visualization")
from qc.plot_preprocessing_qc import plot_bbduk_summary
plot_bbduk_summary(summary_csv="...")
```

Scripts in this folder should be self-contained — they accept paths as CLI
arguments (or via a config file) and do not depend on script-level state.
