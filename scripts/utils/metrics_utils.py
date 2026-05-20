"""
Pipeline-wide metrics writing utility.

Every pipeline step calls write_sample_metrics() once per processed sample.
All steps write to the same long-format TSV schema:

    sample_id | sample_name | step | metric | value

This makes it trivial to concatenate across steps and plot the full read-loss
funnel with visualization/qc/plot_pipeline_funnel.py.

Filename convention (one file per step output directory):
    {output_dir}/per_sample_metrics.tsv
"""

import csv
import logging
from pathlib import Path
from typing import Dict, Union

logger = logging.getLogger(__name__)

# Fixed column order for the metrics TSV
_HEADER = ["sample_id", "sample_name", "step", "metric", "value"]


def write_sample_metrics(
    metrics_path: Union[str, Path],
    sample_id: str,
    sample_name: str,
    step: str,
    metrics: Dict[str, Union[int, float]],
) -> None:
    """
    Append per-sample metrics as rows to a long-format TSV.

    Creates the file with a header row if it does not exist yet.
    Appends without re-reading the file, so it is safe to call once per sample
    inside a loop (sequential). 
    # Not safe for concurrent writes (e.g. SLURM array jobs) without locking or per-sample output files!!!

    Args:
        metrics_path: Path to the TSV file (created if needed).
        sample_id:    Unique sample identifier.
        sample_name:  Human-readable sample name.
        step:         Pipeline step label, e.g. "01_preprocessing".
        metrics:      Dict mapping metric name -> numeric value.
    """
    path = Path(metrics_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    write_header = not path.exists() or path.stat().st_size == 0

    with open(path, "a", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        if write_header:
            writer.writerow(_HEADER)
        for metric, value in metrics.items():
            writer.writerow([sample_id, sample_name, step, metric, value])

    logger.debug("Wrote %d metric(s) for %s_%s to %s", len(metrics), sample_id, sample_name, path)
