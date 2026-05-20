import argparse
import logging
import os
import re
from pathlib import Path
from typing import List, Dict

import pandas as pd
from utils.logging_utils import setup_pipeline_logging

SCRIPT_NAME = "script0_generate_sample_table"
logger = logging.getLogger(__name__)


def collect_fastq_metadata(root: str) -> pd.DataFrame:
    """
    Traverse a directory tree and collect metadata from FASTQ filenames.
    """
    rows: List[Dict[str, str]] = []
    unmatched_files: List[str] = []

    try:
        for dirpath, dirnames, filenames in os.walk(root):
            for f in filenames:
                if not f.endswith(".fastq.gz"):
                    continue

                match = re.match(
                    r"([A-Z0-9-]+)_(\d+)_\d+_(.+?)_S\d+_L(\d+)_R(\d+)_\d+\.fastq\.gz$",
                    f
                )

                full_path = os.path.join(dirpath, f)

                if match:
                    sample_id, _, sample_name, lane, read = match.groups()

                    # extraction of sequencing depth folder (e.g. run1-1000M)
                    run_depth = None
                    for part in dirpath.split(os.sep):
                        if re.match(r"run\d+-\d+M", part):
                            run_depth = part
                            break

                    rows.append({
                        "sample_id": sample_id,
                        "sample_name": sample_name,
                        "run_depth": run_depth,
                        "lane": lane,
                        "read": read,
                        "fastq": full_path,
                    })
                else:
                    unmatched_files.append(os.path.join(dirpath, f))

    except OSError as e:
        logger.error(f"Error accessing files: {e}")
        raise

    if unmatched_files:
        logger.warning(
            f"{len(unmatched_files)} files did not match the expected pattern. Examples:\n"
            + "\n".join(unmatched_files[:5])
        )

    if not rows:
        logger.error("No matching FASTQ files found. Please check your directory or filename pattern.")
        raise ValueError("No matching FASTQ files found.")

    df = pd.DataFrame(rows)

    # sort for downstream pipelines
    df = df.sort_values(["sample_id", "sample_name", "lane", "read"])

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Traverse a raw FASTQ directory and generate a sample sheet TSV."
    )
    parser.add_argument(
        "--root_dir",
        required=True,
        help="Root directory containing raw FASTQ files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path for the output sample sheet TSV (e.g. scripts/P3408_samples.tsv).",
    )
    parser.add_argument(
        "--logs_dir",
        default=None,
        help="Directory for log files. Defaults to NGS_pipeline/logs/script0/.",
    )
    args = parser.parse_args()

    output_path = Path(args.output)
    logs_dir = (
        Path(args.logs_dir)
        if args.logs_dir
        else Path(__file__).parent.parent / "logs" / SCRIPT_NAME
    )

    setup_pipeline_logging(logs_dir=logs_dir, script_name=SCRIPT_NAME, scope="run")

    logger.info("Root directory : %s", args.root_dir)
    logger.info("Output file    : %s", output_path)
    logger.info("Logs directory : %s", logs_dir)

    df_samples: pd.DataFrame = collect_fastq_metadata(args.root_dir)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_samples.to_csv(output_path, sep="\t", index=False)

    logger.info("Wrote %d rows to %s", len(df_samples), output_path)
    logger.info("First few rows:\n%s", df_samples.head().to_string())
