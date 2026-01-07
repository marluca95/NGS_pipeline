import os
import re
import pandas as pd
from typing import List, Dict
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


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
                        "fastq_full_name": f
                    })
                else:
                    unmatched_files.append(os.path.join(dirpath, f))

    except OSError as e:
        logging.error(f"Error accessing files: {e}")
        raise

    if unmatched_files:
        logging.warning(
            f"{len(unmatched_files)} files did not match the expected pattern. Examples:\n"
            + "\n".join(unmatched_files[:5])
        )

    if not rows:
        logging.error("No matching FASTQ files found. Please check your directory or filename pattern.")
        raise ValueError("No matching FASTQ files found.")

    df = pd.DataFrame(rows)

    # sort nicely for downstream pipelines
    df = df.sort_values(["sample_id", "sample_name", "lane", "read"])

    return df


if __name__ == "__main__":
    root_dir: str = "/cluster/project/reddy/marluca/NGS_pipeline/data/raw/P3408_LUCA-TCRA3"
    df_samples: pd.DataFrame = collect_fastq_metadata(root_dir)
    df_samples.to_csv("samples.tsv", sep="\t", index=False)

    logging.info("Metadata extraction complete. First few rows:")
    print(df_samples.head())
