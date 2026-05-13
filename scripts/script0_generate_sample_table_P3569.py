import argparse
import re
from pathlib import Path
from typing import Dict, List

import pandas as pd


DEFAULT_ROOT = "/cluster/project/reddy/katja/data/raw/P3569_LUCA-TCRA3_KH157"
DEFAULT_OUTPUT = (
    "/cluster/project/reddy/katja/NGS_pipeline/data/"
    "P3569_LUCA-TCRA3_KH157/P3569_LUCA-TCRA3_KH157_samples.tsv"
)


def collect_sample_fastqs(root_dir: Path) -> pd.DataFrame:
    """Collect FASTQ metadata for the P3569 sample tree."""
    pattern = re.compile(
        r"([A-Z0-9-]+)_(\d+)_\d+_(.+?)_S\d+_L(\d+)_R(\d+)_\d+\.fastq\.gz$"
    )

    rows: List[Dict[str, str]] = []

    for fastq in sorted(root_dir.rglob("*.fastq.gz")):
        match = pattern.match(fastq.name)
        if not match:
            continue

        sample_id, _, sample_name, lane, read = match.groups()
        rows.append(
            {
                "sample_id": sample_id,
                "sample_name": sample_name,
                "run_depth": "",
                "lane": lane,
                "read": read,
                "fastq": str(fastq),
            }
        )

    if not rows:
        raise ValueError(
            f"No matching FASTQ files found under: {root_dir}"
        )

    return pd.DataFrame(rows).sort_values(
        ["sample_id", "sample_name", "lane", "read"]
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate sample table for P3569_LUCA-TCRA3_KH157 only."
    )
    parser.add_argument("--root_dir", default=DEFAULT_ROOT)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    root_dir = Path(args.root_dir)
    output_path = Path(args.output)

    df = collect_sample_fastqs(root_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)

    print(f"Wrote {len(df)} rows to {output_path}")


if __name__ == "__main__":
    main()
