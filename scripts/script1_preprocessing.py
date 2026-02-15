import argparse
import logging
import os
import re
import subprocess
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd
from utils.config_utils import load_yaml_config
from utils.sample_utils import load_sample_sheet

# -----------------------------------------------------------------------------
# Logging setup
# -----------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler()],
)
logger = logging.getLogger(__name__)


# -----------------------------------------------------------------------------
# Utilities: per-sample logfile
# -----------------------------------------------------------------------------
def add_sample_logfile(sample_label: str, log_dir: Path) -> logging.FileHandler:
    """
    Add a per-sample logfile handler with a timestamped filename.

    Why this exists:
    - We still log to console for the full run.
    - Additionally, each sample gets its own logfile for easier debugging.
    """
    log_dir.mkdir(parents=True, exist_ok=True)

    ts = datetime.now().strftime("%Y%m%d_%H-%M-%S")
    log_path = log_dir / f"{sample_label}_{ts}.log"

    handler = logging.FileHandler(log_path)
    handler.setFormatter(
        logging.Formatter(
            "%(asctime)s %(levelname)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )

    # Important: if this script processes many samples, and you re-run the script,
    # you want to ensure the handler doesn't get duplicated accidentally.
    logger.addHandler(handler)
    logger.info("Per-sample logfile: %s", log_path)

    return handler


def remove_sample_logfile(handler: logging.FileHandler) -> None:
    """Detach and close a per-sample logfile handler."""
    logger.removeHandler(handler)
    handler.close()


# -----------------------------------------------------------------------------
# Config + I/O helpers
# -----------------------------------------------------------------------------
def parse_arguments_from_yaml(yaml_file: str) -> Dict[str, Any]:
    """Load YAML config into a dict."""
    return load_yaml_config(yaml_file)


def safe_name(s: str) -> str:
    """Make a filesystem-safe-ish name (minimal sanitization)."""
    return str(s).strip().replace(" ", "_").replace("/", "_")


def wait_for_file(filepath: Path, timeout: int = 600, interval: int = 5) -> bool:
    """
    Poll until `filepath` exists (or timeout).

    This is useful when an external process returns, but the filesystem might still
    be syncing/writing (rare but can happen on network filesystems).
    """
    start_time = time.time()
    while not filepath.exists():
        if time.time() - start_time > timeout:
            logger.error("Timeout: %s not found within %d seconds.", filepath, timeout)
            return False
        time.sleep(interval)
    return True


# -----------------------------------------------------------------------------
# BBDuk parsing + execution
# -----------------------------------------------------------------------------
def _parse_bbduk_summary(bbduk_lines: List[str]) -> Dict[str, Any]:
    """
    Parse BBDuk stdout lines to extract summary metrics.

    Notes:
    - BBDuk prints these lines to stdout/stderr depending on version.
    - We capture stdout+stderr merged, then parse the merged stream.
    """
    stats: Dict[str, Any] = {}

    # Extract version if present
    for line in bbduk_lines:
        m = re.search(r"\bVersion\s+(\S+)", line)
        if m:
            stats["bbduk_version"] = m.group(1)
            break

    def grab_counts(prefix: str) -> tuple[Optional[int], Optional[int]]:
        """
        For lines like:
          "Input: 123 reads 456 bases ..."
        return (reads, bases) as ints.
        """
        for line in bbduk_lines:
            if line.strip().startswith(prefix):
                m = re.search(r"(\d+)\s+reads.*?(\d+)\s+bases", line)
                if m:
                    return int(m.group(1)), int(m.group(2))
        return None, None

    reads_in, bases_in = grab_counts("Input:")
    reads_out, bases_out = grab_counts("Result:")
    reads_removed, bases_removed = grab_counts("Total Removed:")
    reads_qtrimmed, bases_qtrimmed = grab_counts("QTrimmed:")

    stats.update(
        {
            "reads_in": reads_in,
            "bases_in": bases_in,
            "reads_out": reads_out,
            "bases_out": bases_out,
            "reads_removed": reads_removed,
            "bases_removed": bases_removed,
            "reads_qtrimmed": reads_qtrimmed,
            "bases_qtrimmed": bases_qtrimmed,
        }
    )
    return stats


def _build_bbduk_cmd(
    bbmap_dir: str,
    input_fastq_gz: str,
    output_fastq_gz: str,
    qtrim: str,
    quality_threshold: int,
    min_length: int,
    xmx_gb: int,
) -> List[str]:
    """
    Build the BBDuk command list.

    Keeping this separate makes it easier to inspect/modify parameters later.
    """
    return [
        f"{bbmap_dir}/bbduk.sh",
        f"-Xmx{xmx_gb}g",
        f"in={input_fastq_gz}",
        f"out={output_fastq_gz}",
        f"qtrim={qtrim}",
        f"trimq={quality_threshold}",
        f"minlength={min_length}",
    ]


def trimming_single_end(
    input_fastq_gz: str,
    output_fastq_gz: str,
    bbmap_dir: str,
    qtrim: str,
    quality_threshold: int,
    min_length: int,
    xmx_gb: int,
) -> Dict[str, Any]:
    """
    Run BBDuk for a single-end FASTQ(.gz).

    Returns a dict like:
      {"status": "success", "bbduk_stats": {...}}
    or
      {"status": "error", "message": "..."}
    """
    cmd = _build_bbduk_cmd(
        bbmap_dir=bbmap_dir,
        input_fastq_gz=input_fastq_gz,
        output_fastq_gz=output_fastq_gz,
        qtrim=qtrim,
        quality_threshold=quality_threshold,
        min_length=min_length,
        xmx_gb=xmx_gb,
    )

    logger.info("Running: %s", " ".join(cmd))

    t0 = time.time()
    bbduk_lines: List[str] = []

    try:
        # Merge stderr into stdout so we can log and parse everything in one stream
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip()
            if line:
                bbduk_lines.append(line)
                logger.info("[bbduk] %s", line)

        ret = proc.wait()
        runtime_s = time.time() - t0

        if ret != 0:
            return {
                "status": "error",
                "message": f"bbduk exited with code {ret}",
                "runtime_s": runtime_s,
            }

        parsed = _parse_bbduk_summary(bbduk_lines)
        parsed["runtime_s"] = runtime_s
        return {"status": "success", "bbduk_stats": parsed}

    except Exception as e:
        return {"status": "error", "message": str(e)}


# -----------------------------------------------------------------------------
# FASTQ combination
# -----------------------------------------------------------------------------
def combine_gz_fastqs(input_files: List[str], combined_out: str) -> Dict[str, str]:
    """
    Concatenate multiple .gz FASTQ files into one .gz file.

    Important assumption:
    - gzip streams can be concatenated safely (standard behavior).
    - Downstream tools (zcat/gunzip/readers) will typically read concatenated streams.
    """
    combined_out_path = Path(combined_out)

    logger.info("Combine input files (%d):", len(input_files))
    total_in = 0

    for fp in input_files:
        size = Path(fp).stat().st_size
        total_in += size
        logger.info("  [combine] %s (%d bytes)", fp, size)

    # Stream copy (binary) to avoid any decoding issues
    with open(combined_out_path, "wb") as out_handle:
        for fp in input_files:
            with open(fp, "rb") as in_handle:
                out_handle.write(in_handle.read())

    out_size = combined_out_path.stat().st_size
    logger.info("  [combine] total_in=%d bytes, out=%d bytes", total_in, out_size)

    return {"status": "success"}


def _make_trimmed_name(input_name: str) -> str:
    """
    Create output filename for trimmed FASTQ.

    Keeps your original behavior for .fastq.gz.
    Adds a small improvement: also supports .fq.gz by falling back to a generic suffix.
    """
    if input_name.endswith(".fastq.gz"):
        return input_name.replace(".fastq.gz", ".trimmed.fastq.gz")
    if input_name.endswith(".fq.gz"):
        return input_name.replace(".fq.gz", ".trimmed.fq.gz")
    # Fallback: just append
    return input_name + ".trimmed.fastq.gz"


# -----------------------------------------------------------------------------
# Main pipeline
# -----------------------------------------------------------------------------
def main() -> None:
    """
    Pipeline overview:

    1) Read config from YAML
    2) Read sample sheet (sample_id, sample_name, fastq)
    3) For each sample_id:
        - create per-sample output dir
        - run bbduk on each listed FASTQ
        - write per-sample CSV summarizing each FASTQ
        - concatenate all trimmed FASTQs into one combined file
    4) Write a run-level CSV (all FASTQs across all samples)
    """
    parser = argparse.ArgumentParser(
        description="Trim single-end FASTQs and combine per sample_id."
    )
    parser.add_argument("--yaml_config", required=True)
    parser.add_argument("--sample_sheet", required=True)
    parser.add_argument("--sample_id", default=None)
    args = parser.parse_args()

    config = parse_arguments_from_yaml(args.yaml_config)

    output_dir = Path(config["output_dir"])
    bbmap_dir = str(config["bbmap_dir"])

    logs_base = Path(config["logs_dir"])
    logs_dir = logs_base / "Preprocessing"
    logs_dir.mkdir(parents=True, exist_ok=True)

    # Trimming parameters (with defaults)
    qtrim = config.get("qtrim", "r")
    quality_threshold = int(config.get("quality_threshold", 20))
    min_length = int(config.get("min_length", 190))
    bbduk_xmx_gb = int(config.get("bbduk_xmx_gb", 8))

    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_sample_sheet(
        args.sample_sheet,
        required_cols=("sample_id", "fastq", "sample_name"),
    )

    # Optional: process only one sample_id
    if args.sample_id is not None:
        df = df[df["sample_id"].astype(str) == str(args.sample_id)]
        if df.empty:
            logger.warning("No rows found for sample_id=%s", args.sample_id)
            return

    # Collect ALL per-FASTQ summary rows across the whole run
    run_summary_rows: List[Dict[str, Any]] = []

    # groupby ensures we process all FASTQs belonging to one sample together
    for sample_id, group in df.groupby("sample_id", sort=False):
        sample_id_str = str(sample_id)

        # Enforce 1:1 mapping sample_id -> sample_name
        sample_names = group["sample_name"].unique()
        if len(sample_names) != 1:
            raise ValueError(
                f"sample_id={sample_id_str} has multiple sample_name values: {sample_names}"
            )

        sample_name = safe_name(sample_names[0])
        sample_label = f"{sample_id_str}_{sample_name}"

        logger.info("=== Processing sample: %s (%d FASTQs) ===", sample_label, len(group))
        sample_log_handler = add_sample_logfile(sample_label, logs_dir)

        try:
            sample_out_dir = output_dir / sample_label
            sample_out_dir.mkdir(parents=True, exist_ok=True)

            trimmed_files: List[str] = []
            per_fastq_rows: List[Dict[str, Any]] = []

            # Iterate each FASTQ row for this sample
            for _, row in group.iterrows():
                in_fp = Path(row["fastq"])

                if not in_fp.exists():
                    logger.warning("Missing input FASTQ: %s", in_fp)
                    continue

                out_name = _make_trimmed_name(in_fp.name)
                out_fp = sample_out_dir / out_name

                res = trimming_single_end(
                    input_fastq_gz=str(in_fp),
                    output_fastq_gz=str(out_fp),
                    bbmap_dir=bbmap_dir,
                    qtrim=qtrim,
                    quality_threshold=quality_threshold,
                    min_length=min_length,
                    xmx_gb=bbduk_xmx_gb,
                )

                if res["status"] != "success":
                    logger.error("Trimming failed for %s: %s", in_fp, res.get("message"))
                    continue

                if not wait_for_file(out_fp):
                    continue

                out_size = out_fp.stat().st_size
                logger.info("[trim] output=%s size=%d bytes", out_fp, out_size)

                trimmed_files.append(str(out_fp))

                # Merge bbduk stats into the per-fastq summary row
                bbduk_stats = res.get("bbduk_stats", {})
                per_fastq_rows.append(
                    {
                        "sample_id": sample_id_str,
                        "sample_name": sample_name,
                        "sample_label": sample_label,
                        "input_fastq": str(in_fp),
                        "output_fastq": str(out_fp),
                        "output_size_bytes": out_size,
                        **bbduk_stats,
                    }
                )

            if not trimmed_files:
                logger.error("No trimmed files produced for sample_id=%s", sample_id_str)
                continue

            # Write per-sample bbduk summary CSV (per FASTQ, NOT aggregated)
            summary_csv = sample_out_dir / f"{sample_label}.bbduk_summary.csv"
            pd.DataFrame(per_fastq_rows).to_csv(summary_csv, index=False)
            logger.info("Wrote bbduk summary CSV: %s", summary_csv)

            # Add per-FASTQ rows to run-level summary
            run_summary_rows.extend(per_fastq_rows)

            # Combine all trimmed FASTQs for this sample
            combined_out = sample_out_dir / f"{sample_label}_combined_trimmed.fastq.gz"
            logger.info(
                "Combining %d trimmed FASTQs -> %s", len(trimmed_files), combined_out
            )
            combine_gz_fastqs(trimmed_files, str(combined_out))

            logger.info("Done sample_id=%s", sample_id_str)
            logger.info("Trimmed files in: %s", sample_out_dir)
            logger.info("Combined output : %s", combined_out)

        finally:
            remove_sample_logfile(sample_log_handler)

    # Write one CSV for the whole run
    if run_summary_rows:
        ts = datetime.now().strftime("%Y%m%d_%H-%M-%S")
        run_summary_csv = output_dir / f"bbduk_summary_{ts}.csv"
        pd.DataFrame(run_summary_rows).to_csv(run_summary_csv, index=False)
        logger.info("Wrote RUN-level bbduk summary CSV: %s", run_summary_csv)
    else:
        logger.warning("No bbduk summary rows collected; run-level CSV not written.")


if __name__ == "__main__":
    main()
