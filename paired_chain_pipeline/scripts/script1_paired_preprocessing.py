"""
Paired-Chain Pipeline — Step 01: Paired-End Quality Trimming
=============================================================

Adapts script1_preprocessing.py for paired-end (R1 + R2) sequencing.

Key differences vs. single-end:
  - Groups sample sheet rows by (sample_id, lane) to locate R1 and R2 files.
  - Runs BBDuk in paired mode: in={r1} in2={r2} out={r1_out} out2={r2_out}
    BBDuk in paired mode guarantees R1 and R2 outputs stay in lockstep
    (if one read is dropped, its pair is also dropped).
  - Concatenates all lane R1 trimmed files into {sample_label}_combined_trimmed_R1.fastq.gz
    and lane R2 trimmed files into {sample_label}_combined_trimmed_R2.fastq.gz.

Outputs per sample in output_dir/{sample_label}/:
  {sample_label}_combined_trimmed_R1.fastq.gz
  {sample_label}_combined_trimmed_R2.fastq.gz
  {sample_label}.bbduk_summary.csv
"""

import argparse
import logging
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from utils.config_utils import load_pipeline_config
from utils.logging_utils import (
    add_logger_file_handler,
    remove_logger_handler,
    setup_pipeline_logging,
)
from utils.metrics_utils import write_sample_metrics
from utils.sample_utils import load_sample_sheet, safe_name

SCRIPT_NAME = "script1_paired_preprocessing"
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Per-sample log helpers (identical to script1)
# ---------------------------------------------------------------------------
def add_sample_logfile(sample_label: str, log_dir: Path, run_label: Optional[str] = None) -> logging.FileHandler:
    handler, log_path = add_logger_file_handler(
        logger=logger,
        logs_dir=log_dir,
        script_name=SCRIPT_NAME,
        scope=f"sample_{sample_label}",
        run_label=run_label,
    )
    logger.info("Per-sample logfile: %s", log_path)
    return handler


def remove_sample_logfile(handler: logging.FileHandler) -> None:
    remove_logger_handler(logger, handler)


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
def parse_arguments_from_yaml(yaml_file: str) -> Dict[str, Any]:
    return load_pipeline_config(
        yaml_file,
        required_keys=("output_dir", "bbmap_dir"),
        default_values={
            "qtrim": "r",
            "quality_threshold": 20,
            "min_length": 101,
            "bbduk_xmx_gb": 8,
            "logs_dir": None,
        },
        path_keys=("output_dir", "bbmap_dir", "logs_dir"),
    )


# ---------------------------------------------------------------------------
# BBDuk helpers
# ---------------------------------------------------------------------------
def wait_for_file(filepath: Path, timeout: int = 600, interval: int = 5) -> bool:
    start_time = time.time()
    while not filepath.exists():
        if time.time() - start_time > timeout:
            logger.error("Timeout: %s not found within %d seconds.", filepath, timeout)
            return False
        time.sleep(interval)
    return True


def _parse_bbduk_summary(bbduk_lines: List[str]) -> Dict[str, Any]:
    stats: Dict[str, Any] = {}
    for line in bbduk_lines:
        m = re.search(r"\bVersion\s+(\S+)", line)
        if m:
            stats["bbduk_version"] = m.group(1)
            break

    def grab_counts(prefix: str) -> Tuple[Optional[int], Optional[int]]:
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
    stats.update({
        "reads_in": reads_in, "bases_in": bases_in,
        "reads_out": reads_out, "bases_out": bases_out,
        "reads_removed": reads_removed, "bases_removed": bases_removed,
        "reads_qtrimmed": reads_qtrimmed, "bases_qtrimmed": bases_qtrimmed,
    })
    return stats


def _build_bbduk_paired_cmd(
    bbmap_dir: str,
    r1_in: str, r2_in: str,
    r1_out: str, r2_out: str,
    qtrim: str, quality_threshold: int,
    min_length: int, xmx_gb: int,
) -> List[str]:
    return [
        f"{bbmap_dir}/bbduk.sh",
        f"-Xmx{xmx_gb}g",
        f"in={r1_in}", f"in2={r2_in}",
        f"out={r1_out}", f"out2={r2_out}",
        f"qtrim={qtrim}",
        f"trimq={quality_threshold}",
        f"minlength={min_length}",
    ]


def trimming_paired(
    r1_in: str, r2_in: str,
    r1_out: str, r2_out: str,
    bbmap_dir: str,
    qtrim: str, quality_threshold: int,
    min_length: int, xmx_gb: int,
) -> Dict[str, Any]:
    cmd = _build_bbduk_paired_cmd(
        bbmap_dir=bbmap_dir,
        r1_in=r1_in, r2_in=r2_in,
        r1_out=r1_out, r2_out=r2_out,
        qtrim=qtrim, quality_threshold=quality_threshold,
        min_length=min_length, xmx_gb=xmx_gb,
    )
    logger.info("Running: %s", " ".join(cmd))
    t0 = time.time()
    bbduk_lines: List[str] = []
    try:
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
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
            return {"status": "error", "message": f"bbduk exited with code {ret}", "runtime_s": runtime_s}
        parsed = _parse_bbduk_summary(bbduk_lines)
        parsed["runtime_s"] = runtime_s
        return {"status": "success", "bbduk_stats": parsed}
    except Exception as e:
        return {"status": "error", "message": str(e)}


# ---------------------------------------------------------------------------
# FASTQ combination (identical to script1)
# ---------------------------------------------------------------------------
def combine_gz_fastqs(input_files: List[str], combined_out: str) -> Dict[str, str]:
    combined_out_path = Path(combined_out)
    logger.info("Combining %d files into %s", len(input_files), combined_out_path)
    with open(combined_out_path, "wb") as out_handle:
        for fp in input_files:
            with open(fp, "rb") as in_handle:
                shutil.copyfileobj(in_handle, out_handle, length=1024 * 1024)
    return {"status": "success"}


def _make_trimmed_name(input_name: str, suffix: str) -> str:
    """Create output filename for trimmed paired FASTQ, appending a suffix (e.g. '_R1')."""
    base = input_name
    if base.endswith(".fastq.gz"):
        base = base[:-9]
    elif base.endswith(".fq.gz"):
        base = base[:-6]
    return f"{base}{suffix}.trimmed.fastq.gz"


# ---------------------------------------------------------------------------
# Per-sample processing
# ---------------------------------------------------------------------------
def process_one_sample(sample_id: str, sample_name: str, df_sample: pd.DataFrame, cfg: dict) -> Dict[str, Any]:
    """
    Process all lanes for a single sample in paired-end mode.

    df_sample: rows from sample sheet for this sample_id (contains both read=1 and read=2 rows).
    """
    sample_label = f"{sample_id}_{sample_name}"
    out_dir = Path(cfg["output_dir"]) / sample_label
    out_dir.mkdir(parents=True, exist_ok=True)

    log_dir = Path(cfg.get("logs_dir") or (Path(cfg["output_dir"]) / "_logs"))
    run_label = cfg.get("run_label")
    handler = add_sample_logfile(sample_label, log_dir, run_label)

    bbmap_dir = cfg["bbmap_dir"]
    qtrim = cfg.get("qtrim", "r")
    quality_threshold = int(cfg.get("quality_threshold", 20))
    min_length = int(cfg.get("min_length", 101))
    xmx_gb = int(cfg.get("bbduk_xmx_gb", 8))

    # Group rows by lane, then find R1 (read=1) and R2 (read=2) for each lane
    lanes = df_sample["lane"].unique()
    lane_r1_trimmed: List[str] = []
    lane_r2_trimmed: List[str] = []
    summary_rows: List[Dict[str, Any]] = []

    total_reads_in = 0
    total_reads_out = 0

    for lane in sorted(lanes):
        df_lane = df_sample[df_sample["lane"].astype(str) == str(lane)]
        r1_rows = df_lane[df_lane["read"].astype(str) == "1"]
        r2_rows = df_lane[df_lane["read"].astype(str) == "2"]

        if r1_rows.empty or r2_rows.empty:
            logger.warning("Lane %s for sample %s is missing R1 or R2; skipping.", lane, sample_label)
            continue

        r1_path = str(r1_rows.iloc[0]["fastq"])
        r2_path = str(r2_rows.iloc[0]["fastq"])

        r1_name = Path(r1_path).name
        r2_name = Path(r2_path).name
        r1_out = str(out_dir / _make_trimmed_name(r1_name, "_R1"))
        r2_out = str(out_dir / _make_trimmed_name(r2_name, "_R2"))

        logger.info("Lane %s — R1: %s", lane, r1_path)
        logger.info("Lane %s — R2: %s", lane, r2_path)

        result = trimming_paired(
            r1_in=r1_path, r2_in=r2_path,
            r1_out=r1_out, r2_out=r2_out,
            bbmap_dir=bbmap_dir, qtrim=qtrim,
            quality_threshold=quality_threshold,
            min_length=min_length, xmx_gb=xmx_gb,
        )

        if result["status"] != "success":
            raise RuntimeError(f"BBDuk failed for {sample_label} lane {lane}: {result.get('message')}")

        wait_for_file(Path(r1_out))
        wait_for_file(Path(r2_out))

        lane_r1_trimmed.append(r1_out)
        lane_r2_trimmed.append(r2_out)

        stats = result["bbduk_stats"]
        summary_rows.append({
            "sample_id": sample_id, "sample_name": sample_name,
            "lane": lane, "r1_input": r1_path, "r2_input": r2_path,
            **{f"bbduk_{k}": v for k, v in stats.items()},
        })
        total_reads_in += stats.get("reads_in") or 0
        total_reads_out += stats.get("reads_out") or 0

    if not lane_r1_trimmed:
        raise ValueError(f"No lanes processed for sample {sample_label}")

    # Combine lanes into single R1/R2 files
    combined_r1 = str(out_dir / f"{sample_label}_combined_trimmed_R1.fastq.gz")
    combined_r2 = str(out_dir / f"{sample_label}_combined_trimmed_R2.fastq.gz")

    combine_gz_fastqs(lane_r1_trimmed, combined_r1)
    combine_gz_fastqs(lane_r2_trimmed, combined_r2)

    logger.info("Combined R1: %s", combined_r1)
    logger.info("Combined R2: %s", combined_r2)

    # Write per-sample BBDuk summary CSV
    summary_path = out_dir / f"{sample_label}.bbduk_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)

    remove_sample_logfile(handler)

    return {
        "n_input_fastqs": len(lane_r1_trimmed) * 2,
        "reads_in_total": total_reads_in,
        "reads_out_total": total_reads_out,
        "reads_lost_trimming": total_reads_in - total_reads_out,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="Paired-Chain Pipeline Step 01: Paired-end quality trimming with BBDuk."
    )
    ap.add_argument("--yaml_config", required=True)
    ap.add_argument(
        "--sample_id", default=None,
        help="Process only this sample_id (optional; default processes all samples).",
    )
    args = ap.parse_args()

    cfg = parse_arguments_from_yaml(args.yaml_config)

    logs_dir = Path(cfg.get("logs_dir") or (Path(cfg["output_dir"]) / "_logs"))
    log_file = setup_pipeline_logging(
        logs_dir=logs_dir,
        script_name=SCRIPT_NAME,
        scope="run",
        run_label=cfg.get("run_label"),
    )
    logger.info("Logging to: %s", log_file)
    logger.info("Loaded config: %s", cfg.get("_config_path", args.yaml_config))

    sample_sheet = cfg.get("sample_sheet")
    if not sample_sheet:
        raise ValueError("Config must specify 'sample_sheet'.")

    df = load_sample_sheet(sample_sheet, required_cols=("sample_id", "sample_name", "lane", "read", "fastq"))

    if args.sample_id is not None:
        df = df[df["sample_id"].astype(str) == str(args.sample_id)]
        if df.empty:
            raise ValueError(f"No rows for sample_id={args.sample_id}")

    metrics_path = Path(cfg["output_dir"]) / "per_sample_metrics.tsv"
    run_summary_rows: List[Dict[str, Any]] = []

    for sample_id, g in df.groupby("sample_id", sort=False):
        names = g["sample_name"].dropna().unique()
        if len(names) != 1:
            raise ValueError(f"sample_id={sample_id} has multiple sample_name values: {names}")
        sample_name = safe_name(str(names[0]))

        logger.info("=== Sample: %s_%s ===", sample_id, sample_name)
        stats = process_one_sample(str(sample_id), sample_name, g, cfg)

        write_sample_metrics(
            metrics_path=metrics_path,
            sample_id=str(sample_id),
            sample_name=sample_name,
            step="01_preprocessing",
            metrics=stats,
        )
        run_summary_rows.append({"sample_id": sample_id, "sample_name": sample_name, **stats})

    run_label = cfg.get("run_label", "paired_run")
    run_summary_path = Path(cfg["output_dir"]) / f"bbduk_run_summary_{run_label}.csv"
    pd.DataFrame(run_summary_rows).to_csv(run_summary_path, index=False)
    logger.info("Run summary written: %s", run_summary_path)


if __name__ == "__main__":
    main()
