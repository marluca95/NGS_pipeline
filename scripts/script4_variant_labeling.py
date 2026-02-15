#!/usr/bin/env python3
import argparse
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from utils.config_utils import load_yaml_config, require_config_keys


logger = logging.getLogger(__name__)


def setup_logging(output_dir: Path, logs_subdir: str, log_name: str) -> Path:
    log_dir = output_dir / logs_subdir
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / log_name

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
        force=True,
    )
    return log_path


def parse_yaml(yaml_path: str) -> dict:
    cfg = load_yaml_config(yaml_path)
    require_config_keys(cfg, ["input_dir", "output_dir"])

    cfg.setdefault("input_glob", "*.filtered.PASS.aa.tsv")
    cfg.setdefault("aa_column", "aa_seq")
    cfg.setdefault("pseudocount", 1e-9)
    cfg.setdefault("logs_subdir", "_logs")
    cfg.setdefault("log_name", "script4_variant_labeling.log")

    cfg.setdefault("lib_suffix", "originallib")
    cfg.setdefault("neg_suffix", "2xnegative")
    cfg.setdefault("pos3x_suffix", "3xpositive")
    cfg.setdefault("pos1x_suffix", "31xpositive")
    cfg.setdefault("required_conditions", ["neg", "pos3x", "pos1x"])

    return cfg


def parse_condition_token(
    token: str,
    lib_suffix: str,
    neg_suffix: str,
    pos3x_suffix: str,
    pos1x_suffix: str,
) -> Tuple[Optional[str], Optional[str]]:
    if token == lib_suffix:
        return "LIBRARY", "lib"
    if token.endswith(neg_suffix):
        return token[: -len(neg_suffix)], "neg"
    if token.endswith(pos3x_suffix):
        return token[: -len(pos3x_suffix)], "pos3x"
    if token.endswith(pos1x_suffix):
        return token[: -len(pos1x_suffix)], "pos1x"
    return None, None


def discover_input_files(cfg: dict) -> Tuple[Path, Dict[str, Dict[str, Path]], List[Path]]:
    input_dir = Path(cfg["input_dir"])
    files = sorted(input_dir.glob(cfg["input_glob"]))
    if not files:
        raise FileNotFoundError(f"No input files found in {input_dir} with glob {cfg['input_glob']!r}")

    groups: Dict[str, Dict[str, Path]] = {}
    unknown: List[Path] = []
    library_file: Optional[Path] = None

    for path in files:
        m = re.search(r"_clib(.+)\.filtered\.PASS\.aa\.tsv$", path.name)
        if not m:
            unknown.append(path)
            continue

        token = m.group(1)
        peptide_key, condition = parse_condition_token(
            token=token,
            lib_suffix=cfg["lib_suffix"],
            neg_suffix=cfg["neg_suffix"],
            pos3x_suffix=cfg["pos3x_suffix"],
            pos1x_suffix=cfg["pos1x_suffix"],
        )

        if peptide_key is None or condition is None:
            unknown.append(path)
            continue

        if condition == "lib":
            library_file = path
            continue

        if peptide_key not in groups:
            groups[peptide_key] = {}
        groups[peptide_key][condition] = path

    if library_file is None:
        raise FileNotFoundError(
            f"Could not find library file using suffix '{cfg['lib_suffix']}' in {input_dir}"
        )

    for peptide_key in groups:
        groups[peptide_key]["lib"] = library_file

    return library_file, groups, unknown


def count_variants(path: Path, aa_column: str, label: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if aa_column not in df.columns:
        raise ValueError(f"{path} is missing required column {aa_column!r}. Columns: {list(df.columns)}")

    counts = df.groupby(aa_column).size().reset_index(name=f"count_{label}")
    total = counts[f"count_{label}"].sum()
    if total == 0:
        counts[f"freq_{label}"] = 0.0
    else:
        counts[f"freq_{label}"] = counts[f"count_{label}"] / total
    return counts


def merge_counts(dfs: List[pd.DataFrame], aa_column: str) -> pd.DataFrame:
    out = dfs[0]
    for df in dfs[1:]:
        out = out.merge(df, on=aa_column, how="outer")
    return out.fillna(0)


def add_enrichment_columns(merged: pd.DataFrame, pseudocount: float, numerator: str, denominator: str) -> None:
    ncol = f"freq_{numerator}"
    dcol = f"freq_{denominator}"
    if ncol in merged.columns and dcol in merged.columns:
        out_col = f"enrich_{numerator}_vs_{denominator}"
        merged[out_col] = (merged[ncol] + pseudocount) / (merged[dcol] + pseudocount)


def build_variant_table_for_peptide(peptide_key: str, files: Dict[str, Path], cfg: dict) -> pd.DataFrame:
    aa_column = cfg["aa_column"]
    labels = ["lib", "neg", "pos3x", "pos1x"]
    dfs = []
    for label in labels:
        if label not in files:
            continue
        logger.info("Counting %s for peptide %s from %s", label, peptide_key, files[label].name)
        dfs.append(count_variants(files[label], aa_column=aa_column, label=label))

    merged = merge_counts(dfs, aa_column=aa_column)
    merged["peptide"] = peptide_key

    pseudocount = float(cfg["pseudocount"])
    add_enrichment_columns(merged, pseudocount, numerator="pos3x", denominator="neg")
    add_enrichment_columns(merged, pseudocount, numerator="pos1x", denominator="neg")
    add_enrichment_columns(merged, pseudocount, numerator="pos3x", denominator="lib")
    add_enrichment_columns(merged, pseudocount, numerator="pos1x", denominator="lib")
    return merged


def sanitize_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name).strip())
    return cleaned.strip("_") or "peptide"


def run(cfg: dict) -> None:
    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = setup_logging(output_dir, cfg["logs_subdir"], cfg["log_name"])
    logger.info("Loaded config and started logging to %s", log_path)

    library_file, groups, unknown = discover_input_files(cfg)
    logger.info("Detected shared library file: %s", library_file.name)
    logger.info("Detected %d peptide group(s)", len(groups))
    if unknown:
        logger.warning("Ignored %d file(s) that did not match configured suffix patterns", len(unknown))

    required = set(cfg["required_conditions"])
    summary_rows = []
    processed = 0
    skipped = 0

    for peptide_key in sorted(groups):
        files = groups[peptide_key]
        missing = sorted(required - set(files.keys()))
        if missing:
            logger.warning("Skipping peptide %s because required conditions are missing: %s", peptide_key, missing)
            skipped += 1
            continue

        table = build_variant_table_for_peptide(peptide_key, files, cfg)
        out_file = output_dir / f"{sanitize_name(peptide_key)}.variant_labeling.csv"
        table.sort_values(cfg["aa_column"]).to_csv(out_file, index=False)
        processed += 1

        summary_rows.append(
            {
                "peptide": peptide_key,
                "n_variants": int(table.shape[0]),
                "output_csv": str(out_file),
                "input_lib": str(files["lib"]),
                "input_neg": str(files["neg"]),
                "input_pos3x": str(files["pos3x"]),
                "input_pos1x": str(files["pos1x"]),
            }
        )
        logger.info("Wrote peptide table: %s (%d variants)", out_file, table.shape[0])

    summary_path = output_dir / "script4_variant_labeling_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    logger.info("Wrote summary table: %s", summary_path)
    logger.info("Done. Processed=%d, Skipped=%d", processed, skipped)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate per-peptide variant labeling CSVs from extracted PASS AA tables."
    )
    parser.add_argument(
        "--yaml_file",
        type=str,
        default="/cluster/project/reddy/katja/NGS_pipeline/config/script4_variant_labeling.yaml",
        help="Path to YAML config file.",
    )
    args = parser.parse_args()

    cfg = parse_yaml(args.yaml_file)
    run(cfg)


if __name__ == "__main__":
    main()
