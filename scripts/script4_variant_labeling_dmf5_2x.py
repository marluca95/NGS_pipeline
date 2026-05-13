#!/usr/bin/env python3
import argparse
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from utils.config_utils import load_pipeline_config
from utils.logging_utils import setup_pipeline_logging

SCRIPT_NAME = "script4_variant_labeling_dmf5_2x"
logger = logging.getLogger(__name__)


def parse_yaml(yaml_path: str) -> dict:
    cfg = load_pipeline_config(
        yaml_path,
        required_keys=("input_dir", "output_dir"),
        default_values={
            "input_glob": "*.filtered.PASS.aa.tsv",
            "aa_column": "aa_seq",
            "pseudocount": 1e-9,
            "logs_subdir": "_logs",  # backward-compatible fallback
            "lib_suffix": "orilib",
            "neg_suffix": "2xnegative",
            "pos_suffix": "2xpositive",
            "required_conditions": ["neg", "pos"],
            "logs_dir": None,
            "status_up": 2.0,     # Enriched if >= 2
            "status_down": 0.5,   # Depleted if <= 0.5
            # Token prefix used to parse conditions from filenames.
            # Leave empty ("") to match any text before the condition suffix.
            "file_token_prefix": "",
        },
        path_keys=("input_dir", "output_dir", "logs_dir"),
    )
    if not cfg.get("logs_dir"):
        cfg["logs_dir"] = str(Path(cfg["output_dir"]) / str(cfg.get("logs_subdir", "_logs")))
    return cfg


def parse_condition_token(
    token: str,
    lib_suffix: str,
    neg_suffix: str,
    pos_suffix: str,
) -> Tuple[Optional[str], Optional[str]]:
    if token == lib_suffix:
        return "LIBRARY", "lib"
    if token.endswith(neg_suffix):
        return token[: -len(neg_suffix)], "neg"
    if token.endswith(pos_suffix):
        return token[: -len(pos_suffix)], "pos"
    return None, None


def discover_input_files(cfg: dict) -> Tuple[Path, Dict[str, Dict[str, Path]], List[Path]]:
    input_dir = Path(cfg["input_dir"])
    files = sorted(input_dir.glob(cfg["input_glob"]))
    if not files:
        raise FileNotFoundError(f"No input files found in {input_dir} with glob {cfg['input_glob']!r}")

    groups: Dict[str, Dict[str, Path]] = {}
    unknown: List[Path] = []
    library_file: Optional[Path] = None

    token_prefix = cfg.get("file_token_prefix", "")
    prefix_pattern = re.escape(token_prefix) if token_prefix else ""
    file_re = re.compile(rf"_{prefix_pattern}(.+)\.filtered\.PASS\.aa\.tsv$")

    for path in files:
        m = file_re.search(path.name)
        if not m:
            unknown.append(path)
            continue

        token = m.group(1)
        peptide_key, condition = parse_condition_token(
            token=token,
            lib_suffix=cfg["lib_suffix"],
            neg_suffix=cfg["neg_suffix"],
            pos_suffix=cfg["pos_suffix"],
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


# --- NEW: same rules as variant_analysis.ipynb ---
def flag(val: float, up: float = 2.0, down: float = 0.5) -> str:
    if val >= up:
        return "Enriched"
    if val <= down:
        return "Depleted"
    return "NoChange"


def annotate_row(row: pd.Series) -> int:
    """
    Exactly as in variant_analysis.ipynb:

    if ((pos_vs_neg_status == Enriched OR pos_vs_lib_status == Enriched) AND neg_vs_lib_status != Enriched) -> 1
    elif ((pos_vs_neg_status != Enriched AND pos_vs_lib_status != Enriched) AND neg_vs_lib_status == Enriched) -> 0
    else -> 2
    """
    if (
        (row["pos_vs_neg_status"] == "Enriched" or row["pos_vs_lib_status"] == "Enriched")
        and row["neg_vs_lib_status"] != "Enriched"
    ):
        return 1
    if (
        (row["pos_vs_neg_status"] != "Enriched" and row["pos_vs_lib_status"] != "Enriched")
        and row["neg_vs_lib_status"] == "Enriched"
    ):
        return 0
    return 2

def add_status_and_specificity(merged: pd.DataFrame, cfg: dict) -> None:
    up = float(cfg.get("status_up", 2.0))
    down = float(cfg.get("status_down", 0.5))

    # Always required for the single-positive (2xpositive) workflow:
    base_needed = [
        "enrich_pos_vs_neg",
        "enrich_pos_vs_lib",
        "deplete_neg_lib",
    ]
    missing = [c for c in base_needed if c not in merged.columns]
    if missing:
        raise ValueError(f"Missing required enrichment columns for labeling: {missing}")

    merged["pos_vs_neg_status"] = merged["enrich_pos_vs_neg"].apply(lambda v: flag(float(v), up, down))
    merged["pos_vs_lib_status"] = merged["enrich_pos_vs_lib"].apply(lambda v: flag(float(v), up, down))

    merged["neg_vs_lib_status"] = merged["deplete_neg_lib"].apply(lambda v: flag(float(v), up, down))
    merged["specificity"] = merged.apply(annotate_row, axis=1)

# def add_status_and_specificity(merged: pd.DataFrame, cfg: dict) -> None:
#     """
#     Adds:
#       - deplete_neg_lib (same naming style as your old specificity CSV)
#       - pos_vs_neg_status, pos_vs_lib_status, neg_vs_lib_status
#       - specificity (0/1/2) using the same logic as variant_analysis.ipynb

#     Because you have *two* positive conditions (pos3x and pos1x), we:
#       - compute statuses for each (pos3x vs neg, pos1x vs neg, pos3x vs lib, pos1x vs lib)
#       - then combine them into one pos_vs_neg_status / pos_vs_lib_status
#         using combine_status() (Enriched if either enriched; Depleted if both depleted; else NoChange).
#     """
#     up = float(cfg.get("status_up", 2.0))
#     down = float(cfg.get("status_down", 0.5))

#     # Ensure needed enrichment columns exist
#     needed = [
#         "enrich_pos3x_vs_neg",
#         "enrich_pos1x_vs_neg",
#         "enrich_pos3x_vs_lib",
#         "enrich_pos1x_vs_lib",
#         "deplete_neg_lib",
#     ]
#     missing = [c for c in needed if c not in merged.columns]
#     if missing:
#         raise ValueError(f"Missing required enrichment columns for labeling: {missing}")

#     # Per-condition statuses
#     merged["pos3x_vs_neg_status"] = merged["enrich_pos3x_vs_neg"].apply(lambda v: flag(float(v), up, down))
#     merged["pos1x_vs_neg_status"] = merged["enrich_pos1x_vs_neg"].apply(lambda v: flag(float(v), up, down))
#     merged["pos3x_vs_lib_status"] = merged["enrich_pos3x_vs_lib"].apply(lambda v: flag(float(v), up, down))
#     merged["pos1x_vs_lib_status"] = merged["enrich_pos1x_vs_lib"].apply(lambda v: flag(float(v), up, down))

#     # Combined "pos" statuses (because you have two pos conditions)
#     merged["pos_vs_neg_status"] = [
#         combine_status(a, b) for a, b in zip(merged["pos3x_vs_neg_status"], merged["pos1x_vs_neg_status"])
#     ]
#     merged["pos_vs_lib_status"] = [
#         combine_status(a, b) for a, b in zip(merged["pos3x_vs_lib_status"], merged["pos1x_vs_lib_status"])
#     ]

#     # neg_vs_lib_status uses deplete_neg_lib (neg vs lib fold-change)
#     merged["neg_vs_lib_status"] = merged["deplete_neg_lib"].apply(lambda v: flag(float(v), up, down))

#     # Specificity assignment (0/1/2)
#     merged["specificity"] = merged.apply(annotate_row, axis=1)


def build_variant_table_for_peptide(peptide_key: str, files: Dict[str, Path], cfg: dict) -> pd.DataFrame:
    aa_column = cfg["aa_column"]
    labels = ["lib", "neg", "pos"]
    dfs = []
    for label in labels:
        if label not in files:
            continue
        logger.info("Counting %s for peptide %s from %s", label, peptide_key, files[label].name)
        dfs.append(count_variants(files[label], aa_column=aa_column, label=label))

    merged = merge_counts(dfs, aa_column=aa_column)
    merged["peptide"] = peptide_key

    pseudocount = float(cfg["pseudocount"])
    add_enrichment_columns(merged, pseudocount, numerator="pos", denominator="neg")
    add_enrichment_columns(merged, pseudocount, numerator="pos", denominator="lib")

    # NEW: neg vs lib enrichment (named deplete_neg_lib to match your earlier CSV naming)
    add_enrichment_columns(merged, pseudocount, numerator="neg", denominator="lib")
    if "enrich_neg_vs_lib" in merged.columns:
        merged.rename(columns={"enrich_neg_vs_lib": "deplete_neg_lib"}, inplace=True)

    # NEW: statuses + specificity (0/1/2) using notebook rules
    add_status_and_specificity(merged, cfg)

    return merged


def sanitize_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name).strip())
    return cleaned.strip("_") or "peptide"


def run(cfg: dict) -> None:
    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = setup_pipeline_logging(
        logs_dir=cfg["logs_dir"],
        script_name=SCRIPT_NAME,
        scope="run",
        run_label=cfg.get("run_label"),
    )
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

        # Log specificity distribution
        spec_counts = table["specificity"].value_counts(dropna=False).to_dict()
        logger.info("Specificity distribution for %s: %s", peptide_key, spec_counts)

        out_file = output_dir / f"{sanitize_name(peptide_key)}.variant_labeling.csv"
        table.sort_values(cfg["aa_column"]).to_csv(out_file, index=False)
        processed += 1

        summary_rows.append(
            {
                "peptide": peptide_key,
                "n_variants": int(table.shape[0]),
                "n_spec_1": int((table["specificity"] == 1).sum()),
                "n_spec_0": int((table["specificity"] == 0).sum()),
                "n_spec_2": int((table["specificity"] == 2).sum()),
                "output_csv": str(out_file),
                "input_lib": str(files["lib"]),
                "input_neg": str(files["neg"]),
                "input_pos": str(files.get("pos", "")),
            }
        )
        logger.info("Wrote peptide table: %s (%d variants)", out_file, table.shape[0])

    summary_path = output_dir / "variant_labeling_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    logger.info("Wrote summary table: %s", summary_path)
    logger.info("Done. Processed=%d, Skipped=%d", processed, skipped)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate per-peptide variant labeling CSVs (with specificity 0/1/2) from extracted PASS AA tables "
            "using orilib, 2xnegative, and 2xpositive only."
        )
    )
    parser.add_argument("--yaml_config", type=str, default=None, help="Path to YAML config file.")
    parser.add_argument(
        "--yaml_file",
        dest="yaml_config_legacy",
        type=str,
        default=None,
        help="Legacy alias of --yaml_config",
    )
    args = parser.parse_args()

    yaml_path = args.yaml_config or args.yaml_config_legacy
    if not yaml_path:
        raise ValueError("A config file must be supplied via --yaml_config.")
    cfg = parse_yaml(yaml_path)
    run(cfg)


if __name__ == "__main__":
    main()
