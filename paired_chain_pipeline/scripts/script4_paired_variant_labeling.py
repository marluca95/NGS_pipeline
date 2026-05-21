"""
Paired-Chain Pipeline — Step 04: Paired Variant Labeling
=========================================================

Adapts script4_variant_labeling.py for paired (alpha, beta) variants.

The variant unit here is a (beta_aa, alpha_aa) pair rather than a single
amino acid sequence. A combined key "{beta_aa}|{alpha_aa}" is used for
counting and frequency calculations. Enrichment and specificity logic is
identical to the single-chain pipeline.

Input files: *.filtered.PASS.paired.tsv  (from script3_paired_filtering_extraction.py)
Required columns in each file: beta_aa, alpha_aa

Output per peptide group: {peptide}.variant_labeling.csv
  Columns include: beta_aa, alpha_aa, variant_key, count_*, freq_*, enrich_*, specificity
"""

import argparse
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from utils.config_utils import load_pipeline_config
from utils.logging_utils import setup_pipeline_logging

SCRIPT_NAME = "script4_paired_variant_labeling"
logger = logging.getLogger(__name__)

PAIRED_TSV_SUFFIX = ".filtered.PASS.paired.tsv"
VARIANT_KEY_COL = "variant_key"


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
def parse_yaml(yaml_path: str) -> dict:
    cfg = load_pipeline_config(
        yaml_path,
        required_keys=("input_dir", "output_dir"),
        default_values={
            "input_glob": f"*{PAIRED_TSV_SUFFIX}",
            "beta_aa_column":  "beta_aa",
            "alpha_aa_column": "alpha_aa",
            "pseudocount": 1e-9,
            "logs_subdir": "_logs",
            "lib_suffix": "originallib",
            "neg_suffix": "2xnegative",
            "pos3x_suffix": "3xpositive",
            "pos1x_suffix": "31xpositive",
            "required_conditions": ["neg", "pos3x"],
            "logs_dir": None,
            "status_up": 2.0,
            "status_down": 0.5,
            "file_token_prefix": "",
        },
        path_keys=("input_dir", "output_dir", "logs_dir"),
    )
    if not cfg.get("logs_dir"):
        cfg["logs_dir"] = str(Path(cfg["output_dir"]) / str(cfg.get("logs_subdir", "_logs")))
    return cfg


# ---------------------------------------------------------------------------
# Condition parsing (identical to script4)
# ---------------------------------------------------------------------------
def parse_condition_token(
    token: str,
    lib_suffix: str, neg_suffix: str,
    pos3x_suffix: str, pos1x_suffix: str,
    file_token_prefix: str = "",
) -> Tuple[Optional[str], Optional[str]]:
    if file_token_prefix and token.startswith(file_token_prefix):
        token = token[len(file_token_prefix):].lstrip("-_")
    if token.endswith(lib_suffix):
        return token[: -len(lib_suffix)].rstrip("-_") or "lib", "lib"
    if token.endswith(neg_suffix):
        return token[: -len(neg_suffix)].rstrip("-_") or token, "neg"
    if token.endswith(pos3x_suffix):
        return token[: -len(pos3x_suffix)].rstrip("-_") or token, "pos3x"
    if token.endswith(pos1x_suffix):
        return token[: -len(pos1x_suffix)].rstrip("-_") or token, "pos1x"
    return None, None


def discover_input_files(cfg: dict) -> Tuple[Path, Dict[str, Dict[str, Path]], List[Path]]:
    input_dir = Path(cfg["input_dir"])
    files = sorted(input_dir.glob(cfg["input_glob"]))
    if not files:
        raise FileNotFoundError(
            f"No input files found in {input_dir} with glob {cfg['input_glob']!r}"
        )

    groups: Dict[str, Dict[str, Path]] = {}
    unknown: List[Path] = []
    library_file: Optional[Path] = None

    for path in files:
        if not path.name.endswith(PAIRED_TSV_SUFFIX):
            unknown.append(path)
            continue
        stem = path.name[: -len(PAIRED_TSV_SUFFIX)]
        if "_" not in stem:
            unknown.append(path)
            continue
        _, token = stem.split("_", 1)
        peptide_key, condition = parse_condition_token(
            token=token,
            lib_suffix=cfg["lib_suffix"], neg_suffix=cfg["neg_suffix"],
            pos3x_suffix=cfg["pos3x_suffix"], pos1x_suffix=cfg["pos1x_suffix"],
            file_token_prefix=cfg["file_token_prefix"],
        )
        if peptide_key is None or condition is None:
            unknown.append(path)
            continue
        if condition == "lib":
            library_file = path
            continue
        groups.setdefault(peptide_key, {})[condition] = path

    if library_file is None:
        raise FileNotFoundError(
            f"Could not find library file using suffix '{cfg['lib_suffix']}' in {input_dir}"
        )
    for peptide_key in groups:
        groups[peptide_key]["lib"] = library_file

    return library_file, groups, unknown


# ---------------------------------------------------------------------------
# Variant counting for paired TSV
# ---------------------------------------------------------------------------
def count_paired_variants(path: Path, beta_col: str, alpha_col: str, label: str) -> pd.DataFrame:
    """
    Read a paired TSV, create a combined variant_key = "{beta_aa}|{alpha_aa}",
    count occurrences, and compute frequency.
    """
    df = pd.read_csv(path, sep="\t")
    for col in (beta_col, alpha_col):
        if col not in df.columns:
            raise ValueError(f"{path} is missing required column {col!r}. Columns: {list(df.columns)}")

    df[VARIANT_KEY_COL] = df[beta_col].astype(str) + "|" + df[alpha_col].astype(str)

    # Also keep beta_aa and alpha_aa for output
    key_cols = [VARIANT_KEY_COL, beta_col, alpha_col]
    counts = (
        df.groupby(key_cols, as_index=False)
        .size()
        .rename(columns={"size": f"count_{label}"})
    )
    total = counts[f"count_{label}"].sum()
    counts[f"freq_{label}"] = counts[f"count_{label}"] / total if total > 0 else 0.0
    return counts


def merge_counts(dfs: List[pd.DataFrame]) -> pd.DataFrame:
    out = dfs[0]
    for df in dfs[1:]:
        # Merge on variant_key; beta_aa and alpha_aa are repeated but identical
        merge_cols = [VARIANT_KEY_COL]
        # Add beta_aa / alpha_aa to merge keys if present in both frames
        for c in ("beta_aa", "alpha_aa"):
            if c in out.columns and c in df.columns:
                merge_cols.append(c)
        # Drop duplicate columns from right frame before merging
        df = df.drop(columns=[c for c in merge_cols[1:] if c in df.columns and c in out.columns],
                     errors="ignore")
        out = out.merge(df, on=VARIANT_KEY_COL, how="outer")
    return out.fillna(0)


# ---------------------------------------------------------------------------
# Enrichment + specificity (identical logic to script4)
# ---------------------------------------------------------------------------
def add_enrichment_columns(merged: pd.DataFrame, pseudocount: float,
                            numerator: str, denominator: str) -> None:
    ncol = f"freq_{numerator}"
    dcol = f"freq_{denominator}"
    if ncol in merged.columns and dcol in merged.columns:
        merged[f"enrich_{numerator}_vs_{denominator}"] = (
            (merged[ncol] + pseudocount) / (merged[dcol] + pseudocount)
        )


def flag(val: float, up: float = 2.0, down: float = 0.5) -> str:
    if val >= up: return "Enriched"
    if val <= down: return "Depleted"
    return "NoChange"


def combine_status(a: str, b: str) -> str:
    if a == "Enriched" or b == "Enriched": return "Enriched"
    if a == "Depleted" and b == "Depleted": return "Depleted"
    return "NoChange"


def annotate_row(row: pd.Series) -> int:
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

    base_needed = ["enrich_pos3x_vs_neg", "enrich_pos3x_vs_lib", "deplete_neg_lib"]
    missing = [c for c in base_needed if c not in merged.columns]
    if missing:
        raise ValueError(f"Missing required enrichment columns: {missing}")

    merged["pos3x_vs_neg_status"] = merged["enrich_pos3x_vs_neg"].apply(lambda v: flag(float(v), up, down))
    merged["pos3x_vs_lib_status"] = merged["enrich_pos3x_vs_lib"].apply(lambda v: flag(float(v), up, down))

    has_pos1x = ("enrich_pos1x_vs_neg" in merged.columns) and ("enrich_pos1x_vs_lib" in merged.columns)
    if has_pos1x:
        merged["pos1x_vs_neg_status"] = merged["enrich_pos1x_vs_neg"].apply(lambda v: flag(float(v), up, down))
        merged["pos1x_vs_lib_status"] = merged["enrich_pos1x_vs_lib"].apply(lambda v: flag(float(v), up, down))
        merged["pos_vs_neg_status"] = [
            combine_status(a, b)
            for a, b in zip(merged["pos3x_vs_neg_status"], merged["pos1x_vs_neg_status"])
        ]
        merged["pos_vs_lib_status"] = [
            combine_status(a, b)
            for a, b in zip(merged["pos3x_vs_lib_status"], merged["pos1x_vs_lib_status"])
        ]
    else:
        merged["pos_vs_neg_status"] = merged["pos3x_vs_neg_status"]
        merged["pos_vs_lib_status"] = merged["pos3x_vs_lib_status"]

    merged["neg_vs_lib_status"] = merged["deplete_neg_lib"].apply(lambda v: flag(float(v), up, down))
    merged["specificity"] = merged.apply(annotate_row, axis=1)


# ---------------------------------------------------------------------------
# Per-peptide table builder
# ---------------------------------------------------------------------------
def build_variant_table_for_peptide(
    peptide_key: str, files: Dict[str, Path], cfg: dict,
) -> pd.DataFrame:
    beta_col  = cfg["beta_aa_column"]
    alpha_col = cfg["alpha_aa_column"]
    labels = ["lib", "neg", "pos3x", "pos1x"]

    dfs = []
    for label in labels:
        if label not in files:
            continue
        logger.info("Counting %s for peptide %s from %s", label, peptide_key, files[label].name)
        dfs.append(count_paired_variants(files[label], beta_col, alpha_col, label=label))

    merged = merge_counts(dfs)
    merged["peptide"] = peptide_key

    pseudocount = float(cfg["pseudocount"])
    add_enrichment_columns(merged, pseudocount, "pos3x", "neg")
    add_enrichment_columns(merged, pseudocount, "pos1x", "neg")
    add_enrichment_columns(merged, pseudocount, "pos3x", "lib")
    add_enrichment_columns(merged, pseudocount, "pos1x", "lib")
    add_enrichment_columns(merged, pseudocount, "neg",   "lib")
    if "enrich_neg_vs_lib" in merged.columns:
        merged.rename(columns={"enrich_neg_vs_lib": "deplete_neg_lib"}, inplace=True)

    add_status_and_specificity(merged, cfg)
    return merged


def sanitize_name(name: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name).strip())
    return cleaned.strip("_") or "peptide"


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def run(cfg: dict) -> None:
    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = setup_pipeline_logging(
        logs_dir=cfg["logs_dir"],
        script_name=SCRIPT_NAME,
        scope="run",
        run_label=cfg.get("run_label"),
    )
    logger.info("Logging to: %s", log_path)

    library_file, groups, unknown = discover_input_files(cfg)
    logger.info("Library file: %s", library_file.name)
    logger.info("Peptide groups: %d", len(groups))
    if unknown:
        logger.warning("Ignored %d unrecognised file(s)", len(unknown))

    required = set(cfg["required_conditions"])
    summary_rows = []
    processed = skipped = 0

    for peptide_key in sorted(groups):
        files = groups[peptide_key]
        missing = sorted(required - set(files.keys()))
        if missing:
            logger.warning("Skipping %s — missing conditions: %s", peptide_key, missing)
            skipped += 1
            continue

        table = build_variant_table_for_peptide(peptide_key, files, cfg)
        spec_counts = table["specificity"].value_counts(dropna=False).to_dict()
        logger.info("Specificity for %s: %s", peptide_key, spec_counts)

        out_file = output_dir / f"{sanitize_name(peptide_key)}.variant_labeling.csv"
        table.sort_values(VARIANT_KEY_COL).to_csv(out_file, index=False)
        processed += 1

        summary_rows.append({
            "peptide":     peptide_key,
            "n_variants":  int(table.shape[0]),
            "n_spec_1":    int((table["specificity"] == 1).sum()),
            "n_spec_0":    int((table["specificity"] == 0).sum()),
            "n_spec_2":    int((table["specificity"] == 2).sum()),
            "output_csv":  str(out_file),
        })
        logger.info("Wrote: %s (%d variants)", out_file, table.shape[0])

    summary_path = output_dir / "variant_labeling_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    logger.info("Summary: %s", summary_path)
    logger.info("Done. Processed=%d, Skipped=%d", processed, skipped)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Paired-Chain Pipeline Step 04: Paired variant labeling and enrichment scoring."
    )
    parser.add_argument("--yaml_config", required=True)
    args = parser.parse_args()
    cfg = parse_yaml(args.yaml_config)
    run(cfg)


if __name__ == "__main__":
    main()
