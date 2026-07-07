#!/usr/bin/env python3
import argparse
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from utils.config_utils import load_pipeline_config
from utils.logging_utils import setup_pipeline_logging

SCRIPT_NAME = "script4_variant_labeling"
logger = logging.getLogger(__name__)
VALID_POSITIVE_MODES = {"single_positive", "combined_positive", "monotonic_strict"}


def parse_yaml(yaml_path: str) -> dict:
    cfg = load_pipeline_config(
        yaml_path,
        required_keys=("input_dir", "output_dir"),
        default_values={
            "input_glob": "*.filtered.PASS.aa.tsv",
            "aa_column": "aa_seq",
            "pseudocount": 1e-9,
            "logs_subdir": "_logs",  # backward-compatible fallback
            "lib_suffix": "originallib",
            "neg_suffix": "2xnegative",
            "positive_suffixes": {
                "pos3x": "3xpositive",
                "pos3x1x": "31xpositive",
            },
            "positive_order": ["pos3x", "pos3x1x"],
            "required_conditions": ["neg"],
            "logs_dir": None,
            "status_up": 2.0,     # Enriched if >= 2
            "status_down": 0.5,   # Depleted if <= 0.5
            "positive_mode": "combined_positive",
            "single_positive_source": None,
            # Token prefix used to parse conditions from filenames.
            # Leave empty ("") to match any text before the condition suffix.
            "file_token_prefix": "clib",
        },
        path_keys=("input_dir", "output_dir", "logs_dir"),
    )
    if not cfg.get("logs_dir"):
        cfg["logs_dir"] = str(Path(cfg["output_dir"]) / str(cfg.get("logs_subdir", "_logs")))
    cfg["positive_suffixes"] = normalize_positive_suffixes(cfg.get("positive_suffixes"))
    cfg["positive_order"] = normalize_positive_order(cfg.get("positive_order"), cfg["positive_suffixes"])
    cfg["positive_mode"] = normalize_positive_mode(cfg.get("positive_mode"))
    cfg["single_positive_source"] = normalize_single_positive_source(
        cfg.get("single_positive_source"),
        cfg["positive_suffixes"],
    )
    cfg["required_conditions"] = normalize_required_conditions(
        cfg.get("required_conditions"),
        cfg["positive_suffixes"],
    )
    return cfg


def normalize_positive_suffixes(value: Optional[dict]) -> Dict[str, str]:
    if not isinstance(value, dict) or not value:
        raise ValueError("positive_suffixes must be a non-empty mapping of label -> suffix")

    suffixes: Dict[str, str] = {}
    for raw_label, raw_suffix in value.items():
        label = str(raw_label).strip().lower()
        suffix = str(raw_suffix).strip()
        if not label:
            raise ValueError("positive condition labels cannot be empty")
        if not re.fullmatch(r"[a-z][a-z0-9_]*", label):
            raise ValueError(
                f"Invalid positive label {raw_label!r}. Use lowercase letters/numbers/underscore and start with a letter"
            )
        if not suffix:
            raise ValueError(f"Suffix for positive label {label!r} cannot be empty")
        suffixes[label] = suffix

    suffix_values = list(suffixes.values())
    if len(suffix_values) != len(set(suffix_values)):
        raise ValueError("positive_suffixes values must be unique")
    return suffixes


def normalize_positive_order(value: Optional[List[str]], positive_suffixes: Dict[str, str]) -> List[str]:
    labels = list(positive_suffixes.keys())
    if value is None:
        return labels
    if not isinstance(value, list) or not value:
        raise ValueError("positive_order must be a non-empty list of positive labels")

    order = [str(label).strip().lower() for label in value]
    if len(order) != len(set(order)):
        raise ValueError("positive_order contains duplicate labels")

    unknown = sorted(set(order) - set(labels))
    missing = sorted(set(labels) - set(order))
    if unknown:
        raise ValueError(f"positive_order contains unknown labels: {unknown}")
    if missing:
        raise ValueError(f"positive_order is missing labels from positive_suffixes: {missing}")
    return order


def normalize_required_conditions(value: Optional[List[str]], positive_suffixes: Dict[str, str]) -> List[str]:
    if value is None:
        value = ["neg"] + list(positive_suffixes.keys())
    if not isinstance(value, list) or not value:
        raise ValueError("required_conditions must be a non-empty list")

    normalized = [str(cond).strip().lower() for cond in value]
    allowed = {"lib", "neg", *positive_suffixes.keys()}
    unknown = sorted(set(normalized) - allowed)
    if unknown:
        raise ValueError(
            f"required_conditions contains unknown labels: {unknown}. Allowed: {sorted(allowed)}"
        )
    return normalized


def normalize_positive_mode(value: Optional[str]) -> str:
    mode = str(value or "combined_positive").strip().lower()
    if mode not in VALID_POSITIVE_MODES:
        raise ValueError(f"Unsupported positive_mode {value!r}. Expected one of {sorted(VALID_POSITIVE_MODES)}")
    return mode


def normalize_single_positive_source(value: Optional[str], positive_suffixes: Dict[str, str]) -> str:
    default_source = next(iter(positive_suffixes.keys()))
    source = str(value or default_source).strip().lower()
    if source not in positive_suffixes:
        raise ValueError(
            f"Unsupported single_positive_source {value!r}. Expected one of {sorted(positive_suffixes.keys())}"
        )
    return source


def parse_condition_token(
    token: str,
    lib_suffix: str,
    neg_suffix: str,
    positive_suffixes: Dict[str, str],
    file_token_prefix: str = "",
) -> Tuple[Optional[str], Optional[str]]:
    # Strip known prefix if provided (e.g., "clib")
    if file_token_prefix and token.startswith(file_token_prefix):
        token = token[len(file_token_prefix):].lstrip("-_")
    
    if token.endswith(lib_suffix):
        return token[: -len(lib_suffix)].rstrip("-_") or "lib", "lib"
    if token.endswith(neg_suffix):
        return token[: -len(neg_suffix)].rstrip("-_") or token, "neg"

    # Match longest suffixes first to avoid collisions like 31xpositive vs 1xpositive.
    ordered_suffixes = sorted(positive_suffixes.items(), key=lambda item: len(item[1]), reverse=True)
    for label, suffix in ordered_suffixes:
        if token.endswith(suffix):
            return token[: -len(suffix)].rstrip("-_") or token, label
    return None, None


def discover_input_files(cfg: dict) -> Tuple[Path, Dict[str, Dict[str, Path]], List[Path]]:
    input_dir = Path(cfg["input_dir"])
    files = sorted(input_dir.glob(cfg["input_glob"]))
    if not files:
        raise FileNotFoundError(f"No input files found in {input_dir} with glob {cfg['input_glob']!r}")

    groups: Dict[str, Dict[str, Path]] = {}
    unknown: List[Path] = []
    library_file: Optional[Path] = None

    file_suffix = ".filtered.PASS.aa.tsv"

    for path in files:
        if not path.name.endswith(file_suffix):
            unknown.append(path)
            continue

        stem = path.name[: -len(file_suffix)]
        if "_" not in stem:
            unknown.append(path)
            continue

        _, token = stem.split("_", 1)
        peptide_key, condition = parse_condition_token(
            token=token,
            lib_suffix=cfg["lib_suffix"],
            neg_suffix=cfg["neg_suffix"],
            positive_suffixes=cfg["positive_suffixes"],
            file_token_prefix=cfg["file_token_prefix"],
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


def add_tcr_population_column(merged: pd.DataFrame, positive_order: List[str]) -> None:
    def get_population(row: pd.Series) -> str:
        sources = [
            label
            for label in positive_order
            if f"count_{label}" in row and float(row[f"count_{label}"]) > 0
        ]
        return ";".join(sources) if sources else "none"

    merged["tcr_population"] = merged.apply(get_population, axis=1)


# --- NEW: same rules as variant_analysis.ipynb ---
def flag(val: float, up: float = 2.0, down: float = 0.5) -> str:
    if val >= up:
        return "Enriched"
    if val <= down:
        return "Depleted"
    return "NoChange"


def combine_status(a: str, b: str) -> str:
    """
    Combine two statuses into one:
      - if either is Enriched -> Enriched
      - else if both are Depleted -> Depleted
      - else -> NoChange
    """
    if a == "Enriched" or b == "Enriched":
        return "Enriched"
    if a == "Depleted" and b == "Depleted":
        return "Depleted"
    return "NoChange"


def combine_statuses(statuses: List[str]) -> str:
    if not statuses:
        raise ValueError("At least one status is required to combine")
    combined = statuses[0]
    for status in statuses[1:]:
        combined = combine_status(combined, status)
    return combined


def get_positive_labels_for_mode(files: Dict[str, Path], cfg: dict) -> List[str]:
    mode = cfg["positive_mode"]
    ordered_labels = cfg["positive_order"]
    available = [label for label in ordered_labels if label in files]

    if mode == "single_positive":
        source = cfg["single_positive_source"]
        if source not in files:
            raise ValueError(f"selected positive source {source!r} is missing")
        return [source]

    if mode == "combined_positive":
        if not available:
            raise ValueError("no positive condition files are available")
        return available

    if mode == "monotonic_strict":
        missing = [label for label in ordered_labels if label not in files]
        if missing:
            raise ValueError(f"missing positive condition files: {missing}")
        return ordered_labels

    raise ValueError(f"Unsupported positive_mode {mode!r}")


def monotonic_constraint_satisfied(row: pd.Series, positive_order: List[str]) -> bool:
    values = [float(row["freq_lib"])]
    for label in positive_order:
        freq_col = f"freq_{label}"
        if freq_col not in row:
            return False
        values.append(float(row[freq_col]))
    return all(left < right for left, right in zip(values, values[1:]))


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

def add_status_and_specificity(merged: pd.DataFrame, cfg: dict, positive_labels: List[str]) -> None:
    up = float(cfg.get("status_up", 2.0))
    down = float(cfg.get("status_down", 0.5))
    mode = cfg["positive_mode"]

    base_needed = ["deplete_neg_lib"]
    for label in positive_labels:
        base_needed.extend([f"enrich_{label}_vs_neg", f"enrich_{label}_vs_lib"])
    missing = [c for c in base_needed if c not in merged.columns]
    if missing:
        raise ValueError(f"Missing required enrichment columns for labeling: {missing}")

    for label in positive_labels:
        merged[f"{label}_vs_neg_status"] = merged[f"enrich_{label}_vs_neg"].apply(lambda v: flag(float(v), up, down))
        merged[f"{label}_vs_lib_status"] = merged[f"enrich_{label}_vs_lib"].apply(lambda v: flag(float(v), up, down))

    if mode == "single_positive":
        source = cfg["single_positive_source"]
        merged["pos_vs_neg_status"] = merged[f"{source}_vs_neg_status"]
        merged["pos_vs_lib_status"] = merged[f"{source}_vs_lib_status"]
    elif mode == "combined_positive":
        merged["pos_vs_neg_status"] = [
            combine_statuses([row[f"{label}_vs_neg_status"] for label in positive_labels])
            for _, row in merged.iterrows()
        ]
        merged["pos_vs_lib_status"] = [
            combine_statuses([row[f"{label}_vs_lib_status"] for label in positive_labels])
            for _, row in merged.iterrows()
        ]
    elif mode == "monotonic_strict":
        pos_vs_neg_statuses: List[str] = []
        pos_vs_lib_statuses: List[str] = []
        for _, row in merged.iterrows():
            monotonic_ok = monotonic_constraint_satisfied(row, positive_labels)
            neg_status = combine_statuses([row[f"{label}_vs_neg_status"] for label in positive_labels])
            lib_status = combine_statuses([row[f"{label}_vs_lib_status"] for label in positive_labels])
            if not monotonic_ok:
                if neg_status == "Enriched":
                    neg_status = "NoChange"
                if lib_status == "Enriched":
                    lib_status = "NoChange"
            pos_vs_neg_statuses.append(neg_status)
            pos_vs_lib_statuses.append(lib_status)
        merged["pos_vs_neg_status"] = pos_vs_neg_statuses
        merged["pos_vs_lib_status"] = pos_vs_lib_statuses
    else:
        raise ValueError(f"Unsupported positive_mode {mode!r}")

    merged["neg_vs_lib_status"] = merged["deplete_neg_lib"].apply(lambda v: flag(float(v), up, down))
    merged["specificity"] = merged.apply(annotate_row, axis=1)


def build_variant_table_for_peptide(peptide_key: str, files: Dict[str, Path], cfg: dict, positive_labels: List[str]) -> pd.DataFrame:
    aa_column = cfg["aa_column"]
    labels = ["lib", "neg", *cfg["positive_order"]]
    dfs = []
    for label in labels:
        if label not in files:
            continue
        logger.info("Counting %s for peptide %s from %s", label, peptide_key, files[label].name)
        dfs.append(count_variants(files[label], aa_column=aa_column, label=label))

    merged = merge_counts(dfs, aa_column=aa_column)
    merged["peptide"] = peptide_key
    add_tcr_population_column(merged, cfg["positive_order"])

    pseudocount = float(cfg["pseudocount"])
    for label in cfg["positive_order"]:
        add_enrichment_columns(merged, pseudocount, numerator=label, denominator="neg")
        add_enrichment_columns(merged, pseudocount, numerator=label, denominator="lib")
    add_enrichment_columns(merged, pseudocount, numerator="neg", denominator="lib")
    if "enrich_neg_vs_lib" in merged.columns:
        merged.rename(columns={"enrich_neg_vs_lib": "deplete_neg_lib"}, inplace=True)

    # add statuses + specificity (0/1/2)
    add_status_and_specificity(merged, cfg, positive_labels)

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

        try:
            positive_labels = get_positive_labels_for_mode(files, cfg)
        except ValueError as exc:
            logger.warning("Skipping peptide %s because %s", peptide_key, exc)
            skipped += 1
            continue

        table = build_variant_table_for_peptide(peptide_key, files, cfg, positive_labels)

        # Log specificity distribution
        spec_counts = table["specificity"].value_counts(dropna=False).to_dict()
        logger.info("Specificity distribution for %s: %s", peptide_key, spec_counts)

        out_file = output_dir / f"{sanitize_name(peptide_key)}.variant_labeling.csv"
        table.sort_values(cfg["aa_column"]).to_csv(out_file, index=False)
        processed += 1

        summary_row = {
            "peptide": peptide_key,
            "n_variants": int(table.shape[0]),
            "n_spec_1": int((table["specificity"] == 1).sum()),
            "n_spec_0": int((table["specificity"] == 0).sum()),
            "n_spec_2": int((table["specificity"] == 2).sum()),
            "output_csv": str(out_file),
            "input_lib": str(files["lib"]),
            "input_neg": str(files["neg"]),
        }
        for label in cfg["positive_order"]:
            summary_row[f"input_{label}"] = str(files.get(label, ""))
        summary_rows.append(summary_row)
        logger.info("Wrote peptide table: %s (%d variants)", out_file, table.shape[0])

    summary_path = output_dir / "variant_labeling_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    logger.info("Wrote summary table: %s", summary_path)
    logger.info("Done. Processed=%d, Skipped=%d", processed, skipped)


def main() -> None:

    # Parse yaml file and get all inputs
    parser = argparse.ArgumentParser(
        description="Generate per-peptide variant labeling CSVs (with specificity 0/1/2) from extracted PASS AA tables."
    )
    parser.add_argument("--yaml_config", type=str, default=None, help="Path to YAML config file.")
    args = parser.parse_args()

    yaml_path = args.yaml_config
    if not yaml_path:
        raise ValueError("A config file must be supplied via --yaml_config.")
    
    cfg = parse_yaml(yaml_path)
    run(cfg)


if __name__ == "__main__":
    main()
