#!/usr/bin/env python3
"""Create source tables for thesis dataset-analysis panels."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "results" / "thesis_dataset_analysis"


@dataclass(frozen=True)
class DatasetConfig:
    dataset: str
    path: Path
    wt_tcr: str


DATASETS = (
    DatasetConfig(
        dataset="A3",
        path=Path("/cluster/project/reddy/katja/ml_refactor/data/data_from_NGS_pipline/A3_Q20_3x1x_tcr_peptide_label_all.csv"),
        wt_tcr="ASSPNMADEQY",
    ),
    DatasetConfig(
        dataset="DMF5",
        path=Path("/cluster/project/reddy/katja/ml_refactor/data/data_from_NGS_pipline/DMF5_Q30_2x_tcr_peptide_label_all.csv"),
        wt_tcr="ASSLSFGTEAF",
    ),
)

REQUIRED_COLUMNS = ("tcr", "peptide", "label")
ANALYSIS_TCR_SUBSAMPLE_N = None
ANALYSIS_SUBSAMPLE_SEED = 42
PAIRWISE_SAMPLE_N = 250
N_RANDOM_PAIRS = 20_000
RANDOM_SEED = 42


def hamming_distance_equal_len(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(x != y for x, y in zip(a, b))


def load_dataset(config: DatasetConfig) -> pd.DataFrame:
    if not config.path.exists():
        raise FileNotFoundError(config.path)

    df = pd.read_csv(config.path)
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"[{config.dataset}] Missing required columns: {missing}")

    df = df.copy()
    df["tcr"] = df["tcr"].astype(str).str.strip().str.upper()
    df["peptide"] = df["peptide"].astype(str).str.strip()
    df["label"] = pd.to_numeric(df["label"], errors="coerce").fillna(0).astype(int)
    return df


def analysis_dataset(df: pd.DataFrame) -> pd.DataFrame:
    if ANALYSIS_TCR_SUBSAMPLE_N is None:
        return df.copy()
    sampled_tcrs = (
        df["tcr"]
        .drop_duplicates()
        .sample(
            n=min(ANALYSIS_TCR_SUBSAMPLE_N, df["tcr"].nunique()),
            random_state=ANALYSIS_SUBSAMPLE_SEED,
        )
    )
    return df[df["tcr"].isin(sampled_tcrs)].copy()


def compute_wt_distance_counts(df: pd.DataFrame, wt_tcr: str) -> tuple[pd.DataFrame, dict[str, int | float]]:
    wt_tcr = wt_tcr.strip().upper()
    wt_len = len(wt_tcr)
    tcr_series = df["tcr"].astype(str).str.strip().str.upper()
    valid = tcr_series[tcr_series.str.len() == wt_len]
    row_dists = valid.apply(lambda seq: hamming_distance_equal_len(seq, wt_tcr))
    counts = (
        row_dists.value_counts()
        .sort_index()
        .rename_axis("wt_hamming_distance")
        .reset_index(name="n_tcr_rows")
    )
    stats = {
        "wt_length": wt_len,
        "rows_valid_wt_length": int(len(valid)),
        "rows_invalid_wt_length": int(len(tcr_series) - len(valid)),
        "mean_distance_row_weighted": float(row_dists.mean()) if len(row_dists) else np.nan,
    }
    return counts, stats


def compute_crossreactive_binder_depth(df: pd.DataFrame) -> pd.DataFrame:
    binders = df[df["label"] == 1].copy()
    if binders.empty:
        return pd.DataFrame(columns=["n_bound_peptides", "n_tcrs"])

    depth = binders.groupby("tcr")["peptide"].nunique()
    depth = depth[depth >= 1]
    if depth.empty:
        return pd.DataFrame(columns=["n_bound_peptides", "n_tcrs"])

    return (
        depth.value_counts()
        .sort_index()
        .rename_axis("n_bound_peptides")
        .reset_index(name="n_tcrs")
    )


def compute_single_binder_label_context_by_peptide(
    df: pd.DataFrame,
    all_peptides: list[str],
) -> pd.DataFrame:
    peptide_labels = df.groupby(["tcr", "peptide"], as_index=False)["label"].max()
    positive_peptides = (
        peptide_labels[peptide_labels["label"] == 1]
        .groupby("tcr")["peptide"]
        .agg(list)
        .rename("positive_peptides")
    )
    negative_counts = (
        peptide_labels[peptide_labels["label"] == 0]
        .groupby("tcr")["peptide"]
        .nunique()
        .rename("n_negative_peptides")
    )
    per_tcr = positive_peptides.to_frame().join(negative_counts, how="left")
    per_tcr["n_negative_peptides"] = per_tcr["n_negative_peptides"].fillna(0).astype(int)
    per_tcr["n_bound_peptides"] = per_tcr["positive_peptides"].str.len()
    per_tcr = per_tcr[per_tcr["n_bound_peptides"] == 1].copy()

    context_order = [
        "Positive on one peptide, absent elsewhere",
        "Positive on one peptide, only 0 on other peptides",
    ]
    full_index = pd.MultiIndex.from_product(
        [all_peptides, context_order],
        names=["peptide", "single_binder_context"],
    )
    if per_tcr.empty:
        return full_index.to_frame(index=False).assign(n_tcrs=0)

    per_tcr["peptide"] = per_tcr["positive_peptides"].str[0]
    per_tcr["single_binder_context"] = np.where(
        per_tcr["n_negative_peptides"] > 0,
        "Positive on one peptide, only 0 on other peptides",
        "Positive on one peptide, absent elsewhere",
    )

    return (
        per_tcr.groupby(["peptide", "single_binder_context"])
        .size()
        .reindex(full_index, fill_value=0)
        .rename("n_tcrs")
        .reset_index()
    )


def sample_pairwise_hamming_distances(
    df: pd.DataFrame,
    sample_n: int,
    n_pairs: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, int | float]]:
    unique_tcrs = pd.Index(df["tcr"].dropna().astype(str).str.strip().str.upper().unique())
    if len(unique_tcrs) < 2:
        return pd.DataFrame({"pairwise_hamming_distance": []}), {
            "unique_tcr_total": int(len(unique_tcrs)),
            "sampled_unique_tcrs": int(len(unique_tcrs)),
            "pairs_computed": 0,
        }

    lengths = pd.Series(unique_tcrs).str.len()
    mode_len = int(lengths.mode().iat[0])
    candidate_tcrs = unique_tcrs[(lengths == mode_len).to_numpy()]

    rng = np.random.default_rng(seed)
    sampled = rng.choice(candidate_tcrs.to_numpy(), size=min(sample_n, len(candidate_tcrs)), replace=False)
    arr = np.array([list(seq) for seq in sampled])
    n = len(arr)
    if n < 2:
        return pd.DataFrame({"pairwise_hamming_distance": []}), {
            "unique_tcr_total": int(len(unique_tcrs)),
            "candidate_equal_len_tcrs": int(len(candidate_tcrs)),
            "sampled_unique_tcrs": int(n),
            "pairs_computed": 0,
            "sequence_length_used": mode_len,
        }

    i = rng.integers(0, n, size=n_pairs)
    j = rng.integers(0, n, size=n_pairs)
    mask = i != j
    dists = np.sum(arr[i[mask]] != arr[j[mask]], axis=1)

    stats = {
        "unique_tcr_total": int(len(unique_tcrs)),
        "candidate_equal_len_tcrs": int(len(candidate_tcrs)),
        "sampled_unique_tcrs": int(n),
        "pairs_requested": int(n_pairs),
        "pairs_computed": int(len(dists)),
        "sequence_length_used": int(mode_len),
        "mean_pairwise_distance": float(np.mean(dists)) if len(dists) else np.nan,
        "median_pairwise_distance": float(np.median(dists)) if len(dists) else np.nan,
        "p05_pairwise_distance": float(np.percentile(dists, 5)) if len(dists) else np.nan,
        "p95_pairwise_distance": float(np.percentile(dists, 95)) if len(dists) else np.nan,
    }
    return pd.DataFrame({"pairwise_hamming_distance": dists}), stats


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    wt_frames = []
    crossreact_frames = []
    pairwise_frames = []
    single_binder_frames = []
    validation_rows = []

    for config in DATASETS:
        full_df = load_dataset(config)
        analysis_df = analysis_dataset(full_df)
        peptides = sorted(full_df["peptide"].dropna().unique())

        wt_counts, wt_stats = compute_wt_distance_counts(analysis_df, config.wt_tcr)
        wt_counts["dataset"] = config.dataset
        wt_frames.append(wt_counts[["dataset", "wt_hamming_distance", "n_tcr_rows"]])

        crossreact = compute_crossreactive_binder_depth(analysis_df)
        crossreact["dataset"] = config.dataset
        crossreact_frames.append(crossreact[["dataset", "n_bound_peptides", "n_tcrs"]])

        pairwise, pairwise_stats = sample_pairwise_hamming_distances(
            analysis_df,
            sample_n=PAIRWISE_SAMPLE_N,
            n_pairs=N_RANDOM_PAIRS,
            seed=RANDOM_SEED,
        )
        pairwise["dataset"] = config.dataset
        pairwise_frames.append(pairwise[["dataset", "pairwise_hamming_distance"]])

        single_binder = compute_single_binder_label_context_by_peptide(analysis_df, peptides)
        single_binder["dataset"] = config.dataset
        single_binder_frames.append(single_binder[["dataset", "peptide", "single_binder_context", "n_tcrs"]])

        validation_rows.append(
            {
                "dataset": config.dataset,
                "source_file": str(config.path),
                "rows_total": int(len(full_df)),
                "rows_analyzed": int(len(analysis_df)),
                "unique_tcrs_total": int(full_df["tcr"].nunique()),
                "unique_tcrs_analyzed": int(analysis_df["tcr"].nunique()),
                "unique_peptides_total": int(full_df["peptide"].nunique()),
                "binder_rows_analyzed": int((analysis_df["label"] == 1).sum()),
                "crossreactive_binder_tcrs": int(crossreact.loc[crossreact["n_bound_peptides"] >= 2, "n_tcrs"].sum()) if not crossreact.empty else 0,
                **wt_stats,
                **pairwise_stats,
            }
        )

    pd.concat(wt_frames, ignore_index=True).to_csv(OUTPUT_DIR / "dataset_wt_distance_counts.csv", index=False)
    pd.concat(crossreact_frames, ignore_index=True).to_csv(OUTPUT_DIR / "dataset_crossreactive_binder_depth.csv", index=False)
    pd.concat(pairwise_frames, ignore_index=True).to_csv(OUTPUT_DIR / "dataset_pairwise_hamming_distances.csv", index=False)
    pd.concat(single_binder_frames, ignore_index=True).to_csv(OUTPUT_DIR / "dataset_single_binder_context_by_peptide.csv", index=False)
    pd.DataFrame(validation_rows).to_csv(OUTPUT_DIR / "dataset_analysis_validation.csv", index=False)
    print(f"Wrote dataset-analysis source tables to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
