#!/usr/bin/env python3
"""Create thesis-ready NGS pipeline overview figures."""

from __future__ import annotations

import argparse
import os
import time
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig-katja")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/xdg-cache-katja")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.gridspec import GridSpec


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "results" / "thesis_pipeline_overview"


@dataclass(frozen=True)
class MetricsSource:
    dataset: str
    project: str
    paths: tuple[str, ...]


METRIC_SOURCES = (
    MetricsSource(
        dataset="DMF5",
        project="P3481_LUCA-TCRDMF5",
        paths=(
            "data/P3481_LUCA-TCRDMF5/01_preprocessed/22_05_2026_minLenght110_QScore30/per_sample_metrics.tsv",
            "data/P3481_LUCA-TCRDMF5/02_umi/per_sample_metrics.tsv",
            "data/P3481_LUCA-TCRDMF5/03_extracted/per_sample_metrics.tsv",
        ),
    ),
    MetricsSource(
        dataset="TCRA3",
        project="P3408_LUCA-TCRA3",
        paths=(
            "data/P3408_LUCA-TCRA3/01_preprocessed/10_05_2026_minLenght191_QScore30/per_sample_metrics.tsv",
            "data/P3408_LUCA-TCRA3/02_umi/14_05_2026_FULLRUN/per_sample_metrics.tsv",
            "data/P3408_LUCA-TCRA3/03_extracted/19_05_2026_minlenght191_Q30/per_sample_metrics.tsv",
        ),
    ),
    MetricsSource(
        dataset="TCRA3",
        project="P3569_LUCA-TCRA3_KH157",
        paths=(
            "data/P3569_LUCA-TCRA3_KH157/01_preprocessed/26_04_2026_minLenght191_QScore20/per_sample_metrics.tsv",
            "data/P3569_LUCA-TCRA3_KH157/02_umi/27_04_2026_KH157/per_sample_metrics.tsv",
            "data/P3569_LUCA-TCRA3_KH157/03_extracted/27_04_2026_KH157/per_sample_metrics.tsv",
        ),
    ),
)

VARIANT_SUMMARIES = {
    "DMF5": "data/P3481_LUCA-TCRDMF5/04_variant_labeling/Q30/combined/variant_labeling_summary.csv",
    "TCRA3": "data/P3408_LUCA-TCRA3/04_variant_labeling/19_05_2026_minlenght191_Q30/combined/variant_labeling_summary.csv",
}

STEP_LABELS = {
    "01_preprocessing": "01 preprocessing",
    "02_umi_consensus": "02 UMI consensus",
    "02a_singleton_rescue": "02a singleton rescue",
    "03_extraction": "03 extraction",
}

LINEAR_STAGES = (
    ("Raw reads", "reads"),
    ("After trimming", "reads"),
    ("UMI input", "reads"),
    ("Consensus UMI", "umi_groups"),
    ("Consensus UMI + singleton rescue", "consensus_umi_plus_rescued_singletons"),
    ("Extraction input", "consensus_umi_plus_rescued_singletons"),
    ("Extraction PASS", "reads"),
)

LOSS_LABELS = {
    "reads_lost_trimming": "Trimming",
    "reads_lost_no_anchor": "No anchor",
    "reads_lost_anchor_wrong_pos": "Anchor wrong position",
    "reads_lost_empty_insert": "Empty insert",
    "reads_lost_too_short": "Too short",
    "reads_lost_contains_N": "Contains N",
    "reads_lost_wrong_length": "Wrong length",
    "reads_lost_stop_codon": "Stop codon",
    "reads_lost_validation": "Library validation",
}

DATASET_COLORS = {"DMF5": "#3B6EA8", "TCRA3": "#C4563A"}
LABEL_COLORS = {
    "Specific (label 1)": "#2A9D8F",
    "Non-specific (label 0)": "#6C757D",
    "Ambiguous/other (label 2)": "#E9A03B",
}

SERIF_FONT_STACK = [
    "Times New Roman",
    "Times",
    "Nimbus Roman",
    "Liberation Serif",
    "DejaVu Serif",
]

DEDUP_COLUMNS = [
    "dataset",
    "project",
    "source_file",
    "sample_id",
    "sample_name",
    "step",
    "metric",
    "value",
]


def format_count(value: float) -> str:
    value = float(value)
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if abs(value) >= 1_000:
        return f"{value / 1_000:.1f}k"
    return f"{value:,.0f}"


def configure_plot_theme(font_scale: float = 1.1) -> None:
    sns.set_theme(style="whitegrid", context="paper", font="serif", font_scale=font_scale)
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": SERIF_FONT_STACK,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def load_metric_sources() -> tuple[pd.DataFrame, pd.DataFrame]:
    raw_frames = []
    integrity_rows = []
    for source in METRIC_SOURCES:
        for rel_path in source.paths:
            path = ROOT / rel_path
            if not path.exists():
                raise FileNotFoundError(path)
            frame = pd.read_csv(path, sep="\t")
            frame["dataset"] = source.dataset
            frame["project"] = source.project
            frame["source_file"] = rel_path
            frame["value"] = pd.to_numeric(frame["value"], errors="raise")
            raw_frames.append(frame)

            dup_mask = frame.duplicated(["sample_id", "sample_name", "step", "metric", "value"], keep=False)
            dedup_frame = frame.drop_duplicates(["sample_id", "sample_name", "step", "metric", "value"])
            integrity_rows.append(
                {
                    "dataset": source.dataset,
                    "project": source.project,
                    "source_file": rel_path,
                    "rows_raw": len(frame),
                    "rows_after_dedup": len(dedup_frame),
                    "duplicate_rows": int(dup_mask.sum()),
                    "duplicate_row_pairs_removed": len(frame) - len(dedup_frame),
                    "has_exact_duplicates": bool(dup_mask.any()),
                }
            )

    raw_metrics = pd.concat(raw_frames, ignore_index=True)
    metrics = raw_metrics.drop_duplicates(DEDUP_COLUMNS).copy()
    metrics["step_label"] = metrics["step"].map(STEP_LABELS).fillna(metrics["step"])
    integrity = pd.DataFrame(integrity_rows)
    return metrics, integrity


def metric_total(metrics: pd.DataFrame, dataset: str, step: str, metric: str) -> int:
    mask = (
        (metrics["dataset"] == dataset)
        & (metrics["step"] == step)
        & (metrics["metric"] == metric)
    )
    return int(metrics.loc[mask, "value"].sum())


def metric_total_for_projects(
    metrics: pd.DataFrame,
    dataset: str,
    projects: set[str],
    step: str,
    metric: str,
) -> int:
    if not projects:
        return 0
    mask = (
        (metrics["dataset"] == dataset)
        & (metrics["project"].isin(projects))
        & (metrics["step"] == step)
        & (metrics["metric"] == metric)
    )
    return int(metrics.loc[mask, "value"].sum())


def build_stage_counts(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset in sorted(metrics["dataset"].unique()):
        raw = metric_total(metrics, dataset, "01_preprocessing", "reads_in_total")
        consensus_umi = metric_total(metrics, dataset, "02_umi_consensus", "umis_consensus")
        bridge = build_consensus_bridge(metrics)
        bridge_row = bridge[bridge["dataset"] == dataset].iloc[0]
        consensus_plus_rescue = int(bridge_row["consensus_plus_rescue"])
        stages = {
            "Raw reads": metric_total(metrics, dataset, "01_preprocessing", "reads_in_total"),
            "After trimming": metric_total(metrics, dataset, "01_preprocessing", "reads_out_total"),
            "UMI input": metric_total(metrics, dataset, "02_umi_consensus", "reads_in"),
            "Consensus UMI": consensus_umi,
            "Consensus UMI + singleton rescue": consensus_plus_rescue,
            "Extraction input": metric_total(metrics, dataset, "03_extraction", "reads_in_total"),
            "Extraction PASS": metric_total(metrics, dataset, "03_extraction", "reads_pass"),
        }
        for order, (stage, unit) in enumerate(LINEAR_STAGES):
            count = stages[stage]
            rows.append(
                {
                    "dataset": dataset,
                    "stage": stage,
                    "stage_order": order,
                    "count": count,
                    "count_unit": unit,
                    "retention_from_raw": count / raw if raw else 0.0,
                }
            )
    return pd.DataFrame(rows)


def build_consensus_bridge(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset in sorted(metrics["dataset"].unique()):
        dataset_projects = set(
            metrics.loc[metrics["dataset"] == dataset, "project"].unique()
        )
        rescue_projects = set(
            metrics.loc[
                (metrics["dataset"] == dataset) & (metrics["step"] == "02a_singleton_rescue"),
                "project",
            ].unique()
        )
        non_rescue_projects = dataset_projects - rescue_projects
        consensus_umi = metric_total(metrics, dataset, "02_umi_consensus", "umis_consensus")
        rescued_singletons = metric_total(metrics, dataset, "02a_singleton_rescue", "reads_rescued")
        direct_singletons_without_rescue = metric_total_for_projects(
            metrics,
            dataset,
            non_rescue_projects,
            "02_umi_consensus",
            "reads_in_singletons",
        )
        singleton_addition = rescued_singletons + direct_singletons_without_rescue
        extraction_input = metric_total(metrics, dataset, "03_extraction", "reads_in_total")
        rows.append(
            {
                "dataset": dataset,
                "consensus_umis": consensus_umi,
                "rescued_singletons": rescued_singletons,
                "direct_singletons_without_rescue": direct_singletons_without_rescue,
                "singleton_addition_to_extraction": singleton_addition,
                "consensus_plus_rescue": consensus_umi + singleton_addition,
                "extraction_input": extraction_input,
                "delta": consensus_umi + singleton_addition - extraction_input,
            }
        )
    return pd.DataFrame(rows)


def build_rescue_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    rescue_metrics = {
        "reads_in",
        "reads_lost_no_anchor",
        "reads_lost_too_short",
        "reads_extracted",
        "reads_rescued",
        "variants_total",
        "variants_kept",
    }
    rescue = metrics[
        (metrics["step"] == "02a_singleton_rescue")
        & (metrics["metric"].isin(rescue_metrics))
    ].copy()
    if rescue.empty:
        return pd.DataFrame(
            columns=[
                "dataset",
                "samples_included",
                "reads_in",
                "reads_lost_no_anchor",
                "reads_lost_too_short",
                "reads_extracted",
                "reads_rescued",
                "variants_total",
                "variants_kept",
                "rescue_retention_of_input",
                "rescue_retention_of_extracted",
            ]
        )

    pivot = (
        rescue.pivot_table(
            index=["dataset", "sample_id", "sample_name"],
            columns="metric",
            values="value",
            aggfunc="sum",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )
    summary = (
        pivot.groupby("dataset", as_index=False)
        .agg(
            samples_included=("sample_id", "count"),
            reads_in=("reads_in", "sum"),
            reads_lost_no_anchor=("reads_lost_no_anchor", "sum"),
            reads_lost_too_short=("reads_lost_too_short", "sum"),
            reads_extracted=("reads_extracted", "sum"),
            reads_rescued=("reads_rescued", "sum"),
            variants_total=("variants_total", "sum"),
            variants_kept=("variants_kept", "sum"),
        )
    )
    summary["rescue_retention_of_input"] = summary["reads_rescued"] / summary["reads_in"]
    summary["rescue_retention_of_extracted"] = summary["reads_rescued"] / summary["reads_extracted"]
    return summary


def build_stage_consistency_checks(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset in sorted(metrics["dataset"].unique()):
        trimming_out = metric_total(metrics, dataset, "01_preprocessing", "reads_out_total")
        umi_in = metric_total(metrics, dataset, "02_umi_consensus", "reads_in")
        umi_consensus_reads = metric_total(metrics, dataset, "02_umi_consensus", "reads_in_consensus")
        umi_singletons = metric_total(metrics, dataset, "02_umi_consensus", "reads_in_singletons")
        umi_losses = sum(
            metric_total(metrics, dataset, "02_umi_consensus", metric)
            for metric in (
                "reads_lost_too_short",
                "reads_lost_no_anchor",
                "reads_lost_anchor_wrong_pos",
                "reads_lost_empty_insert",
            )
        )
        dataset_projects = set(
            metrics.loc[metrics["dataset"] == dataset, "project"].unique()
        )
        rescue_projects = set(
            metrics.loc[
                (metrics["dataset"] == dataset) & (metrics["step"] == "02a_singleton_rescue"),
                "project",
            ].unique()
        )
        non_rescue_projects = dataset_projects - rescue_projects
        rescue_in = metric_total(metrics, dataset, "02a_singleton_rescue", "reads_in")
        rescue_extracted = metric_total(metrics, dataset, "02a_singleton_rescue", "reads_extracted")
        rescue_rescued = metric_total(metrics, dataset, "02a_singleton_rescue", "reads_rescued")
        consensus_umis = metric_total(metrics, dataset, "02_umi_consensus", "umis_consensus")
        extraction_input = metric_total(metrics, dataset, "03_extraction", "reads_in_total")
        singletons_rescue_scope = metric_total_for_projects(metrics, dataset, rescue_projects, "02_umi_consensus", "reads_in_singletons")
        direct_singletons_without_rescue = metric_total_for_projects(
            metrics,
            dataset,
            non_rescue_projects,
            "02_umi_consensus",
            "reads_in_singletons",
        )
        singleton_addition = rescue_rescued + direct_singletons_without_rescue

        checks = [
            ("trimming_out_equals_umi_input", trimming_out, umi_in, trimming_out == umi_in),
            ("umi_accounting", umi_in, umi_consensus_reads + umi_singletons + umi_losses, umi_in == umi_consensus_reads + umi_singletons + umi_losses),
            ("singletons_match_rescue_input_on_rescue_projects", singletons_rescue_scope, rescue_in, singletons_rescue_scope == rescue_in),
            ("consensus_plus_singleton_addition_equals_extraction_input", consensus_umis + singleton_addition, extraction_input, consensus_umis + singleton_addition == extraction_input),
            ("rescued_le_extracted", rescue_rescued, rescue_extracted, rescue_rescued <= rescue_extracted),
            ("extracted_le_rescue_input", rescue_extracted, rescue_in, rescue_extracted <= rescue_in),
        ]
        for check_name, left_value, right_value, passes in checks:
            rows.append(
                {
                    "dataset": dataset,
                    "check": check_name,
                    "left_value": int(left_value),
                    "right_value": int(right_value),
                    "delta": int(left_value - right_value),
                    "passes": bool(passes),
                }
            )
    return pd.DataFrame(rows)


def build_dropout_table(metrics: pd.DataFrame) -> pd.DataFrame:
    loss_metrics = set(LOSS_LABELS)
    dropout = metrics[metrics["metric"].isin(loss_metrics)].copy()
    dropout["loss_reason"] = dropout["metric"].map(LOSS_LABELS)
    grouped = (
        dropout.groupby(["dataset", "step", "step_label", "loss_reason"], as_index=False)["value"]
        .sum()
        .rename(columns={"value": "reads_lost"})
    )
    total_inputs = (
        metrics[metrics["metric"].isin(("reads_in_total", "reads_in"))]
        .groupby(["dataset", "step"], as_index=False)["value"]
        .sum()
        .rename(columns={"value": "step_reads_in"})
    )
    grouped = grouped.merge(total_inputs, on=["dataset", "step"], how="left")
    grouped["fraction_of_step_input"] = grouped["reads_lost"] / grouped["step_reads_in"]
    return grouped


def read_variant_summaries() -> pd.DataFrame:
    frames = []
    for dataset, rel_path in VARIANT_SUMMARIES.items():
        path = ROOT / rel_path
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path)
        frame["dataset"] = dataset
        frame["source_file"] = rel_path
        frames.append(frame)
    variants = pd.concat(frames, ignore_index=True)
    long = variants.melt(
        id_vars=["dataset", "peptide", "n_variants", "source_file"],
        value_vars=["n_spec_1", "n_spec_0", "n_spec_2"],
        var_name="label_metric",
        value_name="variant_count",
    )
    long["label_class"] = long["label_metric"].map(
        {
            "n_spec_1": "Specific (label 1)",
            "n_spec_0": "Non-specific (label 0)",
            "n_spec_2": "Ambiguous/other (label 2)",
        }
    )
    long["variant_count"] = pd.to_numeric(long["variant_count"], errors="raise")
    return long


def pivot_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    return (
        metrics.pivot_table(
            index=["dataset", "project", "sample_id", "sample_name", "step"],
            columns="metric",
            values="value",
            aggfunc="sum",
        )
        .reset_index()
        .rename_axis(None, axis=1)
    )


def validate_accounting(metrics: pd.DataFrame) -> pd.DataFrame:
    wide = pivot_metrics(metrics)
    rows = []
    specs = {
        "01_preprocessing": {
            "input": "reads_in_total",
            "outputs": ["reads_out_total", "reads_lost_trimming"],
        },
        "02_umi_consensus": {
            "input": "reads_in",
            "outputs": [
                "reads_in_consensus",
                "reads_in_singletons",
                "reads_lost_too_short",
                "reads_lost_no_anchor",
                "reads_lost_anchor_wrong_pos",
                "reads_lost_empty_insert",
            ],
        },
        "02a_singleton_rescue": {
            "input": "reads_in",
            "outputs": [
                "reads_lost_no_anchor",
                "reads_lost_too_short",
                "reads_extracted",
            ],
        },
        "03_extraction": {
            "input": "reads_in_total",
            "outputs": [
                "reads_pass",
                "reads_lost_no_anchor",
                "reads_lost_too_short",
                "reads_lost_contains_N",
                "reads_lost_wrong_length",
                "reads_lost_stop_codon",
                "reads_lost_validation",
            ],
        },
    }
    for _, row in wide.iterrows():
        spec = specs.get(row["step"])
        if spec is None or spec["input"] not in row:
            continue
        observed_input = row.get(spec["input"], 0)
        output_sum = sum(row.get(col, 0) for col in spec["outputs"])
        rows.append(
            {
                "dataset": row["dataset"],
                "project": row["project"],
                "sample_id": row["sample_id"],
                "sample_name": row["sample_name"],
                "step": row["step"],
                "input_metric": spec["input"],
                "input_count": observed_input,
                "accounted_count": output_sum,
                "delta": observed_input - output_sum,
                "passes_exact_accounting": observed_input == output_sum,
                "reads_rescued": row.get("reads_rescued"),
                "rescued_le_extracted": (
                    row.get("reads_rescued", 0) <= row.get("reads_extracted", 0)
                    if row["step"] == "02a_singleton_rescue"
                    else pd.NA
                ),
            }
        )
    return pd.DataFrame(rows)


def write_outputs(
    metrics: pd.DataFrame,
    stage_counts: pd.DataFrame,
    consensus_bridge: pd.DataFrame,
    duplicate_report: pd.DataFrame,
    rescue_summary: pd.DataFrame,
    dropout: pd.DataFrame,
    variants: pd.DataFrame,
    validation: pd.DataFrame,
    stage_checks: pd.DataFrame,
) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(OUTPUT_DIR / "pipeline_read_metrics_long.csv", index=False)
    stage_counts.to_csv(OUTPUT_DIR / "pipeline_stage_counts.csv", index=False)
    consensus_bridge.to_csv(OUTPUT_DIR / "pipeline_consensus_rescue_bridge.csv", index=False)
    duplicate_report.to_csv(OUTPUT_DIR / "pipeline_source_integrity_checks.csv", index=False)
    rescue_summary.to_csv(OUTPUT_DIR / "pipeline_rescue_summary.csv", index=False)
    dropout.to_csv(OUTPUT_DIR / "pipeline_dropout_reasons.csv", index=False)
    variants.to_csv(OUTPUT_DIR / "variant_label_counts_long.csv", index=False)
    validation.to_csv(OUTPUT_DIR / "pipeline_accounting_validation.csv", index=False)
    stage_checks.to_csv(OUTPUT_DIR / "pipeline_stage_consistency_checks.csv", index=False)


def prepare_thesis_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load all inputs and build the tables used by the thesis figure."""
    metrics, duplicate_report = load_metric_sources()
    stage_counts = build_stage_counts(metrics)
    consensus_bridge = build_consensus_bridge(metrics)
    rescue_summary = build_rescue_summary(metrics)
    dropout = build_dropout_table(metrics)
    variants = read_variant_summaries()
    validation = validate_accounting(metrics)
    stage_checks = build_stage_consistency_checks(metrics)
    return metrics, stage_counts, consensus_bridge, duplicate_report, rescue_summary, dropout, variants, validation, stage_checks


def draw_linear_retention(ax: plt.Axes, stage_counts: pd.DataFrame, consensus_bridge: pd.DataFrame) -> None:
    stage_order = [
        "Raw reads",
        "After trimming",
        "Consensus UMI",
        "Consensus UMI + singleton rescue",
        "Extraction PASS",
    ]
    x_positions = list(range(len(stage_order)))

    for dataset, frame in stage_counts.groupby("dataset"):
        frame = frame.set_index("stage").loc[stage_order].reset_index()
        y = frame["count"].tolist()
        ax.plot(
            x_positions[:2],
            y[:2],
            marker="o",
            linewidth=2.5,
            markersize=6.5,
            color=DATASET_COLORS[dataset],
            label=dataset,
            zorder=3,
        )
        ax.plot(
            x_positions[1:3],
            y[1:3],
            marker="o",
            linewidth=2.0,
            markersize=6.5,
            linestyle="--",
            color=DATASET_COLORS[dataset],
            alpha=0.75,
            zorder=3,
        )
        ax.plot(
            x_positions[2:],
            y[2:],
            marker="o",
            linewidth=2.5,
            markersize=6.5,
            color=DATASET_COLORS[dataset],
            zorder=3,
        )

        for idx, row in frame.iterrows():
            y_value = float(row["count"])
            y_text = y_value * (1.12 if dataset == "DMF5" else 0.9)
            va = "bottom" if y_text >= y_value else "top"
            if row["stage"] == "Consensus UMI + singleton rescue":
                y_text = y_value * (1.2 if dataset == "DMF5" else 0.86)
                va = "bottom" if y_text >= y_value else "top"
            ax.text(
                idx,
                y_text,
                format_count(y_value),
                ha="center",
                va=va,
                fontsize=7.8,
                color=DATASET_COLORS[dataset],
            )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(
        [
            "Raw reads",
            "After trimming",
            "Consensus UMI",
            "Consensus UMI\n+ singleton rescue",
            "Extraction PASS",
        ],
        rotation=18,
        ha="right",
    )
    ax.set_yscale("log")
    visible_counts = stage_counts[stage_counts["stage"].isin(stage_order)]["count"]
    ax.set_ylim(visible_counts.min() * 0.55, visible_counts.max() * 2.2)
    ax.set_xlim(-0.35, len(stage_order) - 0.65)
    ax.axvspan(1.45, 1.62, color="0.88", alpha=0.8, zorder=0)
    ax.text(
        1.535,
        0.98,
        "UNIT switch:\nreads -> variants",
        transform=ax.get_xaxis_transform(),
        ha="center",
        va="top",
        rotation=90,
        fontsize=7.5,
        color="0.25",
    )
    ax.set_ylabel("Count, log scale")
    ax.set_title("B  Linear retention with consensus UMI and singleton rescue", loc="left", fontsize=13)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False)


def _dropout_bar_data(dropout: pd.DataFrame) -> tuple[list[float], list[str], dict[str, list[float]], list[str]]:
    steps = ["01_preprocessing", "02_umi_consensus", "02a_singleton_rescue", "03_extraction"]
    bar_labels = []
    positions = []
    step_keys = []
    x = 0
    for dataset in ("DMF5", "TCRA3"):
        for step in steps:
            if not dropout[(dropout["dataset"] == dataset) & (dropout["step"] == step)].empty:
                positions.append(x)
                label = STEP_LABELS.get(step, step).replace(" ", "\n", 1)
                bar_labels.append(f"{dataset}\n{label}")
                step_keys.append((dataset, step))
                x += 1
        x += 0.45

    reason_totals = dropout.groupby("loss_reason")["reads_lost"].sum()
    reason_order = [LOSS_LABELS[key] for key in LOSS_LABELS if reason_totals.get(LOSS_LABELS[key], 0) > 0]
    heights_by_reason = {}
    for reason in reason_order:
        heights = []
        for dataset, step in step_keys:
            match = dropout[
                (dropout["dataset"] == dataset)
                & (dropout["step"] == step)
                & (dropout["loss_reason"] == reason)
            ]
            heights.append(float(match["reads_lost"].sum()))
        heights_by_reason[reason] = heights
    return positions, bar_labels, heights_by_reason, reason_order


def draw_dropout_broken(
    ax_top: plt.Axes,
    ax_bottom: plt.Axes,
    dropout: pd.DataFrame,
    show_legend: bool = True,
) -> None:
    positions, bar_labels, heights_by_reason, reason_order = _dropout_bar_data(dropout)
    palette = sns.color_palette("tab20", n_colors=len(reason_order))
    color_map = dict(zip(reason_order, palette))

    totals = pd.Series(0.0, index=positions)
    for heights in heights_by_reason.values():
        totals = totals + pd.Series(heights, index=positions)
    max_total = float(totals.max())
    non_preprocessing_totals = totals[
        ["preprocessing" not in label for label in bar_labels]
    ]
    lower_focus = (
        float(non_preprocessing_totals.max())
        if not non_preprocessing_totals.empty
        else float(totals.sort_values(ascending=False).iloc[-1])
    )
    bottom_ylim = max(lower_focus * 1.25, max_total * 0.04)
    clipped_totals = totals[totals > bottom_ylim]
    top_low = max(float(clipped_totals.min()) * 0.72, bottom_ylim * 1.35)

    for ax in (ax_top, ax_bottom):
        bottom = pd.Series(0.0, index=positions)
        for reason in reason_order:
            heights = heights_by_reason[reason]
            ax.bar(positions, heights, bottom=bottom, color=color_map[reason], label=reason, width=0.72)
            bottom = bottom + pd.Series(heights, index=positions)
        ax.grid(axis="y", alpha=0.25)

    ax_top.set_ylim(top_low, max_total * 1.05)
    ax_bottom.set_ylim(0, bottom_ylim)
    ax_top.spines["bottom"].set_visible(False)
    ax_bottom.spines["top"].set_visible(False)
    ax_top.tick_params(labelbottom=False, bottom=False)
    ax_bottom.set_xticks(positions)
    ax_bottom.set_xticklabels(bar_labels, rotation=35, ha="right", fontsize=8)
    ax_bottom.set_ylabel("Reads lost")
    ax_top.set_title("C  Read loss reasons by pipeline step", loc="left", fontsize=13)
    for pos, total in totals.items():
        if total > bottom_ylim:
            ax_top.text(pos, total * 1.015, format_count(total), ha="center", va="bottom", fontsize=7.5)
            ax_bottom.text(pos, bottom_ylim * 0.96, format_count(total), ha="center", va="top", fontsize=7, rotation=90, color="0.25")

    diagonal_kwargs = dict(marker=[(-1, -0.5), (1, 0.5)], markersize=8, linestyle="none", color="k", mec="k", mew=1)
    ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **diagonal_kwargs)
    ax_bottom.plot([0, 1], [1, 1], transform=ax_bottom.transAxes, **diagonal_kwargs)

    if show_legend:
        ax_top.legend(ncol=2, fontsize=7, frameon=True, loc="upper left", framealpha=0.92)


def draw_variant_panel(ax: plt.Axes, variants: pd.DataFrame, dataset: str) -> None:
    frame = variants[variants["dataset"] == dataset].copy()
    order = (
        frame.groupby("peptide")["variant_count"]
        .sum()
        .sort_values(ascending=True)
        .index.tolist()
    )
    y_positions = range(len(order))
    left = pd.Series(0.0, index=order)
    for label_class, color in LABEL_COLORS.items():
        values = (
            frame[frame["label_class"] == label_class]
            .set_index("peptide")
            .reindex(order)["variant_count"]
            .fillna(0)
        )
        ax.barh(y_positions, values, left=left.loc[order], color=color, label=label_class)
        left = left.add(values, fill_value=0)
    ax.set_yticks(list(y_positions))
    ax.set_yticklabels(order)
    ax.set_xlabel("Labeled variants")
    ax.set_title(dataset, fontsize=12)
    ax.grid(axis="x", alpha=0.25)
    ax.xaxis.set_major_formatter(lambda x, _pos: format_count(x))
    for y, peptide in zip(y_positions, order):
        total = left.loc[peptide]
        ax.text(total * 1.01, y, format_count(total), va="center", fontsize=8)


def draw_variant_panels(fig: plt.Figure, left_spec, right_spec, variants: pd.DataFrame) -> tuple[plt.Axes, plt.Axes]:
    ax_d1 = fig.add_subplot(left_spec)
    ax_d2 = fig.add_subplot(right_spec, sharex=ax_d1)
    draw_variant_panel(ax_d1, variants, "DMF5")
    draw_variant_panel(ax_d2, variants, "TCRA3")
    ax_d1.set_title("D  DMF5", loc="left", fontsize=13)
    ax_d2.set_title("TCRA3", fontsize=12)
    return ax_d1, ax_d2


def draw_figure(stage_counts: pd.DataFrame, consensus_bridge: pd.DataFrame, dropout: pd.DataFrame, variants: pd.DataFrame) -> plt.Figure:
    configure_plot_theme(font_scale=1.05)
    fig = plt.figure(figsize=(16, 10.8))
    grid = GridSpec(2, 2, figure=fig, height_ratios=[1.12, 1.08], width_ratios=[1.08, 1.22])

    ax_b = fig.add_subplot(grid[0, 0])
    draw_linear_retention(ax_b, stage_counts, consensus_bridge)

    dropout_grid = grid[0, 1].subgridspec(2, 1, height_ratios=[0.45, 1.0], hspace=0.05)
    ax_c_top = fig.add_subplot(dropout_grid[0])
    ax_c_bottom = fig.add_subplot(dropout_grid[1], sharex=ax_c_top)
    draw_dropout_broken(ax_c_top, ax_c_bottom, dropout)

    ax_d1, ax_d2 = draw_variant_panels(fig, grid[1, 0], grid[1, 1], variants)
    handles, labels = ax_d1.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.01))
    ax_d1.legend_.remove() if ax_d1.legend_ else None
    ax_d2.legend_.remove() if ax_d2.legend_ else None

    fig.suptitle("NGS pipeline read retention and final TCR-peptide label yield", fontsize=16, y=0.985)
    fig.subplots_adjust(left=0.07, right=0.96, top=0.92, bottom=0.1, hspace=0.62, wspace=0.34)
    return fig


def make_panel_b_figure(stage_counts: pd.DataFrame, consensus_bridge: pd.DataFrame) -> plt.Figure:
    configure_plot_theme(font_scale=1.1)
    fig, ax = plt.subplots(figsize=(10.4, 6.2))
    draw_linear_retention(ax, stage_counts, consensus_bridge)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.9, bottom=0.3)
    return fig


def make_panel_c_figure(dropout: pd.DataFrame) -> plt.Figure:
    configure_plot_theme(font_scale=1.1)
    fig = plt.figure(figsize=(9.5, 6.2))
    grid = fig.add_gridspec(2, 1, height_ratios=[0.45, 1.0], hspace=0.05)
    ax_top = fig.add_subplot(grid[0])
    ax_bottom = fig.add_subplot(grid[1], sharex=ax_top)
    draw_dropout_broken(ax_top, ax_bottom, dropout)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.9, bottom=0.25)
    return fig


def make_panel_d_figure(variants: pd.DataFrame) -> plt.Figure:
    configure_plot_theme(font_scale=1.1)
    fig = plt.figure(figsize=(13, 6.5))
    grid = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.2], wspace=0.34)
    ax_d1, ax_d2 = draw_variant_panels(fig, grid[0], grid[1], variants)
    handles, labels = ax_d1.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False)
    ax_d1.legend_.remove() if ax_d1.legend_ else None
    ax_d2.legend_.remove() if ax_d2.legend_ else None
    fig.subplots_adjust(left=0.08, right=0.98, top=0.9, bottom=0.16)
    return fig


def write_validation_summary(
    duplicate_report: pd.DataFrame,
    validation: pd.DataFrame,
    stage_checks: pd.DataFrame,
    variants: pd.DataFrame,
) -> None:
    lines = [
        "NGS thesis pipeline overview validation",
        "",
        "Canonical inputs:",
        "- DMF5: Q30 preprocessing, UMI, extraction, and combined variant labeling.",
        "- TCRA3: Q30 P3408 run plus KH157 Q30 run for read metrics; rescue totals cover the canonical P3408 run only.",
        "",
        "Source integrity checks:",
    ]
    for _, row in duplicate_report.iterrows():
        lines.append(
            f"- {row['dataset']} {row['project']} {row['source_file']}: "
            f"duplicate_rows={int(row['duplicate_rows'])}; rows_removed={int(row['duplicate_row_pairs_removed'])}; "
            f"has_exact_duplicates={bool(row['has_exact_duplicates'])}"
        )

    lines.extend(["", "Read-accounting deltas by step:"])
    summary = (
        validation.assign(abs_delta=validation["delta"].abs())
        .groupby(["dataset", "step"], as_index=False)
        .agg(
            samples=("sample_id", "count"),
            exact_samples=("passes_exact_accounting", "sum"),
            max_abs_delta=("abs_delta", "max"),
            total_delta=("delta", "sum"),
        )
    )
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['dataset']} {row['step']}: {int(row['exact_samples'])}/{int(row['samples'])} "
            f"samples exact; max_abs_delta={int(row['max_abs_delta'])}; total_delta={int(row['total_delta'])}"
        )

    lines.extend(["", "Cross-step consistency checks:"])
    for _, row in stage_checks.iterrows():
        lines.append(
            f"- {row['dataset']} {row['check']}: left={int(row['left_value'])}, "
            f"right={int(row['right_value'])}, delta={int(row['delta'])}, passes={bool(row['passes'])}"
        )

    lines.extend(["", "Variant-label peptide counts:"])
    peptide_counts = variants.groupby("dataset")["peptide"].nunique()
    for dataset, count in peptide_counts.items():
        lines.append(f"- {dataset}: {count} peptides")

    (OUTPUT_DIR / "validation_summary.txt").write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create thesis-ready NGS read-retention and label-yield figure."
    )
    parser.add_argument(
        "--formats",
        nargs="+",
        default=("png",),
        choices=("png", "pdf", "svg"),
        help="Figure formats to export. PNG is the default and fastest on the cluster.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    start = time.perf_counter()
    print("Reading metrics...", flush=True)
    metrics, duplicate_report = load_metric_sources()
    print(f"Read {len(metrics):,} metric rows after deduplication.", flush=True)
    stage_counts = build_stage_counts(metrics)
    consensus_bridge = build_consensus_bridge(metrics)
    rescue_summary = build_rescue_summary(metrics)
    dropout = build_dropout_table(metrics)
    print("Reading variant summaries...", flush=True)
    variants = read_variant_summaries()
    validation = validate_accounting(metrics)
    stage_checks = build_stage_consistency_checks(metrics)

    print("Writing tables and validation summary...", flush=True)
    write_outputs(
        metrics,
        stage_counts,
        consensus_bridge,
        duplicate_report,
        rescue_summary,
        dropout,
        variants,
        validation,
        stage_checks,
    )
    write_validation_summary(duplicate_report, validation, stage_checks, variants)

    print("Drawing figure...", flush=True)
    fig = draw_figure(stage_counts, consensus_bridge, dropout, variants)
    for ext in args.formats:
        print(f"Saving {ext}...", flush=True)
        fig.savefig(OUTPUT_DIR / f"thesis_ngs_pipeline_overview.{ext}", dpi=300)
    plt.close(fig)

    elapsed = time.perf_counter() - start
    print(f"Wrote outputs to {OUTPUT_DIR} in {elapsed:.1f} seconds.", flush=True)


if __name__ == "__main__":
    main()
