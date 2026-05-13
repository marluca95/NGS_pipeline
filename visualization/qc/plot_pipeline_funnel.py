#!/usr/bin/env python3
"""Plot per-sample read-loss funnels across preprocessing, UMI, singleton rescue, and extraction.

Two input modes are supported:
1) Legacy: pass explicit per_sample_metrics.tsv files via --step01/--step02/--step02a/--step03
2) Project autodetect: pass --project_dir and parse current pipeline outputs directly

Examples
--------
Project autodetect (recommended for existing runs):
    python visualization/qc/plot_pipeline_funnel.py \
      --project_dir data/P3481_LUCA-TCRDMF5 \
      --output visualization/qc/P3481_pipeline_funnel.pdf

Legacy per-sample metrics mode:
    python visualization/qc/plot_pipeline_funnel.py \
      --step01 data/P3408_LUCA-TCRA3/01_preprocessed/per_sample_metrics.tsv \
      --step02 data/P3408_LUCA-TCRA3/02_umi_consensus/per_sample_metrics.tsv \
      --step02a data/P3408_LUCA-TCRA3/02_umi_consensus/per_sample_metrics.tsv \
      --step03 data/P3408_LUCA-TCRA3/03_extracted/per_sample_metrics.tsv \
      --output visualization/qc/pipeline_funnel.pdf
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------
# Schema + plot config
# ---------------------------------------------------------------------------
HEADER = ["sample_id", "sample_name", "step", "metric", "value"]

STEPS = [
    ("01_preprocessing", "01\nTrimming"),
    ("02_umi_consensus", "02\nUMI\nconsensus"),
    ("02a_singleton_rescue", "02a\nSingleton\nrescue"),
    ("03_extraction", "03\nExtraction"),
]

READS_IN_METRIC = {
    "01_preprocessing": "reads_in_total",
    "02_umi_consensus": "reads_in",
    "02a_singleton_rescue": "reads_in",
    "03_extraction": "reads_in_total",
}

READS_PASS_METRIC = {
    "01_preprocessing": "reads_out_total",
    "02_umi_consensus": "reads_in_consensus",
    "02a_singleton_rescue": "reads_rescued",
    "03_extraction": "reads_pass",
}

LOSS_REASONS = {
    "01_preprocessing": [
        ("reads_lost_trimming", "Too short / low-quality"),
    ],
    "02_umi_consensus": [
        ("reads_lost_too_short", "Too short (pre-anchor)"),
        ("reads_lost_no_anchor", "No anchor found"),
        ("reads_lost_anchor_wrong_pos", "Anchor wrong position"),
        ("reads_lost_empty_insert", "Empty insert"),
        ("reads_in_singletons", "Singletons (1 read / UMI)"),
    ],
    "02a_singleton_rescue": [
        ("reads_lost_no_anchor", "No anchor in insert"),
        ("reads_lost_too_short", "Too short after anchor"),
        ("reads_extracted_not_rescued", "Extracted but not rescued"),
    ],
    "03_extraction": [
        ("reads_lost_no_anchor", "No anchor"),
        ("reads_lost_too_short", "Region too short"),
        ("reads_lost_contains_N", "Contains N"),
        ("reads_lost_wrong_length", "Wrong length"),
        ("reads_lost_stop_codon", "Stop codon"),
        ("reads_lost_validation", "Fails library check"),
    ],
}

COLOUR_PASS = "#4CAF50"
COLOUR_LOST = "#F44336"
COLOUR_DETAIL = [
    "#EF9A9A",
    "#EF5350",
    "#B71C1C",
    "#FFCC80",
    "#FF9800",
    "#E65100",
    "#80DEEA",
    "#00ACC1",
    "#006064",
]


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------
def add_rows(rows: List[dict], sample_id: str, sample_name: str, step: str, metrics: Dict[str, float]) -> None:
    for metric, value in metrics.items():
        rows.append(
            {
                "sample_id": sample_id,
                "sample_name": sample_name,
                "step": step,
                "metric": metric,
                "value": float(value),
            }
        )


def parse_sample_label(label: str) -> Tuple[str, str]:
    parts = label.split("_", 1)
    if len(parts) == 2:
        return parts[0], parts[1]
    return label, label


def find_existing_dir(project_dir: Path, candidates: Iterable[str]) -> Optional[Path]:
    for name in candidates:
        p = project_dir / name
        if p.exists() and p.is_dir():
            return p
    return None


# ---------------------------------------------------------------------------
# Legacy per_sample_metrics.tsv loading
# ---------------------------------------------------------------------------
def load_per_sample_metrics(paths: Iterable[Optional[Path]]) -> pd.DataFrame:
    frames = []
    for p in paths:
        if p is None or not p.exists():
            continue
        df = pd.read_csv(p, sep="\t")
        missing = [c for c in HEADER if c not in df.columns]
        if missing:
            raise ValueError(f"{p} missing required columns: {missing}")
        frames.append(df[HEADER])

    if not frames:
        raise FileNotFoundError("No metrics files found. Check --step01 … --step03 paths.")

    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Project autodetect parsers
# ---------------------------------------------------------------------------
def load_step01_from_preprocess(pre_dir: Path) -> pd.DataFrame:
    per_sample = pre_dir / "per_sample_metrics.tsv"
    if per_sample.exists():
        return pd.read_csv(per_sample, sep="\t")

    csv_files = sorted(pre_dir.rglob("bbduk_summary*.csv"), key=lambda p: p.stat().st_mtime)
    if not csv_files:
        return pd.DataFrame(columns=HEADER)

    bbduk_path = csv_files[-1]
    df = pd.read_csv(bbduk_path)
    if df.empty:
        return pd.DataFrame(columns=HEADER)

    grouped = (
        df.groupby(["sample_id", "sample_name"], as_index=False)[["reads_in", "reads_out"]]
        .sum()
        .rename(columns={"reads_in": "reads_in_total", "reads_out": "reads_out_total"})
    )
    grouped["reads_lost_trimming"] = grouped["reads_in_total"] - grouped["reads_out_total"]

    rows: List[dict] = []
    for _, r in grouped.iterrows():
        add_rows(
            rows,
            str(r["sample_id"]),
            str(r["sample_name"]),
            "01_preprocessing",
            {
                "reads_in_total": int(r["reads_in_total"]),
                "reads_out_total": int(r["reads_out_total"]),
                "reads_lost_trimming": int(r["reads_lost_trimming"]),
            },
        )
    return pd.DataFrame(rows, columns=HEADER)


def parse_umi_summary_file(path: Path) -> Optional[Dict[str, int]]:
    key_patterns = {
        "reads_in": r"^Total reads processed:\s*(\d+)$",
        "reads_lost_too_short": r"^Reads too short .*\(dropped\):\s*(\d+)$",
        "reads_lost_no_anchor": r"^Reads with no anchor found \(dropped\):\s*(\d+)$",
        "reads_lost_anchor_wrong_pos": r"^Reads with anchor at invalid position \(dropped\):\s*(\d+)$",
        "reads_lost_empty_insert": r"^Reads with empty insert after anchor \(dropped\):\s*(\d+)$",
        "reads_in_consensus": r"^Reads contributing to consensus UMIs:\s*(\d+)$",
        "reads_in_singletons": r"^Reads in singleton UMIs:\s*(\d+)$",
    }

    label = path.name.replace("_umi_summary.txt", "")
    sample_id, sample_name = parse_sample_label(label)

    metrics: Dict[str, int] = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        for key, pat in key_patterns.items():
            m = re.match(pat, line)
            if m:
                metrics[key] = int(m.group(1))

    if not metrics:
        return None

    return {
        "sample_id": sample_id,
        "sample_name": sample_name,
        **metrics,
    }


def load_step02_from_umi(umi_dir: Path) -> pd.DataFrame:
    per_sample = umi_dir / "per_sample_metrics.tsv"
    if per_sample.exists():
        df = pd.read_csv(per_sample, sep="\t")
        return df[df["step"] == "02_umi_consensus"][HEADER]

    rows: List[dict] = []
    for p in sorted(umi_dir.glob("*_umi_summary.txt")):
        parsed = parse_umi_summary_file(p)
        if not parsed:
            continue
        sample_id = parsed.pop("sample_id")
        sample_name = parsed.pop("sample_name")
        add_rows(rows, sample_id, sample_name, "02_umi_consensus", parsed)

    return pd.DataFrame(rows, columns=HEADER)


def load_step02a_from_umi(umi_dir: Path) -> pd.DataFrame:
    per_sample = umi_dir / "per_sample_metrics.tsv"
    if per_sample.exists():
        df = pd.read_csv(per_sample, sep="\t")
        step = df[df["step"] == "02a_singleton_rescue"][HEADER]
        if not step.empty:
            return step

    summary = umi_dir / "singleton_rescue_summary.tsv"
    if not summary.exists():
        return pd.DataFrame(columns=HEADER)

    df = pd.read_csv(summary, sep="\t")
    rows: List[dict] = []
    for _, r in df.iterrows():
        extracted_not_rescued = int(r["extracted_reads"]) - int(r["kept_reads"])
        add_rows(
            rows,
            str(r["sample_id"]),
            str(r["sample_name"]),
            "02a_singleton_rescue",
            {
                "reads_in": int(r["total_singleton_reads"]),
                "reads_lost_no_anchor": int(r["no_anchor_reads"]),
                "reads_lost_too_short": int(r["too_short_after_anchor_reads"]),
                "reads_extracted": int(r["extracted_reads"]),
                "reads_rescued": int(r["kept_reads"]),
                "reads_extracted_not_rescued": max(0, extracted_not_rescued),
            },
        )
    return pd.DataFrame(rows, columns=HEADER)


def parse_extraction_summary(path: Path) -> Optional[Dict[str, int]]:
    wanted = {
        "total_reads": "reads_in_total",
        "pass": "reads_pass",
        "fail_anchor_not_found": "reads_lost_no_anchor",
        "fail_region_too_short": "reads_lost_too_short",
        "fail_contains_N": "reads_lost_contains_N",
        "fail_wrong_length": "reads_lost_wrong_length",
        "fail_stop_codon": "reads_lost_stop_codon",
        "fail_degenerate_mismatch": "reads_lost_validation_1",
        "fail_too_many_mutations": "reads_lost_validation_2",
        "fail_mutation_not_NNK": "reads_lost_validation_3",
    }

    label = path.name.replace(".filtered.summary.tsv", "")
    sample_id, sample_name = parse_sample_label(label)

    metrics: Dict[str, int] = {}
    with open(path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#") or "\t" not in line:
                continue
            key, value = line.split("\t", 1)
            if key not in wanted:
                continue
            try:
                metrics[wanted[key]] = int(value)
            except ValueError:
                continue

    if not metrics:
        return None

    validation = (
        metrics.pop("reads_lost_validation_1", 0)
        + metrics.pop("reads_lost_validation_2", 0)
        + metrics.pop("reads_lost_validation_3", 0)
    )
    metrics["reads_lost_validation"] = validation

    return {
        "sample_id": sample_id,
        "sample_name": sample_name,
        **metrics,
    }


def load_step03_from_extraction(ex_dir: Path) -> pd.DataFrame:
    per_sample = ex_dir / "per_sample_metrics.tsv"
    if per_sample.exists():
        df = pd.read_csv(per_sample, sep="\t")
        return df[df["step"] == "03_extraction"][HEADER]

    rows: List[dict] = []
    for p in sorted(ex_dir.glob("*.filtered.summary.tsv")):
        parsed = parse_extraction_summary(p)
        if not parsed:
            continue
        sample_id = parsed.pop("sample_id")
        sample_name = parsed.pop("sample_name")
        add_rows(rows, sample_id, sample_name, "03_extraction", parsed)

    return pd.DataFrame(rows, columns=HEADER)


def load_metrics_from_project(project_dir: Path) -> pd.DataFrame:
    pre_dir = find_existing_dir(project_dir, ["01_preprocessed", "preprocessed"])
    umi_dir = find_existing_dir(project_dir, ["02_umi_consensus", "umi"])
    ex_dir = find_existing_dir(project_dir, ["03_extracted", "extracted"])

    frames = []
    if pre_dir:
        frames.append(load_step01_from_preprocess(pre_dir))
    if umi_dir:
        frames.append(load_step02_from_umi(umi_dir))
        frames.append(load_step02a_from_umi(umi_dir))
    if ex_dir:
        frames.append(load_step03_from_extraction(ex_dir))

    frames = [f for f in frames if f is not None and not f.empty]
    if not frames:
        raise FileNotFoundError(
            f"No metrics were found in project directory: {project_dir}. "
            "Expected folders like preprocessed/umi/extracted or 01_preprocessed/02_umi_consensus/03_extracted."
        )

    out = pd.concat(frames, ignore_index=True)
    out["value"] = pd.to_numeric(out["value"], errors="coerce").fillna(0)
    return out


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def pivot_sample(df: pd.DataFrame, sample_id: str) -> pd.DataFrame:
    sub = df[df["sample_id"] == sample_id].copy()
    return sub.pivot_table(index="step", columns="metric", values="value", aggfunc="sum")


def plot_funnel_for_sample(
    wide: pd.DataFrame,
    sample_id: str,
    sample_name: str,
    ax_top: plt.Axes,
    ax_bot: plt.Axes,
) -> None:
    present_steps = []
    for step_key, step_label in STEPS:
        if step_key in wide.index:
            present_steps.append((step_key, step_label))

    if not present_steps:
        ax_top.set_title(f"{sample_id}  {sample_name}")
        ax_top.text(0.5, 0.5, "No metrics for this sample", ha="center", va="center")
        ax_top.axis("off")
        ax_bot.axis("off")
        return

    x_map = {step: i for i, (step, _) in enumerate(present_steps)}
    x = list(range(len(present_steps)))
    labels = [label for _, label in present_steps]

    pass_vals = []
    lost_vals = []

    for step_key, _ in present_steps:
        row = wide.loc[step_key]
        r_in = int(row.get(READS_IN_METRIC[step_key], 0) or 0)
        r_pass = int(row.get(READS_PASS_METRIC[step_key], 0) or 0)
        r_lost = max(0, r_in - r_pass)
        pass_vals.append(r_pass)
        lost_vals.append(r_lost)

    ax_top.bar(x, pass_vals, color=COLOUR_PASS, label="Passed")
    ax_top.bar(x, lost_vals, bottom=pass_vals, color=COLOUR_LOST, label="Lost")
    ax_top.set_xticks(x)
    ax_top.set_xticklabels(labels, fontsize=8)
    ax_top.set_ylabel("Read count")
    ax_top.set_title(f"{sample_id}  {sample_name}", fontsize=10, fontweight="bold")
    ax_top.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax_top.legend(fontsize=8, loc="upper right")

    bottom = [0] * len(present_steps)
    colour_idx = 0

    for step_key, _ in present_steps:
        row = wide.loc[step_key]
        x_idx = x_map[step_key]
        for metric_key, reason_label in LOSS_REASONS.get(step_key, []):
            val = int(row.get(metric_key, 0) or 0)
            if val <= 0:
                continue
            colour = COLOUR_DETAIL[colour_idx % len(COLOUR_DETAIL)]
            ax_bot.bar(x_idx, val, bottom=bottom[x_idx], color=colour, label=reason_label)
            bottom[x_idx] += val
            colour_idx += 1

    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(labels, fontsize=8)
    ax_bot.set_ylabel("Reads lost")
    ax_bot.set_title("Loss breakdown", fontsize=9)
    ax_bot.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}"))

    handles, labels = ax_bot.get_legend_handles_labels()
    dedup = dict(zip(labels, handles))
    if dedup:
        ax_bot.legend(dedup.values(), dedup.keys(), fontsize=6, loc="upper right", ncol=2, framealpha=0.7)


def make_funnel_figure(df: pd.DataFrame, sample_id: str) -> plt.Figure:
    sub = df[df["sample_id"] == sample_id]
    if sub.empty:
        raise ValueError(f"sample_id={sample_id!r} not found in data")

    sample_name = str(sub["sample_name"].iloc[0])
    wide = pivot_sample(df, sample_id)

    fig, (ax_top, ax_bot) = plt.subplots(
        2,
        1,
        figsize=(max(6, len(STEPS) * 2), 8),
        gridspec_kw={"height_ratios": [3, 2]},
    )
    fig.subplots_adjust(hspace=0.45)

    plot_funnel_for_sample(wide, sample_id, sample_name, ax_top, ax_bot)
    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pipeline read-loss funnel plot.")

    p.add_argument("--project_dir", type=Path, default=None, help="Project data dir (autodetect mode)")

    p.add_argument("--step01", type=Path, default=None, help="per_sample_metrics.tsv from step 01")
    p.add_argument("--step02", type=Path, default=None, help="per_sample_metrics.tsv from step 02")
    p.add_argument("--step02a", type=Path, default=None, help="per_sample_metrics.tsv from step 02a")
    p.add_argument("--step03", type=Path, default=None, help="per_sample_metrics.tsv from step 03")

    p.add_argument(
        "--output",
        type=Path,
        default=Path("pipeline_funnel.pdf"),
        help="Output file (.pdf or .png). Default: pipeline_funnel.pdf",
    )
    p.add_argument("--sample_id", default=None, help="Plot only this sample_id (default: all samples)")

    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.project_dir is not None:
        df = load_metrics_from_project(args.project_dir)
    else:
        df = load_per_sample_metrics([args.step01, args.step02, args.step02a, args.step03])

    sample_ids = [args.sample_id] if args.sample_id else sorted(df["sample_id"].astype(str).unique())
    if not sample_ids:
        print("No samples found in metrics files.", file=sys.stderr)
        sys.exit(1)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    suffix = args.output.suffix.lower()

    if suffix == ".pdf":
        from matplotlib.backends.backend_pdf import PdfPages

        pages = 0
        with PdfPages(args.output) as pdf:
            for sid in sample_ids:
                if sid not in set(df["sample_id"].astype(str)):
                    print(f"WARNING: sample_id={sid!r} not found in metrics, skipping.")
                    continue
                fig = make_funnel_figure(df, sid)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
                pages += 1
        print(f"Wrote {pages} page(s) to {args.output}")
    else:
        written = 0
        for sid in sample_ids:
            if sid not in set(df["sample_id"].astype(str)):
                continue
            fig = make_funnel_figure(df, sid)
            out = args.output.parent / f"{sid}_{args.output.name}"
            fig.savefig(out, bbox_inches="tight", dpi=150)
            plt.close(fig)
            print(f"Wrote {out}")
            written += 1
        if written == 0:
            print("No output written (no matching samples).", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
