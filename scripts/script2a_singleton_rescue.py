#!/usr/bin/env python3
import argparse
import gzip
import logging
import traceback
from pathlib import Path
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from utils.config_utils import load_pipeline_config
from utils.logging_utils import setup_pipeline_logging
from utils.metrics_utils import write_sample_metrics
from utils.sample_utils import load_sample_sheet

SCRIPT_NAME = "script2a_singleton_rescue"
logger = logging.getLogger(__name__)


# ---------------- config / IO ----------------
def parse_yaml(path: str) -> dict:
    cfg = load_pipeline_config(
        path,
        required_keys=("samples_tsv", "singleton_dir", "output_dir", "wt_anchor"),
        default_values={
            "var_len": 33,
            "min_umis": 2,
            "singleton_suffix": "_singletons.fastq.gz",
            "logs_dir": None,
        },
        path_keys=("samples_tsv", "singleton_dir", "output_dir", "logs_dir"),
    )
    if not cfg.get("logs_dir"):
        cfg["logs_dir"] = str(Path(cfg["output_dir"]).parent / "_logs")
    return cfg


def load_samples(samples_tsv: str) -> pd.DataFrame:
    return load_sample_sheet(
        samples_tsv,
        required_cols=("sample_id", "sample_name"),
        dedupe_on=("sample_id", "sample_name"),
    )


# ---------------- core logic ----------------
def parse_umi(rec_id: str) -> str:
    return rec_id.split("|n=")[0]


def extract_var(insert: str, wt_anchor: str, var_len: int):
    """
    Extract variable region as the next `var_len` bases immediately after the WT anchor.
    Uses left-most match via .find(). If you ever observe multiple hits, consider switching to .rfind().
    """
    p = insert.find(wt_anchor)
    if p == -1:
        return None
    start = p + len(wt_anchor)
    end = start + var_len
    if end > len(insert):
        return None
    return insert[start:end]


def rescue_one_sample(
    singletons_fastq_gz: Path,
    out_fastq_gz: Path,
    out_tsv: Path,
    wt_anchor: str,
    var_len: int,
    min_umis: int,
):
    """
    Reads singleton FASTQ (where sequence = INSERT), extracts 33bp variable region after WT anchor,
    keeps only variants observed in >= min_umis distinct UMIs.
    Writes:
      - filtered singleton FASTQ (only reads whose variant passes threshold)
      - per-sample variants TSV: variant_33bp, n_umis
    Returns stats dict.
    """
    variant_to_umis = defaultdict(set)
    cached = []  # (record, variant)

    n_total = 0
    n_no_anchor = 0
    n_too_short_after_anchor = 0
    n_extracted = 0

    with gzip.open(singletons_fastq_gz, "rt") as handle:
        for rec in SeqIO.parse(handle, "fastq"):
            n_total += 1
            insert = str(rec.seq)

            var = extract_var(insert, wt_anchor, var_len)
            if var is None:
                # distinguish reasons
                if insert.find(wt_anchor) == -1:
                    n_no_anchor += 1
                else:
                    n_too_short_after_anchor += 1
                continue

            n_extracted += 1
            umi = parse_umi(rec.id)
            variant_to_umis[var].add(umi)
            cached.append((rec, var))

    # Filter variants by UMI count
    keep_vars = {v for v, umis in variant_to_umis.items() if len(umis) >= min_umis}

    out_fastq_gz.parent.mkdir(parents=True, exist_ok=True)

    kept_reads = 0
    with gzip.open(out_fastq_gz, "wt") as out:
        for rec, var in cached:
            if var in keep_vars:
                SeqIO.write(rec, out, "fastq")
                kept_reads += 1

    with open(out_tsv, "w") as f:
        f.write("variant_33bp\tn_umis\n")
        for v in sorted(keep_vars, key=lambda x: len(variant_to_umis[x]), reverse=True):
            f.write(f"{v}\t{len(variant_to_umis[v])}\n")

    # Compute summary stats
    extracted_rate = (n_extracted / n_total) if n_total else 0.0
    no_anchor_rate = (n_no_anchor / n_total) if n_total else 0.0
    too_short_rate = (n_too_short_after_anchor / n_total) if n_total else 0.0
    kept_read_rate = (kept_reads / n_extracted) if n_extracted else 0.0

    return {
        "total_singleton_reads": n_total,
        "extracted_reads": n_extracted,
        "no_anchor_reads": n_no_anchor,
        "too_short_after_anchor_reads": n_too_short_after_anchor,
        "n_variants_total": len(variant_to_umis),
        "n_variants_kept_ge_min_umis": len(keep_vars),
        "kept_reads": kept_reads,
        "extracted_rate": extracted_rate,
        "no_anchor_rate": no_anchor_rate,
        "too_short_rate": too_short_rate,
        "kept_read_rate_of_extracted": kept_read_rate,
    }


# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser(description="Singleton rescue: extract 33bp after WT anchor, keep variants seen >=2 UMIs")
    ap.add_argument("--yaml_config", required=True)
    ap.add_argument("--sample_id", default=None, help="Optional: run only one sample_id")
    args = ap.parse_args()

    cfg = parse_yaml(args.yaml_config)

    # logging
    log_file = setup_pipeline_logging(
        logs_dir=cfg["logs_dir"],
        script_name=SCRIPT_NAME,
        scope="run",
        run_label=cfg.get("run_label"),
    )
    logger.info("Logging to: %s", log_file)
    logger.info("Using config: %s", cfg.get("_config_path", args.yaml_config))

    df = load_samples(cfg["samples_tsv"])
    logger.info(f"Loaded {len(df)} unique samples from {cfg['samples_tsv']}")

    if args.sample_id is not None:
        df = df[df["sample_id"].astype(str) == str(args.sample_id)].reset_index(drop=True)
        if df.empty:
            raise ValueError(f"No rows for sample_id={args.sample_id}")
        logger.info(f"Filtered to sample_id={args.sample_id}: {len(df)} sample(s)")

    singleton_dir = Path(cfg["singleton_dir"])
    out_dir = Path(cfg["output_dir"])
    suffix = cfg["singleton_suffix"]

    wt_anchor = str(cfg["wt_anchor"])
    var_len = int(cfg["var_len"])
    min_umis = int(cfg["min_umis"])

    logger.info(f"Singleton dir: {singleton_dir}")
    logger.info(f"Output dir:    {out_dir}")
    logger.info(f"WT anchor:     {wt_anchor}")
    logger.info(f"Var length:    {var_len}")
    logger.info(f"Min UMIs:      {min_umis}")

    all_stats = []
    n_missing = 0

    for _, row in df.iterrows():
        sample_id = str(row["sample_id"])
        sample_name = str(row["sample_name"])
        sample_label = f"{sample_id}_{sample_name}"

        in_fq = singleton_dir / f"{sample_label}{suffix}"
        if not in_fq.exists():
            logger.warning(f"[SKIP] Missing singleton file: {in_fq}")
            n_missing += 1
            continue

        out_fq = out_dir / f"{sample_label}_singletons_rescued.fastq.gz"
        out_tsv = out_dir / f"{sample_label}_singleton_variants_{min_umis}umis.tsv"

        logger.info(f"--- Processing {sample_label} ---")
        logger.info(f"Input : {in_fq}")
        logger.info(f"Outfq : {out_fq}")
        logger.info(f"Outtsv: {out_tsv}")

        try:
            stats = rescue_one_sample(in_fq, out_fq, out_tsv, wt_anchor, var_len, min_umis)
            logger.info(f"Done {sample_label}: {stats}")
            all_stats.append({"sample_id": sample_id, "sample_name": sample_name, **stats})
            write_sample_metrics(
                metrics_path=out_dir / "per_sample_metrics.tsv",
                sample_id=sample_id,
                sample_name=sample_name,
                step="02a_singleton_rescue",
                metrics={
                    "reads_in":              stats["total_singleton_reads"],
                    "reads_lost_no_anchor":  stats["no_anchor_reads"],
                    "reads_lost_too_short":  stats["too_short_after_anchor_reads"],
                    "reads_extracted":       stats["extracted_reads"],
                    "variants_total":        stats["n_variants_total"],
                    "variants_kept":         stats["n_variants_kept_ge_min_umis"],
                    "reads_rescued":         stats["kept_reads"],
                },
            )
        except Exception as e:
            logger.error(f"ERROR processing {sample_label}: {e}")
            logger.error(traceback.format_exc())
            continue

    # write combined summary
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "singleton_rescue_summary.tsv"
    if all_stats:
        pd.DataFrame(all_stats).to_csv(summary_path, sep="\t", index=False)
        logger.info(f"Wrote summary TSV: {summary_path}")
    else:
        logger.warning("No per-sample stats collected (no samples processed successfully).")

    logger.info(f"Finished. Missing singleton files: {n_missing} / {len(df)}")


if __name__ == "__main__":
    main()
