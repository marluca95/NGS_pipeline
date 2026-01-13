#!/usr/bin/env python3
import gzip
import pandas as pd
from Bio import SeqIO
from collections import defaultdict, Counter
from pathlib import Path
from fuzzysearch import find_near_matches
import yaml
import argparse
import traceback
from typing import Dict, List, Tuple

from logger import Logger  # your custom logger


def parse_yaml(yaml_file: str) -> dict:
    with open(yaml_file, "r") as f:
        cfg = yaml.safe_load(f) or {}

    required = ["sample_sheet", "output_dir", "umi_length", "anchor_sequence"]
    missing = [k for k in required if k not in cfg]
    if missing:
        raise ValueError(f"Missing required YAML keys: {missing}")

    # default logs_dir if not provided
    cfg.setdefault("logs_dir", str(Path(cfg["output_dir"]) / "_logs"))

    return cfg



def load_sample_sheet(sample_sheet_path: str) -> pd.DataFrame:
    try:
        df = pd.read_csv(sample_sheet_path, sep="\t")
        if df.shape[1] == 1:
            df = pd.read_csv(sample_sheet_path)
    except Exception:
        df = pd.read_csv(sample_sheet_path)

    required = {"sample_id", "sample_name"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Sample sheet missing required columns: {missing}. Found: {list(df.columns)}")
    return df


def safe_name(s: str) -> str:
    return str(s).strip().replace(" ", "_").replace("/", "_")


def fuzzy_find_anchor(seq: str, anchor: str, max_dist: int = 1) -> int:
    matches = find_near_matches(anchor, seq, max_l_dist=max_dist)
    return matches[0].start if matches else -1


def collapse_umis(reads_per_umi: List[str]) -> str:
    if not reads_per_umi:
        return ""
    read_length = max(len(r) for r in reads_per_umi)
    consensus = []
    for i in range(read_length):
        bases = [r[i] for r in reads_per_umi if len(r) > i]
        if not bases:
            consensus.append("N")
            continue
        consensus.append(Counter(bases).most_common(1)[0][0])
    return "".join(consensus)


def soft_consensus_per_umi(
    umi_groups: Dict[str, List[str]],
    top_k: int = 3,
    cumulative_frac: float = 0.9,
    min_count_for_report: int = 2,
) -> Dict[str, List[Tuple[str, int]]]:
    out = {}
    for umi, seqs in umi_groups.items():
        counts = Counter(seqs)
        items = [(s, c) for s, c in counts.items() if c >= min_count_for_report]
        items.sort(key=lambda x: x[1], reverse=True)
        if not items:
            out[umi] = []
            continue
        total = sum(c for _, c in items)
        cum = 0
        keep = []
        for s, c in items:
            keep.append((s, c))
            cum += c
            if cumulative_frac and (cum / total) >= cumulative_frac:
                break
            if top_k and len(keep) >= top_k:
                break
        out[umi] = keep
    return out


def process_one_sample(sample_id: str, sample_name: str, cfg: dict, logs_dir: Path):
    sample_label = f"{sample_id}_{safe_name(sample_name)}"
    out_base = Path(cfg["output_dir"])
    sample_dir = out_base / sample_label

    combined_suffix = cfg.get("combined_suffix", "_combined_trimmed.fastq.gz")
    matches = list(sample_dir.glob(f"*{combined_suffix}"))
    if not matches:
        raise FileNotFoundError(f"No *{combined_suffix} found in {sample_dir}")
    input_fastq = matches[0]

    umi_dir = sample_dir / "umi"
    umi_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    logger = Logger(str(logs_dir), f"{sample_label}_umi.log")
    logger.log(f"Processing sample: {sample_label}")
    logger.log(f"Input FASTQ: {input_fastq}")
    logger.log(f"UMI output dir: {umi_dir}")

    umi_len = int(cfg["umi_length"])
    anchor_seq = cfg["anchor_sequence"]
    max_anchor_dist = int(cfg.get("max_mismatches", 1))
    min_reads_per_umi = int(cfg.get("min_reads_per_umi", 2))

    top_k = int(cfg.get("soft_top_k", 3))
    cumulative_frac = float(cfg.get("soft_cumulative_frac", 0.9))
    min_count_for_report = int(cfg.get("soft_min_count", 2))

    umi_groups = defaultdict(list)
    total_reads, no_anchor_found, short_umi = 0, 0, 0

    try:
        with gzip.open(input_fastq, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                total_reads += 1
                seq = str(record.seq)

                if len(seq) < umi_len:
                    short_umi += 1
                    continue

                # anchor must exist somewhere after the UMI
                anchor_pos = fuzzy_find_anchor(seq[umi_len:], anchor_seq, max_dist=max_anchor_dist)
                if anchor_pos == -1:
                    no_anchor_found += 1
                    continue

                umi = seq[:umi_len]

                # Current behavior (matches your old script):
                extracted_seq = seq[umi_len:]

                umi_groups[umi].append(extracted_seq)

    except Exception as e:
        logger.log(f"ERROR reading FASTQ: {e}")
        logger.log(traceback.format_exc())
        raise

    # HARD consensus
    hard_out = umi_dir / f"{sample_label}_consensus_hard.fastq.gz"
    with gzip.open(hard_out, "wt") as out_handle:
        for umi, seqs in umi_groups.items():
            if len(seqs) < min_reads_per_umi:
                continue
            cons = collapse_umis(seqs)
            out_handle.write(f"@{umi}\n{cons}\n+\n{'I'*len(cons)}\n")
    logger.log(f"Hard consensus written: {hard_out}")

    # SOFT consensus variants
    soft_out = umi_dir / f"{sample_label}_consensus_soft.fastq.gz"
    soft_results = soft_consensus_per_umi(
        umi_groups,
        top_k=top_k,
        cumulative_frac=cumulative_frac,
        min_count_for_report=min_count_for_report,
    )
    with gzip.open(soft_out, "wt") as out_handle:
        for umi, kept in soft_results.items():
            for idx, (seq, count) in enumerate(kept, start=1):
                out_handle.write(f"@{umi}_var{idx}_count{count}\n{seq}\n+\n{'I'*len(seq)}\n")
    logger.log(f"Soft consensus written: {soft_out}")

    # Summary
    summary = umi_dir / f"{sample_label}_umi_summary.txt"
    counts = [len(v) for v in umi_groups.values()]
    with open(summary, "w") as s:
        s.write(f"Sample: {sample_label}\n")
        s.write(f"Total reads: {total_reads}\n")
        s.write(f"Total UMIs detected: {len(umi_groups)}\n")
        s.write(f"Short reads (<UMI length): {short_umi}\n")
        s.write(f"Reads with no anchor: {no_anchor_found}\n")
        s.write("Reads per UMI distribution:\n")
        for c in sorted(set(counts)):
            s.write(f"  {c} reads: {counts.count(c)} UMIs\n")
    logger.log(f"Summary written: {summary}")
    logger.log("UMI processing complete.")


def main():
    ap = argparse.ArgumentParser(description="UMI consensus")
    ap.add_argument("--yaml_config", required=True)
    ap.add_argument("--sample_id", default=None, help="Process only one sample_id")
    args = ap.parse_args()

    cfg = parse_yaml(args.yaml_config)
    df = load_sample_sheet(cfg["sample_sheet"])

    if args.sample_id is not None:
        df = df[df["sample_id"].astype(str) == str(args.sample_id)]
        if df.empty:
            raise ValueError(f"No rows for sample_id={args.sample_id}")

    # Ensure one sample_name per sample_id
    grouped = df.groupby("sample_id", sort=False)

    logs_dir = Path(cfg["logs_dir"])

    for sample_id, g in grouped:
        names = g["sample_name"].dropna().unique()
        if len(names) != 1:
            raise ValueError(f"sample_id={sample_id} has multiple sample_name values: {names}")
        process_one_sample(str(sample_id), str(names[0]), cfg, logs_dir)


if __name__ == "__main__":
    main()
