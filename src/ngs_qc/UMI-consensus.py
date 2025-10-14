#!/usr/bin/env python3
import gzip
from Bio import SeqIO
from collections import defaultdict, Counter
from pathlib import Path
from fuzzysearch import find_near_matches
import yaml
import sys
import traceback
from typing import Dict, List, Tuple

from logger import Logger   # <-- your custom logger


def parse_config(config_file):
    """Load YAML config into dictionary and validate required fields."""
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    required = ["sample_name", "umi_length", "anchor_sequence"]
    for key in required:
        if key not in config:
            raise ValueError(f"Missing required config parameter: {key}")

    return config


def fuzzy_find_anchor(seq, anchor, max_dist=1):
    """Find approximate match of anchor in a sequence."""
    matches = find_near_matches(anchor, seq, max_l_dist=max_dist)
    return matches[0].start if matches else -1


def collapse_umis(reads_per_umi):
    """Return consensus sequence using majority rule."""
    if not reads_per_umi:
        return ""
    read_length = len(reads_per_umi[0])
    consensus_seq = []
    for i in range(read_length):
        bases = [read[i] for read in reads_per_umi if len(read) > i]
        if not bases:
            consensus_seq.append("N")
            continue
        most_common_base, _ = Counter(bases).most_common(1)[0]
        consensus_seq.append(most_common_base)
    return "".join(consensus_seq)

def soft_consensus_per_umi(umi_groups: Dict[str, Counter],
                           top_k: int = 3,
                           cumulative_frac: float = 0.9,
                           min_count_for_report: int = 1
                          ) -> Dict[str, List[Tuple[str,int]]]:
    """
    For each UMI, return a list (ordered) of (sequence, count) that are kept.
    Works on umi_groups = {umi: [seq1, seq2, ...]}.
    """
    out = {}
    for umi, seqs in umi_groups.items():
        counts = Counter(seqs)
        # filter by count
        filtered = {s: c for s, c in counts.items() if c >= min_count_for_report}
        if not filtered:
            out[umi] = []
            continue
        items = sorted(filtered.items(), key=lambda x: x[1], reverse=True)
        total = sum(c for _, c in items)
        cum = 0
        keep = []
        for idx, (s, cnt) in enumerate(items):
            keep.append((s, cnt))
            cum += cnt
            if cumulative_frac and (cum / total) >= cumulative_frac:
                break
            if top_k and len(keep) >= top_k:
                break
        out[umi] = keep
    return out

def process_sample(config_file):
    config = parse_config(config_file)
    base_sample_name = config["sample_name"]  # e.g. GFB-59053

    # Define directories
    base_dir = Path(f"/cluster/project/reddy/marluca/NGS_pipeline/data/processed/{config['pool']}")
    sample_dir = base_dir / base_sample_name
    matches = list(sample_dir.glob("*_combined_final.fastq.gz"))
    if not matches:
        raise FileNotFoundError(f"No *_combined_final.fastq.gz file found in {sample_dir}")
    input_fastq = matches[0]

    # Extract descriptive prefix from filename
    descriptive_prefix = input_fastq.stem.replace("_combined_final.fastq", "")
    sample_name = f"{base_sample_name}_{descriptive_prefix}"

    output_dir = sample_dir / "umi"
    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir = Path("/cluster/project/reddy/marluca/NGS_pipeline/logs")
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / f"{sample_name}_umi.log"
    logger = Logger(str(log_dir), log_file.name)

    logger.log(f"Processing sample {sample_name}")
    logger.log(f"Config file: {config_file}")
    logger.log(f"Input FASTQ: {input_fastq}")
    logger.log(f"Output dir: {output_dir}")

    # Parameters
    umi_len = int(config["umi_length"])
    anchor_seq = config["anchor_sequence"]
    max_anchor_dist = int(config.get("max_mismatches", 1))
    min_reads_per_umi = int(config.get("min_reads_per_umi", 2))

    # soft collapse parameters
    top_k = int(config.get("soft_top_k", 3))
    cumulative_frac = float(config.get("soft_cumulative_frac", 0.9))
    min_count_for_report = int(config.get("soft_min_count", 2))

    umi_groups = defaultdict(list)
    total_reads, discarded_short_umi, no_anchor_found = 0, 0, 0

    try:
        with gzip.open(input_fastq, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                total_reads += 1
                seq = str(record.seq)

                # Check for anchor downstream of UMI
                anchor_pos = fuzzy_find_anchor(seq[umi_len:], anchor_seq, max_dist=max_anchor_dist)
                if anchor_pos == -1:
                    no_anchor_found += 1
                    continue

                umi = seq[:umi_len]
                extracted_seq = seq[umi_len:]

                if len(umi) < umi_len:
                    discarded_short_umi += 1
                    continue

                umi_groups[umi].append(extracted_seq)
    except Exception as e:
        logger.log(f"ERROR while reading FASTQ: {e}")
        logger.log(traceback.format_exc())
        raise

    # ---- HARD CONSENSUS ----
    hard_out = output_dir / f"{sample_name}_consensus_hard.fastq.gz"
    try:
        with gzip.open(hard_out, "wt") as out_handle:
            for umi, seqs in umi_groups.items():
                if len(seqs) < min_reads_per_umi:
                    continue
                consensus = collapse_umis(seqs)
                out_handle.write(f"@{umi}\n{consensus}\n+\n{'I'*len(consensus)}\n")
        logger.log(f"Hard consensus FASTQ written to {hard_out}")
    except Exception as e:
        logger.log(f"ERROR writing hard consensus: {e}")
        logger.log(traceback.format_exc())
        raise

    # ---- SOFT CONSENSUS ----
    # Convert umi_groups -> Counter per UMI
    umi_counters = {umi: Counter(seqs) for umi, seqs in umi_groups.items()}
    soft_results = soft_consensus_per_umi(
        umi_counters,
        top_k=top_k,
        cumulative_frac=cumulative_frac,
        min_count_for_report=min_count_for_report
    )

    soft_out = output_dir / f"{sample_name}_consensus_soft.fastq.gz"
    try:
        with gzip.open(soft_out, "wt") as out_handle:
            for umi, kept_list in soft_results.items():
                for idx, (seq, count) in enumerate(kept_list, start=1):
                    umi_id = f"{umi}_var{idx}_count{count}"
                    out_handle.write(f"@{umi_id}\n{seq}\n+\n{'I'*len(seq)}\n")
        logger.log(f"Soft consensus FASTQ written to {soft_out}")
    except Exception as e:
        logger.log(f"ERROR writing soft consensus: {e}")
        logger.log(traceback.format_exc())
        raise

    # ---- SUMMARY ----
    summary_file = output_dir / f"{sample_name}_umi_summary.txt"
    try:
        with open(summary_file, "w") as s:
            s.write(f"Sample: {sample_name}\n")
            s.write(f"Total reads: {total_reads}\n")
            s.write(f"Total UMIs detected: {len(umi_groups)}\n")
            s.write(f"Discarded short UMIs: {discarded_short_umi}\n")
            s.write(f"Reads with no anchor: {no_anchor_found}\n")
            s.write("Reads per UMI distribution:\n")
            counts = [len(v) for v in umi_groups.values()]
            for c in sorted(set(counts)):
                s.write(f"  {c} reads: {counts.count(c)} UMIs\n")
        logger.log(f"Summary written to {summary_file}")
    except Exception as e:
        logger.log(f"ERROR writing summary: {e}")
        logger.log(traceback.format_exc())
        raise

    logger.log("UMI processing complete.")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="UMI consensus script")
    parser.add_argument("--config", required=True, help="Path to YAML config")
    args = parser.parse_args()
    process_sample(args.config)