"""
Paired-Chain Pipeline — Step 02a: Paired Singleton Rescue
==========================================================

Rescues singleton UMI pairs from script2's output.

Strategy:
  1. Load singleton FASTQ files (_singletons_beta.fastq.gz and _singletons_alpha.fastq.gz).
     Headers have format: @{UMI_beta}|{UMI_alpha}|n=1 — beta and alpha files are
     already in lockstep (same header in both files at the same position).
  2. For each singleton pair, extract the variable region from each chain using
     the AA-extraction anchor (not the UMI-extraction anchor).
  3. Build a dict: (beta_var_nt, alpha_var_nt) → set of UMI keys (UMI_beta|UMI_alpha).
  4. Keep only (beta_var, alpha_var) pairs that appear across >= min_umis distinct
     UMI pairs — these are likely real variants seen multiple times on different molecules.
  5. For each rescued pair: pick one representative read (the first encountered) and
     write it to rescued beta and alpha output FASTQ files, still in lockstep.

The rescued outputs use the same header format as consensus reads:
  @{UMI_beta}|{UMI_alpha}|n={n_umis}
but with n= reflecting the count of distinct UMI pairs supporting the variant.

Outputs per sample:
  {sample_label}_singletons_rescued_beta.fastq.gz
  {sample_label}_singletons_rescued_alpha.fastq.gz
  {sample_label}_singleton_rescue_summary.txt
"""

import argparse
import gzip
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Set, Tuple

from Bio import SeqIO
from fuzzysearch import find_near_matches
from utils.config_utils import load_pipeline_config
from utils.logging_utils import setup_pipeline_logging
from utils.metrics_utils import write_sample_metrics
from utils.sample_utils import load_sample_sheet, safe_name

SCRIPT_NAME = "script2a_paired_singleton_rescue"
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
def parse_yaml(yaml_file: str) -> dict:
    cfg = load_pipeline_config(
        yaml_file,
        required_keys=(
            "samples_tsv", "singleton_dir", "output_dir",
            "anchor_beta", "anchor_alpha",
            "var_len_beta", "var_len_alpha",
        ),
        default_values={
            "logs_dir": None,
            "min_umis": 3,
            "singleton_suffix_beta": "_singletons_beta.fastq.gz",
            "singleton_suffix_alpha": "_singletons_alpha.fastq.gz",
        },
        path_keys=("samples_tsv", "singleton_dir", "output_dir", "logs_dir"),
    )
    if not cfg.get("logs_dir"):
        cfg["logs_dir"] = str(Path(cfg["output_dir"]) / "_logs")
    return cfg


# ---------------------------------------------------------------------------
# Anchor search (exact or fuzzy) for variable region extraction
# ---------------------------------------------------------------------------
def find_var_region(seq: str, anchor: str, var_len: int) -> Optional[str]:
    """
    Find anchor in seq (exact match preferred, no fuzzy).
    Returns the var_len nucleotides immediately following the anchor, or None.
    """
    seq = seq.upper()
    anchor = anchor.upper()
    idx = seq.find(anchor)
    if idx == -1:
        return None
    start = idx + len(anchor)
    region = seq[start: start + var_len]
    if len(region) != var_len:
        return None
    return region


# ---------------------------------------------------------------------------
# Per-sample processing
# ---------------------------------------------------------------------------
def process_one_sample(sample_id: str, sample_name: str, cfg: dict) -> Dict[str, int]:
    sample_label = f"{sample_id}_{sample_name}"
    sing_dir = Path(cfg["singleton_dir"])
    out_dir = Path(cfg["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    suf_beta  = cfg.get("singleton_suffix_beta",  "_singletons_beta.fastq.gz")
    suf_alpha = cfg.get("singleton_suffix_alpha", "_singletons_alpha.fastq.gz")

    sing_beta_path  = sing_dir / f"{sample_label}{suf_beta}"
    sing_alpha_path = sing_dir / f"{sample_label}{suf_alpha}"

    for p in (sing_beta_path, sing_alpha_path):
        if not p.exists():
            raise FileNotFoundError(f"Singleton file not found: {p}")

    anchor_beta  = str(cfg["anchor_beta"])
    anchor_alpha = str(cfg["anchor_alpha"])
    var_len_beta  = int(cfg["var_len_beta"])
    var_len_alpha = int(cfg["var_len_alpha"])
    min_umis = int(cfg.get("min_umis", 3))

    logger.info("Sample: %s", sample_label)
    logger.info("Singleton beta:  %s", sing_beta_path)
    logger.info("Singleton alpha: %s", sing_alpha_path)

    reads_in = 0
    lost_beta_no_anchor = 0
    lost_alpha_no_anchor = 0

    # variant_pair → set of UMI keys
    # also store representative inserts for writing output
    variant_to_umis: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    variant_to_repr: Dict[Tuple[str, str], Tuple[str, str, str, str]] = {}
    # repr tuple: (umi_key, beta_insert, alpha_insert, alpha_insert_qual)

    with gzip.open(sing_beta_path, "rt") as bh, gzip.open(sing_alpha_path, "rt") as ah:
        for beta_rec, alpha_rec in zip(SeqIO.parse(bh, "fastq"), SeqIO.parse(ah, "fastq")):
            reads_in += 1

            # Parse UMI key from header (same in both files)
            # Header format: {UMI_beta}|{UMI_alpha}|n=1
            umi_key = beta_rec.id  # BioPython strips after space

            beta_seq  = str(beta_rec.seq)
            alpha_seq = str(alpha_rec.seq)

            beta_var = find_var_region(beta_seq, anchor_beta, var_len_beta)
            if beta_var is None:
                lost_beta_no_anchor += 1
                continue

            alpha_var = find_var_region(alpha_seq, anchor_alpha, var_len_alpha)
            if alpha_var is None:
                lost_alpha_no_anchor += 1
                continue

            pair_key = (beta_var, alpha_var)
            variant_to_umis[pair_key].add(umi_key)

            if pair_key not in variant_to_repr:
                # Store full insert sequences as representative
                variant_to_repr[pair_key] = (
                    umi_key,
                    beta_seq,
                    str(beta_rec.letter_annotations.get("phred_quality", [])),
                    alpha_seq,
                    str(alpha_rec.letter_annotations.get("phred_quality", [])),
                )

    variants_total = len(variant_to_umis)
    variants_kept = sum(1 for umis in variant_to_umis.values() if len(umis) >= min_umis)
    reads_rescued = 0

    rescued_beta_path  = out_dir / f"{sample_label}_singletons_rescued_beta.fastq.gz"
    rescued_alpha_path = out_dir / f"{sample_label}_singletons_rescued_alpha.fastq.gz"

    with gzip.open(rescued_beta_path, "wt") as rb_h, gzip.open(rescued_alpha_path, "wt") as ra_h:
        for pair_key, umiset in variant_to_umis.items():
            if len(umiset) < min_umis:
                continue

            umi_key, beta_seq, beta_qual_str, alpha_seq, alpha_qual_str = variant_to_repr[pair_key]
            # Build a new header reflecting the rescue count
            n_supporting = len(umiset)
            # Use representative UMI key for traceability, with n=count
            header = f"{umi_key.rsplit('|n=', 1)[0]}|n={n_supporting}"

            rb_h.write(f"@{header}\n{beta_seq}\n+\n{'I' * len(beta_seq)}\n")
            ra_h.write(f"@{header}\n{alpha_seq}\n+\n{'I' * len(alpha_seq)}\n")
            reads_rescued += 1

    # Write summary
    summary_path = out_dir / f"{sample_label}_singleton_rescue_summary.txt"
    with open(summary_path, "w") as s:
        s.write(f"Sample: {sample_label}\n")
        s.write(f"min_umis threshold: {min_umis}\n\n")
        s.write(f"Singleton pairs read in:               {reads_in}\n")
        s.write(f"Pairs lost (beta anchor not found):    {lost_beta_no_anchor}\n")
        s.write(f"Pairs lost (alpha anchor not found):   {lost_alpha_no_anchor}\n")
        s.write(f"Distinct (beta_var, alpha_var) pairs:  {variants_total}\n")
        s.write(f"Pairs with >= {min_umis} distinct UMIs:        {variants_kept}\n")
        s.write(f"Reads written to rescued output:       {reads_rescued}\n")

    logger.info("Rescued beta:  %s (%d pairs)", rescued_beta_path, reads_rescued)
    logger.info("Rescued alpha: %s", rescued_alpha_path)
    logger.info("Summary: %s", summary_path)

    return {
        "reads_in":                   reads_in,
        "reads_lost_beta_no_anchor":  lost_beta_no_anchor,
        "reads_lost_alpha_no_anchor": lost_alpha_no_anchor,
        "variants_total":             variants_total,
        "variants_kept":              variants_kept,
        "reads_rescued":              reads_rescued,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="Paired-Chain Pipeline Step 02a: Paired singleton rescue."
    )
    ap.add_argument("--yaml_config", required=True)
    args = ap.parse_args()

    cfg = parse_yaml(args.yaml_config)

    log_file = setup_pipeline_logging(
        logs_dir=cfg["logs_dir"],
        script_name=SCRIPT_NAME,
        scope="run",
        run_label=cfg.get("run_label"),
    )
    logger.info("Logging to: %s", log_file)

    df = load_sample_sheet(cfg["samples_tsv"], required_cols=("sample_id", "sample_name"))
    metrics_path = Path(cfg["output_dir"]) / "per_sample_metrics.tsv"

    for sample_id, g in df.groupby("sample_id", sort=False):
        names = g["sample_name"].dropna().unique()
        if len(names) != 1:
            raise ValueError(f"sample_id={sample_id} has multiple sample_name values: {names}")
        sample_name = safe_name(str(names[0]))
        stats = process_one_sample(str(sample_id), sample_name, cfg)
        write_sample_metrics(
            metrics_path=metrics_path,
            sample_id=str(sample_id),
            sample_name=sample_name,
            step="02a_singleton_rescue",
            metrics=stats,
        )


if __name__ == "__main__":
    main()
