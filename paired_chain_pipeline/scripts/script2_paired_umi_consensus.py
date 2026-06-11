"""
Paired-Chain Pipeline — Step 02: Paired-End UMI Consensus
==========================================================

Processes R1 (beta chain) and R2 (alpha chain) reads simultaneously.

Read structure expected:
  R1 (beta):  [spacer 0–7 nt] [UMI 12 nt] [anchor_beta] [beta_insert]
  R2 (alpha): [spacer 0–7 nt] [UMI 12 nt] [anchor_alpha] [alpha_insert]

Algorithm:
  1. Read combined_trimmed_R1 and combined_trimmed_R2 files in LOCKSTEP.
     Illumina paired FASTQ files are always in the same order; read IDs are
     verified to match at each position as a sanity check.
  2. Extract UMI_beta from R1, UMI_alpha from R2, and their respective inserts.
  3. Both chains must pass anchor search; pairs where either fails are dropped.
  4. Group pairs by compound key "{UMI_beta}_{UMI_alpha}" using a bucketed
     temp-file approach (same bucketing strategy as script2_umi_consensus.py
     to keep memory usage bounded).
  5. For UMI groups with >= min_reads_per_umi reads:
     - Call posterior quality-aware consensus on R1 inserts → cons_beta
     - Call posterior quality-aware consensus on R2 inserts → cons_alpha
  6. Write TWO PARALLEL FASTQ files with IDENTICAL headers so the paired
     relationship is preserved by position:
       Header: @{UMI_beta}|{UMI_alpha}|n={n}
  7. Singletons (n=1) are written to separate beta/alpha FASTQ files with the
     same cross-referencing header format.

Outputs per sample:
  {sample_label}_consensus_beta.fastq.gz
  {sample_label}_consensus_alpha.fastq.gz
  {sample_label}_singletons_beta.fastq.gz
  {sample_label}_singletons_alpha.fastq.gz
  {sample_label}_paired_umi_summary.txt
  {sample_label}_consensus_qc.tsv
"""

import argparse
import gzip
import logging
import math
import shutil
import traceback
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from Bio import SeqIO
from fuzzysearch import find_near_matches
from utils.config_utils import load_pipeline_config
from utils.logging_utils import setup_pipeline_logging
from utils.metrics_utils import write_sample_metrics
from utils.sample_utils import load_sample_sheet, safe_name

SCRIPT_NAME = "script2_paired_umi_consensus"
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
def parse_yaml(yaml_file: str) -> dict:
    cfg = load_pipeline_config(
        yaml_file,
        required_keys=("sample_sheet", "input_dir", "output_dir",
                       "umi_length", "anchor_beta", "anchor_alpha"),
        default_values={
            "logs_dir": None,
            "combined_suffix_r1": "_combined_trimmed_R1.fastq.gz",
            "combined_suffix_r2": "_combined_trimmed_R2.fastq.gz",
            "max_mismatches": 1,
            "min_reads_per_umi": 2,
            "umi_bucket_count": 64,
            "spacer_min": 0,
            "spacer_max": 7,
            "posterior_min_log_delta": 2.0,
            "consensus_placeholder_qual_char": "I",
        },
        path_keys=("sample_sheet", "input_dir", "output_dir", "logs_dir"),
    )
    if not cfg.get("logs_dir"):
        cfg["logs_dir"] = str(Path(cfg["output_dir"]) / "_logs")
    return cfg


# ---------------------------------------------------------------------------
# Quality / consensus helpers (adapted from script2_umi_consensus.py)
# ---------------------------------------------------------------------------
def phred_to_fastq_ascii(phred: List[int]) -> str:
    return "".join(chr(min(93, max(0, q)) + 33) for q in phred)


def placeholder_qual(length: int, char: str = "I") -> str:
    return char * max(0, int(length))


def phred_to_error_prob(q: int) -> float:
    return 10 ** (-q / 10.0)


def posterior_base_call(
    bases: List[str], quals: List[int], min_log_delta: float = 2.0,
) -> Tuple[str, float]:
    obs = [(b, q) for b, q in zip(bases, quals) if b in "ACGT"]
    if not obs:
        return "N", 0.0
    logL: Dict[str, float] = {}
    for candidate in "ACGT":
        ll = 0.0
        for b, q in obs:
            p_err = phred_to_error_prob(int(q))
            p_err = min(max(p_err, 1e-12), 1 - 1e-12)
            ll += math.log(1.0 - p_err) if b == candidate else math.log(p_err / 3.0)
        logL[candidate] = ll
    vals_sorted = sorted(logL.items(), key=lambda kv: kv[1], reverse=True)
    best_base, best_ll = vals_sorted[0]
    second_ll = vals_sorted[1][1]
    log_delta = best_ll - second_ll
    if log_delta < min_log_delta:
        return "N", log_delta
    return best_base, log_delta


def collapse_umis_posterior(
    seqs: List[str], quals: List[List[int]], min_log_delta: float = 2.0,
) -> Tuple[str, Dict[str, float]]:
    if not seqs:
        return "", {"len": 0.0, "nN": 0.0, "min_log_delta": 0.0,
                    "median_log_delta": 0.0, "p10_log_delta": 0.0}
    read_length = max(len(s) for s in seqs)
    consensus_bases: List[str] = []
    log_deltas: List[float] = []
    nN = 0
    for i in range(read_length):
        bases_i = [s[i] for s, _ in zip(seqs, quals) if i < len(s)]
        quals_i = [q[i] if i < len(q) else 0 for s, q in zip(seqs, quals) if i < len(s)]
        if not bases_i:
            consensus_bases.append("N"); nN += 1; continue
        call, log_delta = posterior_base_call(bases_i, quals_i, min_log_delta=min_log_delta)
        consensus_bases.append(call)
        log_deltas.append(float(log_delta))
        if call == "N":
            nN += 1
    cons = "".join(consensus_bases)

    def _median(xs: List[float]) -> float:
        if not xs: return 0.0
        ys = sorted(xs); m = len(ys) // 2
        return ys[m] if len(ys) % 2 == 1 else 0.5 * (ys[m - 1] + ys[m])

    def _p10(xs: List[float]) -> float:
        if not xs: return 0.0
        ys = sorted(xs); idx = max(0, int(0.10 * (len(ys) - 1)))
        return ys[idx]

    qc = {
        "len": float(len(cons)), "nN": float(nN),
        "min_log_delta": float(min(log_deltas)) if log_deltas else 0.0,
        "median_log_delta": float(_median(log_deltas)),
        "p10_log_delta": float(_p10(log_deltas)),
    }
    return cons, qc


# ---------------------------------------------------------------------------
# Anchor search (adapted from script2_umi_consensus.py)
# ---------------------------------------------------------------------------
def find_anchor_start_valid(
    seq: str, anchor: str, umi_len: int, max_dist: int,
    spacer_min: int = 0, spacer_max: int = 7,
) -> Tuple[Optional[Tuple[int, int]], bool]:
    matches = find_near_matches(anchor, seq, max_l_dist=max_dist)
    if not matches:
        return None, False
    lo = umi_len + spacer_min
    hi = umi_len + spacer_max
    for m in sorted(matches, key=lambda x: x.start):
        if lo <= m.start <= hi:
            return (m.start, m.end), True
    return None, True


# ---------------------------------------------------------------------------
# Main per-sample processing
# ---------------------------------------------------------------------------
def process_one_sample(
    sample_id: str, sample_name: str, cfg: dict,
) -> Dict[str, Union[int, float]]:
    sample_label = f"{sample_id}_{sample_name}"
    in_base = Path(cfg["input_dir"])
    sample_in_dir = in_base / sample_label
    out_dir = Path(cfg["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    suf_r1 = cfg.get("combined_suffix_r1", "_combined_trimmed_R1.fastq.gz")
    suf_r2 = cfg.get("combined_suffix_r2", "_combined_trimmed_R2.fastq.gz")

    r1_fastq = sample_in_dir / f"{sample_label}{suf_r1}"
    r2_fastq = sample_in_dir / f"{sample_label}{suf_r2}"

    for p in (r1_fastq, r2_fastq):
        if not p.exists():
            raise FileNotFoundError(f"Expected input not found: {p}")

    logger.info("Sample: %s", sample_label)
    logger.info("R1 (beta):  %s", r1_fastq)
    logger.info("R2 (alpha): %s", r2_fastq)

    umi_len = int(cfg["umi_length"])
    anchor_beta = str(cfg["anchor_beta"])
    anchor_alpha = str(cfg["anchor_alpha"])
    max_dist = int(cfg.get("max_mismatches", 1))
    min_reads_per_umi = int(cfg.get("min_reads_per_umi", 2))
    umi_bucket_count = int(cfg.get("umi_bucket_count", 64))
    spacer_min = int(cfg.get("spacer_min", 0))
    spacer_max = int(cfg.get("spacer_max", 7))
    posterior_min_log_delta = float(cfg.get("posterior_min_log_delta", 2.0))
    placeholder_char = str(cfg.get("consensus_placeholder_qual_char", "I"))
    if len(placeholder_char) != 1:
        raise ValueError("consensus_placeholder_qual_char must be a single ASCII character.")

    min_prefix = spacer_max + umi_len + max(len(anchor_beta), len(anchor_alpha))

    # Counters
    total_pairs = 0
    r1_too_short = 0; r1_no_anchor = 0; r1_anchor_wrong_pos = 0
    r2_too_short = 0; r2_no_anchor = 0; r2_anchor_wrong_pos = 0
    r1_empty_insert = 0; r2_empty_insert = 0
    id_mismatch_count = 0

    # Bucketed temp files: each bucket stores tab-separated lines
    # format: {umi_key}\t{beta_insert}\t{beta_q_hex}\t{alpha_insert}\t{alpha_q_hex}
    bucket_temp_dir = out_dir / f".{sample_label}_umi_buckets"
    bucket_handles: Dict[int, object] = {}
    shutil.rmtree(bucket_temp_dir, ignore_errors=True)

    def _bucket_path(bid: int) -> Path:
        return bucket_temp_dir / f"bucket_{bid:04d}.tsv.gz"

    def _bucket_id(umi_key: str) -> int:
        return hash(umi_key) % umi_bucket_count

    def _get_bucket(bid: int):
        if bid not in bucket_handles:
            bucket_temp_dir.mkdir(parents=True, exist_ok=True)
            bucket_handles[bid] = gzip.open(_bucket_path(bid), "wt")
        return bucket_handles[bid]

    try:
        # --- Pass 1: parse paired reads, extract UMIs, write to buckets ---
        try:
            with gzip.open(r1_fastq, "rt") as r1_h, gzip.open(r2_fastq, "rt") as r2_h:
                for r1_rec, r2_rec in zip(SeqIO.parse(r1_h, "fastq"), SeqIO.parse(r2_h, "fastq")):
                    total_pairs += 1

                    # Verify pairing by read ID
                    if r1_rec.id != r2_rec.id:
                        id_mismatch_count += 1
                        if id_mismatch_count <= 5:
                            logger.warning(
                                "Read ID mismatch at pair %d: R1=%s R2=%s",
                                total_pairs, r1_rec.id, r2_rec.id,
                            )

                    seq_r1 = str(r1_rec.seq)
                    seq_r2 = str(r2_rec.seq)
                    qual_r1 = r1_rec.letter_annotations.get("phred_quality", [])
                    qual_r2 = r2_rec.letter_annotations.get("phred_quality", [])

                    # --- R1 (beta) anchor search ---
                    if len(seq_r1) < min_prefix:
                        r1_too_short += 1; continue
                    hit_r1, has_r1 = find_anchor_start_valid(
                        seq_r1, anchor_beta, umi_len, max_dist, spacer_min, spacer_max,
                    )
                    if hit_r1 is None:
                        if has_r1: r1_anchor_wrong_pos += 1
                        else: r1_no_anchor += 1
                        continue

                    # --- R2 (alpha) anchor search ---
                    if len(seq_r2) < min_prefix:
                        r2_too_short += 1; continue
                    hit_r2, has_r2 = find_anchor_start_valid(
                        seq_r2, anchor_alpha, umi_len, max_dist, spacer_min, spacer_max,
                    )
                    if hit_r2 is None:
                        if has_r2: r2_anchor_wrong_pos += 1
                        else: r2_no_anchor += 1
                        continue

                    a_start_r1, a_end_r1 = hit_r1
                    a_start_r2, a_end_r2 = hit_r2

                    if a_start_r1 < umi_len or a_start_r2 < umi_len:
                        r1_too_short += 1; continue

                    umi_beta  = seq_r1[a_start_r1 - umi_len: a_start_r1]
                    umi_alpha = seq_r2[a_start_r2 - umi_len: a_start_r2]

                    insert_beta  = seq_r1[a_end_r1:]
                    insert_alpha = seq_r2[a_end_r2:]

                    if not insert_beta:
                        r1_empty_insert += 1; continue
                    if not insert_alpha:
                        r2_empty_insert += 1; continue

                    iq_beta  = qual_r1[a_end_r1:][:len(insert_beta)]
                    iq_alpha = qual_r2[a_end_r2:][:len(insert_alpha)]

                    umi_key = f"{umi_beta}_{umi_alpha}"
                    bid = _bucket_id(umi_key)
                    _get_bucket(bid).write(
                        f"{umi_key}\t{insert_beta}\t{bytes(iq_beta).hex()}"
                        f"\t{insert_alpha}\t{bytes(iq_alpha).hex()}\n"
                    )

        except Exception as e:
            logger.error("ERROR reading FASTQs for %s: %s", sample_label, e)
            logger.error(traceback.format_exc())
            raise

        for h in bucket_handles.values():
            h.close()
        bucket_handles.clear()

        if id_mismatch_count > 0:
            logger.warning("Total read ID mismatches: %d / %d pairs", id_mismatch_count, total_pairs)

        # --- Pass 2: call consensus per UMI group, write outputs ---
        cons_beta_out  = out_dir / f"{sample_label}_consensus_beta.fastq.gz"
        cons_alpha_out = out_dir / f"{sample_label}_consensus_alpha.fastq.gz"
        sing_beta_out  = out_dir / f"{sample_label}_singletons_beta.fastq.gz"
        sing_alpha_out = out_dir / f"{sample_label}_singletons_alpha.fastq.gz"
        qc_out         = out_dir / f"{sample_label}_consensus_qc.tsv"

        n_consensus_umis = 0; n_singleton_umis = 0; n_other_umis = 0
        reads_in_consensus = 0; reads_in_singletons = 0
        total_umis = 0
        counts: List[int] = []

        with (gzip.open(cons_beta_out, "wt") as cb_h,
              gzip.open(cons_alpha_out, "wt") as ca_h,
              gzip.open(sing_beta_out, "wt") as sb_h,
              gzip.open(sing_alpha_out, "wt") as sa_h,
              open(qc_out, "w") as qc_h):

            qc_h.write("umi_key\tn_reads\tbeta_cons_len\tbeta_nN\talpha_cons_len\talpha_nN"
                       "\tmin_log_delta_beta\tmedian_log_delta_beta"
                       "\tmin_log_delta_alpha\tmedian_log_delta_alpha\n")

            for bid in range(umi_bucket_count):
                bp = _bucket_path(bid)
                if not bp.exists():
                    continue

                umi_groups: Dict[str, List[Tuple]] = defaultdict(list)
                with gzip.open(bp, "rt") as bh:
                    for line in bh:
                        parts = line.rstrip("\n").split("\t")
                        umi_key, ib, iq_hex_b, ia, iq_hex_a = parts
                        umi_groups[umi_key].append((
                            ib, list(bytes.fromhex(iq_hex_b)),
                            ia, list(bytes.fromhex(iq_hex_a)),
                        ))

                for umi_key, records in umi_groups.items():
                    n = len(records)
                    total_umis += 1
                    counts.append(n)

                    umi_beta, umi_alpha = umi_key.rsplit("_", 1)
                    header = f"{umi_beta}|{umi_alpha}|n={n}"

                    if n >= min_reads_per_umi:
                        beta_seqs  = [r[0] for r in records]
                        beta_quals = [r[1] for r in records]
                        alpha_seqs  = [r[2] for r in records]
                        alpha_quals = [r[3] for r in records]

                        cons_b, qc_b = collapse_umis_posterior(
                            beta_seqs, beta_quals, posterior_min_log_delta)
                        cons_a, qc_a = collapse_umis_posterior(
                            alpha_seqs, alpha_quals, posterior_min_log_delta)

                        cb_h.write(f"@{header}\n{cons_b}\n+\n{placeholder_qual(len(cons_b), placeholder_char)}\n")
                        ca_h.write(f"@{header}\n{cons_a}\n+\n{placeholder_qual(len(cons_a), placeholder_char)}\n")

                        qc_h.write(
                            f"{umi_key}\t{n}\t{int(qc_b['len'])}\t{int(qc_b['nN'])}"
                            f"\t{int(qc_a['len'])}\t{int(qc_a['nN'])}"
                            f"\t{qc_b['min_log_delta']:.3f}\t{qc_b['median_log_delta']:.3f}"
                            f"\t{qc_a['min_log_delta']:.3f}\t{qc_a['median_log_delta']:.3f}\n"
                        )

                        n_consensus_umis += 1
                        reads_in_consensus += n

                    elif n == 1:
                        ib, iq_b, ia, iq_a = records[0]
                        sb_h.write(f"@{header}\n{ib}\n+\n{phred_to_fastq_ascii(iq_b)}\n")
                        sa_h.write(f"@{header}\n{ia}\n+\n{phred_to_fastq_ascii(iq_a)}\n")
                        n_singleton_umis += 1
                        reads_in_singletons += 1
                    else:
                        n_other_umis += 1

        # --- Summary ---
        summary_path = out_dir / f"{sample_label}_paired_umi_summary.txt"
        dist = Counter(counts)
        with open(summary_path, "w") as s:
            s.write(f"Sample: {sample_label}\n")
            s.write(f"R1 (beta):  {r1_fastq}\n")
            s.write(f"R2 (alpha): {r2_fastq}\n\n")
            s.write(f"Total read pairs processed: {total_pairs}\n\n")
            s.write("--- R1 (beta) dropout ---\n")
            s.write(f"  R1 too short:              {r1_too_short}\n")
            s.write(f"  R1 no anchor:              {r1_no_anchor}\n")
            s.write(f"  R1 anchor wrong position:  {r1_anchor_wrong_pos}\n")
            s.write(f"  R1 empty insert:           {r1_empty_insert}\n")
            s.write("--- R2 (alpha) dropout ---\n")
            s.write(f"  R2 too short:              {r2_too_short}\n")
            s.write(f"  R2 no anchor:              {r2_no_anchor}\n")
            s.write(f"  R2 anchor wrong position:  {r2_anchor_wrong_pos}\n")
            s.write(f"  R2 empty insert:           {r2_empty_insert}\n\n")
            s.write(f"Total paired UMIs detected:          {total_umis}\n")
            s.write(f"UMIs with >= {min_reads_per_umi} reads (consensus): {n_consensus_umis}\n")
            s.write(f"UMIs with 1 read  (singletons):      {n_singleton_umis}\n")
            s.write(f"UMIs with 0 reads:                   {n_other_umis}\n\n")
            s.write(f"Reads in consensus UMIs:  {reads_in_consensus}\n")
            s.write(f"Reads in singleton UMIs:  {reads_in_singletons}\n\n")
            s.write("Reads-per-UMI distribution:\n")
            for c in sorted(dist.keys()):
                s.write(f"  {c} reads: {dist[c]} UMIs\n")

        logger.info("Consensus beta:  %s", cons_beta_out)
        logger.info("Consensus alpha: %s", cons_alpha_out)
        logger.info("Singletons beta: %s", sing_beta_out)
        logger.info("Singletons alpha: %s", sing_alpha_out)
        logger.info("Summary: %s", summary_path)

        return {
            "pairs_in":                        total_pairs,
            "pairs_lost_r1_too_short":         r1_too_short,
            "pairs_lost_r1_no_anchor":         r1_no_anchor,
            "pairs_lost_r1_anchor_wrong_pos":  r1_anchor_wrong_pos,
            "pairs_lost_r2_too_short":         r2_too_short,
            "pairs_lost_r2_no_anchor":         r2_no_anchor,
            "pairs_lost_r2_anchor_wrong_pos":  r2_anchor_wrong_pos,
            "reads_in_consensus":              reads_in_consensus,
            "reads_in_singletons":             reads_in_singletons,
            "umis_total":                      total_umis,
            "umis_consensus":                  n_consensus_umis,
            "umis_singleton":                  n_singleton_umis,
        }

    finally:
        for h in bucket_handles.values():
            h.close()
        shutil.rmtree(bucket_temp_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="Paired-Chain Pipeline Step 02: Paired-end UMI consensus."
    )
    ap.add_argument("--yaml_config", required=True)
    ap.add_argument("--sample_id", default=None)
    args = ap.parse_args()

    cfg = parse_yaml(args.yaml_config)

    log_file = setup_pipeline_logging(
        logs_dir=cfg["logs_dir"],
        script_name=SCRIPT_NAME,
        scope="run",
        run_label=cfg.get("run_label"),
    )
    logger.info("Logging to: %s", log_file)
    logger.info("Loaded config: %s", cfg.get("_config_path", args.yaml_config))

    df = load_sample_sheet(cfg["sample_sheet"], required_cols=("sample_id", "sample_name"))

    if args.sample_id is not None:
        df = df[df["sample_id"].astype(str) == str(args.sample_id)]
        if df.empty:
            raise ValueError(f"No rows for sample_id={args.sample_id}")

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
            step="02_paired_umi_consensus",
            metrics=stats,
        )


if __name__ == "__main__":
    main()
