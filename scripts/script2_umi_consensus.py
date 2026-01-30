"""
UMI consensus generation (posterior/quality-aware hard consensus + keep singletons),
adapted for library design:

  [spacer (0–3 nt)] [UMI (umi_length nt)] [anchor_sequence] [insert]

Per read:
  - Find anchor (fuzzy, up to max_mismatches)
  - Accept only if anchor_start in [umi_length + spacer_min, umi_length + spacer_max]
  - Extract UMI = umi_length bases immediately before anchor
  - Extract insert = everything after anchor (anchor removed)
  - Keep insert qualities as well

Outputs per sample:
  output_dir/umi/<sample>_consensus.fastq.gz
  output_dir/umi/<sample>_singletons.fastq.gz
  output_dir/umi/<sample>_umi_summary.txt
  output_dir/umi/<sample>_consensus_qc.tsv

FASTQ headers:
  @{UMI}|n={read_count}

Consensus calling:
  - For each position, compute log-likelihood for candidate A/C/G/T under a simple
    sequencing error model: match => log(1-p_err), mismatch => log(p_err/3)
  - Pick best base only if (best - second_best) >= min_log_delta; else output 'N'
  - Produce a consensus confidence measure derived from best-vs-second evidence.

NOTE:
  - FASTQ qualities for consensus reads are replaced with a PLACEHOLDER character
    (default 'I') to avoid implying they are real sequencer-derived base qualities.
  - Real diagnostics are written to *_consensus_qc.tsv.

This method assumes errors are independent across reads and substitutions are uniform.
"""

import argparse
import gzip
import logging
import math
import traceback
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import pandas as pd
import yaml
from Bio import SeqIO
from fuzzysearch import find_near_matches


# -------- logging --------
def setup_logging(logs_dir: Path, log_name: str = "umi_pipeline.log") -> Path:
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_file = logs_dir / log_name

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
        force=True,
    )
    return log_file


logger = logging.getLogger(__name__)


# -------- config / IO helpers --------
def parse_yaml(yaml_file: str) -> dict:
    with open(yaml_file, "r") as f:
        cfg = yaml.safe_load(f) or {}

    required = ["sample_sheet", "input_dir", "output_dir", "umi_length", "anchor_sequence"]
    missing = [k for k in required if k not in cfg]
    if missing:
        raise ValueError(f"Missing required YAML keys: {missing}")

    cfg.setdefault("logs_dir", str(Path(cfg["output_dir"]) / "_logs"))

    cfg.setdefault("combined_suffix", "_combined_trimmed.fastq.gz")
    cfg.setdefault("max_mismatches", 1)
    cfg.setdefault("min_reads_per_umi", 2)
    cfg.setdefault("write_readcount_in_header", True)

    # library layout
    cfg.setdefault("spacer_min", 0)
    cfg.setdefault("spacer_max", 3)

    # posterior consensus params
    cfg.setdefault("posterior_min_log_delta", 2.0)

    # placeholder qualities for consensus FASTQ
    # 'I' => Q40 (still interpretable, but constant); '!' => Q0 (strong "ignore" signal)
    cfg.setdefault("consensus_placeholder_qual_char", "I")

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
        raise ValueError(
            f"Sample sheet missing required columns: {missing}. Found: {str(list(df.columns))}"
        )
    return df


def phred_to_fastq_ascii(phred: List[int]) -> str:
    # Cap at 93 to keep ASCII printable after +33
    return "".join(chr(min(93, max(0, q)) + 33) for q in phred)


def placeholder_qual(length: int, char: str = "I") -> str:
    # FASTQ quality string must match sequence length
    return char * max(0, int(length))


def find_anchor_start_valid(
    seq: str,
    anchor: str,
    umi_len: int,
    max_dist: int,
    spacer_min: int = 0,
    spacer_max: int = 3,
) -> Tuple[Optional[Tuple[int, int]], bool]:
    """
    Fuzzy anchor search + position validation.

    Returns:
      (hit, has_any_anchor)

    - hit: (anchor_start, anchor_end) if a valid-position match exists, else None
    - has_any_anchor: True if ANY fuzzy match exists anywhere in the read, else False

    This preserves the original fuzzysearch behavior for what counts as a "match"
    (substitutions + indels up to max_dist), but lets you separate:
      - no anchor found
      - anchor found but not at valid position
    """
    matches = find_near_matches(anchor, seq, max_l_dist=max_dist)
    if not matches:
        return None, False  # no anchor anywhere

    lo = umi_len + spacer_min
    hi = umi_len + spacer_max

    # Keep original behavior: choose the left-most valid-position match.
    for m in sorted(matches, key=lambda x: x.start):
        if lo <= m.start <= hi:
            return (m.start, m.end), True

    # anchor exists, but only at invalid position(s)
    return None, True


# -------- posterior consensus logic --------
def phred_to_error_prob(q: int) -> float:
    return 10 ** (-q / 10.0)


def posterior_base_call(
    bases: List[str],
    quals: List[int],
    min_log_delta: float = 2.0,
) -> Tuple[str, float]:
    """
    Returns:
      (base_call, log_delta)
    where log_delta = best_logL - second_best_logL.

    Ignores any base not in A/C/G/T.
    """
    obs = [(b, q) for b, q in zip(bases, quals) if b in "ACGT"]
    if not obs:
        return "N", 0.0

    logL: Dict[str, float] = {}
    for candidate in "ACGT":
        ll = 0.0
        for b, q in obs:
            p_err = phred_to_error_prob(int(q))
            p_err = min(max(p_err, 1e-12), 1 - 1e-12)
            if b == candidate:
                ll += math.log(1.0 - p_err)
            else:
                ll += math.log(p_err / 3.0)
        logL[candidate] = ll

    vals_sorted = sorted(logL.items(), key=lambda kv: kv[1], reverse=True)
    best_base, best_ll = vals_sorted[0]
    second_ll = vals_sorted[1][1]
    log_delta = best_ll - second_ll

    if log_delta < min_log_delta:
        return "N", log_delta

    return best_base, log_delta


def collapse_umis_posterior(
    seqs: List[str],
    quals: List[List[int]],
    min_log_delta: float = 2.0,
) -> Tuple[str, List[int], Dict[str, float]]:
    if not seqs:
        return "", [], {"len": 0.0, "nN": 0.0, "min_log_delta": 0.0, "median_log_delta": 0.0, "p10_log_delta": 0.0}

    read_length = max(len(s) for s in seqs)
    consensus_bases: List[str] = []

    log_deltas: List[float] = []
    nN = 0

    for i in range(read_length):
        bases_i: List[str] = []
        quals_i: List[int] = []

        for s, q in zip(seqs, quals):
            if i >= len(s):
                continue
            bases_i.append(s[i])
            quals_i.append(q[i] if i < len(q) else 0)

        if not bases_i:
            consensus_bases.append("N")
            nN += 1
            continue

        call, log_delta = posterior_base_call(bases_i, quals_i, min_log_delta=min_log_delta)
        consensus_bases.append(call)
        log_deltas.append(float(log_delta))

        if call == "N":
            nN += 1

    cons = "".join(consensus_bases)

    def _median(xs: List[float]) -> float:
        if not xs:
            return 0.0
        ys = sorted(xs)
        m = len(ys) // 2
        return ys[m] if len(ys) % 2 == 1 else 0.5 * (ys[m - 1] + ys[m])

    def _p10(xs: List[float]) -> float:
        if not xs:
            return 0.0
        ys = sorted(xs)
        idx = max(0, int(0.10 * (len(ys) - 1)))
        return ys[idx]

    qc = {
        "len": float(len(cons)),
        "nN": float(nN),
        "min_log_delta": float(min(log_deltas)) if log_deltas else 0.0,
        "median_log_delta": float(_median(log_deltas)),
        "p10_log_delta": float(_p10(log_deltas)),
    }

    return cons, qc


# -------- main per-sample processing --------
def process_one_sample(sample_id: str, sample_name: str, cfg: dict) -> None:
    
    # get sample paths and fastaQ files
    sample_label = f"{sample_id}_{(sample_name)}"
    in_base = Path(cfg["input_dir"])
    input_sample_dir = in_base / sample_label
    
    output_sample_dir = Path(cfg["output_dir"])
    output_sample_dir.mkdir(parents=True, exist_ok=True)
    

    combined_suffix = cfg.get("combined_suffix", "_combined_trimmed.fastq.gz")
    matches = list(input_sample_dir.glob(f"*{combined_suffix}"))
    if not matches:
        raise FileNotFoundError(f"No *{combined_suffix} found in input sample dir: {input_sample_dir}")
    input_fastq = matches[0]  

    # logging
    logger.info(f"Processing sample: {sample_label}")
    logger.info(f"Input sample dir: {input_sample_dir}")
    logger.info(f"Input FASTQ: {input_fastq}")
    logger.info(f"UMI output dir: {output_sample_dir}")

    # get papameters 
    umi_len = int(cfg["umi_length"])
    anchor_seq = str(cfg["anchor_sequence"])
    max_anchor_dist = int(cfg.get("max_mismatches", 1))
    min_reads_per_umi = int(cfg.get("min_reads_per_umi", 2))

    spacer_min = int(cfg.get("spacer_min", 0))
    spacer_max = int(cfg.get("spacer_max", 3))

    posterior_min_log_delta = float(cfg.get("posterior_min_log_delta", 2.0))

    consensus_placeholder_char = str(cfg.get("consensus_placeholder_qual_char", "I"))
    if len(consensus_placeholder_char) != 1:
        raise ValueError("consensus_placeholder_qual_char must be a single ASCII character (e.g. 'I' or '!').")

    # set variables and Dictionary
    umi_groups: Dict[str, List[Tuple[str, List[int]]]] = defaultdict(list)

    total_reads = 0
    too_short = 0
    no_anchor = 0
    anchor_wrong_pos = 0
    empty_insert = 0

    # min lenght of read to fin anchor
    min_prefix = spacer_max + umi_len + len(anchor_seq)

    # Parse one fastQ file, record by record
    try:
        with gzip.open(input_fastq, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                total_reads += 1
                seq = str(record.seq)

                if len(seq) < min_prefix:
                    too_short += 1
                    continue

                hit, has_any_anchor = find_anchor_start_valid(
                    seq=seq,
                    anchor=anchor_seq,
                    umi_len=umi_len,
                    max_dist=max_anchor_dist,
                    spacer_min=spacer_min,
                    spacer_max=spacer_max,
                )
                if hit is None:
                    if has_any_anchor:
                        anchor_wrong_pos += 1
                    else:
                        no_anchor += 1
                    continue

                anchor_start, anchor_end = hit

                if anchor_start < umi_len:
                    too_short += 1
                    continue

                umi = seq[anchor_start - umi_len: anchor_start]

                insert = seq[anchor_end:]
                phred_full = record.letter_annotations.get("phred_quality", [])
                insert_q = phred_full[anchor_end:]

                if not insert:
                    empty_insert += 1
                    continue


                if len(insert_q) != len(insert):
                    insert_q = insert_q[:len(insert)]

                umi_groups[umi].append((insert, bytes(insert_q)))


    except Exception as e:
        logger.error(f"ERROR reading FASTQ for sample {sample_label}: {e}")
        logger.error(traceback.format_exc())
        raise

    # Write Outputs
    hard_out = output_sample_dir / f"{sample_label}_consensus.fastq.gz"
    singletons_out = output_sample_dir / f"{sample_label}_singletons.fastq.gz"
    qc_out = output_sample_dir / f"{sample_label}_consensus_qc.tsv"

    n_consensus_umis = 0
    n_singleton_umis = 0
    n_other_umis = 0
    reads_in_consensus = 0
    reads_in_singletons = 0

    with gzip.open(hard_out, "wt") as hard_handle, gzip.open(singletons_out, "wt") as sing_handle, open(qc_out, "w") as qc_handle:
        qc_handle.write("umi\tn_reads\tcons_len\tnN\tmin_log_delta\tmedian_log_delta\tp10_log_delta\n")

        for umi, seq_qual_pairs in umi_groups.items():
            n = len(seq_qual_pairs)

            if n >= min_reads_per_umi:
                seqs = [s for s, _ in seq_qual_pairs]
                quals = [list(q) for _, q in seq_qual_pairs]  # convert bytes back to ints


                # UMI Consensus generation
                cons, qc = collapse_umis_posterior(
                    seqs,
                    quals,
                    min_log_delta=posterior_min_log_delta,
                )

                header = f"{umi}|n={n}"

                hard_handle.write(
                    f"@{header}\n{cons}\n+\n{placeholder_qual(len(cons), char=consensus_placeholder_char)}\n"
                )

                qc_handle.write(
                    f"{umi}\t{n}\t{int(qc['len'])}\t{int(qc['nN'])}\t"
                    f"{qc['min_log_delta']:.3f}\t{qc['median_log_delta']:.3f}\t{qc['p10_log_delta']:.3f}\n"
                )

                n_consensus_umis += 1
                reads_in_consensus += n

            elif n == 1:
                s, q = seq_qual_pairs[0]
                header = f"{umi}|n={n}"
                sing_handle.write(f"@{header}\n{s}\n+\n{phred_to_fastq_ascii(q)}\n")

                n_singleton_umis += 1
                reads_in_singletons += 1

            else:
                n_other_umis += 1

    logger.info(f"Hard consensus written: {hard_out} (UMIs>={min_reads_per_umi}: {n_consensus_umis})")
    logger.info(f"Singletons written: {singletons_out} (UMIs==1: {n_singleton_umis})")
    logger.info(f"Consensus QC written: {qc_out}")

    # Summary
    summary = output_sample_dir / f"{sample_label}_umi_summary.txt"
    counts = [len(v) for v in umi_groups.values()]
    dist = Counter(counts)
    total_umis = len(umi_groups)

    with open(summary, "w") as s:
        s.write(f"Sample: {sample_label}\n")
        s.write(f"Input FASTQ: {input_fastq}\n")
        s.write("\n")
        s.write(f"Total reads processed: {total_reads}\n")
        s.write(f"Reads too short for spacer+UMI+anchor (dropped): {too_short}\n")
        s.write(f"Reads with no anchor found (dropped): {no_anchor}\n")
        s.write(f"Reads with anchor at invalid position (dropped): {anchor_wrong_pos}\n")
        s.write(f"Reads with empty insert after anchor (dropped): {empty_insert}\n")
        s.write("\n")
        s.write(f"Total UMIs detected: {total_umis}\n")
        s.write(f"UMIs with >= {min_reads_per_umi} reads (consensus): {n_consensus_umis}\n")
        s.write(f"UMIs with 1 read (singletons kept): {n_singleton_umis}\n")
        s.write(f"UMIs with 0 reads: {n_other_umis}\n")
        s.write("\n")
        s.write(f"Reads contributing to consensus UMIs: {reads_in_consensus}\n")
        s.write(f"Reads in singleton UMIs: {reads_in_singletons}\n")
        s.write("\n")
        s.write("Reads per UMI distribution:\n")
        for c in sorted(dist.keys()):
            s.write(f"  {c} reads: {dist[c]} UMIs\n")
        s.write("\n")
        s.write(f"Consensus QC TSV: {qc_out}\n")
        s.write(f"Consensus FASTQ placeholder qual char: '{consensus_placeholder_char}'\n")

    logger.info(f"Summary written: {summary}")
    logger.info(f"UMI processing complete for {sample_label}.")


def main() -> None:
    ap = argparse.ArgumentParser(description="UMI consensus (posterior quality-aware, keep singletons)")
    ap.add_argument("--yaml_config", required=True)
    ap.add_argument("--sample_id", default=None, help="Process only one sample_id")
    args = ap.parse_args()

    cfg = parse_yaml(args.yaml_config)

    logs_dir = Path(cfg["logs_dir"])
    log_file = setup_logging(logs_dir)
    logger.info(f"Logging to: {log_file}")

    df = load_sample_sheet(cfg["sample_sheet"])

    if args.sample_id is not None:
        df = df[df["sample_id"].astype(str) == str(args.sample_id)]
        if df.empty:
            raise ValueError(f"No rows for sample_id={args.sample_id}")

    grouped = df.groupby("sample_id", sort=False)

    for sample_id, g in grouped:
        names = g["sample_name"].dropna().unique()
        if len(names) != 1:
            raise ValueError(f"sample_id={sample_id} has multiple sample_name values: {names}")
        process_one_sample(str(sample_id), str(names[0]), cfg)


if __name__ == "__main__":
    main()
