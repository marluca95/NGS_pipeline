"""
Paired-Chain Pipeline — Step 03: Paired Filtering and Extraction
================================================================

Implements the new 'ab-1xNNK-1xNNK' library mode for paired alpha-beta TCR
mutagenesis: exactly one codon (at most) may be mutated in each chain, and
any mutated codon must be of NNK pattern (no stop codons allowed).

Read structure coming in (from script2 output):
  beta FASTQ:  consensus sequences for the beta chain (CDR3β region)
  alpha FASTQ: consensus sequences for the alpha chain (CDR3α region)
  Both files are in LOCKSTEP — record N in beta = record N in alpha.
  Headers match: @{UMI_beta}|{UMI_alpha}|n={n}

Processing per pair:
  1. Find anchor_beta in beta read → extract 33 nt region → reverse complement.
  2. Find anchor_alpha in alpha read → extract 30 nt region → reverse complement.
  3. Validate (beta_region, alpha_region) against AbOneXNNKStrategy:
     - Each chain: max 1 mutated codon vs WT, must be NNK, no stop codons.
  4. PASS: translate both, write to {prefix}.filtered.PASS.paired.tsv
  5. FAIL: count by failure reason.

Inputs (per sample, resolved from samples_tsv + input_dir):
  {sample_label}_consensus_beta.fastq.gz
  {sample_label}_consensus_alpha.fastq.gz
  {sample_label}_singletons_rescued_beta.fastq.gz   (optional)
  {sample_label}_singletons_rescued_alpha.fastq.gz  (optional)

Outputs:
  {sample_label}.filtered.PASS.paired.tsv
  {sample_label}.filtered.FAIL.paired.tsv  (if write_fail_fastq=true)
  {sample_label}.filtered.summary.tsv
"""

import argparse
import csv
import gzip
import logging
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import IUPACData
from fuzzysearch import find_near_matches
from utils.config_utils import load_pipeline_config
from utils.logging_utils import setup_pipeline_logging
from utils.metrics_utils import write_sample_metrics
from utils.sample_utils import load_sample_sheet

SCRIPT_NAME = "script3_paired_filtering_extraction"
logger = logging.getLogger(__name__)

ALLOWED_LIBRARY_MODES = {"ab-1xNNK-1xNNK"}

# IUPAC ambiguous DNA values for NNK codon validation
IUPAC: Dict[str, set] = {k: set(v) for k, v in IUPACData.ambiguous_dna_values.items()}


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
def parse_yaml(path: str) -> dict:
    cfg = load_pipeline_config(
        path,
        required_keys=(
            "output_dir", "library_mode",
            "anchor_beta", "anchor_alpha",
            "anchor_max_mismatches",
            "region_length_beta_nt", "region_length_alpha_nt",
            "drop_if_contains_N",
            "write_fail_fastq",
            "WT_nt_beta", "WT_nt_alpha",
            "samples_tsv", "input_dir",
            "input_suffixes_consensus_beta",
            "input_suffixes_consensus_alpha",
        ),
        default_values={
            "max_mutated_codons_per_chain": 1,
            "write_aa_tsv": True,
            "logs_dir": None,
            "input_suffixes_singletons_beta": "",
            "input_suffixes_singletons_alpha": "",
            # Beta (R1) is anti-sense → RC needed; alpha (R2) is sense → no RC.
            # Override in config if your library orientation differs.
            "apply_rc_beta": True,
            "apply_rc_alpha": False,
        },
        path_keys=("output_dir", "samples_tsv", "input_dir", "logs_dir"),
    )

    mode = str(cfg["library_mode"]).strip()
    if mode not in ALLOWED_LIBRARY_MODES:
        raise ValueError(f"Invalid library_mode='{mode}'. Must be one of: {sorted(ALLOWED_LIBRARY_MODES)}")

    if not cfg.get("logs_dir"):
        cfg["logs_dir"] = str(Path(cfg["output_dir"]) / "_logs")

    wt_beta  = str(cfg["WT_nt_beta"]).upper().strip()
    wt_alpha = str(cfg["WT_nt_alpha"]).upper().strip()
    rlen_beta  = int(cfg["region_length_beta_nt"])
    rlen_alpha = int(cfg["region_length_alpha_nt"])

    if len(wt_beta) != rlen_beta:
        raise ValueError(
            f"WT_nt_beta length {len(wt_beta)} != region_length_beta_nt {rlen_beta}"
        )
    if len(wt_alpha) != rlen_alpha:
        raise ValueError(
            f"WT_nt_alpha length {len(wt_alpha)} != region_length_alpha_nt {rlen_alpha}"
        )

    return cfg


# ---------------------------------------------------------------------------
# Utility functions (adapted from script3_library_filtering_and_extraction.py)
# ---------------------------------------------------------------------------
def find_anchor(seq: str, anchor: str, max_mm: int) -> Optional[Tuple[int, int]]:
    """Find rightmost anchor in seq. Returns (start, end) or None."""
    seq = seq.upper()
    anchor = anchor.upper()
    if max_mm == 0:
        i = seq.rfind(anchor)
        return None if i == -1 else (i, i + len(anchor))
    matches = find_near_matches(anchor, seq, max_substitutions=max_mm,
                                max_insertions=0, max_deletions=0)
    if not matches:
        return None
    m = sorted(matches, key=lambda x: x.start)[-1]  # rightmost
    return (m.start, m.end)


def matches_degenerate_codon(codon: str, pattern: str) -> bool:
    codon = codon.upper()
    pattern = pattern.upper()
    if len(codon) != 3 or len(pattern) != 3:
        return False
    return all(codon[i] in IUPAC[pattern[i]] for i in range(3))


def translate_codon(codon: str) -> str:
    codon = codon.upper()
    if len(codon) != 3:
        return "X"
    try:
        aa = str(Seq(codon).translate(table=1, to_stop=False))
        return aa if aa else "X"
    except Exception:
        return "X"


def translate_region_nt(region_nt: str) -> str:
    region_nt = region_nt.upper()
    region_nt = region_nt[: (len(region_nt) // 3) * 3]
    try:
        return str(Seq(region_nt).translate(table=1, to_stop=False))
    except Exception:
        return "".join(translate_codon(region_nt[i:i + 3]) for i in range(0, len(region_nt), 3))


def revcomp(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


# ---------------------------------------------------------------------------
# Library strategy: ab-1xNNK-1xNNK
# ---------------------------------------------------------------------------
class AbOneXNNKStrategy:
    """
    Validates a (beta_nt, alpha_nt) pair against the ab-1xNNK-1xNNK library design:
      - Beta chain:  exactly 0 or 1 codon mutated vs WT_nt_beta (33 nt, 11 codons).
                     Any mutated codon must be NNK; no stop codons.
      - Alpha chain: exactly 0 or 1 codon mutated vs WT_nt_alpha (30 nt, 10 codons).
                     Same rules.
    """

    def __init__(self, wt_beta: str, wt_alpha: str, max_mut: int = 1):
        self.wt_beta  = wt_beta.upper()
        self.wt_alpha = wt_alpha.upper()
        self.max_mut  = max_mut
        self.n_codons_beta  = len(self.wt_beta) // 3
        self.n_codons_alpha = len(self.wt_alpha) // 3

    def _validate_chain(
        self, region_nt: str, wt_nt: str, n_codons: int, chain_name: str,
    ) -> Tuple[bool, str, List[int]]:
        """
        Validate one chain's region against its WT.
        Returns (ok, fail_reason, mutated_positions_1based).
        """
        if len(region_nt) != len(wt_nt):
            return False, f"{chain_name}_wrong_length", []

        mutated: List[int] = []
        for pos in range(n_codons):
            cod = region_nt[pos * 3: (pos + 1) * 3]
            ref = wt_nt[pos * 3: (pos + 1) * 3]
            if cod != ref:
                mutated.append(pos + 1)

        if len(mutated) > self.max_mut:
            return False, f"{chain_name}_too_many_mutations", mutated

        for pos in mutated:
            cod = region_nt[(pos - 1) * 3: pos * 3]
            if translate_codon(cod) == "*":
                return False, f"{chain_name}_stop_codon", mutated
            if not matches_degenerate_codon(cod, "NNK"):
                return False, f"{chain_name}_mutation_not_NNK", mutated

        return True, "ok", mutated

    def validate_pair(
        self, beta_nt: str, alpha_nt: str,
    ) -> Tuple[bool, str, Dict[str, Any]]:
        """
        Validate (beta_nt, alpha_nt) pair.
        Returns (ok, fail_reason, details).
        details keys: beta_mutated_positions, alpha_mutated_positions
        """
        ok_b, reason_b, muts_b = self._validate_chain(
            beta_nt, self.wt_beta, self.n_codons_beta, "beta")
        if not ok_b:
            return False, reason_b, {"beta_mutated_positions": muts_b,
                                      "alpha_mutated_positions": []}

        ok_a, reason_a, muts_a = self._validate_chain(
            alpha_nt, self.wt_alpha, self.n_codons_alpha, "alpha")
        if not ok_a:
            return False, reason_a, {"beta_mutated_positions": muts_b,
                                      "alpha_mutated_positions": muts_a}

        return True, "ok", {"beta_mutated_positions": muts_b,
                             "alpha_mutated_positions": muts_a}

    def update_stats(self, ok: bool, reason: str, details: Dict, stats: Counter) -> None:
        muts_b = details.get("beta_mutated_positions", [])
        muts_a = details.get("alpha_mutated_positions", [])
        if ok:
            stats[f"beta_mutations_{len(muts_b)}"] += 1
            stats[f"alpha_mutations_{len(muts_a)}"] += 1
        else:
            stats[f"fail_reason_{reason}"] += 1

    def write_strategy_summary(self, stats: Counter, out_handle) -> None:
        out_handle.write("\n# strategy=ab-1xNNK-1xNNK\n")
        out_handle.write("metric\tvalue\n")
        for k in sorted(stats):
            out_handle.write(f"{k}\t{stats[k]}\n")


# ---------------------------------------------------------------------------
# Input resolution
# ---------------------------------------------------------------------------
def resolve_inputs_for_sample(
    cfg: dict, sample_id: str,
) -> Tuple[str, List[Tuple[Path, Path]]]:
    """
    Returns (prefix_base, [(beta_path, alpha_path), ...])
    Each tuple is a pair of FASTQ files to process in lockstep.
    """
    df = load_sample_sheet(Path(cfg["samples_tsv"]), required_cols=("sample_id", "sample_name"))
    matched = df[df["sample_id"].astype(str).str.strip() == str(sample_id).strip()]
    if matched.empty:
        raise ValueError(f"sample_id '{sample_id}' not found in {cfg['samples_tsv']}")
    sample_name = str(matched.iloc[0]["sample_name"]).strip()
    prefix_base = f"{sample_id}_{sample_name}"

    in_dir = Path(cfg["input_dir"])
    pairs: List[Tuple[Path, Path]] = []

    cons_beta_suf  = str(cfg.get("input_suffixes_consensus_beta",  "") or "").strip()
    cons_alpha_suf = str(cfg.get("input_suffixes_consensus_alpha", "") or "").strip()
    sing_beta_suf  = str(cfg.get("input_suffixes_singletons_beta",  "") or "").strip()
    sing_alpha_suf = str(cfg.get("input_suffixes_singletons_alpha", "") or "").strip()

    if cons_beta_suf and cons_alpha_suf:
        bp = in_dir / f"{prefix_base}{cons_beta_suf}"
        ap = in_dir / f"{prefix_base}{cons_alpha_suf}"
        for p in (bp, ap):
            if not p.is_file():
                raise FileNotFoundError(f"Expected input not found: {p}")
        pairs.append((bp, ap))
    else:
        logger.warning("Consensus suffixes not set; skipping consensus inputs.")

    if sing_beta_suf and sing_alpha_suf:
        bp = in_dir / f"{prefix_base}{sing_beta_suf}"
        ap = in_dir / f"{prefix_base}{sing_alpha_suf}"
        if bp.is_file() and ap.is_file():
            pairs.append((bp, ap))
        elif bp.is_file() or ap.is_file():
            logger.warning(
                "Only one of the rescued singleton files exists for %s; skipping singletons.",
                prefix_base,
            )

    if not pairs:
        raise ValueError(f"No input FASTQ pairs resolved for sample_id={sample_id}")

    return prefix_base, pairs


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------
def process_sample(
    prefix: str,
    input_pairs: List[Tuple[Path, Path]],
    out_dir: Path,
    strategy: AbOneXNNKStrategy,
    anchor_beta: str, anchor_alpha: str,
    anchor_mm: int,
    region_len_beta: int, region_len_alpha: int,
    dropN: bool,
    write_fail: bool,
    apply_rc_beta: bool = True,
    apply_rc_alpha: bool = False,
) -> Counter:
    out_dir.mkdir(parents=True, exist_ok=True)

    pass_path    = out_dir / f"{prefix}.filtered.PASS.paired.tsv"
    fail_path    = out_dir / f"{prefix}.filtered.FAIL.paired.tsv"
    summary_path = out_dir / f"{prefix}.filtered.summary.tsv"

    counts = Counter()
    strategy_stats = Counter()

    pass_fh = open(pass_path, "w", newline="")
    pass_writer = csv.writer(pass_fh, delimiter="\t")
    pass_writer.writerow([
        "read_id",
        "beta_nt", "alpha_nt",
        "beta_aa", "alpha_aa",
        "beta_mut_pos", "alpha_mut_pos",
    ])

    fail_fh = open(fail_path, "w", newline="") if write_fail else None
    if fail_fh:
        fail_writer = csv.writer(fail_fh, delimiter="\t")
        fail_writer.writerow(["read_id", "fail_reason", "beta_seq", "alpha_seq"])

    for beta_path, alpha_path in input_pairs:
        logger.info("Processing pair: %s | %s", beta_path.name, alpha_path.name)

        with gzip.open(beta_path, "rt") as bh, gzip.open(alpha_path, "rt") as ah:
            for beta_rec, alpha_rec in zip(SeqIO.parse(bh, "fastq"), SeqIO.parse(ah, "fastq")):
                counts["total_reads"] += 1

                beta_seq  = str(beta_rec.seq).upper()
                alpha_seq = str(alpha_rec.seq).upper()
                read_id = beta_rec.id  # same as alpha_rec.id

                # --- Beta region extraction ---
                hit_b = find_anchor(beta_seq, anchor_beta, anchor_mm)
                if hit_b is None:
                    counts["fail_beta_anchor_not_found"] += 1
                    if fail_fh:
                        fail_writer.writerow([read_id, "beta_anchor_not_found", beta_seq, alpha_seq])
                    continue
                _, b_end = hit_b
                beta_region = beta_seq[b_end: b_end + region_len_beta]
                if len(beta_region) != region_len_beta:
                    counts["fail_beta_region_too_short"] += 1
                    if fail_fh:
                        fail_writer.writerow([read_id, "beta_region_too_short", beta_seq, alpha_seq])
                    continue
                if apply_rc_beta:
                    beta_region = revcomp(beta_region)

                # --- Alpha region extraction ---
                hit_a = find_anchor(alpha_seq, anchor_alpha, anchor_mm)
                if hit_a is None:
                    counts["fail_alpha_anchor_not_found"] += 1
                    if fail_fh:
                        fail_writer.writerow([read_id, "alpha_anchor_not_found", beta_seq, alpha_seq])
                    continue
                _, a_end = hit_a
                alpha_region = alpha_seq[a_end: a_end + region_len_alpha]
                if len(alpha_region) != region_len_alpha:
                    counts["fail_alpha_region_too_short"] += 1
                    if fail_fh:
                        fail_writer.writerow([read_id, "alpha_region_too_short", beta_seq, alpha_seq])
                    continue
                if apply_rc_alpha:
                    alpha_region = revcomp(alpha_region)

                # --- N-content filter ---
                if dropN and ("N" in beta_region):
                    counts["fail_beta_contains_N"] += 1
                    if fail_fh:
                        fail_writer.writerow([read_id, "beta_contains_N", beta_seq, alpha_seq])
                    continue
                if dropN and ("N" in alpha_region):
                    counts["fail_alpha_contains_N"] += 1
                    if fail_fh:
                        fail_writer.writerow([read_id, "alpha_contains_N", beta_seq, alpha_seq])
                    continue

                # --- Pair validation ---
                ok, reason, details = strategy.validate_pair(beta_region, alpha_region)
                strategy.update_stats(ok, reason, details, strategy_stats)

                if ok:
                    counts["pass"] += 1
                    beta_aa  = translate_region_nt(beta_region)
                    alpha_aa = translate_region_nt(alpha_region)
                    muts_b = details.get("beta_mutated_positions", [])
                    muts_a = details.get("alpha_mutated_positions", [])
                    muts_b_str = ",".join(str(x) for x in muts_b)
                    muts_a_str = ",".join(str(x) for x in muts_a)
                    pass_writer.writerow([
                        read_id,
                        beta_region, alpha_region,
                        beta_aa, alpha_aa,
                        muts_b_str, muts_a_str,
                    ])
                else:
                    counts[f"fail_{reason}"] += 1
                    if fail_fh:
                        fail_writer.writerow([read_id, reason, beta_seq, alpha_seq])

    pass_fh.close()
    if fail_fh:
        fail_fh.close()

    # --- Summary ---
    with open(summary_path, "w") as s:
        s.write("# inputs\n")
        for bp, ap in input_pairs:
            s.write(f"# beta:  {bp}\n")
            s.write(f"# alpha: {ap}\n")
        s.write("\nmetric\tvalue\n")
        for k in sorted(counts.keys()):
            s.write(f"{k}\t{counts[k]}\n")
        strategy.write_strategy_summary(strategy_stats, s)

    logger.info("PASS: %s", pass_path)
    if write_fail:
        logger.info("FAIL: %s", fail_path)
    logger.info("Summary: %s", summary_path)

    return counts


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="Paired-Chain Pipeline Step 03: Filter and extract paired alpha-beta regions."
    )
    ap.add_argument("--yaml_config", required=True)
    ap.add_argument("--sample_id", required=True)
    args = ap.parse_args()

    cfg = parse_yaml(args.yaml_config)

    log_path = setup_pipeline_logging(
        logs_dir=cfg["logs_dir"],
        script_name=SCRIPT_NAME,
        scope="run",
        run_label=cfg.get("run_label"),
    )
    logger.info("Logging to: %s", log_path)
    logger.info("Loaded config: %s", cfg.get("_config_path", args.yaml_config))

    prefix, input_pairs = resolve_inputs_for_sample(cfg, sample_id=str(args.sample_id))
    logger.info("sample_id=%s → prefix=%s, %d input pair(s)", args.sample_id, prefix, len(input_pairs))

    strategy = AbOneXNNKStrategy(
        wt_beta=str(cfg["WT_nt_beta"]).upper(),
        wt_alpha=str(cfg["WT_nt_alpha"]).upper(),
        max_mut=int(cfg.get("max_mutated_codons_per_chain", 1)),
    )

    counts = process_sample(
        prefix=prefix,
        input_pairs=input_pairs,
        out_dir=Path(cfg["output_dir"]),
        strategy=strategy,
        anchor_beta=str(cfg["anchor_beta"]),
        anchor_alpha=str(cfg["anchor_alpha"]),
        anchor_mm=int(cfg["anchor_max_mismatches"]),
        region_len_beta=int(cfg["region_length_beta_nt"]),
        region_len_alpha=int(cfg["region_length_alpha_nt"]),
        dropN=bool(cfg["drop_if_contains_N"]),
        write_fail=bool(cfg["write_fail_fastq"]),
        apply_rc_beta=bool(cfg.get("apply_rc_beta", True)),
        apply_rc_alpha=bool(cfg.get("apply_rc_alpha", False)),
    )

    _parts = prefix.split("_", 1)
    _sample_id   = _parts[0]
    _sample_name = _parts[1] if len(_parts) > 1 else prefix

    write_sample_metrics(
        metrics_path=Path(cfg["output_dir"]) / "per_sample_metrics.tsv",
        sample_id=_sample_id,
        sample_name=_sample_name,
        step="03_paired_extraction",
        metrics={
            "reads_in_total":                     counts.get("total_reads", 0),
            "reads_pass":                         counts.get("pass", 0),
            "reads_lost_beta_no_anchor":          counts.get("fail_beta_anchor_not_found", 0),
            "reads_lost_alpha_no_anchor":         counts.get("fail_alpha_anchor_not_found", 0),
            "reads_lost_beta_too_short":          counts.get("fail_beta_region_too_short", 0),
            "reads_lost_alpha_too_short":         counts.get("fail_alpha_region_too_short", 0),
            "reads_lost_contains_N":              counts.get("fail_beta_contains_N", 0)
                                                  + counts.get("fail_alpha_contains_N", 0),
            "reads_lost_beta_too_many_mutations": counts.get("fail_beta_too_many_mutations", 0),
            "reads_lost_alpha_too_many_mutations":counts.get("fail_alpha_too_many_mutations", 0),
            "reads_lost_stop_codon":              counts.get("fail_beta_stop_codon", 0)
                                                  + counts.get("fail_alpha_stop_codon", 0),
            "reads_lost_not_NNK":                 counts.get("fail_beta_mutation_not_NNK", 0)
                                                  + counts.get("fail_alpha_mutation_not_NNK", 0),
        },
    )


if __name__ == "__main__":
    main()
