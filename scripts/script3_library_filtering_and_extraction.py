#!/usr/bin/env python3
import argparse
import gzip
import logging
import csv
from collections import Counter
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import IUPACData
from fuzzysearch import find_near_matches
from utils.config_utils import load_pipeline_config
from utils.logging_utils import setup_pipeline_logging
from utils.metrics_utils import write_sample_metrics
from utils.sample_utils import load_sample_sheet


# ---------------- IUPAC + translation via Biopython ----------------
# IUPAC ambiguous DNA values, e.g. "R" -> "AG"
IUPAC: Dict[str, set] = {k: set(v) for k, v in IUPACData.ambiguous_dna_values.items()}

# ---- IMPORTANT: design-orientation codons from the paper/figure ----
DEFAULT_DEG_CODONS = ["SCN", "RGY", "RGY", "VNB", "ARY", "VNB", "GYN", "RRY", "GRR", "CWR", "TWY"]
ALLOWED_LIBRARY_MODES = {"combinatorial", "3xNNK"}
SCRIPT_NAME = "script3_library_filtering_and_extraction"
logger = logging.getLogger(__name__)


# ---------------- utility ----------------
def setup_logging(cfg: dict) -> Path:
    logs_dir = cfg.get("logs_dir") or (Path(cfg["output_dir"]) / "_logs")
    return setup_pipeline_logging(
        logs_dir=logs_dir,
        script_name=SCRIPT_NAME,
        scope="run",
        run_label=cfg.get("run_label"),
    )


def parse_yaml(path: str) -> dict:
    cfg = load_pipeline_config(
        path,
        required_keys=(
            "output_dir",
            "anchor_sequence",
            "anchor_max_mismatches",
            "region_length_nt",
            "drop_if_contains_N",
            "write_fail_fastq",
            "library_mode",
            "samples_tsv",
            "input_dir",
            "input_suffixes_consensus",
            "input_suffixes_singletons",
        ),
        default_values={
            "max_mutated_codons": 3,
            "write_aa_tsv": True,
            "logs_dir": None,
        },
        path_keys=("output_dir", "samples_tsv", "input_dir", "logs_dir"),
    )

    mode = str(cfg["library_mode"]).strip()
    if mode not in ALLOWED_LIBRARY_MODES:
        raise ValueError(f"Invalid library_mode='{mode}'. Must be one of: {sorted(ALLOWED_LIBRARY_MODES)}")

    cfg.setdefault("degenerate_codons", DEFAULT_DEG_CODONS)
    if not cfg.get("logs_dir"):
        cfg["logs_dir"] = str(Path(cfg["output_dir"]) / "_logs")

    # mode-specific checks
    if mode == "combinatorial":
        if "degenerate_codons" not in cfg:
            raise ValueError("library_mode=combinatorial requires 'degenerate_codons'.")
        if len(cfg["degenerate_codons"]) != 11:
            raise ValueError(f"'degenerate_codons' must have length 11, got {len(cfg['degenerate_codons'])}.")
    elif mode == "3xNNK":
        if "WT_nt_33" not in cfg:
            raise ValueError("library_mode=3xNNK requires 'WT_nt_33' (exact 33 nt template).")

    return cfg


def find_anchor(seq: str, anchor: str, max_mm: int) -> Optional[Tuple[int, int]]:
    """
    Find anchor in seq. Uses the RIGHTMOST hit to avoid grabbing an upstream anchor
    if the motif appears more than once in the read.
    """
    seq = seq.upper()
    anchor = anchor.upper()

    if max_mm == 0:
        i = seq.rfind(anchor)  # RIGHTMOST exact hit
        if i == -1:
            return None
        return (i, i + len(anchor))

    matches = find_near_matches(anchor, seq, max_substitutions=max_mm, max_insertions=0, max_deletions=0)
    if not matches:
        return None

    m = sorted(matches, key=lambda x: x.start)[-1]  # rightmost hit
    return (m.start, m.end)


def matches_degenerate_codon(codon: str, pattern: str) -> bool:
    codon = codon.upper()
    pattern = pattern.upper()
    if len(codon) != 3 or len(pattern) != 3:
        return False
    return all(codon[i] in IUPAC[pattern[i]] for i in range(3))


def translate_codon(codon: str) -> str:
    """
    Translate a single codon using Biopython.
    Returns "*" for stop, "X" for invalid/ambiguous.
    """
    codon = codon.upper()
    if len(codon) != 3:
        return "X"
    try:
        aa = str(Seq(codon).translate(table=1, to_stop=False))
        return aa if aa else "X"
    except Exception:
        return "X"


def translate_region_nt(region_nt: str) -> str:
    """Translate an nt region in-frame codons (len truncated to multiple of 3)."""
    region_nt = region_nt.upper()
    region_nt = region_nt[: (len(region_nt) // 3) * 3]
    try:
        return str(Seq(region_nt).translate(table=1, to_stop=False))
    except Exception:
        return "".join(translate_codon(region_nt[i:i + 3]) for i in range(0, len(region_nt), 3))


def revcomp(seq: str) -> str:
    """Reverse complement using Biopython (handles IUPAC ambiguous bases)."""
    return str(Seq(seq).reverse_complement())


def lookup_sample_name(samples_tsv: Path, sample_id: str) -> str:
    """
    samples.tsv expected columns (tab-separated):
      sample_id  sample_name  run_depth  lane  read  fastq
    Returns the first sample_name matching sample_id.
    """
    sample_id = sample_id.strip()
    df = load_sample_sheet(samples_tsv, required_cols=("sample_id", "sample_name"))
    matched = df[df["sample_id"].astype(str).str.strip() == sample_id]
    if not matched.empty:
        return str(matched.iloc[0]["sample_name"]).strip()
    raise ValueError(f"sample_id '{sample_id}' not found in {samples_tsv}")


def resolve_inputs_for_sample(cfg: dict, sample_id: str) -> Tuple[str, List[Path]]:
    """
    No mode selection.

    It will look for BOTH (consensus + singletons) using FLAT YAML suffix keys:
      input_suffixes_consensus
      input_suffixes_singletons

    If a suffix value is empty/None -> it is skipped.
    If both are empty -> raise error.
    If a suffix is provided but the resulting file does not exist -> raise error.

    Output prefix: {sample_id}_{sample_name}
    """
    samples_tsv = Path(cfg["samples_tsv"])
    input_dir = Path(cfg["input_dir"])

    sample_name = lookup_sample_name(samples_tsv, sample_id)
    prefix_base = f"{sample_id}_{sample_name}"

    consensus_suf = str(cfg.get("input_suffixes_consensus") or "").strip()
    singletons_suf = str(cfg.get("input_suffixes_singletons") or "").strip()

    if not consensus_suf and not singletons_suf:
        raise ValueError(
            "Both input_suffixes_consensus and input_suffixes_singletons are empty. Nothing to process."
        )

    fqs: List[Path] = []
    if consensus_suf:
        fqs.append(input_dir / f"{prefix_base}{consensus_suf}")
    else:
        logging.warning("input_suffixes_consensus is empty -> skipping consensus input.")

    if singletons_suf:
        fqs.append(input_dir / f"{prefix_base}{singletons_suf}")
    else:
        logging.warning("input_suffixes_singletons is empty -> skipping singletons input.")

    for p in fqs:
        if not p.is_file():
            raise FileNotFoundError(f"Expected input FASTQ not found: {p}")

    return prefix_base, fqs


# ---------------- Strategy base + implementations ----------------
class BaseLibraryStrategy:
    name = "base"

    def validate(self, region_nt: str) -> Tuple[bool, str, Dict[str, Any]]:
        raise NotImplementedError

    def update_stats(self, ok: bool, reason: str, details: Dict[str, Any], stats: Counter) -> None:
        return

    def write_strategy_summary(self, stats: Counter, out_handle) -> None:
        if not stats:
            return
        out_handle.write(f"\n# strategy={self.name}\n")
        out_handle.write("strategy_metric\tvalue\n")
        for k in sorted(stats):
            out_handle.write(f"{k}\t{stats[k]}\n")


class CombinatorialStrategy(BaseLibraryStrategy):
    name = "combinatorial"

    def __init__(self, degenerate_codons: List[str]):
        if len(degenerate_codons) != 11:
            raise ValueError(f"degenerate_codons must have length 11, got {len(degenerate_codons)}")
        self.degenerate_codons = [c.upper() for c in degenerate_codons]

    def validate(self, region_nt: str) -> Tuple[bool, str, Dict[str, Any]]:
        if len(region_nt) != 33:
            return False, "wrong_length", {}

        for pos in range(11):
            codon = region_nt[pos * 3:(pos + 1) * 3]

            if translate_codon(codon) == "*":
                return False, "stop_codon", {"fail_pos": pos + 1}

            pattern = self.degenerate_codons[pos]
            if not matches_degenerate_codon(codon, pattern):
                return False, "degenerate_mismatch", {"fail_pos": pos + 1}

        return True, "ok", {}

    def update_stats(self, ok: bool, reason: str, details: Dict[str, Any], stats: Counter) -> None:
        if ok:
            return
        fail_pos = details.get("fail_pos")
        if fail_pos is not None and reason in ("degenerate_mismatch", "stop_codon"):
            stats[f"fail_pos_{fail_pos:02d}"] += 1

    def write_strategy_summary(self, stats: Counter, out_handle) -> None:
        out_handle.write(f"\n# strategy={self.name}\n")
        out_handle.write("fail_position_1based\tcount\n")
        for p in range(1, 12):
            out_handle.write(f"{p}\t{stats.get(f'fail_pos_{p:02d}', 0)}\n")


class ThreeXNNKStrategy(BaseLibraryStrategy):
    name = "3xNNK"

    def __init__(self, wt33: str, max_mut: int):
        wt33 = wt33.upper()
        if len(wt33) != 33:
            raise ValueError("WT_nt_33 must be exactly 33 nt.")
        self.wt33 = wt33
        self.max_mut = int(max_mut)

    def validate(self, region_nt: str) -> Tuple[bool, str, Dict[str, Any]]:
        if len(region_nt) != 33:
            return False, "wrong_length", {"mutated_positions": []}

        mutated_positions: List[int] = []
        for pos in range(11):
            cod = region_nt[pos * 3:(pos + 1) * 3]
            ref = self.wt33[pos * 3:(pos + 1) * 3]
            if cod != ref:
                mutated_positions.append(pos + 1)

        if len(mutated_positions) > self.max_mut:
            return False, "too_many_mutations", {"mutated_positions": mutated_positions}

        for pos in mutated_positions:
            cod = region_nt[(pos - 1) * 3:pos * 3]
            if translate_codon(cod) == "*":
                return False, "stop_codon", {"mutated_positions": mutated_positions}
            if not matches_degenerate_codon(cod, "NNK"):
                return False, "mutation_not_NNK", {"mutated_positions": mutated_positions}

        return True, "ok", {"mutated_positions": mutated_positions}

    def update_stats(self, ok: bool, reason: str, details: Dict[str, Any], stats: Counter) -> None:
        mutated_positions = details.get("mutated_positions", [])
        if not isinstance(mutated_positions, list):
            return

        stats[f"mutations_{len(mutated_positions)}"] += 1
        for p in mutated_positions:
            stats[f"mut_pos_{p:02d}"] += 1

    def write_strategy_summary(self, stats: Counter, out_handle) -> None:
        out_handle.write(f"\n# strategy={self.name}\n")

        out_handle.write("mutations_per_read\tcount\n")
        max_seen = 0
        for k in stats:
            if k.startswith("mutations_"):
                try:
                    max_seen = max(max_seen, int(k.split("_")[1]))
                except Exception:
                    pass

        for m in range(0, max_seen + 1):
            out_handle.write(f"{m}\t{stats.get(f'mutations_{m}', 0)}\n")

        out_handle.write("\nmutated_position_1based\tcount\n")
        for p in range(1, 12):
            out_handle.write(f"{p}\t{stats.get(f'mut_pos_{p:02d}', 0)}\n")


def build_strategy(cfg: dict) -> BaseLibraryStrategy:
    mode = str(cfg["library_mode"]).strip()

    if mode == "combinatorial":
        deg_codons = list(cfg.get("degenerate_codons", DEFAULT_DEG_CODONS))
        return CombinatorialStrategy(deg_codons)

    if mode == "3xNNK":
        wt33 = str(cfg["WT_nt_33"]).upper()
        max_mut = int(cfg.get("max_mutated_codons", 3))
        return ThreeXNNKStrategy(wt33=wt33, max_mut=max_mut)

    raise ValueError(f"Unsupported library_mode: {mode}")


# ---------------- merged-output processing ----------------
def process_fastqs(
    input_paths: List[Path],
    out_dir: Path,
    prefix: str,
    anchor: str,
    anchor_mm: int,
    region_len: int,
    dropN: bool,
    write_fail: bool,
    write_aa_tsv: bool,
    strategy: BaseLibraryStrategy,
) -> Counter[str]:
    """
    Process one or more FASTQs but write a SINGLE output set:
      {prefix}.filtered.PASS.fastq.gz
      {prefix}.filtered.FAIL.fastq.gz  (optional)
      {prefix}.filtered.PASS.aa.tsv
      {prefix}.filtered.summary.tsv

    IMPORTANT: extracted region is reverse-complemented before validate/translate
    to match DESIGN orientation (paper codon list).
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    pass_path = out_dir / f"{prefix}.filtered.PASS.fastq.gz"
    fail_path = out_dir / f"{prefix}.filtered.FAIL.fastq.gz"
    summary_path = out_dir / f"{prefix}.filtered.summary.tsv"
    aa_path = out_dir / f"{prefix}.filtered.PASS.aa.tsv"

    counts = Counter()
    strategy_stats = Counter()

    if region_len % 3 != 0:
        logging.warning(
            f"region_length_nt={region_len} is not divisible by 3. Translation will truncate to full codons."
        )

    aa_fh = None
    aa_writer = None
    if write_aa_tsv:
        aa_fh = open(aa_path, "w", newline="")
        aa_writer = csv.writer(aa_fh, delimiter="\t")
        aa_writer.writerow(["read_id", "nt_seq", "aa_seq", "mutated_positions"])

    with gzip.open(pass_path, "wt") as pass_h:
        fail_h = gzip.open(fail_path, "wt") if write_fail else None

        for input_path in input_paths:
            logging.info(f"Processing input: {input_path}")

            with gzip.open(input_path, "rt") as handle:
                for rec in SeqIO.parse(handle, "fastq"):
                    counts["total_reads"] += 1

                    seq = str(rec.seq).upper()
                    qual = rec.letter_annotations.get("phred_quality", [])

                    hit = find_anchor(seq, anchor, anchor_mm)
                    if hit is None:
                        counts["fail_anchor_not_found"] += 1
                        if fail_h:
                            SeqIO.write(rec, fail_h, "fastq")
                        continue

                    _, a_end = hit
                    region = seq[a_end:a_end + region_len]
                    region_q = qual[a_end:a_end + region_len]

                    if len(region) != region_len:
                        counts["fail_region_too_short"] += 1
                        if fail_h:
                            SeqIO.write(rec, fail_h, "fastq")
                        continue

                    # normalize to DESIGN orientation
                    region = revcomp(region)
                    if region_q:
                        region_q = list(reversed(region_q))

                    if dropN and ("N" in region):
                        counts["fail_contains_N"] += 1
                        if fail_h:
                            SeqIO.write(rec, fail_h, "fastq")
                        continue

                    ok, reason, details = strategy.validate(region)

                    if ok:
                        counts["pass"] += 1

                        # Write PASS FASTQ as region-only (DESIGN orientation)
                        rec.letter_annotations = {}
                        rec.seq = rec.seq.__class__(region)
                        rec.letter_annotations["phred_quality"] = list(region_q) if region_q else [40] * len(region)
                        SeqIO.write(rec, pass_h, "fastq")

                        if aa_writer is not None:
                            aa_seq = translate_region_nt(region)
                            muts = details.get("mutated_positions", [])
                            muts_str = ",".join(str(x) for x in muts) if isinstance(muts, list) else ""
                            aa_writer.writerow([rec.id, region, aa_seq, muts_str])

                    else:
                        counts[f"fail_{reason}"] += 1
                        if fail_h:
                            SeqIO.write(rec, fail_h, "fastq")

                    strategy.update_stats(ok, reason, details, strategy_stats)

        if fail_h:
            fail_h.close()

    if aa_fh is not None:
        aa_fh.close()

    with open(summary_path, "w") as s:
        s.write("# inputs\n")
        for p in input_paths:
            s.write(f"# {p}\n")
        s.write("\n")

        s.write("metric\tvalue\n")
        for k in sorted(counts.keys()):
            s.write(f"{k}\t{counts[k]}\n")

        strategy.write_strategy_summary(strategy_stats, s)

    logging.info("Done.")
    logging.info(f"  PASS: {pass_path}")
    if write_fail:
        logging.info(f"  FAIL: {fail_path}")
    if write_aa_tsv:
        logging.info(f"  AA TSV: {aa_path}")
    logging.info(f"  Summary: {summary_path}")

    return counts


def main():
    ap = argparse.ArgumentParser(
        description="Filter reads by library rules after anchor + translate passing reads."
    )
    ap.add_argument("--yaml_config", required=True)
    ap.add_argument("--sample_id", required=True, help="Run a single sample by sample_id from samples_tsv")
    args = ap.parse_args()

    cfg = parse_yaml(args.yaml_config)

    log_path = setup_logging(cfg)
    logger.info("Logging to: %s", log_path)
    logger.info("Loaded config: %s", cfg.get("_config_path", args.yaml_config))

    prefix, input_fastqs = resolve_inputs_for_sample(cfg, sample_id=str(args.sample_id))
    logging.info(f"Resolved sample_id={args.sample_id} -> output_prefix={prefix}")
    for p in input_fastqs:
        logging.info(f"  input_fastq: {p}")

    anchor = str(cfg["anchor_sequence"])
    anchor_mm = int(cfg["anchor_max_mismatches"])
    region_len = int(cfg["region_length_nt"])
    dropN = bool(cfg["drop_if_contains_N"])
    write_fail = bool(cfg["write_fail_fastq"])
    write_aa_tsv = bool(cfg.get("write_aa_tsv", True))

    strategy = build_strategy(cfg)
    out_dir = Path(cfg["output_dir"])

    # prefix is "{sample_id}_{sample_name}" — split for metrics
    _parts = prefix.split("_", 1)
    _sample_id   = _parts[0]
    _sample_name = _parts[1] if len(_parts) > 1 else prefix

    counts = process_fastqs(
        input_paths=input_fastqs,
        out_dir=out_dir,
        prefix=prefix,
        anchor=anchor,
        anchor_mm=anchor_mm,
        region_len=region_len,
        dropN=dropN,
        write_fail=write_fail,
        write_aa_tsv=write_aa_tsv,
        strategy=strategy,
    )

    write_sample_metrics(
        metrics_path=out_dir / "per_sample_metrics.tsv",
        sample_id=_sample_id,
        sample_name=_sample_name,
        step="03_extraction",
        metrics={
            "reads_in_total":          counts.get("total_reads", 0),
            "reads_pass":              counts.get("pass", 0),
            "reads_lost_no_anchor":    counts.get("fail_anchor_not_found", 0),
            "reads_lost_too_short":    counts.get("fail_region_too_short", 0),
            "reads_lost_contains_N":   counts.get("fail_contains_N", 0),
            "reads_lost_wrong_length": counts.get("fail_wrong_length", 0),
            "reads_lost_stop_codon":   counts.get("fail_stop_codon", 0),
            "reads_lost_validation":   counts.get("fail_degenerate_mismatch", 0)
                                     + counts.get("fail_too_many_mutations", 0)
                                     + counts.get("fail_mutation_not_NNK", 0),
        },
    )


if __name__ == "__main__":
    main()
