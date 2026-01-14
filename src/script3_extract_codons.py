#!/usr/bin/env python3
import gzip
import csv
import yaml
from pathlib import Path
from Bio import SeqIO, Seq
from fuzzysearch import find_near_matches
from ngs_qc.logger import Logger  # custom logger


def fuzzy_find(seq: str, pattern: str, max_dist: int = 1) -> int:
    """Return start of first approximate match to pattern in seq, else -1."""
    matches = find_near_matches(pattern, seq, max_l_dist=max_dist)
    return matches[0].start if matches else -1


def process_sample(
    consensus_fastq: Path,
    output_csv: Path,
    start_anchor: str,
    end_anchor: str,
    codons: int = 13,
    max_dist: int = 1,
    logger=None,
    debug_examples: int = 5,
) -> None:
    """Extract codon region between anchors, translate, and write to CSV with debug stats."""
    results = []

    total_reads = 0
    extracted_reads = 0

    # Debug counters
    start_found = 0
    end_found = 0
    both_found = 0
    end_before_start = 0
    region_too_short = 0
    empty_region = 0  # start and end found but region length == 0
    exceptions = 0

    # Optional: distribution tracking
    region_lengths = []
    start_positions = []
    end_positions = []

    # Keep a few examples of failures for logging
    example_logs = []

    with gzip.open(consensus_fastq, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            total_reads += 1
            seq = str(record.seq)

            try:
                start_pos = fuzzy_find(seq, start_anchor, max_dist=max_dist)
                end_pos = fuzzy_find(seq, end_anchor, max_dist=max_dist)

                if start_pos != -1:
                    start_found += 1
                    start_positions.append(start_pos)
                if end_pos != -1:
                    end_found += 1
                    end_positions.append(end_pos)

                if start_pos == -1 or end_pos == -1:
                    # record a couple examples to help debugging
                    if len(example_logs) < debug_examples:
                        example_logs.append(
                            f"[NO_ANCHOR] {record.id}: start_pos={start_pos}, end_pos={end_pos}"
                        )
                    continue

                both_found += 1

                # Important sanity check: end must be after start anchor
                start_end = start_pos + len(start_anchor)
                if end_pos <= start_end:
                    end_before_start += 1
                    if len(example_logs) < debug_examples:
                        snippet = seq[max(0, start_pos - 10): min(len(seq), end_pos + len(end_anchor) + 10)]
                        example_logs.append(
                            f"[END_BEFORE_START] {record.id}: start_pos={start_pos}, end_pos={end_pos}, "
                            f"start_end={start_end}, snippet='{snippet}'"
                        )
                    continue

                region = seq[start_end:end_pos]

                if len(region) == 0:
                    empty_region += 1
                    if len(example_logs) < debug_examples:
                        example_logs.append(
                            f"[EMPTY_REGION] {record.id}: start_pos={start_pos}, end_pos={end_pos}"
                        )
                    continue

                region_lengths.append(len(region))

                if len(region) < codons * 3:
                    region_too_short += 1
                    if len(example_logs) < debug_examples:
                        example_logs.append(
                            f"[REGION_TOO_SHORT] {record.id}: region_len={len(region)}, required={codons*3}, "
                            f"start_pos={start_pos}, end_pos={end_pos}"
                        )
                    continue

                region = region[: codons * 3]
                nt_seq = Seq.Seq(region)
                aa_seq = str(nt_seq.translate())

                results.append((record.id, str(nt_seq), aa_seq))
                extracted_reads += 1

            except Exception as e:
                exceptions += 1
                if len(example_logs) < debug_examples:
                    example_logs.append(f"[EXCEPTION] {record.id}: {repr(e)}")
                continue

    # Write CSV
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["umi_id", "nt_seq", "aa_seq"])
        writer.writerows(results)

    # Logging summary
    if logger:
        logger.log(f"Processed {consensus_fastq}")
        logger.log(f"Total reads: {total_reads}")
        logger.log(f"Start anchor found: {start_found}")
        logger.log(f"End anchor found: {end_found}")
        logger.log(f"Both anchors found: {both_found}")
        logger.log(f"End before start_end (invalid order/overlap): {end_before_start}")
        logger.log(f"Empty region (start/end adjacent): {empty_region}")
        logger.log(f"Region too short (<{codons*3} bp): {region_too_short}")
        logger.log(f"Extracted: {extracted_reads}")
        logger.log(f"Exceptions: {exceptions}")
        logger.log(f"Output written: {output_csv}")
        logger.log(f"Extraction rate: {extracted_reads/total_reads:.3%}")


        # Some quick distribution stats (only if we have data)
        if region_lengths:
            region_lengths_sorted = sorted(region_lengths)
            mid = len(region_lengths_sorted) // 2
            median = region_lengths_sorted[mid]
            logger.log(
                f"Region length stats: min={region_lengths_sorted[0]}, median={median}, max={region_lengths_sorted[-1]}"
            )

        if start_positions:
            sp = sorted(start_positions)
            logger.log(f"Start pos stats: min={sp[0]}, median={sp[len(sp)//2]}, max={sp[-1]}")
        if end_positions:
            ep = sorted(end_positions)
            logger.log(f"End pos stats: min={ep[0]}, median={ep[len(ep)//2]}, max={ep[-1]}")

        # Log a few examples
        if example_logs:
            logger.log("Example failures (up to a few):")
            for line in example_logs:
                logger.log(line)



def read_samples_tsv(samples_tsv: Path):
    """
    Read TSV sample sheet.
    Requires columns: sample_id, sample_name
    Returns unique (sample_id, sample_name) pairs (deduplicated across lanes/runs).
    """
    seen = set()
    with open(samples_tsv, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Sample sheet TSV is empty or missing header")

        required = {"sample_id", "sample_name"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Sample sheet TSV missing required columns: {sorted(missing)}")

        for row in reader:
            sample_id = (row.get("sample_id") or "").strip()
            sample_name = (row.get("sample_name") or "").strip()
            if not sample_id or not sample_name:
                continue

            key = (sample_id, sample_name)
            if key in seen:
                continue
            seen.add(key)
            yield sample_id, sample_name


def find_consensus_fastq(sample_dir: Path, consensus_mode: str, config: dict) -> Path | None:
    """
    Find the FASTQ for a sample based on consensus_mode using patterns from YAML.

    YAML keys used:
      no_consensus: "*_combined_final.fastq.gz"
      consensus: "umi/*_consensus_{consensus_mode}.fastq.gz"
    """
    if consensus_mode == "no_consensus":
        pattern = config["no_consensus"]
    else:
        pattern = config["consensus"].format(consensus_mode=consensus_mode)

    matches = list(sample_dir.glob(pattern))
    return matches[0] if matches else None


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Extract codon sequences and translate to amino acids (multi-sample via TSV)"
    )
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    parser.add_argument(
        "--consensus_mode",
        required=False,
        default="no_consensus",
        help="Specify consensus mode (topk_x, min_x, hard, or no_consensus)",
    )
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    processed_base_dir = Path(config["processed_base_dir"])
    extracted_base_dir = Path(config["extracted_base_dir"])
    samples_tsv = Path(config["sample_sheet"])
    logs_dir = Path(config["logs_dir"])

    extracted_base_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    start_anchor = config["start_anchor"]
    end_anchor = config["end_anchor"]
    codons = int(config.get("codons", 13))
    max_dist = int(config.get("max_mismatches", 1))

    for sample_id, sample_name in read_samples_tsv(samples_tsv):
        sample_folder = f"{sample_id}_{sample_name}"
        sample_dir = processed_base_dir / sample_folder

        if not sample_dir.exists():
            print(f"[WARN] Skipping {sample_folder}: sample directory not found: {sample_dir}")
            continue

        consensus_fastq = find_consensus_fastq(sample_dir, args.consensus_mode, config)
        if consensus_fastq is None or not consensus_fastq.exists():
            print(f"[WARN] Skipping {sample_folder}: no FASTQ found for mode={args.consensus_mode} in {sample_dir}")
            continue

        output_csv = extracted_base_dir / f"{sample_folder}_{args.consensus_mode}_aa.csv"
        log_file = logs_dir / f"{sample_folder}_{args.consensus_mode}_extract.log"
        logger = Logger(str(logs_dir), log_file.name)

        logger.log(f"Extracting codons for sample {sample_folder} ({args.consensus_mode})")

        process_sample(
            consensus_fastq=consensus_fastq,
            output_csv=output_csv,
            start_anchor=start_anchor,
            end_anchor=end_anchor,
            codons=codons,
            max_dist=max_dist,
            logger=logger,
        )


if __name__ == "__main__":
    main()
