#!/usr/bin/env python3
import gzip
from Bio import SeqIO, Seq
from pathlib import Path
import csv
import yaml
from fuzzysearch import find_near_matches

from ngs_qc.logger import Logger   # <-- your custom logger


def fuzzy_find(seq, pattern, max_dist=1):
    """Return start of first approximate match to pattern in seq, else -1."""
    matches = find_near_matches(pattern, seq, max_l_dist=max_dist)
    return matches[0].start if matches else -1


def process_sample(consensus_fastq, output_csv, start_anchor, end_anchor,
                   codons=13, max_dist=1, logger=None):
    results = []
    total_reads, extracted_reads = 0, 0

    with gzip.open(consensus_fastq, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            total_reads += 1
            seq = str(record.seq)

            # Fuzzy match anchors
            start_pos = fuzzy_find(seq, start_anchor, max_dist=max_dist)
            end_pos = fuzzy_find(seq, end_anchor, max_dist=max_dist)

            if start_pos == -1 or end_pos == -1:
                continue

            # Extract region between anchors
            region = seq[start_pos + len(start_anchor):end_pos]
            if len(region) < codons * 3:
                continue

            # Take exactly codons*3 bp
            region = region[:codons * 3]

            # Translate
            nt_seq = Seq.Seq(region)
            aa_seq = str(nt_seq.translate())

            results.append((record.id, str(nt_seq), aa_seq))
            extracted_reads += 1

    # Write CSV
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["umi_id", "nt_seq", "aa_seq"])
        writer.writerows(results)

    if logger:
        logger.log(f"Processed {consensus_fastq}")
        logger.log(f"Total reads: {total_reads}, Extracted: {extracted_reads}")
        logger.log(f"Output written: {output_csv}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Extract codon sequences and translate to amino acids")
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    parser.add_argument("--consensus_mode", required=False, help="specify if using hard or soft consensus", default='hard')
    args = parser.parse_args()

    # Load config
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    pool = config["pool"]
    start_anchor = config["start_anchor"]
    end_anchor = config["end_anchor"]
    codons = int(config.get("codons", 13))
    max_dist = int(config.get("max_mismatches", 1))

    base_dir = Path(f"/cluster/project/reddy/marluca/NGS_pipeline/data/processed/{pool}")
    extracted_dir = Path(f"/cluster/project/reddy/marluca/NGS_pipeline/data/extracted/{pool}")
    extracted_dir.mkdir(parents=True, exist_ok=True)

    log_dir = Path("/cluster/project/reddy/marluca/NGS_pipeline/logs")
    log_dir.mkdir(parents=True, exist_ok=True)

    for sample_dir in base_dir.glob("GFB-*"):
        base_sample_name = sample_dir.name
        umi_dir = sample_dir / "umi"
        consensus_files = list(umi_dir.glob(f"*_consensus_{args.consensus_mode}.fastq.gz"))
        if not consensus_files:
            continue

        consensus_fastq = consensus_files[0]

        # safer way to strip the trailing known suffix:
        fname = consensus_fastq.name
        suffix = f"_consensus_{args.consensus_mode}.fastq.gz"
        if fname.endswith(suffix):
            sample_name = fname[:-len(suffix)]
        else:
            # fallback: split at first occurrence of "_consensus"
            sample_name = fname.split("_consensus", 1)[0]

        # sample_name is now e.g. "GFB-59053_A3t4xNNK"
        output_csv = extracted_dir / f"{sample_name}_{args.consensus_mode}_aa.csv"
        log_file = log_dir / f"{sample_name}_{args.consensus_mode}_extract.log"
        logger = Logger(str(log_dir), log_file.name)

        logger.log(f"Extracting codons for sample {sample_name}")
        process_sample(consensus_fastq, output_csv,
                       start_anchor, end_anchor,
                       codons=codons, max_dist=max_dist,
                       logger=logger)

if __name__ == "__main__":
    main()