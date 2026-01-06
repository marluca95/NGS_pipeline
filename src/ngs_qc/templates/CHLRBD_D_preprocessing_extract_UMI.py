import gzip
import pandas as pd
import os
import re
import subprocess
from collections import Counter
import Levenshtein as lev
from Bio import SeqIO
from Bio.Seq import Seq
import shutil
from pathlib import Path
import argparse
import yaml
import time
import sys
import logger as l #Need file logger.py
#from concurrent.futures import ProcessPoolExecutor
import gc

def parse_arguments_from_yaml(yaml_file):
    # Load YAML configuration file
    with open(yaml_file, 'r') as file:
        config = yaml.safe_load(file)  # This should return a dictionary
    return config

def read_fastq_gz(filepath, umi_length=15, config= "test", max_primer_mismatch=1):
    records = []

    with gzip.open(filepath, 'rt') as f:  # 'rt' = read text (auto-decodes)
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            if max_primer_mismatch == 1:
                extraction, extr_qual = extract_seq(seq, qual, fwd_adapter = config['fwd_adapter'], rc_adapter = config['rc_adapter'])
            else:
                extraction, extr_qual = extract_approx_seq(seq, qual, fwd_adapter = config['fwd_adapter'], rc_adapter = config['rc_adapter'], max_primer_mismatch = max_primer_mismatch)
            # if len(extraction) == config['expected_length']:
            #     translation = Seq((extraction)).translate(to_stop=True)
            #     if len(translation) != (config['expected_length']/3):
            #         translation = None
            # else:
            #     translation = None
            len_extraction = len(extraction)

            umi = seq[:umi_length]
            records.append({
                "header": header,
                "umi": umi,
                "sequence": seq,
                "quality": qual,
                "extraction": extraction,
                "extracted_qual": extr_qual,
                #"translation": translation
                "len_extraction": len_extraction
            })

    return records

def extract_seq(seq, qual, fwd_adapter, rc_adapter):
    sequence = str(seq)  # Store the sequence
    quality = str(qual)
    # Initialize extracted_region to None to handle cases where adapters are not found
    extracted_region = "NNNN"
    extracted_qual = "XXXX"

    # Find the position of the forward (fwd_adapter) and the reverse (rc_adapter) primer
    fwd_pos = sequence.find(fwd_adapter)
    rc_pos = sequence.find(rc_adapter)

    # If both primers are found in the sequence, extract the region between them
    if fwd_pos != -1 and rc_pos != -1 and fwd_pos < rc_pos:
        # Extract the region between the primers
        extracted_region = sequence[fwd_pos + len(fwd_adapter):rc_pos]
        extracted_qual = quality[fwd_pos + len(fwd_adapter):rc_pos]
    return extracted_region, extracted_qual

def find_approx(sequence: str, pattern: str, max_mismatch: int) -> int:
    """
    Return the index of the first occurrence of `pattern` in `sequence`
    allowing up to `max_mismatch` mismatches, or -1 if none found.
    """
    plen, slen = len(pattern), len(sequence)
    for i in range(slen - plen + 1):
        window = sequence[i:i+plen]
        # count mismatches
        mismatches = sum(1 for a, b in zip(window, pattern) if a != b)
        if mismatches <= max_mismatch:
            return i
    return -1

def extract_approx_seq(seq: str,
                qual: str,
                fwd_adapter: str,
                rc_adapter: str,
                max_mismatch: int = 0) -> str:
    """
    Extract the subsequence between fwd_adapter and rc_adapter,
    each allowed up to `max_mismatch` mismatches.
    Returns 'NNNN' if adapters not found (or in wrong order).
    """
    sequence = str(seq)
    quality = str(qual)
    # default if not found
    extracted_region = "NNNN"
    extracted_qual = "XXXX"

    # fuzzy-find both adapters
    fwd_pos = find_approx(sequence, fwd_adapter, max_mismatch)
    rc_pos  = find_approx(sequence, rc_adapter,  max_mismatch)

    # ensure both found and forward is before reverse
    if 0 <= fwd_pos < rc_pos:
        start = fwd_pos + len(fwd_adapter)
        extracted_region = sequence[start:rc_pos]
        extracted_qual = quality[start:rc_pos]

    return extracted_region, extracted_qual

def main():
    # Parse the YAML file path as an argument
    parser = argparse.ArgumentParser(description="Process and analyze FASTQ files.")
    parser.add_argument('--yaml_config', type=str, required=True, help="Path to the YAML configuration file")
    parser.add_argument('--samplename', type=str, required=True, help="Name of the sample")
    parser.add_argument('--fastq', type=str, required=True, help="Path to the fastq")
    parser.add_argument('--max_primer_mismatch', type=str, required=True, help="Path to the fastq")

    args = parser.parse_args()

    # Parse the YAML configuration file
    config = parse_arguments_from_yaml(args.yaml_config)
    print("Script has started processing...")
    print(f"Input directory: {config['input_dir']}")
    print(f"Output directory: {config['output_dir']}")
    print(f"BBMap directory: {config['bbmap_dir']}")
    print(f"Quality threshold: {config['quality_threshold']}")
    print(f"Trimming option: {config['qtrim']}")
    print(f"Minimum length: {config['min_length']}")
    print(f"Expected sequence length: {config['expected_length']}")
    print(f"Sequence threshold: {config['seq_threshold']}")
    print(f"Forward adapter: {config['fwd_adapter']}")
    print(f"Reverse adapter: {config['rc_adapter']}")
    print(f"Wild-type sequence: {config['wt']}")

    records = read_fastq_gz(args.fastq, umi_length=15, config = config, max_primer_mismatch = args.max_primer_mismatch)
    df = pd.DataFrame(records)
    df.to_csv(f"UMI_extracted_{args.samplename}.csv", index=False)

if __name__ == "__main__":
    main()