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


import pandas as pd

# Parse the YAML file path as an argument
parser = argparse.ArgumentParser(description="Process and analyze FASTQ files.")
parser.add_argument('--yaml_config', type=str, required=True, help="Path to the YAML configuration file")
parser.add_argument('--samplename', type=str, required=True, help="Name of the sample")
parser.add_argument('--csv', type=str, required=True, help="Path to the csv")

args = parser.parse_args()

def parse_arguments_from_yaml(yaml_file):
    # Load YAML configuration file
    with open(yaml_file, 'r') as file:
        config = yaml.safe_load(file)  # This should return a dictionary
    return config

config = parse_arguments_from_yaml(args.yaml_config)
print("Reading YAML file...")
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


# Load the CSV file
df = pd.read_csv(args.csv)

# Filter for sequences only which have expected length

df = df[df["len_extraction"] == config['expected_length']]

# Keep only rows where 'translation' is not None and not NaN
#df = df[df["translation"].notna() & (df["translation"] != "None")]

print(f"✅ Filtered {len(df)} rows with valid translations.")

# Function to compute consensus sequence
def consensus(sequences):
    if not sequences:
        return ""
    length = len(sequences[0])
    consensus_seq = ""
    for i in range(length):
        col = [s[i] for s in sequences if len(s) > i]
        most_common = Counter(col).most_common(1)
        consensus_seq += most_common[0][0] if most_common else "X"
    return consensus_seq

# Function to compute average Hamming distance to the consensus
def average_divergence(sequences, consensus_seq):
    if not sequences or not consensus_seq:
        return 0.0
    def hamming(s1, s2):
        return sum(a != b for a, b in zip(s1, s2))
    distances = [hamming(seq, consensus_seq) for seq in sequences if len(seq) == len(consensus_seq)]
    return sum(distances) / len(distances) if distances else 0.0

# Group by UMI and compute consensus, read count, divergence, and list original sequences
umi_summary = []
for umi, group in df.groupby("umi"):
    extractions = group["extraction"].tolist()
    extracted_qual = group["extracted_qual"].tolist()
    read_count = len(extractions)
    cons = consensus(extractions)
    div = average_divergence(extractions, cons)
    cons_aa = Seq((cons)).translate(to_stop=True)
    all_sequences = ";".join(group["sequence"])
    umi_summary.append({
        "umi": umi,
        "cons_aa": cons_aa,
        "cons_seq": cons,
        "read_count": read_count,
        "avg_hamming_divergence": div,
        "original_sequences": all_sequences
    })

umi_summary = pd.DataFrame(umi_summary)
# Save the result
umi_summary.to_csv(f"umi_consensus_summary_{args.samplename}.csv", index=False)

# Count unique UMIs and reads per translation
grouped = umi_summary.groupby("cons_aa").agg({
    "umi": lambda x: ";".join(x),
    "read_count": "sum",
    "avg_hamming_divergence": lambda x: ";".join(map(str, x)),
    "original_sequences": lambda x: ";".join(x),
    "cons_seq": lambda x: ";".join(x)
}).reset_index()

# Add UMI count
grouped["umi_count"] = grouped["umi"].apply(lambda x: len(x.split(";")))

# Reorder columns
grouped = grouped[[
    "cons_aa", "umi_count", "read_count",
    "umi", "cons_seq", "original_sequences", "avg_hamming_divergence"
]]

# Filter out Translations with wrong length

grouped["len_cons"] = grouped["cons_aa"].apply(lambda x: len(x))
grouped = grouped[grouped["len_cons"] == (config['expected_length']/3)]

# Save result
output_file = f"umi_translation_{args.samplename}.csv"
grouped.to_csv(output_file, index=False)


print(f"✅ UMI consensus computation complete. Output saved to 'umi_translation_{args.samplename}.csv'.")

