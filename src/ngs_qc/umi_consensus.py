import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict, Counter
import yaml
from Levenshtein import distance as lev
from pathlib import Path
from fuzzysearch import find_near_matches

def parse_config(config_file):
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)

def fuzzy_find_anchor(seq, anchor, max_dist=1):
    """
    Return the start position of the first approximate match to anchor (allow max_dist mismatches).
    Returns -1 if no match found.
    """
    matches = find_near_matches(anchor, seq, max_l_dist=max_dist)
    if matches:
        return matches[0].start
    return -1

def collapse_umis(reads_per_umi):
    """Return consensus sequence for a group of reads with the same UMI."""
    consensus_seq = []
    read_length = len(reads_per_umi[0])
    for i in range(read_length):
        bases = [read[i] for read in reads_per_umi]
        most_common_base, count = Counter(bases).most_common(1)[0]
        consensus_seq.append(most_common_base)
    return "".join(consensus_seq)

def main(config_file):
    config = parse_config(config_file)
    input_fastq = config["input_fastq"]
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)
    
    umi_len = config["umi_length"]
    anchor_seq = config["anchor_seq"]
    max_anchor_dist = config["max_anchor_mismatch"]
    min_reads_per_umi = config["min_reads_per_umi"]
    
    umi_groups = defaultdict(list)
    total_reads = 0
    discarded_short_umi = 0
    
    with gzip.open(input_fastq, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            total_reads += 1
            seq = str(record.seq)
            anchor_pos = fuzzy_find_anchor(seq[umi_len:], anchor_seq, max_dist=max_anchor_dist)
            if anchor_pos == -1:
                continue  # anchor not found
            umi = seq[:umi_len]
            extracted_seq = seq[umi_len:umi_len+anchor_pos]
            if len(umi) < umi_len:
                discarded_short_umi += 1
                continue
            umi_groups[umi].append(extracted_seq)
    
    # Write consensus sequences to FASTQ
    output_fastq = output_dir / "consensus_reads.fastq.gz"
    with gzip.open(output_fastq, "wt") as out_handle:
        for umi, seqs in umi_groups.items():
            if len(seqs) < min_reads_per_umi:
                continue
            consensus = collapse_umis(seqs)
            out_handle.write(f"@{umi}\n{consensus}\n+\n{'I'*len(consensus)}\n")
    
    # Log summary
    summary_file = output_dir / "umi_summary.txt"
    with open(summary_file, "w") as s:
        s.write(f"Total reads: {total_reads}\n")
        s.write(f"Total UMIs detected: {len(umi_groups)}\n")
        s.write(f"Discarded short UMIs: {discarded_short_umi}\n")
        s.write(f"Reads per UMI distribution:\n")
        counts = [len(v) for v in umi_groups.values()]
        for c in sorted(set(counts)):
            s.write(f"{c} reads: {counts.count(c)} UMIs\n")
    
    print(f"UMI consensus complete. Output FASTQ: {output_fastq}")
    print(f"Summary written to: {summary_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="UMI consensus script")
    parser.add_argument("--config", required=True, help="Path to YAML config")
    args = parser.parse_args()
    main(args.config)
