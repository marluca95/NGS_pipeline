import os
import re
import gzip
import subprocess
import pandas as pd
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


def read_fastq_gz_in_chunks(file_path, chunk_size=100000):
    """
    Reads a compressed FASTQ file in chunks and yields batches of sequences.
    """
    with gzip.open(file_path, "rt") as handle:
        batch = []
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
            batch.append(record)
            if (i + 1) % chunk_size == 0:
                yield batch
                batch = []
        if batch:
            yield batch


def find_r1_files(input_dir, identifier):
    # Regex: look for the identifier AND '_R1' before extension
    pattern = re.compile(rf"{re.escape(identifier)}.*_R1.*\.fastq\.gz$")
    
    r1_files = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if pattern.search(file):
                r1_files.append(os.path.join(root, file))
    
    return r1_files

# Function to process a FASTQ file pair (R1 and R2)
def process_fastq_file(output_dir, R1, R2, config):
    # Extract sample name by cleaning the filename (removing parts of the name indicating lane, read direction, etc.)
    # sample_name = re.sub(r'(_L[0-9]+_R[12](_[0-9]+)?(_.*)?\.fastq\.gz$)', '', os.path.basename(R1))   # Clean filename
    # sample_name = re.sub(r'_S\d+_L\d+_R[12]_001\.fastq\.gz$', '', os.path.basename(R1))
    # Alternatively, remove the _1 or _2 (before DS190B1) from the name for grouping similar samples
    # sample_name = re.sub(r'(_[12])', '', sample_name)  # Remove _1 or _2
    # sample_name = re.sub(r'_[12]_(?=[A-Za-z0-9]+$)', '_', sample_name)
    sample_name = re.sub(r'^.*?_([A-Za-z0-9]+)_S\d+_L\d+_R[12]_001\.fastq.gz$', r'\1', os.path.basename(R1))


    # Extract lane information from the filename (e.g., 'L001' for lane 1)
    lane = re.search(r'L[0-9]+', os.path.basename(R1))
    lane = lane.group() if lane else "unknown_lane"

    print(f"Processing sample: {sample_name}, lane: {lane}")

    # Check if the input files exist (R1 and R2). If either is missing, print an error and return None
    if not os.path.exists(R1) or not os.path.exists(R2):
        print(f"Error: Missing input files {R1} or {R2}")
        return None

    # Create an empty dictionary to hold the output filenames
    output_files = {}

    # Conditionally add output files based on the arguments passed
    if config['repair']:  # If the repair step is requested
        output_files['repaired_R1'] = os.path.join(output_dir, f"{sample_name}_{lane}_repaired_R1.fastq.gz")
        output_files['repaired_R2'] = os.path.join(output_dir, f"{sample_name}_{lane}_repaired_R2.fastq.gz")

    if config['trim']:  # If the trimming step is requested
        output_files['trimmed_R1'] = os.path.join(output_dir, f"{sample_name}_{lane}_trimmed_R1.fastq.gz")
        output_files['trimmed_R2'] = os.path.join(output_dir, f"{sample_name}_{lane}_trimmed_R2.fastq.gz")

    output_files['merged_out'] = os.path.join(output_dir, f"{sample_name}_{lane}_merged.fastq.gz")

    if config['filter']:  # If the filtering step is requested
        output_files['filtered_merged_out'] = os.path.join(output_dir, f"{sample_name}_{lane}_filtered_merged.fastq.gz")

    # Always generate a final combined output file (reads from all the lanes for each sample)
    output_files['combined'] = os.path.join(output_dir, f"{sample_name}_combined_final.fastq.gz")
    # Define the output filename for the final extracted aa sequences
    output_files['extracted'] = os.path.join(output_dir, f"{sample_name}_extracted.fastq.gz")

    return R1, R2, output_files


# Function to repair a FASTQ file pair (R1 and R2)
def repairing(R1, R2, output_files, bbmap_dir):
    # Repair command
    repair_cmd = f"{bbmap_dir}/repair.sh -Xmx64g in1={R1} in2={R2} out1={output_files['repaired_R1']} out2={output_files['repaired_R2']} repair"
    try:
        # Run the repair command and check if it completes without errors
        subprocess.run(repair_cmd, check=True, shell=True)
        # Return a success message if repair is completed without errors
        return {'status': 'success', 'message': 'Repair completed successfully'}
    except subprocess.CalledProcessError as e:
        # Catch errors during the subprocess execution and return an error message
        return {'status': 'error', 'message': f"Repair failed: {e}"}


# Function to trim the reads from a FASTQ file pair (R1 and R2)
def trimming(R1, R2, output_files, bbmap_dir, qtrim, quality_threshold, min_length, repair_done):
    # Decide which files to use based on the repair_done status
    if repair_done['status'] == 'success':
        # If repair is successful, use the repaired files (if they exist)
        input_R1 = output_files.get('repaired_R1')
        input_R2 = output_files.get('repaired_R2')
    elif repair_done['status'] == 'skipped':
        # If repair is skipped, use the repaired files (if they exist)
        input_R1 = output_files.get('repaired_R1')
        input_R2 = output_files.get('repaired_R2')
    else:
        # Raw files otherwise
        input_R1, input_R2 = R1, R2

    print(f"Using files: input_R1={input_R1}, input_R2={input_R2}")

    # Check if input_R1 and input_R2 are correctly assigned
    if input_R1 is None or input_R2 is None:
        raise ValueError("Input files for trimming are missing. Ensure that repair or raw files are available.")

    # Trimming command
    trimming_cmd = f"{bbmap_dir}/bbduk.sh -Xmx64g in1={input_R1} in2={input_R2} out1={output_files['trimmed_R1']} out2={output_files['trimmed_R2']} qtrim={qtrim} trimq={quality_threshold} minlength={min_length}"

    try:
        # Run the trimming command and check if it completes without errors
        subprocess.run(trimming_cmd, check=True, shell=True)
        # Return a success message if trimming is completed without errors
        return {'status': 'success', 'message': 'Trimming completed successfully'}
    # Catch errors during the subprocess execution and return an error message
    except subprocess.CalledProcessError as e:
        return {'status': 'error', 'message': f"Trimming failed: {e}"}


# Function to merge paired-end reads (R1 and R2)
def merging(R1, R2, output_files, bbmap_dir, trimming_done, repair_done):
    # Decide which files to use based on the repair_done and trimming_done status
    if trimming_done['status'] == 'success':
        # If trimming was successful or skipped, use the trimmed files
        input_R1 = output_files.get('trimmed_R1')
        input_R2 = output_files.get('trimmed_R2')
    elif repair_done['status'] == 'success':
        # If repair was successful or skipped, use the repaired files
        input_R1 = output_files.get('repaired_R1')
        input_R2 = output_files.get('repaired_R2')
    else:
        # If neither repair nor trimming was done, use the raw R1 and R2 files
        input_R1, input_R2 = R1, R2  # Use raw files if repair or trimming is not done

    # Prepare the output file path for merged output from the dictionary
    merged_out = output_files['merged_out']

    print(f"🔄 Running merging step for {R1} and {R2}...")

    # Merge command
    #merge_command = f"{bbmap_dir}/bbmerge.sh -Xmx64g in1={input_R1} in2={input_R2} out={merged_out}"
    merge_command = f"{bbmap_dir}/bbmerge.sh -Xmx64g in1={input_R2} in2={input_R1} out={merged_out}" # RK for AVITI swapped R1 and R2

    print("🔄 Running 'sync' to flush file system...")
    subprocess.run(["sync"])  # Ensure all writes are flushed

    try:
        # Run the merging command and capture its output
        subprocess.run(merge_command, check=True, shell=True, capture_output=True, text=True)
        # Print and return success message if merging completes without errors
        print(f"Merge completed for {input_R1} and {input_R2}. Merged reads saved to {merged_out}.")
        return {'status': 'success', 'message': f'Merge completed, output saved to {merged_out}'}
    except subprocess.CalledProcessError as e:
        # Print and return error message if the merging process fails
        print(f"Error during merging: {e}")
        return {'status': 'error', 'message': f"Merge failed: {e}"}


# Function to filter reads based on minimum length
def filtering(merged_out, output_files, min_length, bbmap_dir):
    # Access the filtered output file path
    filtered_merged_out = output_files['filtered_merged_out']

    # Filtering command
    filter_command = f"{bbmap_dir}/bbduk.sh -Xmx64g in1={merged_out} out={filtered_merged_out} minlength={min_length}"
    print(f"Running command: {filter_command}")

    try:
        # Run the filtering command and capture its output
        subprocess.run(filter_command, check=True, shell=True)
        # Print and return success message if filtering completes without errors
        print(f"Filtering completed for {merged_out}. Filtered reads saved to {filtered_merged_out}.")
        return {'status': 'success', 'message': f'Filtering completed, output saved to {filtered_merged_out}'}
    except subprocess.CalledProcessError as e:
        # Print and return error message if the merging process fails
        print(f"Error during filtering: {e}")
        return {'status': 'error', 'message': f"Filtering failed: {e}"}


# Combine filtered or merged files from different lanes into one single FASTQ file per sample
def combine_filtered_merged_files(output_dir, all_merged_files, output_files, filtering_done):
    print("Combining lanes for each sample into a single FASTQ")

    # List to store paths of final merged FASTQ files for each sample
    final_fastq_paths = []

    # Initialize a list to store sample file names
    sample_files = []

    # Determine whether filtering on merged files was applied or skipped
    if filtering_done['status'] == 'success':
        # If filtering was applied successfully, use the filtered files
        print("Filtering was applied successfully.")
        sample_files = sorted(set(
            [re.sub(r'_L[0-9]+_filtered_merged.fastq.gz', '', os.path.basename(f))
             for f in all_merged_files if f.endswith("_filtered_merged.fastq.gz")]
        ))
        print(f"Filtered sample files detected: {sample_files}")

    elif filtering_done['status'] == 'skipped':
        # If filtering was skipped, use the merged files without filtering
        print("Filtering was skipped.")
        sample_files = sorted(set(
            [re.sub(r'_L[0-9]+_merged.fastq.gz', '', os.path.basename(f))
             for f in all_merged_files if f.endswith("_merged.fastq.gz")]
        ))
        print(f"Filtered sample files detected: {sample_files}")

    else:

        sample_files = sorted(set(
            [re.sub(r'_L[0-9]+_merged.fastq.gz', '', os.path.basename(f))
             for f in all_merged_files if f.endswith("_merged.fastq.gz")]
        ))
        print(f"Merged sample files detected: {sample_files}")

    # Check if no files were found for the sample combination process
    if not sample_files:
        print("No sample files found for combining!")
        return {'status': 'error', 'message': "No sample files found for combining."}

    # Combine lane-specific files into one final file per sample
    for sample in sample_files:
        final_merged_fastq = os.path.join(output_dir, f"{sample}_combined_final.fastq.gz")
        lane_files = []  # List to hold lane-specific files for the current sample

        # Collect the lane-specific files for the current sample
        for f in os.listdir(output_dir):
            if filtering_done['status'] == 'success':
                # If filtering was successful, match files with '_L' followed by the lane number
                match = re.search(r'_L(\d+)_filtered_merged.fastq.gz', f)
            else:
                # If filtering wasn't done, match merged files based on lane numbers
                match = re.search(r'_L(\d+)_merged.fastq.gz', f)

            # If the file belongs to the current sample and a lane, add it to the list
            if match and f.startswith(f"{sample}_L"):
                lane_files.append(os.path.join(output_dir, f))

        # Ensure that lane-specific files exist before attempting to combine
        if lane_files:
            try:
                # Open the final combined file in binary write mode
                with open(final_merged_fastq, 'wb') as outfile:
                    # Iterate through the lane files and write their content into the final file
                    for fname in lane_files:
                        with open(fname, 'rb') as infile:
                            outfile.write(infile.read())

                print(f"Created combined file for {sample} at {final_merged_fastq}")

                # Remove the intermediate lane-specific files after combining
                for fname in lane_files:
                    os.remove(fname)
                    print(f"Removed {fname}")

                # Add the path of the final merged file to the list
                final_fastq_paths.append(final_merged_fastq)

            except Exception as e:
                # If there is an error during the combining process, return the error message
                print(f"Error during combining: {e}")
                return {'status': 'error', 'message': f"Error during combining: {e}"}

    # If combined files were created, update output_files and return success status
    if final_fastq_paths:
        output_files['combined'] = final_fastq_paths  # Add the combined files to output_files
        return {'status': 'success', 'output_files': output_files}
    else:
        # If no combined files were created, return error status
        return {'status': 'error', 'message': "No lane-specific files found for merging."}


# Function to calculate amino acid counts/frequencies and aa sequence
def calculate_aa_info(sequences):
    # Initialize a counter to keep track of amino acid counts
    aa_counts = Counter()
    # List to store the translated amino acid sequences
    aa_sequences = []

    try:
        for seq in sequences:
            # Translate DNA sequence to amino acid sequence
            aa_sequence = Seq(seq).translate(to_stop=True)  # Translate DNA to AA sequence
            aa_counts.update(aa_sequence)  # Update the amino acid counts
            aa_sequences.append(str(aa_sequence))  # Append the amino acid sequence as a string to the list

        # Calculate the total count of amino acids
        total_count = sum(aa_counts.values())
        # Calculate the frequency of each amino acid by dividing its count by the total count
        aa_frequencies = {aa: count / total_count for aa, count in aa_counts.items()}
        # Return counts, frequencies, and sequences inside a status object
        return {'status': 'success', 'aa_counts': aa_counts, 'aa_frequencies': aa_frequencies,
                'aa_sequences': aa_sequences}

    except Exception as e:
        # If an error occurs, catch it and return an error message with status
        return {'status': 'error', 'message': f"Error during AA calculation: {e}"}


# Function to calculate Levenshtein distance to a wild type sequence for amino acid sequences
def calculate_levenshtein_distance(aa_sequences, wt_aa_sequence):
    lev_distances = []  # List to store the Levenshtein distances for each sequence

    try:
        for aa_seq in aa_sequences:
            # Calculate Levenshtein distance between the amino acid sequence and the wild type amino acid sequence
            dist = lev.distance(aa_seq, wt_aa_sequence)
            lev_distances.append(dist)

        # Return the status and the list of Levenshtein distances
        return {'status': 'success', 'lev_distances': lev_distances}

    except Exception as e:
        # Return error message if anything goes wrong
        return {'status': 'error', 'message': f"Error during Levenshtein distance calculation: {e}"}


def extract_seq_process(final_merged_fastq_paths, fwd_adapter, rc_adapter, wt, seq_threshold, output_files):
    sequence_counts = {}  # Dictionary to hold counts of each unique extracted DNA sequence

    total_sequences_before_trimming = 0  # Initialize total sequence processing variables
    unique_sequences = set()  # Set to track unique sequences
    filtered_sequences_before_after = []  # List to track sequences before and after trimming
    wt_count = 0 # counts occurences of WT seqeunce (RK)

    # Process each file in the list (or single file)
    file_paths = final_merged_fastq_paths if isinstance(final_merged_fastq_paths, list) else [final_merged_fastq_paths]

    try:
        # Extract the sample name from the file path
        sample_name = file_paths[0].split("/")[-1].replace("_combined_final.fastq.gz",
                                                           "")  # Extract sample name from file path
        # Extract only the part of the sample name you want (e.g., Sample5_NonBinding)
        # match = re.search(r'(_(HC[A-Za-z0-9]+|DS\d+[A-Za-z0-9]*))', sample_name)
        # match = re.search(r'(_(HC[A-Za-z0-9]+|DS[A-Za-z0-9]+))', sample_name)

        # match = re.search(r'_([0-9]{2}[A-Z])$', sample_name)
        match = re.search(r'_(DS(?:Pool)?[A-Za-z0-9]+|HC(?:Pool)?[A-Za-z0-9]+)', sample_name, re.IGNORECASE)


        if match:
            # short_sample_name = match.group(0).strip('_')
            short_sample_name = match.group(1)
        else:
            # short_sample_name = "Unknown_Sample"  # Default name if regex doesn't match
            short_sample_name = sample_name # RK addition, if there is not DS or Pool then use sample name as short sample name

        print(f"Processing files for sample: {sample_name}")

        # Define output file name for extracted sequences
        extracted_fastq_file = output_files['extracted']
        with gzip.open(extracted_fastq_file, 'wt') as output_handle:  # Open in write-text mode for FASTQ
            # Process each input file
            for file_path in file_paths:
                print(f"Processing file: {file_path} (Sample: {sample_name})")
                with gzip.open(file_path, 'rt') as f:  # Open input file in read-text mode
                    reads = SeqIO.parse(f, 'fastq')  # Parse the FASTQ file

                    # Iterate through each sequence in the file
                    for record in reads:
                        sequence = str(record.seq)  # Store the sequence
                        total_sequences_before_trimming += 1  # Increment total sequence counter
                        unique_sequences.add(sequence)  # Add sequence to unique sequences set

                        # Initialize extracted_region to None to handle cases where adapters are not found
                        extracted_region = None

                        # Find the position of the forward (fwd_adapter) and the reverse (rc_adapter) primer
                        fwd_pos = sequence.find(fwd_adapter)
                        rc_pos = sequence.find(rc_adapter)

                        # If both primers are found in the sequence, extract the region between them
                        if fwd_pos != -1 and rc_pos != -1 and fwd_pos < rc_pos:
                            # Extract the region between the primers
                            extracted_region = sequence[fwd_pos + len(fwd_adapter):rc_pos]

                        # Only proceed if extracted_region is not None and matches wt length (or expected length)
                        if extracted_region and len(extracted_region) == len(wt) and extracted_region != wt:
                            # Translate the extracted region to amino acids
                            aa_seq = str(Seq(extracted_region).translate(to_stop=True))
                            print(f"AA Seq: {aa_seq}")
                            # Filter out sequences with stop codons length mismatch
                            if "*" not in aa_seq:
                                # Count the occurrences of the extracted region (DNA sequence)
                                sequence_counts[extracted_region] = sequence_counts.get(extracted_region, 0) + 1

                                # Save the before and after trimming sequences in a dictionary or list for later use
                                filtered_sequences_before_after.append({
                                    'before': sequence,
                                    'after': extracted_region,
                                    'before_length': len(sequence),
                                    'after_length': len(extracted_region)
                                })
                        if extracted_region == wt: # RK
                            wt_count = wt_count + 1


        # Counters for sequences that match and don't match the expected length
        matched_length_count = 0
        mismatched_length_count = 0

        # Check each extracted region in sequence_counts
        for extracted_region in sequence_counts:
            if len(extracted_region) != len(wt):  # Length check
                mismatched_length_count += 1  # Increment counter for mismatches
            else:
                matched_length_count += 1  # Increment counter for matches

        # Print the counts
        print(f"Sequences that matched the expected length: {matched_length_count}")
        print(f"Sequences that did not match the expected length: {mismatched_length_count}")
        print(f"Number of WT sequences found: {wt_count}") #RK

        # Print the number of sequences and counts before applying seq_threshold filtering
        print(f"Number of sequences before applying seq_threshold: {len(sequence_counts)}")
        print(f"Counts before filtering:")
        # for seq, count in sequence_counts.items():  # Use the unfiltered sequence_counts
        # print(f"Sequence: {seq}, Count: {count}")

        filtered_sequences = []
        filtered_counts = []  # List to hold counts for the filtered sequences
        for seq, count in sequence_counts.items():
            if count >= seq_threshold:  # Only keep sequences that appear more than the threshold
                filtered_sequences.append(seq)
                filtered_counts.append(count)

        # Print how many sequences remain after applying seq_threshold
        print(f"Number of sequences after applying seq_threshold filtering: {len(filtered_sequences)}")

        # Convert the filtered DNA sequences to amino acid sequences
        aa_sequences = set(str(Seq(seq).translate(to_stop=True)) for seq in filtered_sequences)

        # Print how many unique amino acid sequences remain after conversion
        print(f"Number of unique amino acid sequences after conversion: {len(aa_sequences)}")

        # Get the length of the wild-type amino acid sequence
        wt_aa_length = len(str(Seq(wt).translate(to_stop=True)))

        # Step 2: Counters for sequences that match and don't match the wild-type amino acid length
        matched_aa_count = 0
        mismatched_aa_count = 0
        mismatched_aa_sequences = []  # Store sequences that don't match the expected AA length

        # Filter the unique amino acid sequences to keep only those with the same length as the wild-type
        aa_sequences_same_length_as_wt = []
        aa_counts_same_length_as_wt = []
        for seq, count in zip(aa_sequences, filtered_counts):
            if len(seq) == wt_aa_length:
                aa_sequences_same_length_as_wt.append(seq)
                aa_counts_same_length_as_wt.append(count)
            else:
                mismatched_aa_count += 1  # Increment mismatch count
                mismatched_aa_sequences.append(seq)  # Store the mismatched sequence for further investigation

        # Print how many sequences match the expected amino acid length
        print(f"Number of sequences that match the wild-type amino acid length: {len(aa_sequences_same_length_as_wt)}")
        print(f"Number of sequences that do not match the wild-type amino acid length: {mismatched_aa_count}")

        # Before calculating Levenshtein distance, print the number of sequences being processed
        print(f"Number of sequences being processed for Levenshtein distance: {len(aa_sequences_same_length_as_wt)}")

        # Calculate Levenshtein distances for each filtered amino acid sequence against the wild-type (wt) amino acid sequence
        lev_distances = calculate_levenshtein_distance(aa_sequences_same_length_as_wt,
                                                       str(Seq(wt).translate(to_stop=True)))

        # Print information about total sequences and unique sequences before and after filtering
        print(f"Total sequences before trimming: {total_sequences_before_trimming}")
        print(f"Unique sequences before trimming: {len(unique_sequences)}")  # Raw unique sequences
        print(f"Unique sequences after filtering: {len(filtered_sequences)}")  # Unique after threshold filtering

        # Calculate the number of reads that pass the threshold
        filtered_reads = sum(count for seq, count in sequence_counts.items() if count >= seq_threshold)
        print(f"Total reads after trimming and filtering: {filtered_reads}")

        # Create a DataFrame with the extracted and filtered sequences and relevant information
        info_df = pd.DataFrame({
            'sample_ID': [short_sample_name] * len(aa_sequences_same_length_as_wt),  # Sample name for each sequence
            'aa_sequence': aa_sequences_same_length_as_wt,  # Unique amino acid sequences that match the wt length
            'seq_counts': aa_counts_same_length_as_wt,  # Counts for each sequence with matching length
            'LV': lev_distances['lev_distances'],
            # Levenshtein distance for  each amino acid sequence with matching length
        })

        # Print the first few rows of the resulting DataFrame
        print(info_df.head())

        # Set Pandas display options to show all columns and rows
        pd.set_option('display.max_columns', None)  # Show all columns
        pd.set_option('display.max_rows', None)  # Show all rows (be careful with large DataFrames)
        pd.set_option('display.width', None)  # Prevent truncation of rows in the terminal
        pd.set_option('display.max_colwidth', None)  # Show full content in cells

        # Return a dictionary with the processing results
        return {'status': 'success',
                'info_df': info_df,
                'extracted_fastq_file': extracted_fastq_file,
                'total_sequences': filtered_reads,
                'short_sample_name': short_sample_name,
                }

    except Exception as e:
        # If an error occurs, return a status of error and the exception message
        return {'status': 'error', 'message': f"Error during sequence processing: {e}"}


def organize_files_at_end(output_dir, output_files, sample_id, info_df):
    """
    At the end of the processing, move all files to their final locations.
    - Intermediate files go to lane folders.
    - Final files (including CSVs) go to the sample folder.
    """
    # Create sample folder if it doesn't exist
    sample_folder = Path(output_dir) / sample_id
    sample_folder.mkdir(parents=True, exist_ok=True)

    # Collect lane numbers from all files in the output directory
    lane_numbers = set()
    for file_name in os.listdir(output_dir):
        match = re.search(r'_L(\d+)_', file_name)  # This regex will capture lane numbers like L001, L002
        if match:
            lane_numbers.add(match.group(1))  # Add lane number to the set

    # Move final FASTQ files like 'combined' and 'extracted' to the sample folder
    for key in ['combined', 'extracted']:
        if key in output_files:
            final_output_file = output_files[key]

            # If it's a list, iterate through all files in the list
            if isinstance(final_output_file, list):
                for file in final_output_file:
                    final_output_path = sample_folder / f"{sample_id}_{key}.fastq.gz"
                    if os.path.exists(file):
                        shutil.move(file, final_output_path)
                        print(f"Moved {file} to {final_output_path}")
            else:
                final_output_path = sample_folder / f"{sample_id}_{key}.fastq.gz"
                if os.path.exists(final_output_file):
                    shutil.move(final_output_file, final_output_path)
                    print(f"Moved {final_output_file} to {final_output_path}")

    # Move CSV files to the sample folder (info_df, read_df)
    if info_df is not None:
        info_csv_path = sample_folder / f"{sample_id}_info_df.csv"
        info_df.to_csv(info_csv_path, index=False)
        print(f"Moved {sample_id}_info_df.csv to {info_csv_path}")

    # Move all files in the output directory, check for lane-specific files
    for file_name in os.listdir(output_dir):
        file_path = os.path.join(output_dir, file_name)

        # Wait for file to appear before attempting to move
        if not wait_for_file(file_path, timeout=60, interval=2):
            print(f"⚠ Warning: {file_path} did not appear within the expected time, skipping.")
            continue

        # Check if the file contains lane number in the name
        match = re.search(r'_L(\d+)_', file_name)

        if match:
            lane_number_from_output = match.group(1)

            # If the file has a lane number, move it to the lane folder
            if lane_number_from_output in lane_numbers:
                lane_folder = sample_folder / f"lane_{lane_number_from_output}"
                lane_folder.mkdir(parents=True, exist_ok=True)

                lane_output_path = lane_folder / f"{sample_id}_{file_name}"
                shutil.move(file_path, lane_output_path)
                print(f"Moved {file_path} to {lane_output_path}")


def wait_for_file(filepath, timeout=600, interval=5):
    """Wait for a file (or multiple files if a list is given) to appear within a timeout period."""
    import time
    import os

    if isinstance(filepath, list):  # Handle case where filepath is a list
        for file in filepath:
            wait_for_file(file, timeout, interval)  # Recursively wait for each file
        return  # Exit function after all files are checked

    start_time = time.time()
    while not os.path.exists(filepath):
        if time.time() - start_time > timeout:
            print(f"❌ Timeout: {filepath} not found within {timeout} seconds.")
            return False
        time.sleep(interval)
    return True


def main():
    # Parse the YAML file path as an argument
    parser = argparse.ArgumentParser(description="Process and analyze FASTQ files.")
    parser.add_argument('--yaml_config', type=str, required=True, help="Path to the YAML configuration file")
    parser.add_argument('--sample', type=str, required=True, help="Sample ID to process")

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
    logger = l.Logger('./logs/', f'logger_array_id_{args.sample}.txt')
    # Define directories based on the parsed config
    input_dir = config['input_dir']
    output_dir = config['output_dir']
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists
    bbmap_dir = config['bbmap_dir']

    # After parsing YAML config
    sample_id = args.sample  # This is e.g., '59055' from SLURM array
    sample_output_dir = Path(config['output_dir']) / sample_id
    sample_output_dir.mkdir(parents=True, exist_ok=True)


    # Initialize a list to store the read counts for each file
    all_read_counts = []
    # List to store merged/filtered files for combining later
    all_merged_files = []
    # Initialize the output_files dictionary before processing
    output_files = {}

    # Initialise the total variables to accumulate the counts
    total_initial_R1_reads = 0
    total_repaired_R1_reads = 0
    total_trimmed_R1_reads = 0
    total_merged_reads = 0
    total_filtered_reads = 0

    # Extract the R1 and R2 files corresponding to the sample ID
    print(args.sample)
    R1_files = find_r1_files(input_dir, args.sample)

    samples = {}

    # Process each R1 file and find the corresponding R2
    logger.log('Starting iteration over R1 files')
    for R1 in R1_files:
        logger.log(f'{R1}')
        # Extract sample_id and lane_number using regex from R1 filename
        # match = re.search(r'_([A-Za-z0-9]+)_S\d+_L(\d+)_R1_001.fastq.gz', R1)
        match = re.search(r'^(.*)_S\d+_L(\d+)_R1_001\.fastq\.gz$', R1)
        if match:
            # sample_id = match.group(1)  # Extracted sample identifier (e.g., DS190B1)
            # sample_id = re.sub(r'_[12]_(?=[A-Za-z0-9]+$)', '_', match.group(1))
            raw_sample_id = re.sub(r'_[12]_(?=[A-Za-z0-9]+$)', '_', match.group(1))
            sample_id = re.sub(r'^(GFB-\d+_)?(HTV3VDRX5|000000000-DTH4R)_?', '', raw_sample_id)

            lane_number = match.group(2)  # Extracted lane number (e.g., L001)
            logger.log(f"Detected sample: {sample_id}, Lane: {lane_number}")

            # Construct full path for R1
            R1_path = os.path.join(input_dir, R1)

            # Find the corresponding R2 file based on the same sample_id and lane_number
            R2_file = R1.replace('_R1', '_R2')  # Replace '_R1' with '_R2' to get the corresponding R2 file
            R2_path = os.path.join(input_dir, R2_file)

            # Ensure the R2 file exists
            if os.path.exists(R2_path):
                # If the sample_id and lane_number do not exist in the dictionary, create them
                if sample_id not in samples:
                    samples[sample_id] = {}
                if lane_number not in samples[sample_id]:
                    samples[sample_id][lane_number] = {'R1': [], 'R2': []}

                # Add the full paths of R1 and R2 to the dictionary
                samples[sample_id][lane_number]['R1'].append(R1_path)
                samples[sample_id][lane_number]['R2'].append(R2_path)
                logger.log(f"Matched R1: {R1_path}, R2: {R2_path}")
            else:
                print(f"Warning: No corresponding R2 file found for {R1} in lane {lane_number}")
                logger.log(f"Warning: No corresponding R2 file found for {R1} in lane {lane_number}")

    logger.log('Starting iteration over samples')
    # Process each sample's R1 and R2 files by lane
    for sample_id, lanes in samples.items():
        logger.log(f"Processing sample: {sample_id}")
        print(f"\nProcessing sample: {sample_id}")

        filter_done = {}

        # Loop over each lane for this sample
        for lane_number, files in lanes.items():
            print(f"\nProcessing lane {lane_number} for sample {sample_id}")
            logger.log(f"Processing lane {lane_number} for sample {sample_id}")

            # Loop over all R1 and R2 files for this lane
            for R1, R2 in zip(files['R1'], files['R2']):
                print(f"Processing R1: {R1} and R2: {R2} for lane {lane_number}")

                # Initialize read counts for this sample
                initial_R1_reads = 0
                repaired_R1_reads = 0
                trimmed_R1_reads = 0
                merged_reads = 0
                filtered_reads = 0

                # Process the R1 and R2 files
                R1, R2, output_files = process_fastq_file(sample_output_dir, R1, R2, config)

                # Verify that R1 and R2 exist after processing
                if not os.path.exists(R1) or not os.path.exists(R2):
                    logger.log(f"⚠ ERROR: Processed R1 or R2 file missing for {R1}, skipping this lane.")
                    continue

                # Verify that R1 and R2 are not empty
                if os.path.getsize(R1) == 0 or os.path.getsize(R2) == 0:
                    logger.log(f"⚠ ERROR: Processed R1 or R2 file is empty for {R1}, skipping this lane.")
                    continue

                if output_files:
                    trim_done = {}

                    # Count initial reads
                    initial_R1_reads = sum(len(batch) for batch in read_fastq_gz_in_chunks(R1))
                    all_read_counts.append({'file': R1, 'initial_reads': initial_R1_reads})
                    total_initial_R1_reads += initial_R1_reads

                    logger.log('Starting repair')
                    logger.log(f"Starting repair for {R1} and {R2}")
                    # Repair step (optional based on arguments)
                    if config['repair']:
                        print("Running repair...")
                        repair_done = repairing(R1, R2, output_files, bbmap_dir)
                        if repair_done['status'] != 'success':
                            print(f"Repair failed for {R1}. Error: {repair_done['message']}")
                            continue  # Skip the current file if repair fails
                        print(f"Repair done for {R1}")

                        wait_for_file(output_files.get('repaired_R1', ''))
                        wait_for_file(output_files.get('repaired_R2', ''))

                    else:
                        repair_done = {'status': 'skipped', 'message': 'Repair skipped'}

                    # Check if repaired files exist before proceeding to trimming
                    if 'repaired_R1' not in output_files or not os.path.exists(output_files['repaired_R1']) or \
                            'repaired_R2' not in output_files or not os.path.exists(output_files['repaired_R2']):
                        logger.log(f"⚠ ERROR: Repaired files missing for {R1}, skipping trimming for this lane.")
                        continue
                    if os.path.getsize(output_files['repaired_R1']) == 0 or os.path.getsize(
                            output_files['repaired_R2']) == 0:
                        logger.log(f"⚠ ERROR: Repaired FASTQ files are empty for {R1}, skipping trimming.")
                        continue
                    logger.log(
                        f"✅ Repaired files found, proceeding to trimming: {output_files['repaired_R1']} and {output_files['repaired_R2']}")

                    repaired_R1_path = output_files.get('repaired_R1', '')
                    # Count reads after repair (only if repair is done)
                    if repair_done['status'] == 'success':
                        repaired_R1_reads = sum(len(batch) for batch in read_fastq_gz_in_chunks(repaired_R1_path))
                        all_read_counts.append({'file': repaired_R1_path, 'repaired_reads': repaired_R1_reads})
                        total_repaired_R1_reads += repaired_R1_reads

                    else:
                        repaired_R1_reads = '--'

                    # Immediately delete repaired files after counting
                    del (repaired_R1_path)

                    logger.log(
                        f"Starting trim for repaired files: {output_files.get('repaired_R1', '')} and {output_files.get('repaired_R2', '')}")
                    # Trimming step (optional based on arguments)
                    if config['trim']:
                        print("Running trimming...")
                        if repair_done['status'] == 'success' or repair_done['message'] == 'Repair skipped':
                            trim_done = trimming(R1, R2, output_files, bbmap_dir, qtrim=config['qtrim'],
                                                 quality_threshold=config['quality_threshold'],
                                                 min_length=config['min_length'],
                                                 repair_done=repair_done)

                        if trim_done['status'] != 'success':
                            print(f"Trimming failed for {R1}. Error: {trim_done['message']}")
                            continue  # Skip the current file if trimming fails
                        wait_for_file(output_files.get('trimmed_R1', ''))
                        wait_for_file(output_files.get('trimmed_R2', ''))

                        print(f"Trimming done for {R1}")  # Success case

                        wait_for_file(output_files.get('trimmed_R1', ''))
                        wait_for_file(output_files.get('trimmed_R2', ''))

                        if not os.path.exists(output_files['trimmed_R1']) or not os.path.exists(
                                output_files['trimmed_R2']):
                            logger.log(f"⚠ Warning: Trimmed FASTQ files are missing for {R1}, skipping merging")
                            continue
                        if os.path.getsize(output_files['trimmed_R1']) == 0 or os.path.getsize(
                                output_files['trimmed_R2']) == 0:
                            logger.log(f"⚠ Warning: Trimmed FASTQ files are empty for {R1}, skipping merging")
                            continue
                        if os.path.getsize(output_files['trimmed_R1']) == 0 or os.path.getsize(
                                output_files['trimmed_R2']) == 0:
                            logger.log(f"⚠ Warning: Trimmed FASTQ files are empty for {R1}, skipping merging")
                            continue

                    else:
                        trim_done = {'status': 'skipped', 'message': 'Trimming skipped'}

                    trimmed_R1_path = output_files.get('trimmed_R1', '')
                    # Step 3: Count reads after trimming (only if trimming is done)
                    if trim_done['status'] == 'success':
                        trimmed_R1_reads = sum(len(batch) for batch in read_fastq_gz_in_chunks(trimmed_R1_path))
                        all_read_counts.append({'file': trimmed_R1_path, 'trimmed_reads': trimmed_R1_reads})
                        total_trimmed_R1_reads += trimmed_R1_reads
                    else:
                        trimmed_R1_reads = '--'

                    # Immediately delete repaired files after counting
                    del (trimmed_R1_path)

                    logger.log(
                        f"Starting merge for trimming files: {output_files.get('trimmed_R1', '')} and {output_files.get('trimmed_R2', '')}")

                    # Merging step (combine R1 and R2 files)
                    print("Running merging...")
                    merge_done = merging(R1, R2, output_files, bbmap_dir, trimming_done=trim_done,
                                         repair_done=repair_done)
                    logger.log(f"Merging result: {merge_done}")
                    if merge_done['status'] != 'success':
                        print(f"Merging failed for {R1}. Error: {merge_done['message']}")
                        continue  # Skip the current file if merging fails
                    print(f"Successfully merged the trimmed files: {merge_done['message']}")  # Success case

                    # Verify if the merged file exists immediately after merge
                    merged_file = output_files.get('merged_out', '')

                    if not merged_file or not os.path.exists(merged_file):
                        logger.log(f"❌ ERROR: Merge step completed, but merged file is missing: {merged_file}")
                        continue

                    if os.path.getsize(merged_file) == 0:
                        logger.log(f"❌ ERROR: Merged file is empty: {merged_file}. Skipping this lane.")
                        continue

                    logger.log(f"✅ Merged file successfully created: {merged_file}, proceeding to filtering.")

                    # Wait for the file, but now we're sure it's actually being written
                    wait_for_file(merged_file)

                    if 'merged_out' not in output_files or not os.path.exists(output_files['merged_out']):
                        logger.log(f"⚠ Warning: Merged file missing for {R1}, skipping filtering")
                        continue

                    merged_reads = sum(len(batch) for batch in read_fastq_gz_in_chunks(merged_file))
                    all_read_counts.append({'file': merged_file, 'merged_reads': merged_reads})
                    total_merged_reads += merged_reads

                    logger.log(f"Starting filtering for merged files: {output_files.get('merged_out', '')}")
                    # Filtering step (optional based on arguments)
                    if config['filter']:
                        print("Running filtering...")
                        filter_done = filtering(output_files['merged_out'], output_files, config['min_length'],
                                                bbmap_dir)
                        if filter_done['status'] != 'success':
                            print(f"Filtering failed for {R1}. Error: {filter_done['message']}")
                            continue  # Skip the current file if filtering fails
                        print(f"Filtering completed for {R1}")  # Success case
                        # Wait for filtered file before proceeding
                        wait_for_file(output_files.get('filtered_merged_out', ''))
                    else:
                        filter_done = {'status': 'skipped', 'message': 'Filtering skipped'}

                    filtered_file = output_files.get('filtered_merged_out', '')

                    # Step 5: Count reads after filtering (only if filtering is done)
                    if 'filtered_merged_out' in output_files and os.path.exists(output_files['filtered_merged_out']):
                        filtered_reads = sum(len(batch) for batch in read_fastq_gz_in_chunks(filtered_file))
                        all_read_counts.append({'file': filtered_file, 'filtered_reads': filtered_reads})
                        total_filtered_reads += filtered_reads

                    else:
                        filtered_reads = "--"
                        print("Filtering step was skipped or failed, filtered merged file not found.")

                        # Immediately delete repaired files after counting
                        del (filtered_file)

                    logger.log('Appending')
                    # After all samples are processed, combine the files for final merging
                    all_merged_files.append(
                        output_files['filtered_merged_out'] if 'filtered_merged_out' in output_files else output_files[
                            'merged_out'])

        logger.log(f"Starting combine for sample {sample_id}")
        # Now combine all processed files
        print("Running final combining of all processed files...")
        combine_result = combine_filtered_merged_files(sample_output_dir, all_merged_files, output_files,
                                                       filtering_done=filter_done)

        if combine_result['status'] != 'success':
            print(f"Error during combining: {combine_result['message']}")
            logger.log(f"combine done for sample {sample_id}")
        else:
            print(f"Successfully combined all files.")  # success case

        wait_for_file(output_files['combined'])

        # Track the reads after combining
        combined_reads = 0
        combined_files = combine_result.get('output_files', {}).get('combined', [])
        for file in combined_files:
            if os.path.exists(file):
                for batch in read_fastq_gz_in_chunks(file):
                    combined_reads += len(batch)

        print(f"Total reads in combined files: {combined_reads}")

        logger.log(f"combined files found: {output_files['combined']}")

        # Perform sequence extraction using the combined files
        if combined_files:
            logger.log("Starting extraction")
            extraction_result = extract_seq_process(combined_files, config['fwd_adapter'], config['rc_adapter'],
                                                    config['wt'],
                                                    config['seq_threshold'], output_files=output_files)

            if extraction_result['status'] != 'success':
                print(f"Error during sequence extraction: {extraction_result['message']}")
                logger.log("Extraction okk")
            else:
                info_df = extraction_result['info_df']
                final_filtered_count = extraction_result['total_sequences']
                short_sample_name = extraction_result['short_sample_name']

                # Save the info_df for this sample to a CSV file
                info_df.to_csv(f"{sample_output_dir}/{short_sample_name}_info_df.csv", index=False)
                print(f"Info DataFrame saved for {short_sample_name} as {short_sample_name}_info_df.csv")

                # Sum of the 'seq_counts' column
                sum_seq_counts = info_df['seq_counts'].sum()
                print(f"Total extracted sequences: {sum_seq_counts}")

                # Save the read counts for this sample to a separate CSV file
                file_read_counts = [
                    short_sample_name,  # sample_ID
                    total_initial_R1_reads,  # initial_reads
                    total_repaired_R1_reads,  # reads_after_repair
                    total_trimmed_R1_reads,  # reads_after_trimming
                    total_merged_reads,  # reads_after_merge
                    total_filtered_reads,  # reads_after_filter
                    combined_reads,  # reads_after_combibing
                    final_filtered_count  # reads_after_final_filter
                ]
                all_read_counts.append(file_read_counts)

                # Create a DataFrame with one row per file and 6 columns
                read_columns = [
                    "sample_ID",
                    'initial_reads',
                    'reads_after_repair',
                    'reads_after_trimming',
                    'reads_after_merge',
                    'reads_after_filter',
                    'reads_after_combining',
                    'reads_after_final_filter'
                ]

                read_df = pd.DataFrame([file_read_counts], columns=read_columns)
                print(read_df.head())

                # Add the read_df to the final CSV
                output_filename = os.path.join(sample_output_dir, "final_read_counts.csv")

                # Append the read counts for this sample to the final CSV
                read_df.to_csv(output_filename, mode='a', header=(not os.path.exists(output_filename)), index=False)
                print(f"Appended read counts for {sample_id} to {output_filename}")

                # Debugging output for read counts
                # print(f"Number of reads initially: {total_initial_R1_reads}")
                # print(f"Number of reads after repair: {total_repaired_R1_reads}")
                # print(f"Number of reads after trimming: {total_trimmed_R1_reads}")
                # print(f"Number of reads after merging: {total_merged_reads}")
                # print(f"Number of reads after filtering: {total_filtered_reads}")
                # print(f"Number of reads after combining: {combined_reads}")
                # print(f"Number of reads after final filtering: {final_filtered_count}")

                organize_files_at_end(sample_output_dir, output_files, sample_id, info_df)

                del extraction_result, combined_files
                # del filtered_reads, initial_R1_reads, repaired_R1_reads, trimmed_R1_reads, merged_reads, final_filtered_count
                del repair_done, trim_done, merge_done, filter_done
                del short_sample_name  # Keep sample_name if still needed
                del info_df
                # del all_merged_files  # Clear large lists
                del R1, R2

                # Force garbage collection
                gc.collect()

                # Additional logic for final output and cleanup
                print("Processing complete.")


if __name__ == "__main__":
    main()