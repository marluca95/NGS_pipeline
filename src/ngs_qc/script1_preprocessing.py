import os
import gzip
import subprocess
import pandas as pd
from pathlib import Path
import argparse
import yaml
import time


def parse_arguments_from_yaml(yaml_file: str) -> dict:
    with open(yaml_file, "r") as f:
        return yaml.safe_load(f)


def wait_for_file(filepath, timeout=600, interval=5):
    start_time = time.time()
    while not os.path.exists(filepath):
        if time.time() - start_time > timeout:
            print(f"❌ Timeout: {filepath} not found within {timeout} seconds.")
            return False
        time.sleep(interval)
    return True


def trimming_single_end(input_fastq_gz: str, output_fastq_gz: str,
                       bbmap_dir: str, qtrim: str, quality_threshold: int, min_length: int):
    """
    Single-end trimming with bbduk.
    """
    cmd = (
        f"{bbmap_dir}/bbduk.sh -Xmx64g "
        f"in={input_fastq_gz} out={output_fastq_gz} "
        f"qtrim={qtrim} trimq={quality_threshold} minlength={min_length}"
    )
    print(f"Running: {cmd}")
    try:
        subprocess.run(cmd, check=True, shell=True)
        return {"status": "success"}
    except subprocess.CalledProcessError as e:
        return {"status": "error", "message": str(e)}


def combine_gz_fastqs(input_files, combined_out):
    """
    Concatenate gzipped FASTQ files safely by binary concatenation.
    This is valid for gzip: concatenated gzip members are readable by gzip tools and SeqIO.
    """
    combined_out = str(combined_out)
    with open(combined_out, "wb") as out_handle:
        for fp in input_files:
            with open(fp, "rb") as in_handle:
                out_handle.write(in_handle.read())
    return {"status": "success", "combined": combined_out}


def load_sample_sheet(sample_sheet_path: str) -> pd.DataFrame:
    """
    Auto-detect delimiter (tab or comma). Your pasted format is TSV.
    """
    # Try tab first, fallback to comma
    try:
        df = pd.read_csv(sample_sheet_path, sep="\t")
        if df.shape[1] == 1:
            df = pd.read_csv(sample_sheet_path)  # comma
    except Exception:
        df = pd.read_csv(sample_sheet_path)
    required = {"sample_id", "fastq", "sample_name"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Sample sheet missing required columns: {missing}. Found: {list(df.columns)}")
    return df

def safe_name(s: str) -> str:
    """
    Make a string safe for folder names.
    """
    return (
        str(s)
        .strip()
        .replace(" ", "_")
        .replace("/", "_")
    )



def main():
    parser = argparse.ArgumentParser(description="Trim single-end FASTQs and combine per sample_id.")
    parser.add_argument("--yaml_config", required=True, help="Path to global YAML config")
    parser.add_argument("--sample_sheet", required=True, help="Path to sample sheet (TSV/CSV)")
    parser.add_argument("--sample_id", default=None,
                        help="Optional: process only one sample_id (useful for array jobs)")
    args = parser.parse_args()

    config = parse_arguments_from_yaml(args.yaml_config)

    output_dir = Path(config["output_dir"])
    bbmap_dir = str(config["bbmap_dir"])

    trim_enabled = bool(config.get("trim", True))
    qtrim = config.get("qtrim", "rl")
    quality_threshold = int(config.get("quality_threshold", 15))
    min_length = int(config.get("min_length", 50))

    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_sample_sheet(args.sample_sheet)

    # Optionally filter to one sample_id (HPC array use)
    if args.sample_id is not None:
        df = df[df["sample_id"].astype(str) == str(args.sample_id)]
        if df.empty:
            print(f"No rows found for sample_id={args.sample_id}")
            return

    # Group by sample_id
    for sample_id, g in df.groupby("sample_id", sort=False):
        sample_id = str(sample_id)

        # get unique sample_name for this sample_id
        sample_names = g["sample_name"].unique()
        if len(sample_names) != 1:
            raise ValueError(
                f"sample_id={sample_id} has multiple sample_name values: {sample_names}"
            )

        sample_name = safe_name(sample_names[0])

        sample_label = f"{sample_id}_{sample_name}"

        print(f"\n=== Processing sample: {sample_label} ({len(g)} FASTQs) ===")

        sample_out = output_dir / sample_label
        sample_out.mkdir(parents=True, exist_ok=True)

        trimmed_files = []

        # Trim every listed FASTQ for this sample_id
        for _, row in g.iterrows():
            in_fp = Path(str(row["fastq"]))
            if not in_fp.exists():
                print(f"⚠ Missing input FASTQ: {in_fp}  (skipping)")
                continue

            # output name keeps original filename but adds ".trimmed" before .fastq.gz
            name = in_fp.name
            if name.endswith(".fastq.gz"):
                out_name = name.replace(".fastq.gz", ".trimmed.fastq.gz")
            elif name.endswith(".fq.gz"):
                out_name = name.replace(".fq.gz", ".trimmed.fq.gz")
            else:
                out_name = name + ".trimmed.fastq.gz"

            out_fp = sample_out / out_name

            if trim_enabled:
                res = trimming_single_end(
                    input_fastq_gz=str(in_fp),
                    output_fastq_gz=str(out_fp),
                    bbmap_dir=bbmap_dir,
                    qtrim=qtrim,
                    quality_threshold=quality_threshold,
                    min_length=min_length,
                )
                if res["status"] != "success":
                    print(f"❌ Trimming failed for {in_fp}: {res.get('message', '')}")
                    continue

                if not wait_for_file(str(out_fp), timeout=600, interval=5):
                    print(f"❌ Trimmed file never appeared: {out_fp}")
                    continue

                if out_fp.stat().st_size == 0:
                    print(f"❌ Trimmed file is empty: {out_fp}")
                    continue

                trimmed_files.append(str(out_fp))
            else:
                # If trimming disabled, just pass through
                trimmed_files.append(str(in_fp))

        if not trimmed_files:
            print(f"❌ No trimmed files produced for sample_id={sample_id}. Skipping combine.")
            continue

        # Combine ALL trimmed files for this sample_id into one file
        combined_out = sample_out / f"{sample_label}_combined_trimmed.fastq.gz"
        print(f"Combining {len(trimmed_files)} trimmed FASTQs -> {combined_out}")

        comb = combine_gz_fastqs(trimmed_files, combined_out)
        if comb["status"] != "success":
            print(f"❌ Combine failed for sample_id={sample_id}: {comb}")
            continue

        if not wait_for_file(str(combined_out), timeout=600, interval=5):
            print(f"❌ Combined file never appeared: {combined_out}")
            continue

        print(f"✅ Done sample_id={sample_id}")
        print(f"   Trimmed files in: {sample_out}")
        print(f"   Combined output : {combined_out}")


if __name__ == "__main__":
    main()
