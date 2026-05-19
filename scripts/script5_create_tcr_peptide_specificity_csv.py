import pandas as pd
from pathlib import Path

# Default paths
DATA_DIR = Path("/cluster/project/reddy/katja/NGS_pipeline/data/P3408_LUCA-TCRA3/variant_labeling")
MAPPING_PATH = Path("/cluster/project/reddy/katja/NGS_pipeline/data/peptide_mapping.csv")
OUTPUT_PATH = Path("/cluster/project/reddy/katja/NGS_pipeline/data/tcr_peptide_specificity_A3_all.csv")


def main() -> None:
    # Load mapping
    mapping_df = pd.read_csv(MAPPING_PATH)
    peptide_map = dict(zip(mapping_df["Name"], mapping_df["AA Sequence"]))

    # Read all CSV files in directory
    dfs = []
    for csv_file in DATA_DIR.glob("*.csv"):
        dfs.append(pd.read_csv(csv_file))

    if not dfs:
        raise ValueError(f"No CSV files found in {DATA_DIR}")

    combined_df = pd.concat(dfs, ignore_index=True)

    # Keep only specificity values 0 and 1
    filtered_df = combined_df[combined_df["specificity"].isin([0, 1])].copy()

    # Build output schema
    filtered_df["tcr"] = filtered_df["aa_seq"]
    filtered_df["peptide"] = filtered_df["peptide"].map(peptide_map)

    result_df = filtered_df[["tcr", "peptide", "specificity"]]
    result_df = result_df.rename(columns={"specificity": "label"})

    # Save
    result_df.to_csv(OUTPUT_PATH, index=False)

    missing_mappings = int(result_df["peptide"].isna().sum())

    print(f"Read {len(dfs)} files from: {DATA_DIR}")
    print(f"Combined rows: {len(combined_df)}")
    print(f"Rows kept (specificity in [0, 1]): {len(result_df)}")
    print(f"Rows with missing peptide mapping: {missing_mappings}")
    print(f"Saved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()