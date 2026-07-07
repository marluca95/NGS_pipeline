import argparse
import pandas as pd
from pathlib import Path

# Default paths
DATA_DIR = Path(
    "/cluster/project/reddy/katja/NGS_pipeline/data/P3408_LUCA-TCRA3/04_variant_labeling/"
    "19_05_2026_minlenght191_Q30/combined_with_tcr_population/"
)
MAPPING_PATH = Path("/cluster/project/reddy/katja/NGS_pipeline/data/peptide_mapping.csv")
OUTPUT_PATH = Path(
    "/cluster/project/reddy/katja/NGS_pipeline/data/"
    "A3_Q30_combined_tcr_peptide_label_all_with_tcr_population.csv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a final TCR/peptide/label CSV from script4 variant-labeling outputs."
    )
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR, help="Directory containing *.variant_labeling.csv files.")
    parser.add_argument("--mapping-path", type=Path, default=MAPPING_PATH, help="Peptide name-to-sequence mapping CSV.")
    parser.add_argument("--output-path", type=Path, default=OUTPUT_PATH, help="Output final CSV path.")
    return parser.parse_args()


def create_tcr_peptide_specificity_csv(data_dir: Path, mapping_path: Path, output_path: Path) -> None:
    # Load mapping
    mapping_df = pd.read_csv(mapping_path)
    peptide_map = dict(zip(mapping_df["Name"], mapping_df["AA Sequence"]))

    # Read all CSV files in directory
    dfs = []
    for csv_file in data_dir.glob("*.variant_labeling.csv"):
        dfs.append(pd.read_csv(csv_file))

    if not dfs:
        raise ValueError(f"No *.variant_labeling.csv files found in {data_dir}")

    combined_df = pd.concat(dfs, ignore_index=True)
    if "tcr_population" not in combined_df.columns:
        raise ValueError(
            "Input variant-labeling CSVs are missing required column 'tcr_population'. "
            "Rerun script4_variant_labeling.py with the updated code before creating the final CSV."
        )

    # Keep only specificity values 0 and 1
    filtered_df = combined_df[combined_df["specificity"].isin([0, 1])].copy()

    # Build output schema
    filtered_df["tcr"] = filtered_df["aa_seq"]
    filtered_df["peptide"] = filtered_df["peptide"].map(peptide_map)

    result_df = filtered_df[["tcr", "peptide", "specificity", "tcr_population"]]
    result_df = result_df.rename(columns={"specificity": "label"})

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_path, index=False)

    missing_mappings = int(result_df["peptide"].isna().sum())

    print(f"Read {len(dfs)} files from: {data_dir}")
    print(f"Combined rows: {len(combined_df)}")
    print(f"Rows kept (specificity in [0, 1]): {len(result_df)}")
    print(f"Rows with missing peptide mapping: {missing_mappings}")
    print(f"Saved: {output_path}")


def main() -> None:
    args = parse_args()
    create_tcr_peptide_specificity_csv(
        data_dir=args.data_dir,
        mapping_path=args.mapping_path,
        output_path=args.output_path,
    )


if __name__ == "__main__":
    main()
