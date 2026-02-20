from pathlib import Path
from typing import Iterable, Optional, Union

import pandas as pd


def read_table_with_fallback(path: Union[str, Path]) -> pd.DataFrame:
    """
    Read a table, preferring TSV and falling back to CSV.
    Keeps behavior used across pipeline scripts.
    """
    try:
        df = pd.read_csv(path, sep="\t")
        if df.shape[1] == 1:
            df = pd.read_csv(path)
    except Exception:
        df = pd.read_csv(path)
    return df


def load_sample_sheet(
    path: Union[str, Path],
    required_cols: Iterable[str],
    dedupe_on: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Load sample sheet and validate required columns."""
    df = read_table_with_fallback(path)
    required = set(required_cols)
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Sample sheet missing required columns: {missing}. Found: {list(df.columns)}"
        )

    if dedupe_on is not None:
        df = df.drop_duplicates(subset=list(dedupe_on)).reset_index(drop=True)

    return df

