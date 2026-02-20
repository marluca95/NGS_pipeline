#!/usr/bin/env python3
import re
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

N_RE = re.compile(r"\|n=(\d+)")

def parse_n(read_id: str) -> int:
    m = N_RE.search(str(read_id))
    if not m:
        raise ValueError(f"Could not parse n from read_id: {read_id}")
    return int(m.group(1))

def weighted_counts_matrix(seqs, weights) -> pd.DataFrame:
    seqs = pd.Series(seqs).astype(str)
    weights = pd.Series(weights).astype(int)

    lengths = seqs.str.len().unique()
    if len(lengths) != 1:
        raise ValueError(f"Variable aa_seq lengths found: {sorted(map(int, lengths))}")
    L = int(lengths[0])

    alphabet = sorted(set("".join(seqs.values)))
    a2i = {c: i for i, c in enumerate(alphabet)}

    counts = np.zeros((L, len(alphabet)), dtype=np.int64)
    for seq, w in zip(seqs.values, weights.values):
        for pos, ch in enumerate(seq):
            idx = a2i.get(ch)
            if idx is not None:
                counts[pos, idx] += w

    return pd.DataFrame(counts, columns=alphabet)

def plot_logo(counts_df: pd.DataFrame, out_png: Path, title: str):
    probs = counts_df.div(counts_df.sum(axis=1), axis=0).fillna(0.0)

    L = probs.shape[0]
    fig_w = max(10, min(30, L * 0.35))
    plt.figure(figsize=(fig_w, 4))

    logo = logomaker.Logo(probs)

    # Apply classic (WebLogo-style) protein colour scheme
    logo.style_glyphs(color_scheme='classic')

    # Optional: cleaner look
    logo.ax.set_ylim(0, 1)
    logo.ax.spines[['top', 'right']].set_visible(False)

    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", required=True, help="Folder containing *.aa.tsv files")
    ap.add_argument("--out_dir", required=True, help="Output folder for PNG logos")
    args = ap.parse_args()

    input_dir = Path(args.input_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    aa_files = sorted(input_dir.glob("*.aa.tsv"))
    if not aa_files:
        raise SystemExit(f"No *.aa.tsv files found in: {input_dir}")

    print(f"Input:  {input_dir}")
    print(f"Output: {out_dir}")
    print(f"Found {len(aa_files)} files")

    for f in aa_files:
        df = pd.read_csv(f, sep="\t", dtype=str)

        # tolerate minor header casing differences
        cols = {c.lower(): c for c in df.columns}
        if "read_id" not in cols or "aa_seq" not in cols:
            print(f"[SKIP] {f.name}: missing read_id or aa_seq. Columns={list(df.columns)}")
            continue

        df["n"] = df[cols["read_id"]].apply(parse_n)
        d = df[[cols["aa_seq"], "n"]].dropna()
        d = d[d[cols["aa_seq"]].astype(str).str.len() > 0]
        if d.empty:
            print(f"[SKIP] {f.name}: no sequences after filtering")
            continue

        counts = weighted_counts_matrix(d[cols["aa_seq"]], d["n"])

        sample_name = f.name.replace(".filtered.PASS.aa.tsv", "").replace(".aa.tsv", "")
        out_png = out_dir / f"{sample_name}.aa_logo.png"

        plot_logo(counts, out_png, title=f"{sample_name} (AA logo)")
        print(f"[OK] {f.name} -> {out_png}")

if __name__ == "__main__":
    main()
