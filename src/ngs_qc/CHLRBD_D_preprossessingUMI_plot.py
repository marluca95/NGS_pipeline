import matplotlib.pyplot as plt
import pandas as pd
import logomaker
import argparse

def plot_aa_sequence_logo(df: pd.DataFrame,
                          seq_col: str,
                          color_scheme: str = "chemistry",
                          figsize: tuple = (12, 4),
                          save_path: str = None):
    """
    Generate an amino-acid sequence logo plot from sequences in a specified DataFrame column,
    colored by biophysical properties.

    Parameters:
    - df: pandas DataFrame containing your sequences.
    - seq_col: name of the column with fixed-length amino acid sequences.
    - color_scheme: one of logomaker's built-in schemes, e.g. 
        'chemistry', 'hydrophobicity', 'charge', 'clustal', etc.
      You can also pass a dict mapping each amino acid to a color.
    - figsize: tuple specifying the size of the figure.
    - save_path: optional filepath to save the plot (e.g., 'aa_logo.png').
    """
    seqs = df[seq_col].dropna().astype(str).tolist()
    lens = {len(s) for s in seqs}
    if len(lens) != 1:
        raise ValueError("All sequences must be the same length to plot a logo.")
    L = lens.pop()

    # Build count and frequency DataFrame
    aas = [
        'A','C','D','E','F','G','H','I','K','L',
        'M','N','P','Q','R','S','T','V','W','Y'
    ]
    counts = []
    for pos in range(L):
        col = {aa: 0 for aa in aas}
        for seq in seqs:
            aa = seq[pos]
            if aa in col:
                col[aa] += 1
        counts.append(col)
    count_df = pd.DataFrame(counts)
    freq_df  = count_df.div(count_df.sum(axis=1), axis=0)

    # Create the logo, with coloring by your chosen scheme
    plt.figure(figsize=figsize)
    logo = logomaker.Logo(
        freq_df,
        color_scheme=color_scheme,   # e.g. "chemistry", "hydrophobicity", "charge", etc.
        font_name='ArialRoundedMTBold'  # optional: nicer font for amino acids
    )
    plt.xlabel("Position")
    plt.ylabel("Frequency")
    plt.title("Amino Acid Sequence Logo\ncolored by {}".format(color_scheme))
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
    plt.show()


# Parse the YAML file path as an argument
parser = argparse.ArgumentParser(description="Process and analyze FASTQ files.")
parser.add_argument('--samplename', type=str, required=True, help="Name of the sample")
parser.add_argument('--csv', type=str, required=True, help="Path to the csv")
parser.add_argument('--outpath', type=str, required=True, help="Path to the output")

args = parser.parse_args()

df = pd.read_csv(args.csv)

plot_aa_sequence_logo(df, 'cons_aa', color_scheme = "chemistry", figsize = (14,4),  save_path= f"./{args.outpath}/{args.samplename}_logoplot_consAA.png" )

