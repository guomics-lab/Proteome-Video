# preprocess_data.py — Convert downsampled TSV to PyTorch tensors for Neural ODE

import pandas as pd
import torch
import os
import re
import argparse


def preprocess(input_tsv, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(input_tsv, sep='\t', index_col=0)
    print(f"Loaded: {df.shape[0]} proteins x {df.shape[1]} samples")

    col_pattern = re.compile(r'([CAH])(\d+)')
    groups = {'C': [], 'A': [], 'H': []}

    for col in df.columns:
        m = col_pattern.match(col)
        if m:
            groups[m.group(1)].append((col, int(m.group(2))))

    for k in groups:
        groups[k].sort(key=lambda x: x[1])

    # Save time points (from Control columns, same grid for all conditions)
    time_indices = [t for _, t in groups['C']]
    torch.save(torch.tensor(time_indices, dtype=torch.float32),
               os.path.join(output_dir, "Fourier_time_points.pt"))

    for cond, cols in [('C', groups['C']), ('A', groups['A']), ('H', groups['H'])]:
        col_names = [c[0] for c in cols]
        tensor = torch.tensor(df[col_names].values, dtype=torch.float32)
        path = os.path.join(output_dir, f"Fourier_timeseries_{cond}.pt")
        torch.save(tensor, path)

    print("Done.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='../../Step02_Clustering/output/downsampled_median.tsv',
                        help="Input TSV (default: ../../Step02_Clustering/output/downsampled_median.tsv)")
    parser.add_argument('--output', type=str, default='../output', help="Output directory for .pt files")
    args = parser.parse_args()
    preprocess(args.input, args.output)
