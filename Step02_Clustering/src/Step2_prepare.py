# Step2_prepare.py — Rolling-mean smoothing and median downsampling
import pandas as pd
import numpy as np
import re, os, sys
from collections import defaultdict


def load_and_parse_data(filename):
    df = pd.read_csv(filename, sep='\t', index_col=0)
    treatments = defaultdict(list)
    for col in df.columns:
        match = re.match(r'([A-Za-z]+)(\d+)', col)
        if match:
            t, tp = match.group(1), int(match.group(2))
            treatments[t].append((tp, col))
    for t in treatments:
        treatments[t].sort(key=lambda x: x[0])
    return df, treatments


def apply_moving_average_to_df(df, treatments, window_size=6):
    df_smoothed = df.copy()
    for treatment, time_cols in treatments.items():
        columns = [c[1] for c in time_cols]
        for gene in df.index:
            df_smoothed.loc[gene, columns] = df.loc[gene, columns].rolling(
                window=window_size, center=True, min_periods=1).mean().values
    return df_smoothed


def downsample_data_median(data, bin_size=6):
    if len(data) < bin_size:
        return data.copy()
    num_bins = len(data) // bin_size
    truncated = data[:num_bins * bin_size]
    reshaped = truncated.values.reshape(num_bins, bin_size)
    downsampled = np.median(reshaped, axis=1)
    new_indices = [np.median(np.arange(len(data))[i*bin_size:(i+1)*bin_size])
                   for i in range(num_bins)]
    new_index = [data.index[int(idx)] for idx in new_indices]
    return pd.Series(downsampled, index=new_index)


def apply_downsampling_to_df(df, treatments, bin_size=6):
    downsampled_data = {}
    for treatment, time_cols in treatments.items():
        columns = [c[1] for c in time_cols]
        time_points = [c[0] for c in time_cols]
        num_bins = len(time_points) // bin_size
        new_cols = []
        for i in range(num_bins):
            avg_time = np.median(time_points[i*bin_size:(i+1)*bin_size])
            new_cols.append(f"{treatment}{int(avg_time)}")
        for gene in df.index:
            gene_data = df.loc[gene, columns].dropna()
            if len(gene_data) >= bin_size:
                ds = downsample_data_median(gene_data, bin_size)
                for j, nc in enumerate(new_cols):
                    if j < len(ds):
                        downsampled_data.setdefault(nc, {})[gene] = ds.iloc[j]
    df_out = pd.DataFrame(downsampled_data).reindex(df.index)
    return df_out


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    input_file = os.path.join(script_dir, "..", "..", "Step00_Preprocess",
                              "output", "imputed_temporal_matrix.tsv")
    if not os.path.exists(input_file):
        print(f"Error: {input_file} not found"); sys.exit(1)
    out_dir = os.path.join(script_dir, "..", "output")
    os.makedirs(out_dir, exist_ok=True)
    df, treatments = load_and_parse_data(input_file)
    print(f"Loaded: {df.shape[0]} proteins x {df.shape[1]} samples")

    df_sm = apply_moving_average_to_df(df, treatments, window_size=6)
    df_sm.to_csv(os.path.join(out_dir, "smoothed_average.tsv"), sep='\t')
    print(f"  -> smoothed_average.tsv")

    ds = apply_downsampling_to_df(df, treatments, bin_size=6)
    ds.to_csv(os.path.join(out_dir, "downsampled_median.tsv"), sep='\t')
    print(f"  -> downsampled_median.tsv")


if __name__ == "__main__":
    main()
