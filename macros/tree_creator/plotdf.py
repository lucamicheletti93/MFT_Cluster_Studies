import argparse
import polars as pl
import os
import matplotlib.pyplot as plt
from hipe4ml import plot_utils


def plot_keys_from_parquet(parquet_file, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Read the .gzip parquet file into a Polars DataFrame
    df = pl.read_parquet(parquet_file)
    pd_df = df.to_pandas()
    columns = list(pd_df.columns.values)

    plot_utils.plot_distr(
        pd_df,
        columns,
        100,
        "mip",
        figsize=(12, 7),
        alpha=0.3,
        log=True,
        grid=False,
        density=True,
    )
    plt.subplots_adjust(
        left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55
    )
    plt.savefig(f"{output_dir}/DistributionsAll.png")
    plt.close("all")


    corrmatrix_fig = plot_utils.plot_corr(
        [pd_df],
        columns,
        ["mip"]
    )
    plt.subplots_adjust(left=0.2, bottom=0.25, right=0.95, top=0.9)
    corrmatrix_fig.savefig(f"{output_dir}/CorrAll.png")
    plt.close("all")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Upload gzip.parquet file and plot every key in a single canvas"
    )
    parser.add_argument("parquet_file", type=str, help="Path to the .gzip.parquet file")
    parser.add_argument(
        "output_dir", type=str, help="Directory to save the output plots"
    )

    args = parser.parse_args()
    plot_keys_from_parquet(args.parquet_file, args.output_dir)
