#!/usr/bin/env python3

import glob
import os

import polars as pl


def make_final_summary(hybracter_summary, complete_summary_dir, incomplete_summary_dir):
    """
    reads all individual hybracter summaries and combines them all
    """

    # Use glob to find files with the .score extension in the directory
    complete_file_list = glob.glob(os.path.join(complete_summary_dir, "*summary.tsv"))
    incomplete_file_list = glob.glob(
        os.path.join(incomplete_summary_dir, "*summary.tsv")
    )

    # Read each per-sample summary as strings so the values (ints, floats, and the
    # "Unknown" Number_circular_plasmids for incomplete assemblies) are carried into
    # the combined summary exactly as written - this is the FINAL_OUTPUT summary.
    summary_dfs = [
        pl.read_csv(file_path, separator="\t", infer_schema_length=0)
        for file_path in complete_file_list + incomplete_file_list
    ]

    # Concatenate all DataFrames into a single DataFrame
    combined_df = pl.concat(summary_dfs, how="vertical_relaxed")

    combined_df.write_csv(hybracter_summary, separator="\t")


if "snakemake" in globals():  # only runs under Snakemake; lets pytest import this module
    make_final_summary(
        snakemake.output.hybracter_summary,
        snakemake.params.complete_summaries_dir,
        snakemake.params.incomplete_summaries_dir,
    )
