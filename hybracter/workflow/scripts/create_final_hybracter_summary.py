#!/usr/bin/env python3

import glob
import os

import pandas as pd


def make_final_summary(hybracter_summary, complete_summary_dir, incomplete_summary_dir):
    """
    reads all individual hybracter summaries and combines them all
    """

    # Use glob to find files with the .score extension in the directory
    complete_file_list = glob.glob(os.path.join(complete_summary_dir, "*summary.tsv"))
    incomplete_file_list = glob.glob(
        os.path.join(incomplete_summary_dir, "*summary.tsv")
    )

    # list
    summary_dfs = []

    # complete
    for file_path in complete_file_list:
        # Read the DataFrame from the file and append it to the list
        df = pd.read_csv(file_path, delimiter="\t")  # Assuming tab-separated files
        summary_dfs.append(df)

    # incomplete
    for file_path in incomplete_file_list:
        # Read the DataFrame from the file and append it to the list
        df = pd.read_csv(file_path, delimiter="\t")  # Assuming tab-separated files
        summary_dfs.append(df)

    # Concatenate all DataFrames into a single DataFrame
    combined_df = pd.concat(summary_dfs, ignore_index=True)

    combined_df.to_csv(hybracter_summary, index=False, sep="\t")


make_final_summary(
    snakemake.output.hybracter_summary,
    snakemake.params.complete_summaries_dir,
    snakemake.params.incomplete_summaries_dir,
)
