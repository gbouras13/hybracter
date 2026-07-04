#!/usr/bin/env python3

import glob
import os

import polars as pl


def combine_sample_plassembler(summary_dir, output):
    # read into list

    # Specify a pattern to match files (e.g., all .txt files)
    pattern = "*.tsv"

    # Get a list of files that match the pattern in the directory
    summary_list = glob.glob(os.path.join(summary_dir, pattern))

    # write all the summary dfs to a list
    summaries = []

    # counts the number of summaries
    complete_sample_no = len(summary_list)

    # counts the number of samples with plasmids
    samples_with_plasmids = 0

    if complete_sample_no > 0:
        for a in summary_list:
            # only if > 0
            if os.path.getsize(a) > 0:
                # read as strings to preserve values exactly across the concat
                summaries.append(pl.read_csv(a, separator="\t", infer_schema_length=0))
                # appends to the number of plasmids
                samples_with_plasmids += 1

    # save as tsv
    if samples_with_plasmids > 1:
        # make into combined dataframe
        total_summary_df = pl.concat(summaries, how="vertical_relaxed")
        total_summary_df.write_csv(output, separator="\t")
    elif samples_with_plasmids == 1:  # just print the first
        summaries[0].write_csv(output, separator="\t")
    else:  # touch the output
        with open(output, "a"):
            os.utime(output, None)


if "snakemake" in globals():  # only runs under Snakemake; lets pytest import this module
    combine_sample_plassembler(snakemake.params.summary_dir, snakemake.output.out)
