#!/usr/bin/env python3

import glob
import os

import pandas as pd


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
                tmp_summary = pd.read_csv(a, delimiter="\t", index_col=False, header=0)
                summaries.append(tmp_summary)
                # appends to the number of plasmids
                samples_with_plasmids += 1

    # save as tsv
    if samples_with_plasmids > 1:
        # make into combined dataframe
        total_summary_df = pd.concat(summaries, ignore_index=True)
        total_summary_df.to_csv(output, sep="\t", index=False)
    elif samples_with_plasmids == 1:  # just print the first
        summaries[0].to_csv(output, sep="\t", index=False)
    else:  # touch the output
        with open(output, "a"):
            os.utime(output, None)


combine_sample_plassembler(snakemake.params.summary_dir, snakemake.output.out)
