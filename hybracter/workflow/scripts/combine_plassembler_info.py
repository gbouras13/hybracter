#!/usr/bin/env python3

import pandas as pd
import os


def combine_sample_plassembler(summary_list, output):
    # read into list
    summaries = []
    l = summary_list

    for a in l:
        # only if > 0
        if os.path.getsize(a) > 0:
            tmp_summary = pd.read_csv(a, delimiter="\t", index_col=False, header=0)
            summaries.append(tmp_summary)

    # make into combined dataframe
    total_summary_df = pd.concat(summaries, ignore_index=True)
    total_summary_df.to_csv(output, sep=",", index=False)


combine_sample_plassembler(snakemake.input.summaries, snakemake.output.out)
