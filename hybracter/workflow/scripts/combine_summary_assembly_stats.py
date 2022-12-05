#!/usr/bin/env python3

import pandas as pd

# function to read the csv
def read_csv(csv):
    colnames=['Sample', 'total_contigs', 'max_contig_size', 'complete_assembly', 'plasmid_count', 'coverage_biggest_contig'] 
    df = pd.read_csv(csv, delimiter= ',', index_col=False, header=None, names=colnames)
    return df


def summarise_contigs(summary_list, output):
    # read into list       
    summaries = []
    l =summary_list

    for a in l:
        tmp_summary = read_csv(a)
        # remove first row (from the file)
        tmp_summary = tmp_summary.iloc[1: , :]
        summaries.append(tmp_summary)

    # make into combined dataframe
    total_summary_df = pd.concat(summaries,  ignore_index=True)
    total_summary_df.to_csv(output, sep=",", index=False)
        

summarise_contigs(snakemake.input.summaries, snakemake.output.out)




