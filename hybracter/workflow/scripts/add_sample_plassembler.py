#!/usr/bin/env python3

import os

import pandas as pd


def add_sample_plassembler(input, output, sample):
    if os.path.getsize(input) > 0:
        tmp_summary = pd.read_csv(input, delimiter="\t", index_col=False, header=0)
        tmp_summary["Sample"] = sample

        # shift column 'Sample' to first position
        first_column = tmp_summary.pop("Sample")

        # insert column using insert(position,column_name,
        # first_column) function
        tmp_summary.insert(0, "Sample", first_column)
        tmp_summary.to_csv(output, sep="\t", index=False)

    # touch a file if nothing
    else:
        with open(output, "a"):
            os.utime(output, None)


add_sample_plassembler(
    snakemake.input.inp, snakemake.output.out, snakemake.wildcards.sample
)
