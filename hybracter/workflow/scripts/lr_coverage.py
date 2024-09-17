#!/usr/bin/env python3


import pandas as pd
import sys


def get_sum_len(file_path):
    """
    get sum_len from seqkit output
    """

    # Define the column names based on the file format

    # Read the file into a pandas DataFrame
    df = pd.read_csv(file_path, sep="\t")
    sum_len_value = df["sum_len"].values[0]

    return sum_len_value


def get_quick_coverage_estimate_long(
    long_seqkit, min_chrom_length, min_bases, coverage_output_path
):
    """
    takes seqkit input
    """

    long_bases = get_sum_len(long_seqkit)

    # number of bases is over the min_bases
    if long_bases > int(min_bases):
        coverage = round(long_bases / int(min_chrom_length))

        # write to file
        with open(coverage_output_path, "w") as file:
            file.write(str(coverage))
    else: # number of bases is less than min bases
        sys.exit(f"The total number of bases of long-read sequencing is {long_bases}, which is lower than the required minimum value {min_bases} based on the chromosome size multiplied by the --min_depth parameter.\n Please check --min_depth and/or get more long-read sequencing data! Hybracter will stop here.")


get_quick_coverage_estimate_long(
    snakemake.input.long_bases,
    snakemake.params.min_chrom_length,
    snakemake.params.min_bases,
    snakemake.output.lr_coverage,
)
