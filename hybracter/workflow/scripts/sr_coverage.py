#!/usr/bin/env python3


import pandas as pd


def get_sum_len(file_path):
    """
    get sum_len from seqkit output
    """

    # Define the column names based on the file format

    # Read the file into a pandas DataFrame
    df = pd.read_csv(file_path, sep="\t")
    sum_len_value = df["sum_len"].values[0]

    return sum_len_value


def get_quick_coverage_estimate(
    r1_seqkit, r2_seqkit, min_chrom_length, coverage_output_path
):
    """
    takes seqkit input
    """

    r1_bases = get_sum_len(r1_seqkit)
    r2_bases = get_sum_len(r2_seqkit)

    coverage = round((r1_bases + r2_bases) / int(min_chrom_length))

    # write to file
    with open(coverage_output_path, "w") as file:
        file.write(str(coverage))


get_quick_coverage_estimate(
    snakemake.input.r1,
    snakemake.input.r2,
    snakemake.params.chromlen,
    snakemake.output.sr_coverage,
)
