#!/usr/bin/env python3


def get_sum_len(file_path):
    """
    get sum_len from seqkit stats -T output
    """
    with open(file_path) as fh:
        cols = fh.readline().strip().split("\t")
        vals = fh.readline().strip().split("\t")
    return int(vals[cols.index("sum_len")])


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


if "snakemake" in globals():  # only runs under Snakemake; lets pytest import this module
    get_quick_coverage_estimate(
        snakemake.input.r1,
        snakemake.input.r2,
        snakemake.params.chromlen,
        snakemake.output.sr_coverage,
    )
