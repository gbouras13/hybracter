#!/usr/bin/env python3

import os

import polars as pl


def add_sample_plassembler(input, output, sample):
    if os.path.getsize(input) > 0:
        # read as strings to preserve plassembler's own value formatting exactly
        tmp_summary = pl.read_csv(input, separator="\t", infer_schema_length=0)
        # add the Sample column at the first position (overwriting any existing one)
        tmp_summary = tmp_summary.select(
            pl.lit(sample).alias("Sample"), pl.exclude("Sample")
        )
        tmp_summary.write_csv(output, separator="\t")

    # touch a file if nothing
    else:
        with open(output, "a"):
            os.utime(output, None)


if "snakemake" in globals():  # only runs under Snakemake; lets pytest import this module
    add_sample_plassembler(
        snakemake.input.inp, snakemake.output.out, snakemake.wildcards.sample
    )
