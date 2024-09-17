#!/usr/bin/env python3

import re


def write_chrom_size(
    kmc_log_path, output_file
):

    with open(kmc_log_path, "r") as file:
        for line in file:
            if "No. of unique counted k-mers" in line:
                # keep 80% of kmers as lower bound for chromosome
                chrom_size = int(float(re.search(r"No. of unique counted k-mers\s*:\s*([\d\.eE+-]+)", line).group(1)) * 0.8)

    # write to file
    with open(output_file, "w") as file:
        file.write(str(chrom_size))


write_chrom_size(
    snakemake.input.kmcLOG,
    snakemake.output.chrom_size,
)
