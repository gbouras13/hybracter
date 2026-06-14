#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO

# get chromosome if length > 2.5 million


def get_completeness(input_fasta, completeness_check, min_chrom_length, info_file_path, circular_chromosome):
    # Read the flye info file
    info_df = pd.read_csv(
        info_file_path,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "seq_name",
            "length",
            "cov",
            "circ",
            "repeat",
            "mult",
            "alt_group",
            "graph_path",
        ],
    )

    # flag for if the assembly is complete
    comp = False
    for dna_record in SeqIO.parse(input_fasta, "fasta"):
        circ_value = info_df[info_df["seq_name"] == dna_record.id]["circ"].values[0]
        if circular_chromosome:
            # --circular_chromosome: complete requires length AND Flye circularity
            if len(dna_record.seq) > int(min_chrom_length) and circ_value == "Y":
                comp = True
        else:
            # default: complete means length above threshold only
            if len(dna_record.seq) > int(min_chrom_length):
                comp = True

    with open(completeness_check, "w") as f:
        if comp == True:
            f.write("C")
        else:
            f.write("I")


if "snakemake" in globals():  # only runs under Snakemake; lets pytest import this module
    get_completeness(
        snakemake.input.fasta,
        snakemake.output.completeness_check,
        snakemake.params.min_chrom_length,
        snakemake.input.info,
        snakemake.params.circular_chromosome,
    )
