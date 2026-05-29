#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import os

"""
Extract chromosomes for the complete assembly path.
Default: extract all contigs above min_chrom_length regardless of circularity.
Non-circular contigs are added to ignore_list so dnaapler skips reorientation for them.
With --circular_chromosome: only extract contigs that are BOTH above min_chrom_length AND
marked circular by Flye.
"""

# touches an empty file
def touch_file(path):
    with open(path, "a"):
        os.utime(path, None)



def get_chromosome_plasmids(
    input_fasta, chromosome_fasta, ignore_list, min_chrom_length, info_file_path, polypolish_flag, circular_chromosome
):
    # read in the fasta

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

    # to make sure it gets made regardless for snakemake
    touch_file(ignore_list)

    with open(chromosome_fasta, "w") as fa:
        for dna_record in SeqIO.parse(input_fasta, "fasta"):

            seq_id = dna_record.id

            # get all the values
            circ_values = info_df.loc[info_df["seq_name"] == seq_id, "circ"].values

            if len(circ_values) == 1:
                circ_value = circ_values[0]
            # this should never every happen as there won't be dupes in the Flye info sheet but just in case
            elif len(circ_values) > 1:
                circ_value = circ_values[0]
            # will be 0 (so under this case) if:
                # chromosome didn't circularise
                # plasmids - not part of the chromosome
                # so 'N' to make sure it doesn't get extracted with --circular_chromosome
            else:
                circ_value = "N"

            if circular_chromosome:
                # --circular_chromosome: only extract contigs that are long AND circular
                if len(dna_record.seq) > int(min_chrom_length) and circ_value == "Y":
                    SeqIO.write(dna_record, fa, "fasta")
            else:
                # default: extract all long contigs; add non-circular to ignore_list
                # so dnaapler skips reorientation for them
                if len(dna_record.seq) > int(min_chrom_length):
                    SeqIO.write(dna_record, fa, "fasta")
                    if circ_value == "N":
                        with open(ignore_list, "a") as file:
                            file.write(f"{seq_id}\n")


get_chromosome_plasmids(
    snakemake.input.fasta,
    snakemake.output.fasta,
    snakemake.output.ignore_list,
    snakemake.params.min_chrom_length,
    snakemake.input.info,
    snakemake.params.polypolish_flag,
    snakemake.params.circular_chromosome,
)
