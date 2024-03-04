#!/usr/bin/env python3

import pandas as pd
from Bio import SeqIO
import os

"""
get chromosome if length > min chrom length 
add ignore_list if the contig is not marked as circular by Flye - to make sure Dnaapler doesn't rotate it
"""

# touches an empty file
def touch_file(path):
    with open(path, "a"):
        os.utime(path, None)



def get_chromosome_plasmids(
    input_fasta, chromosome_fasta, ignore_list, min_chrom_length, info_file_path, polypolish_flag
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
            # this should never every happen as there won't be dupes in the Flye info shee but just in case
            elif len(circ_values) > 1:
                circ_value = circ_values[0]
            # will be 0 (so under this case) if:
                # chromosome didn't circularise
                # plasmids - not part of the chromosome 
                # so 'N' to make sure it doesn't get extracted
            else:
                circ_value = "N"
            # will actually extract multiple chromosomes if above the desired chrom length and circularised
            if len(dna_record.seq) > int(min_chrom_length):
                SeqIO.write(dna_record, fa, "fasta")
                if circ_value == "N":
                    with open(ignore_list, 'a') as file:
                        file.write(f'{seq_id}\n')

                    
            

get_chromosome_plasmids(
    snakemake.input.fasta,
    snakemake.output.fasta,
    snakemake.output.ignore_list,
    snakemake.params.min_chrom_length,
    snakemake.input.info,
    snakemake.params.polypolish_flag,
)
