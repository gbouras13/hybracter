#!/usr/bin/env python3

from Bio import SeqIO

# get chromosome if length > 2.5 million


def get_completeness(input_fasta, completeness_check, min_chrom_length):
    # flag for if the assembly is complete
    comp = False
    for dna_record in SeqIO.parse(input_fasta, "fasta"):
        if len(dna_record.seq) > int(min_chrom_length):
            comp = True

    with open(completeness_check, "w") as f:
        if comp == True:
            f.write("C")
        else:
            f.write("I")


get_completeness(
    snakemake.input.fasta,
    snakemake.output.completeness_check,
    snakemake.params.min_chrom_length,
)
