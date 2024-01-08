#!/usr/bin/env python3

from Bio import SeqIO

# gets all contigs for incomplete assemblies


def get_incomplete(input_fasta, incomplete_fasta, min_chrom_length):
    # read in the fasta and save all contigs
    with open(incomplete_fasta, "w") as fa:
        for dna_record in SeqIO.parse(input_fasta, "fasta"):
            SeqIO.write(dna_record, fa, "fasta")


get_incomplete(
    snakemake.input.fasta, snakemake.output.fasta, snakemake.params.min_chrom_length
)
