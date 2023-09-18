#!/usr/bin/env python3

from Bio import SeqIO

# get chromosome if length > min chrom length


def get_chromosome_plasmids(input_fasta, chromosome_fasta, min_chrom_length):
    # read in the fasta
    with open(chromosome_fasta, "w") as fa:
        for dna_record in SeqIO.parse(input_fasta, "fasta"):
            # will actually extract multiple chromosomes if above the desired chrom length
            if len(dna_record.seq) > int(min_chrom_length):
                SeqIO.write(dna_record, fa, "fasta")


get_chromosome_plasmids(
    snakemake.input.fasta, snakemake.output.fasta, snakemake.params.min_chrom_length
)
