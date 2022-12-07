#!/usr/bin/env python3

from Bio import SeqIO

# get chromosome if length > 2.5 million

def get_incomplete( input_fasta, incomplete_fasta,   min_chrom_length):
    # read in the fasta
    with open(incomplete_fasta, 'w') as fa:
        for dna_record in SeqIO.parse(input_fasta, "fasta"):
            if len(dna_record.seq) < int(min_chrom_length):
                SeqIO.write(dna_record, fa, 'fasta')

            
get_incomplete(snakemake.input[0], snakemake.output[0],  snakemake.params[0])



        
