#!/usr/bin/env python3

from Bio import SeqIO
import pandas as pd

# get chromosome if length > min chrom length and circ by Flye
# need to trim '_polypolish' on the end of the seq_name for polypolish and pypolca extractions as polypolish introduces this

def get_chromosome_plasmids(input_fasta, chromosome_fasta, min_chrom_length, info_file_path, polypolish_flag):
    # read in the fasta

    # Read the flye info file
    info_df = pd.read_csv(info_file_path, sep='\t', comment='#', header=None, names=['seq_name', 'length', 'cov', 'circ', 'repeat', 'mult', 'alt_group', 'graph_path'])    

    with open(chromosome_fasta, "w") as fa:
        for dna_record in SeqIO.parse(input_fasta, "fasta"):
            # get the circ value
            # if 
            if polypolish_flag is False:
                circ_value = info_df[info_df['seq_name'] == dna_record.id]['circ'].values[0]
            else:
                # trim the last 11 characters _polypolish
                trim_seq_id = dna_record.id[:-11]
                circ_value = info_df[info_df['seq_name'] == trim_seq_id]['circ'].values[0]
            # will actually extract multiple chromosomes if above the desired chrom length
            if len(dna_record.seq) > int(min_chrom_length) and circ_value == 'Y':
                SeqIO.write(dna_record, fa, "fasta")


get_chromosome_plasmids(
    snakemake.input.fasta, snakemake.output.fasta, snakemake.params.min_chrom_length, snakemake.input.info,
    snakemake.params.polypolish_flag
)
