#!/usr/bin/env python3

import pandas as pd
import re
from Bio import SeqIO

def parse_fastas( is_finder_df,fasta, output):
  
    # read in gff

    colnames=['contig', 'evidence', 'type', 'start', 'end', 'score', 'strand', 'n', 'description', 'locus_tag', 'IS_name', 'product'] 

    gff = pd.read_csv(is_finder_df, delimiter= ',', index_col=False, header=None, names=colnames)

    # read in the fasta
    with open(output, 'w') as aa_fa:
        for dna_record in SeqIO.parse(fasta, "fasta"):
            for locus in gff['locus_tag']:
                if re.search(locus, dna_record.id):
                    SeqIO.write(dna_record, aa_fa, 'fasta')
                
 
parse_fastas(snakemake.input[0],snakemake.input[1], snakemake.output[0])



        
