#!/usr/bin/env python3

import pandas as pd

def get_is_rows(sample, output):
  # read in cleaned df of ISfinder element
  colnames=['contig', 'evidence', 'type', 'start', 'end', 'score', 'strand', 'n', 'description'] 
  abundance_df = pd.read_csv(sample, delimiter= '\t', index_col=False, header=None, names=colnames)
  abundance_df = abundance_df[abundance_df['description'].str.contains('ISfinder')]
  abundance_df.to_csv(output, sep='\t', index=False, header=False)
 
get_is_rows(snakemake.input[0], snakemake.output[0])




