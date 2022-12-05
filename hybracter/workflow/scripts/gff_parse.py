#!/usr/bin/env python3

import pandas as pd

def extract_is_seq(sample, output):
  # read in cleaned df of ISfinder element
  colnames=['contig', 'evidence', 'type', 'start', 'end', 'score', 'strand', 'n', 'description'] 
  # if file is empty
  try:
    abundance_df = pd.read_csv(sample, delimiter= '\t', index_col=False, header=None, names=colnames)
  except pd.errors.EmptyDataError:
    print('No IS.')

  ### if empty
  empty = False
  if abundance_df.empty:
    empty = True

  # function to get the locus tag (between ID= and the first semi colon)

  def find_between( s, first, last ):
      try:
          start = s.index( first ) + len( first )
          end = s.index( last, start )
          return s[start:end]
      except ValueError:
          return ""

  if empty == False:
    # get is_name and locus_tage
    abundance_df['locus_tag'] = abundance_df['description'].apply(lambda x: find_between(x,"ID=", ";"  ) )
    abundance_df['IS_name'] = abundance_df['description'].apply(lambda x: find_between(x,"ISfinder:", ";locus_tag"  ) )
    # get the product
    # https://stackoverflow.com/questions/37333299/splitting-a-pandas-dataframe-column-by-delimiter
    abundance_df[['description','product']] = abundance_df['description'].str.split('product=',expand=True)
    

    # write to csv
    abundance_df.to_csv(output, sep=",", index=False, header=False)
  else:
    colnames=['contig', 'evidence', 'type', 'start', 'end', 'score', 'strand', 'n', 'description', 'locus_tag', 'IS_name', 'product'] 
    abundance_df = pd.DataFrame(columns=colnames)
    abundance_df.to_csv(output, sep=",", index=False, header=False)
 
extract_is_seq(snakemake.input[0], snakemake.output[0])




