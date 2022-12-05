#!/usr/bin/env python3

import pandas as pd

def collate_sample( is_finder_df, cluster_df, sample, merged_out, count_out):
    # read gff        

    colnames=['contig', 'evidence', 'type', 'start', 'end', 'score', 'strand', 'n', 'description', 'locus_tag', 'IS_name', 'product'] 
    gff = pd.read_csv(is_finder_df, delimiter= ',', index_col=False, header=None, names=colnames)
    # read mmseqs output

    colnames_mmseqs=['representative_fasta', 'locus_tag']
    mmseqs = pd.read_csv(cluster_df, delimiter= '\t', index_col=False, header=None, names=colnames_mmseqs) 

    merged = pd.merge(gff, mmseqs)
    merged = merged.drop(columns=['evidence', 'score', 'type', 'n', 'description'])
    merged['sample'] = sample

    # move sample to front 

    cols = merged.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    merged = merged[cols]

    merged.to_csv(merged_out, sep=",", index=False)

    count_df = merged['representative_fasta'].value_counts().rename_axis('representative_fasta').reset_index(name='count')

    merged_subset = merged[["representative_fasta", "IS_name", "product"]].drop_duplicates()

    count_df = pd.merge(count_df, merged_subset)

    count_df.merge(merged_subset, how='left', on='representative_fasta')
    # add total
    total_is = count_df['count'].sum()
    total_df = pd.DataFrame([['Total Insertion Sequences', total_is, '', '']], columns=["representative_fasta", "count", "IS_name", "product"], index=['total'])
    count_df = count_df.append(total_df)

    count_df.to_csv(count_out, sep=",", index=False)

           
 
collate_sample(snakemake.input[0],snakemake.input[1], snakemake.wildcards.sample,  snakemake.output[0],snakemake.output[1])




