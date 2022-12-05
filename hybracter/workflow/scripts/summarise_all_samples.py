#!/usr/bin/env python3

import pandas as pd
import numpy as np

def read_gff_csv(gff_csv):
    colnames=['sample', 'contig',  'start', 'end',  'strand', 'locus_tag', 'IS_name', 'product', 'representative_fasta'] 
    gff = pd.read_csv(gff_csv, delimiter= ',', index_col=False, header=None, names=colnames)
    return gff


def collate_gffs( is_finder_df_list, cluster_df, count_out_per_sample, counts_matrix):

    gffs = []
    l =is_finder_df_list

    for a in l:
        tmp_gff = read_gff_csv(a)
        # remove first row (from the file)
        tmp_gff = tmp_gff.iloc[1: , :]
        gffs.append(tmp_gff)

    # make into combined dataframe
    total_gff = pd.concat(gffs,  ignore_index=True)

    # drop representative_fasta 
    total_gff = total_gff.drop(columns=['representative_fasta'])

    # read in cluster
    colnames_mmseqs=['representative_fasta', 'locus_tag']
    mmseqs = pd.read_csv(cluster_df, delimiter= '\t', index_col=False, header=None, names=colnames_mmseqs) 

    merged = pd.merge(total_gff, mmseqs)

    merged['count'] = merged.groupby(['sample','representative_fasta'])['representative_fasta'].transform('count')

    # final summaries

    count_df = merged[["sample", "representative_fasta","IS_name", "product", "count"]].drop_duplicates()

    count_df.to_csv(count_out_per_sample, sep=",", index=False)

    count_df_p = merged[["sample", "representative_fasta", "count"]].drop_duplicates()
    count_df_p = count_df_p.pivot(index="sample", columns="representative_fasta", values="count").fillna(0)
    count_df_p['Total Insertion Sequences'] = np.sum(count_df_p,axis=1)
    count_df_p.to_csv(counts_matrix, sep=",")






collate_gffs(snakemake.input.finals,snakemake.input.cluster, snakemake.output[0],snakemake.output[1])




