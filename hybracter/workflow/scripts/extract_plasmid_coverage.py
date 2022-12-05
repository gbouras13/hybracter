#!/usr/bin/env python3

import pandas as pd

def calc_coverage(assembly_info, sample,  summary_out ):
    # read gff        

    colnames=['seq_name', 'length', 'cov', 'circ', 'repeat', 'mult', 'alt_group', 'graph_path'] 
    assembly_df = pd.read_csv(assembly_info, delimiter= '\t', index_col=False, header=None, names=colnames)
    #assembly_df = pd.read_csv('assembly_info.txt', delimiter= '\t', index_col=False, header=None, names=colnames)

    # remove first row (from the file)
    assembly_df = assembly_df.iloc[1: , :]

    # add sample
    assembly_df['sample'] = sample

    # Convert copv and length to int
    assembly_df['length'] = assembly_df['length'].astype('int')
    assembly_df['cov'] = assembly_df['cov'].astype('int')

    # count contigs
    total_contigs = len(assembly_df.index)

    # get max contig size
    # need to convert to integer first
    max_contig = assembly_df["length"].max()
    #print(max_contig)

    # covnert to int
    max_contig = int(max_contig)

    # get coverage of largest assembly 
    max_contig_cov = assembly_df[assembly_df["length"] == max_contig].iloc[0]['cov']

    assembly_df["copy_number"] = (assembly_df["cov"] / max_contig_cov).round(2)

    assembly_df = assembly_df.drop(assembly_df[assembly_df["copy_number"] == 1].index)

    assembly_df.to_csv(summary_out, sep=",", index=False)


calc_coverage(snakemake.input[0],  snakemake.wildcards.sample,  snakemake.output[0])




