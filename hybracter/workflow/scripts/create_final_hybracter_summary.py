#!/usr/bin/env python3

import pandas as pd
import glob
import os 
import sys
from Bio import SeqIO



def make_final_summary(hybracter_summary, complete_summary_dir, incomplete_summary_dir ):
    """
    reads all individual hybracter summaries and combines them all
    """

    # Use glob to find files with the .score extension in the directory
    complete_file_list = glob.glob(os.path.join(complete_summary_dir, '*.tsv'))
    incomplete_file_list = glob.glob(os.path.join(incomplete_summary_dir, '*.tsv'))

    # list
    summary_dfs = []

    # complete
    for file_path in complete_file_list:
        # Read the DataFrame from the file and append it to the list
        df = pd.read_csv(file_path, delimiter='\t')  # Assuming tab-separated files
        summary_dfs.append(df)
    
    # incomplete
    for file_path in incomplete_file_list:
        # Read the DataFrame from the file and append it to the list
        df = pd.read_csv(file_path, delimiter='\t')  # Assuming tab-separated files
        summary_dfs.append(df)

    # Concatenate all DataFrames into a single DataFrame
    combined_df = pd.concat(summary_dfs, ignore_index=True)

    combined_df.to_csv(hybracter_summary, index=False, sep = "\t")

select_best_chromosome_assembly_complete(snakemake.output.hybracter_summary, snakemake.params.complete_summaries_dir, snakemake.input.incomplete_summaries_dir  )



