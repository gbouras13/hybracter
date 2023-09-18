#!/usr/bin/env python3

import sys
import pandas as pd
from util import calculate_mean_CDS_length
from Bio import SeqIO


def select_best_chromosome_assembly_long_incomplete(
    hybracter_summary,
    pyrodigal_summary,
    output_fasta,
    pre_polish_fasta,
    medaka_fasta,
    sample,
):
    """
    get prodigal mean length for each assembly
    Then creates summary tsv
    statistics similar to unicycler
    instead of 1,2,3 etc we will use 'contig00001', 'contig00002' etc as more parsable
    Then creates hybracter_summary tsv
    """

    # get mean CDS lengths
    pre_polish_mean_cds = calculate_mean_CDS_length(pre_polish_fasta)
    medaka_mean_cds = calculate_mean_CDS_length(medaka_fasta)

    # create dict with mean CDS lengths
    dict = {
        "Sample": sample,
        "pre_polish_mean_cds_length": pre_polish_mean_cds,
        "medaka_mean_cds_length": medaka_mean_cds,
    }

    # Convert the dictionary to a DataFrame
    summary_df = pd.DataFrame([dict])

    # Write to a CSV file
    summary_df.to_csv(pyrodigal_summary, index=False)

    # determine the best assembly
    best_assembly = medaka_fasta
    best_round = "medaka"

    if pre_polish_fasta > medaka_fasta:
        best_assembly = pre_polish_fasta
        best_round = "pre_polish"

    # count contigs
    number_of_contigs = 1

    # instantiate longest contig length
    longest_contig_length = 0

    total_assembly_length = 0

    # Open the output file in write mode
    with open(output_fasta, "w") as output_handle:
        # Iterate through the records in the best assembly FASTA file and write them to the output file
        for record in SeqIO.parse(best_assembly, "fasta"):
            # to match the 00001 output favoured generally for parsing
            # usually there will be 1 chromosome of course!
            record.id = f"contig{number_of_contigs:05}"

            # Calculate the length of the sequence
            sequence_length = len(record.seq)

            # total assembly length
            total_assembly_length += sequence_length

            # to get longest contig
            if number_of_contigs == 1:
                longest_contig_length = sequence_length
            else:
                if sequence_length > longest_contig_length:
                    longest_contig_length = sequence_length

            # Update the description (header) with the length information
            record.description = f"len={sequence_length}"

            # Write the modified record to the output file
            SeqIO.write(record, output_handle, "fasta")

            number_of_contigs += 1

    # to get the summary df
    summary_dict = {
        "Sample": sample,
        "Complete": "False",
        "Total_assembly_length": total_assembly_length,
        "Number_of_contigs": number_of_contigs,
        "Most_accurate_polishing_round": best_round,
        "Longest_contig_length": longest_contig_length,
        "Number_circular_plasmids": "Unknown",
    }

    # Create a DataFrame from the dictionary
    summary_df = pd.DataFrame([summary_dict])
    summary_df.to_csv(hybracter_summary, index=False, sep="\t")


select_best_chromosome_assembly_long_incomplete(
    snakemake.output.hybracter_summary,
    snakemake.output.pyrodigal_summary,
    snakemake.output.total_fasta,
    snakemake.params.pre_polish_fasta,
    snakemake.params.medaka_fasta,
    snakemake.wildcards.sample,
)
