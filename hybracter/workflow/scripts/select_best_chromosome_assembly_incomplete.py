#!/usr/bin/env python3

import glob
import os

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def select_best_chromosome_assembly_incomplete(
    hybracter_summary,
    per_conting_summary,
    ale_dir,
    output_fasta,
    ale_summary,
    pre_polish_fasta,
    medaka_fasta,
    polypolish_fasta,
    polca_fasta,
    sample,
    flye_info,
    logic,
    no_pypolca,
):
    """
    reads all the .score files in teh ale directory, picks the best one (closest to zero) and then takes that chromosome fasta and writes it to file with length
    statistics similar to unicycler
    instead of 1,2,3 etc we will use 'contig00001', 'contig00002' etc as more parsable (lessons from dragonflye)
    Then creates hybracter_summary tsv
    """

    # Use glob to find files with the .score extension in the directory
    file_list = glob.glob(os.path.join(ale_dir, "*.score"))

    # Create an empty dictionary to store the results
    score_dict = {}

    for file_path in file_list:
        # Strip the ".score" extension and use it as the dictionary key
        file_name = os.path.splitext(os.path.basename(file_path))[0]

        # Initialize the score as None (in case the file doesn't contain a valid score)
        score = None

        # Read the first line of the file
        with open(file_path, "r") as file:
            first_line = file.readline().strip()

        # Check if the first line is a valid float
        try:
            score = float(first_line)
        except ValueError:
            pass  # If it's not a valid float, score remains None

        # Store the score (None if it wasn't a valid float)
        score_dict[file_name] = score

    # Filter out None values from the score_dict
    filtered_score_dict = {k: v for k, v in score_dict.items() if v is not None}

    # Find the earliest polishing round associated with the best score
    best_round = None

    # Check if there are any valid scores left
    if filtered_score_dict:
        # Find the key with the score closest to 0
        # this will be the best score from ALE
        closest_to_zero_key = min(
            filtered_score_dict, key=lambda k: abs(filtered_score_dict[k] - 0)
        )
        best_score = filtered_score_dict[closest_to_zero_key]

    # output df with scores and round
    scores_df = pd.DataFrame(list(score_dict.items()), columns=["Key", "Score"])
    # sorts ascending - worst top, best bottom
    scores_df.sort_values(by="Score", ascending=True, inplace=True)
    scores_df.to_csv(ale_summary, index=False, sep="\t")

    # by default the best assembly is the polca or polypolish
    # check that the best assembly wasn't something else
    if no_pypolca is False:
        best_assembly = polca_fasta
        best_round = "pypolca"
    else:
        best_assembly = polypolish_fasta
        best_round = "polypolish"

    if logic == "best":
        if "incomp_pre_polish" in closest_to_zero_key:
            best_assembly = pre_polish_fasta
            best_round = "pre_polish"
        elif "medaka" in closest_to_zero_key:
            best_assembly = medaka_fasta
            best_round = "medaka"
        elif "polypolish" in closest_to_zero_key:
            best_assembly = polypolish_fasta
            best_round = "polypolish"
        else:  # polca
            best_assembly = polca_fasta
            best_round = "pypolca"

    stats_dict = {}

    # count contigs
    number_of_contigs = 0

    # instantiate longest contig length
    longest_contig_length = 0

    total_assembly_length = 0

    # Open the output file in write mode
    with open(output_fasta, "w") as output_handle:
        # Iterate through the records in the best assembly FASTA file and write them to the output file
        for record in SeqIO.parse(best_assembly, "fasta"):
            # to match the 00001 output favoured generally for parsing

            number_of_contigs += 1

            record.id = f"contig{number_of_contigs:05}"

            # Calculate the length of the sequence
            sequence_length = len(record.seq)

            # total assembly length
            total_assembly_length += sequence_length

            # gc
            gc_content = round(gc_fraction(record.seq) * 100, 2)

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

            # append for stats dict

            stats_dict[record.id] = {}
            stats_dict[record.id]["length"] = sequence_length
            stats_dict[record.id]["gc"] = gc_content

    # read in the flye info and extract longest contig
    flye_df = pd.read_csv(flye_info, sep="\t")

    # Find the row with the largest length.
    longest_contig_row = flye_df[flye_df["length"] == flye_df["length"].max()]

    # Extract the coverage value from the longest contig row.
    longest_contig_coverage = longest_contig_row["cov."].values[0]

    # to get the summary df
    summary_dict = {
        "Sample": sample,
        "Complete": "False",
        "Total_assembly_length": total_assembly_length,
        "Number_of_contigs": number_of_contigs,
        "Most_accurate_polishing_round": best_round,
        "Longest_contig_length": longest_contig_length,
        "Longest_contig_coverage": longest_contig_coverage,
        "Number_circular_plasmids": "Unknown",
    }

    # Create a DataFrame from the dictionary
    summary_df = pd.DataFrame([summary_dict])
    summary_df.to_csv(hybracter_summary, index=False, sep="\t")

    # stats dict
    stats_df = pd.DataFrame.from_dict(stats_dict, orient="index")
    stats_df["contig_name"] = stats_df.index
    # Reorder the columns with 'contig_name' as the first column
    stats_df = stats_df[
        ["contig_name"] + [col for col in stats_df.columns if col != "contig_name"]
    ]

    stats_df.to_csv(per_conting_summary, index=False, sep="\t")


select_best_chromosome_assembly_incomplete(
    snakemake.output.hybracter_summary,
    snakemake.output.per_conting_summary,
    snakemake.params.ale_dir,
    snakemake.output.fasta,
    snakemake.output.ale_summary,
    snakemake.params.pre_polish_fasta,
    snakemake.params.medaka_fasta,
    snakemake.params.polypolish_fasta,
    snakemake.params.polca_fasta,
    snakemake.wildcards.sample,
    snakemake.input.flye_info,
    snakemake.params.logic,
    snakemake.params.no_pypolca,
)
