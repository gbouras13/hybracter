#!/usr/bin/env python3

import glob
import os
import sys

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


# determines whether a file is empty
def is_file_empty(file):
    """
    Determines if file is empty
    :param file: file path
    :return: empty Boolean
    """
    empty = False
    if os.stat(file).st_size == 0:
        empty = True
    return empty


# touches an empty file
def touch_file(path):
    with open(path, "a"):
        os.utime(path, None)


def select_best_chromosome_assembly_complete(
    hybracter_summary,
    per_conting_summary,
    ale_dir,
    input_plassember_fasta,
    output_chromosome_fasta,
    output_plasmid_fasta,
    overall_output_fasta,
    ale_summary,
    chrom_pre_polish_fasta,
    polypolish_fasta,
    polca_fasta,
    sample,
    flye_info,
    logic,
    no_pypolca,
    ignore_list
):
    """
    reads all the .score files in teh ale directory, picks the best one (closest to zero) and then takes that chromosome fasta and writes it to file with length
    statistics similar to unicycler
    instead of 1,2,3 etc we will use 'chromosome00001', 'chromosome00002' etc (for edge cases of multiple chroms/megaplasmids chromids etc)
    Then it reads the plassembler output
    Checks if not empty
    And if it isn't, then adds the plassembler contigs as 'plasmid00001' etc
    Otherwise the plasmid output is empty

    Then creates hybracter_summary tsv
    """

    # Use glob to find files with the .score extension in the directory
    file_list = glob.glob(os.path.join(ale_dir, "*.score"))

    # Create an empty dictionary to store the summary info
    summary_dict = {}

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

    # Find the key associated with the best score
    best_round = None

    # Check if there are any valid scores left
    if filtered_score_dict:
        # Find the key with the score closest to 0
        # this will be the best score from ALE
        closest_to_zero_key = min(
            filtered_score_dict, key=lambda k: abs(filtered_score_dict[k] - 0)
        )
        best_score = filtered_score_dict[closest_to_zero_key]

    # df with scores and files
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

    # if best is chosen then run the logic
    if logic == "best":
        if "chrom_pre_polish" in closest_to_zero_key:
            best_assembly = chrom_pre_polish_fasta
            best_round = "pre_polish"
        elif "polypolish" in closest_to_zero_key:
            best_assembly = polypolish_fasta
            best_round = "polypolish"
        else:  # polca
            best_assembly = polca_fasta
            best_round = "pypolca"

    stats_dict = {}

    # write the chromosome(s)
    # usually should be 1!
    number_of_chromosomes = sum(1 for _ in SeqIO.parse(best_assembly, "fasta"))
    if number_of_chromosomes <= 0:
        sys.exit(f"The assembly FASTA {best_assembly} is empty.")

    # in case there is multiple - counter
    chromosomes = 1

    # instantiate longest contig length
    longest_contig_length = 0

    # total assembly length
    total_assembly_length = 0

    # ignore list 
    # strip '\n'
    with open(ignore_list, 'r') as file:
        non_circular_chromosome_contig_list = [line.strip() for line in file.readlines()]

    # Open the output file in write mode
    with open(output_chromosome_fasta, "w") as output_handle:
        with open(overall_output_fasta, "w") as output_handle_overall:
            # Iterate through the records in the best assembly FASTA file and write them to the output file
            for record in SeqIO.parse(best_assembly, "fasta"):

                # for the circularity match at the end
                contig_id = record.id

                # to match the 00001 output favoured generally for parsing
                # usually there will be 1 chromosome of course!
                record.id = f"chromosome{chromosomes:05}"

                # Calculate the length of the sequence
                sequence_length = len(record.seq)

                # gc
                gc_content = round(gc_fraction(record.seq) * 100, 2)

                # total assembly length
                total_assembly_length += sequence_length

                # to get longest contig
                if number_of_chromosomes == 1:
                    longest_contig_length = sequence_length
                else:
                    if sequence_length > longest_contig_length:
                        longest_contig_length = sequence_length

                # Update the description (header) with the length information
                if contig_id in non_circular_chromosome_contig_list:
                    record.description = f"len={sequence_length}"
                else:
                    record.description = f"len={sequence_length} circular=true"

                # Write the modified record to the output file
                SeqIO.write(record, output_handle, "fasta")
                SeqIO.write(record, output_handle_overall, "fasta")

                # append for stats dict

                stats_dict[record.id] = {}
                stats_dict[record.id]["contig_type"] = "chromosome"
                stats_dict[record.id]["length"] = sequence_length
                stats_dict[record.id]["gc"] = gc_content
                if contig_id in non_circular_chromosome_contig_list:
                    stats_dict[record.id]["circular"] = "False"
                else:
                    stats_dict[record.id]["circular"] = "True"

                chromosomes += 1

    #######################
    # plasmid
    #######################

    # set counter to 0 for number of plasmids
    plasmids = 0
    circular_plasmids = 0

    if (
        is_file_empty(input_plassember_fasta) is False
    ):  # if the plassembler output is not empty
        # Open the output file in write mode
        with open(output_plasmid_fasta, "w") as output_handle:
            with open(
                overall_output_fasta, "a"
            ) as output_handle_overall:  # needs to be append
                # Iterate through the records in the best assembly FASTA file and write them to the output file
                for record in SeqIO.parse(input_plassember_fasta, "fasta"):
                    plasmids += 1

                    # to match the 00001 output favoured generally for parsing
                    # usually there will be 1 chromosome of course!
                    record.id = f"plasmid{plasmids:05}"

                    sequence_length = len(record.seq)

                    # add record length
                    total_assembly_length += sequence_length

                    # gc
                    gc_content = round(gc_fraction(record.seq) * 100, 2)

                    # get rid off the contig id (1, 2, 3) etc from plassembler
                    record.description = record.description.split(" ", 1)[1]

                    completeness_flag = False
                    if "circular=true" in record.description:
                        completeness_flag = True
                        circular_plasmids += 1

                    # take description from plassembler (length and copy number)
                    # Write the modified record to the output file
                    SeqIO.write(record, output_handle, "fasta")
                    SeqIO.write(record, output_handle_overall, "fasta")

                    stats_dict[record.id] = {}
                    stats_dict[record.id]["contig_type"] = "plasmid"
                    stats_dict[record.id]["length"] = sequence_length
                    stats_dict[record.id]["gc"] = gc_content
                    stats_dict[record.id]["circular"] = str(completeness_flag)

    else:
        touch_file(output_plasmid_fasta)

    number_of_contigs = chromosomes + plasmids - 1

    # read in the flye info and extract longest contig
    # Read the TSV file into a Pandas DataFrame.
    flye_df = pd.read_csv(flye_info, sep="\t")

    # Find the row with the largest length.
    longest_contig_row = flye_df[flye_df["length"] == flye_df["length"].max()]

    # Extract the coverage value from the longest contig row.
    longest_contig_coverage = longest_contig_row["cov."].values[0]

    # to get the summary df
    summary_dict = {
        "Sample": sample,
        "Complete": "True",
        "Total_assembly_length": total_assembly_length,
        "Number_of_contigs": number_of_contigs,
        "Most_accurate_polishing_round": best_round,
        "Longest_contig_length": longest_contig_length,
        "Longest_contig_coverage": longest_contig_coverage,
        "Number_circular_plasmids": int(circular_plasmids),
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


select_best_chromosome_assembly_complete(
    snakemake.output.hybracter_summary,
    snakemake.output.per_conting_summary,
    snakemake.params.ale_dir,
    snakemake.input.plassembler_fasta,
    snakemake.output.chromosome_fasta,
    snakemake.output.plasmid_fasta,
    snakemake.output.total_fasta,
    snakemake.output.ale_summary,
    snakemake.input.chrom_pre_polish_fasta,
    snakemake.params.polypolish_fasta,
    snakemake.params.polca_fasta,
    snakemake.wildcards.sample,
    snakemake.input.flye_info,
    snakemake.params.logic,
    snakemake.params.no_pypolca,
    snakemake.input.ignore_list
)
