#!/usr/bin/env python3

import pandas as pd
import glob
import os 
import sys
from Bio import SeqIO

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


def select_best_chromosome_assembly_complete(ale_dir, input_plassember_fasta, output_chromosome_fasta, output_plasmid_fasta, overall_output_fasta, ale_summary, chrom_pre_polish_fasta, medaka_rd_1_fasta, medaka_rd_2_fasta, polypolish_fasta, polca_fasta):
    """
    reads all the .score files in teh ale directory, picks the best one (closest to zero) and then takes that chromosome fasta and writes it to file with length 
    statistics similar to unicycler 
    instead of 1,2,3 etc we will use 'chromosome00001', 'chromosome00002' etc (for edge cases of multiple chroms/megaplasmids chromids etc)
    Then it reads the plassembler output
    Checks if not empty
    And if it isn't, then adds the plassembler contigs as 'plasmid00001' etc
    Otherwise the plasmid output is empty
    """

    # Use glob to find files with the .score extension in the directory
    file_list = glob.glob(os.path.join(ale_dir, '*.score'))

    # Create an empty dictionary to store the results
    score_dict = {}

    for file_path in file_list:
        # Strip the ".score" extension and use it as the dictionary key
        file_name = os.path.splitext(os.path.basename(file_path))[0]

        # Initialize the score as None (in case the file doesn't contain a valid score)
        score = None
        
        # Read the first line of the file
        with open(file_path, 'r') as file:
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

    # Check if there are any valid scores left
    if filtered_score_dict:
        # Find the key with the score closest to 0
        # this will be the best score from ALE
        closest_to_zero_key = min(filtered_score_dict, key=lambda k: abs(filtered_score_dict[k] - 0))
        closest_score = filtered_score_dict[closest_to_zero_key]

    # df with scores and files
    scores_df = pd.DataFrame(list(score_dict.items()), columns=['Key', 'Score'])
    scores_df.to_csv(ale_summary, index=False)


    # by default the best assembly is the polca fasta
    # check that the best assembly wasn't something else
    best_assembly = polca_fasta
    if "chrom_pre_polish" in closest_to_zero_key:
        best_assembly = chrom_pre_polish_fasta
    elif "medaka_rd_1" in closest_to_zero_key:
        best_assembly = medaka_rd_1_fasta
    elif "medaka_rd_2" in closest_to_zero_key:
        best_assembly = medaka_rd_2_fasta
    elif "polypolish" in closest_to_zero_key:
        best_assembly = polypolish_fasta
    else: # polca 
        best_assembly = polca_fasta


    # write the chromosome(s)
    # usually should be 1!
    number_of_chromosomes = sum(1 for _ in SeqIO.parse(best_assembly, "fasta"))
    if number_of_chromosomes <= 0:
        sys.exit(f"The assembly FASTA {best_assembly} is empyu")

    # in case there is multiple - counter
    i = 1

    # Open the output file in write mode
    with open(output_chromosome_fasta, "w") as output_handle:
        with open(overall_output_fasta, "w") as output_handle_overall:
            # Iterate through the records in the best assembly FASTA file and write them to the output file
            for record in SeqIO.parse(best_assembly, "fasta"):

                # to match the 00001 output favoured generally for parsing
                # usually there will be 1 chromosome of course!
                record.id = f"chromosome{i:05}" 
                
                # Calculate the length of the sequence
                sequence_length = len(record.seq)
                
                # Update the description (header) with the length information
                record.description = f"len={sequence_length}"
                
                # Write the modified record to the output file
                SeqIO.write(record, output_handle, "fasta")
                SeqIO.write(record, output_handle_overall, "fasta")

    #######################
    # plasmid
    #######################

    # reset counter to 1
    i = 1

    if is_file_empty(input_plassember_fasta) is False:  # if the plassembler output is not empty    
        # Open the output file in write mode
        with open(output_plasmid_fasta, "w") as output_handle:
            with open(overall_output_fasta, "a") as output_handle_overall: # needs to be append
                # Iterate through the records in the best assembly FASTA file and write them to the output file
                for record in SeqIO.parse(input_plassember_fasta, "fasta"):

                    # to match the 00001 output favoured generally for parsing
                    # usually there will be 1 chromosome of course!
                    record.id = f"plasmid{i:05}" 
                    
                    # take description from plassembler (length and copy number)
                    # Write the modified record to the output file
                    SeqIO.write(record, output_handle, "fasta")
                    SeqIO.write(record, output_handle_overall, "fasta")
    else:
        touch_file(output_plasmid_fasta)


select_best_chromosome_assembly_complete(snakemake.params.ale_dir, snakemake.input.plassembler_fasta, snakemake.output.chromosome_fasta, snakemake.output.plasmid_fasta, snakemake.output.total_fasta, snakemake.output.ale_summary, snakemake.params.chrom_pre_polish_fasta, snakemake.params.medaka_rd_1_fasta, snakemake.params.medaka_rd_2_fasta, snakemake.params.polypolish_fasta, snakemake.params.polca_fasta   )



