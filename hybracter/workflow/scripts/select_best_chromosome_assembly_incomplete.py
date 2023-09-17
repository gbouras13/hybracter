#!/usr/bin/env python3

import pandas as pd
import glob
import os 
import sys
from Bio import SeqIO


def select_best_chromosome_assembly_incomplete(ale_dir, output_fasta, ale_summary, pre_polish_fasta, medaka_fasta, polypolish_fasta, polca_fasta):
    """
    reads all the .score files in teh ale directory, picks the best one (closest to zero) and then takes that chromosome fasta and writes it to file with length 
    statistics similar to unicycler 
    instead of 1,2,3 etc we will use 'contig00001', 'contig00002' etc as more parsable (lessons from dragonflye)
    """

    # Use glob to find files with the .score extension in the directory
    print(ale_dir)
    file_list = glob.glob(os.path.join(ale_dir, '*.score'))
    print(file_list)

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

    
    #print(score_dict)
    # Filter out None values from the score_dict
    filtered_score_dict = {k: v for k, v in score_dict.items() if v is not None}

    #print(filtered_score_dict)

    # Check if there are any valid scores left
    if filtered_score_dict:
        # Find the key with the score closest to 0
        # this will be the best score from ALE
        closest_to_zero_key = min(filtered_score_dict, key=lambda k: abs(filtered_score_dict[k] - 0))
        closest_score = filtered_score_dict[closest_to_zero_key]

    print(closest_to_zero_key)

    # df with scores and files
    scores_df = pd.DataFrame(list(score_dict.items()), columns=['Key', 'Score'])
    scores_df.to_csv(ale_summary, index=False)


    # by default the best assembly is the polca fasta
    # check that the best assembly wasn't something else
    best_assembly = polca_fasta
    if "incomp_pre_polish" in closest_to_zero_key:
        best_assembly = pre_polish_fasta
    elif "medaka" in closest_to_zero_key:
        best_assembly = medaka_fasta
    elif "polypolish" in closest_to_zero_key:
        best_assembly = polypolish_fasta
    else: # polca 
        best_assembly = polca_fasta



    # in case there is multiple - counter
    i = 1


    # Open the output file in write mode
    with open(output_fasta, "w") as output_handle:
        # Iterate through the records in the best assembly FASTA file and write them to the output file
        for record in SeqIO.parse(best_assembly, "fasta"):

            # to match the 00001 output favoured generally for parsing
            # usually there will be 1 chromosome of course!
            record.id = f"contig{i:05}" 
            
            # Calculate the length of the sequence
            sequence_length = len(record.seq)
            
            # Update the description (header) with the length information
            record.description = f"len={sequence_length}"
            
            # Write the modified record to the output file
            SeqIO.write(record, output_handle, "fasta")


select_best_chromosome_assembly_incomplete(snakemake.params.ale_dir, snakemake.output.fasta, snakemake.output.ale_summary, snakemake.params.pre_polish_fasta, snakemake.params.medaka_fasta, snakemake.params.polypolish_fasta, snakemake.params.polca_fasta   )



