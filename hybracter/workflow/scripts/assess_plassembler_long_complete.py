#!/usr/bin/env python3

import pandas as pd
import pyrodigal_helpers
from Bio import SeqIO


def determine_best_plassembler_assembly(
    plassembler_fasta,
    medaka_fasta,
    final_plasmid_fasta,
    plassembler_prodigal_summary,
    sample,
):
    if pyrodigal_helpers.is_file_empty(plassembler_fasta) is True:  # if empty, just create empty outputs
        pyrodigal_helpers.touch_file(plassembler_prodigal_summary)
        pyrodigal_helpers.touch_file(final_plasmid_fasta)

    else:  # plassembler not empty
        best_assembly = medaka_fasta

        plassembler_mean_cds = pyrodigal_helpers.calculate_mean_CDS_length(plassembler_fasta)
        medaka_mean_cds = pyrodigal_helpers.calculate_mean_CDS_length(medaka_fasta)

        dict = {
            "Sample": sample,
            "plassembler_raw_mean_cds_length": plassembler_mean_cds,
            "plassembler_polished_mean_cds_length": medaka_mean_cds,
        }

        # Convert the dictionary to a DataFrame
        summary_df = pd.DataFrame([dict])

        # Write to a CSV file
        summary_df.to_csv(plassembler_prodigal_summary, index=False)

        # determine the best assembly
        if plassembler_mean_cds > medaka_mean_cds:
            best_assembly = plassembler_fasta

        # set counter to 0 for number of plasmids
        plasmids = 0

        # Open the output file in write mode
        with open(final_plasmid_fasta, "w") as output_handle:
            # Iterate through the records in the best assembly FASTA file and write them to the output file
            for record in SeqIO.parse(best_assembly, "fasta"):
                plasmids += 1

                # to match the 00001 output favoured generally for parsing
                # usually there will be 1 chromosome of course!
                record.id = f"plasmid{plasmids:05}"

                # get rid off the contig id (1, 2, 3) and len from plassembler
                # extra)description will keep the circularity and copy number info
                extra_description = record.description.split(" ", 2)[2]

                # Calculate the length of the sequence
                sequence_length = len(record.seq)

                # Update the description (header) with the length information
                record.description = f"len={sequence_length} {extra_description}"

                # Write the modified record to the output file
                SeqIO.write(record, output_handle, "fasta")


determine_best_plassembler_assembly(
    snakemake.input.plassembler_fasta,
    snakemake.input.medaka_fasta,
    snakemake.output.final_plasmid_fasta,
    snakemake.output.plassembler_prodigal_summary,
    snakemake.wildcards.sample,
)
