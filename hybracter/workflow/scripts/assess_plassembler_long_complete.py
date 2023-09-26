#!/usr/bin/env python3

import os
import sys

import pandas as pd
import pyrodigal
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


def calculate_mean_CDS_length(filepath_in):
    """
    calculates the mean CDS length in a fasta file and returns it as a float

    :param filepath_in: The path to the input FASTA file.
    :return: The mean length of sequences in the FASTA file as a float.
    """

    prodigal_metamode = False
    coding_table = 11

    # get total length of plasmids
    total_length = 0

    with open(filepath_in, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            total_length += len(record.seq)

    # if under 20000, pyrodigal will only work in meta mode
    # https://github.com/hyattpd/prodigal/wiki/Advice-by-Input-Type#plasmids-phages-viruses-and-other-short-sequences
    # https://github.com/hyattpd/Prodigal/issues/51
    # so make sure of this

    if total_length < 20001:
        prodigal_metamode = True
        orf_finder = pyrodigal.GeneFinder(meta=prodigal_metamode)
    else:  # plasmids over 20000 in length
        # try:
        # for training if you want different coding table
        seqs = [bytes(record.seq) for record in SeqIO.parse(filepath_in, "fasta")]
        record = SeqIO.parse(filepath_in, "fasta")

        # train pyrodigal
        orf_finder = pyrodigal.GeneFinder(meta=prodigal_metamode)
        trainings_info = orf_finder.train(*seqs, translation_table=int(coding_table))
        orf_finder = pyrodigal.GeneFinder(trainings_info, meta=prodigal_metamode)

    # run pyrodigal

    total_genes = 0
    total_length = 0

    for i, record in enumerate(SeqIO.parse(filepath_in, "fasta")):
        genes = orf_finder.find_genes(str(record.seq))
        for gene in genes:
            total_genes += 1
            total_length += len(gene.sequence())  # add the length of the gene called

    mean_cds_len = float(total_length / total_genes)
    mean_cds_len = float("{:.2f}".format(mean_cds_len))

    return mean_cds_len


def determine_best_plassembler_assembly(
    plassembler_fasta,
    medaka_fasta,
    final_plasmid_fasta,
    plassembler_prodigal_summary,
    sample,
):
    if is_file_empty(plassembler_fasta) is True:  # if empty, just create empty outputs
        touch_file(plassembler_prodigal_summary)
        touch_file(final_plasmid_fasta)

    else:  # plassembler not empty
        best_assembly = medaka_fasta

        plassembler_mean_cds = calculate_mean_CDS_length(plassembler_fasta)
        medaka_mean_cds = calculate_mean_CDS_length(medaka_fasta)

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

        # dict to store headers for the original plassembler assembly to get copy number info
        header_mapping = {}
        # read in headers in plassembler_fasta
        i = 1
        for record in SeqIO.parse(plassembler_fasta, "fasta"):
            header_mapping[i] = record.description
            i = +1

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
                extra_description = header_mapping[i].split(" ", 2)[2]

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
