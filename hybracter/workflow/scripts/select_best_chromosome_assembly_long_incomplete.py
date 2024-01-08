#!/usr/bin/env python3

import pandas as pd
import pyrodigal
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


def calculate_mean_CDS_length(filepath_in):
    """
    calculates the mean CDS length in a fasta file and returns it as a float

    :param filepath_in: The path to the input FASTA file.
    :return: The mean length of sequences in the FASTA file as a float.
    """

    prodigal_metamode = False
    coding_table = 11

    # for training if you want different coding table
    seqs = [bytes(record.seq) for record in SeqIO.parse(filepath_in, "fasta")]
    record = SeqIO.parse(filepath_in, "fasta")

    # train pyrodigal
    orf_finder = pyrodigal.GeneFinder(meta=prodigal_metamode)

    too_short = False

    try:
        trainings_info = orf_finder.train(*seqs, translation_table=int(coding_table))
    except ValueError:
        # less than 20000bp - prodigal train will break
        too_short = True
        print("The assembly is less than 20000bp. Please get some more sequencing :)")
        mean_cds_len = 0.00

    if too_short is False:
        orf_finder = pyrodigal.GeneFinder(trainings_info, meta=prodigal_metamode)
        # run pyrodigal
        total_genes = 0
        total_length = 0
        for i, record in enumerate(SeqIO.parse(filepath_in, "fasta")):
            genes = orf_finder.find_genes(str(record.seq))
            for gene in genes:
                total_genes += 1
                total_length += len(
                    gene.sequence()
                )  # add the length of the gene called

        mean_cds_len = float(total_length / total_genes)
        mean_cds_len = float("{:.2f}".format(mean_cds_len))

    return mean_cds_len


def select_best_chromosome_assembly_long_incomplete(
    hybracter_summary,
    per_conting_summary,
    pyrodigal_summary,
    output_fasta,
    pre_polish_fasta,
    medaka_fasta,
    sample,
    flye_info,
    logic,
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

    # Write to a tsv file
    summary_df.to_csv(pyrodigal_summary, index=False, sep="\t")

    # determine the best assembly
    # plausible medaka makes it worse
    best_assembly = medaka_fasta
    best_round = "medaka"

    if logic == "best":
        if pre_polish_fasta > medaka_fasta:
            best_assembly = pre_polish_fasta
            best_round = "pre_polish"

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
            number_of_contigs += 1

            # to match the 00001 output favoured generally for parsing
            # usually there will be 1 chromosome of course!
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


select_best_chromosome_assembly_long_incomplete(
    snakemake.output.hybracter_summary,
    snakemake.output.per_conting_summary,
    snakemake.output.pyrodigal_summary,
    snakemake.output.total_fasta,
    snakemake.params.pre_polish_fasta,
    snakemake.params.medaka_fasta,
    snakemake.wildcards.sample,
    snakemake.input.flye_info,
    snakemake.params.logic,
)
