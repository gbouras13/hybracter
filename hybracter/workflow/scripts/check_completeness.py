#!/usr/bin/env python3

from Bio import SeqIO

# get chromosome if length > 2.5 million


def circ_by_seq_name(info_file_path):
    """Map each Flye contig seq_name -> its circ flag ('Y'/'N') from assembly_info.txt.

    Replaces a pandas read of the tab-separated Flye info table
    (cols: seq_name, length, cov, circ, ...); '#' comment lines are skipped.
    """
    circ = {}
    with open(info_file_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            circ[fields[0]] = fields[3]
    return circ


def get_completeness(input_fasta, completeness_check, min_chrom_length, info_file_path, circular_chromosome):
    # Read the flye info file (seq_name -> circ flag)
    circ_by_name = circ_by_seq_name(info_file_path)

    # flag for if the assembly is complete
    comp = False
    for dna_record in SeqIO.parse(input_fasta, "fasta"):
        circ_value = circ_by_name[dna_record.id]
        if circular_chromosome:
            # --circular_chromosome: complete requires length AND Flye circularity
            if len(dna_record.seq) > int(min_chrom_length) and circ_value == "Y":
                comp = True
        else:
            # default: complete means length above threshold only
            if len(dna_record.seq) > int(min_chrom_length):
                comp = True

    with open(completeness_check, "w") as f:
        if comp == True:
            f.write("C")
        else:
            f.write("I")


if "snakemake" in globals():  # only runs under Snakemake; lets pytest import this module
    get_completeness(
        snakemake.input.fasta,
        snakemake.output.completeness_check,
        snakemake.params.min_chrom_length,
        snakemake.input.info,
        snakemake.params.circular_chromosome,
    )
