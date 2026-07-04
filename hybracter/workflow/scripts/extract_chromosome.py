#!/usr/bin/env python3

from Bio import SeqIO
import os

"""
Extract chromosomes for the complete assembly path.
Default: extract all contigs above min_chrom_length regardless of circularity.
Non-circular contigs are added to ignore_list so dnaapler skips reorientation for them.
With --circular_chromosome: only extract contigs that are BOTH above min_chrom_length AND
marked circular by Flye.
"""

# touches an empty file
def touch_file(path):
    with open(path, "a"):
        os.utime(path, None)



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


def get_chromosome_plasmids(
    input_fasta, chromosome_fasta, ignore_list, min_chrom_length, info_file_path, polypolish_flag, circular_chromosome
):
    # read in the fasta

    # Read the flye info file (seq_name -> circ flag)
    circ_by_name = circ_by_seq_name(info_file_path)

    # to make sure it gets made regardless for snakemake
    touch_file(ignore_list)

    with open(chromosome_fasta, "w") as fa:
        for dna_record in SeqIO.parse(input_fasta, "fasta"):

            seq_id = dna_record.id

            # look up the Flye circularity; absent (chromosome didn't circularise,
            # or a plasmid not part of the chromosome set) defaults to 'N' so it
            # isn't extracted under --circular_chromosome
            circ_value = circ_by_name.get(seq_id, "N")

            if circular_chromosome:
                # --circular_chromosome: only extract contigs that are long AND circular
                if len(dna_record.seq) > int(min_chrom_length) and circ_value == "Y":
                    SeqIO.write(dna_record, fa, "fasta")
            else:
                # default: extract all long contigs; add non-circular to ignore_list
                # so dnaapler skips reorientation for them
                if len(dna_record.seq) > int(min_chrom_length):
                    SeqIO.write(dna_record, fa, "fasta")
                    if circ_value == "N":
                        with open(ignore_list, "a") as file:
                            file.write(f"{seq_id}\n")


if "snakemake" in globals():  # only runs under Snakemake; lets pytest import this module
    get_chromosome_plasmids(
        snakemake.input.fasta,
        snakemake.output.fasta,
        snakemake.output.ignore_list,
        snakemake.params.min_chrom_length,
        snakemake.input.info,
        snakemake.params.polypolish_flag,
        snakemake.params.circular_chromosome,
    )
