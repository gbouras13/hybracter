#!/usr/bin/env python3

import os

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

    try:
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
            total_length += len(gene.seq)  # add the length of the gene called

        mean_cds_len = float(total_length / total_genes)

    except Exception:  # in case there are no genes called or some other error
        mean_cds_len = 0

    return mean_cds_len
