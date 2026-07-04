#!/usr/bin/env python3

import os
import sys

# snakemake only puts this script's own directory (no_medaka/) on sys.path, so
# add the parent scripts/ dir to import the shared library.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import select_best_lib as sbl


def select_best_chromosome_assembly_long_incomplete(
    hybracter_summary,
    per_contig_summary,
    pyrodigal_summary,
    output_fasta,
    pre_polish_fasta,
    sample,
    flye_info,
):
    """
    Long-read-only, incomplete assembly, no-medaka path.

    With no medaka rounds the pre-polish assembly is the answer. Shared IO lives
    in select_best_lib.
    """

    pre_polish_mean_cds = sbl.calculate_mean_CDS_length(pre_polish_fasta)

    sbl.write_single_row_tsv(
        pyrodigal_summary,
        {
            "Sample": sample,
            "pre_polish_mean_cds_length": pre_polish_mean_cds,
        },
    )

    best_assembly = pre_polish_fasta
    best_round = "pre_polish"

    stats_dict = {}

    number_of_contigs, longest_contig_length, total_assembly_length = sbl.write_contigs(
        best_assembly, output_fasta, stats_dict
    )

    sbl.write_summary(
        hybracter_summary,
        sample,
        False,
        total_assembly_length,
        number_of_contigs,
        best_round,
        longest_contig_length,
        sbl.longest_contig_coverage(flye_info),
        "Unknown",
    )

    sbl.write_per_contig_stats(per_contig_summary, stats_dict)


if "snakemake" in globals():  # only runs under Snakemake; lets pytest import this module
    select_best_chromosome_assembly_long_incomplete(
        snakemake.output.hybracter_summary,
        snakemake.output.per_contig_summary,
        snakemake.output.pyrodigal_summary,
        snakemake.output.total_fasta,
        snakemake.params.pre_polish_fasta,
        snakemake.wildcards.sample,
        snakemake.input.flye_info,
    )
