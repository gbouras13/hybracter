#!/usr/bin/env python3

import os
import sys

# snakemake only puts this script's own directory (no_medaka/) on sys.path, so
# add the parent scripts/ dir to import the shared library.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import select_best_lib as sbl


def select_best_chromosome_assembly_long_complete(
    hybracter_summary,
    per_conting_summary,
    pyrodigal_summary,
    plassembler_plasmid_fasta,
    output_chromosome_fasta,
    output_plasmid_fasta,
    overall_output_fasta,
    chrom_pre_polish_fasta,
    sample,
    flye_info,
    ignore_list,
):
    """
    Long-read-only, complete assembly, no-medaka path.

    With no medaka rounds there is nothing to select between: the pre-polish
    chromosome is the answer. Shared IO lives in select_best_lib.
    """

    chrom_pre_polish_mean_cds = sbl.calculate_mean_CDS_length(chrom_pre_polish_fasta)

    sbl.write_single_row_tsv(
        pyrodigal_summary,
        {
            "Sample": sample,
            "pre_polish_mean_cds_length": chrom_pre_polish_mean_cds,
        },
    )

    best_assembly = chrom_pre_polish_fasta
    best_round = "pre_polish"

    stats_dict = {}

    chromosomes, longest_contig_length, total_assembly_length = sbl.write_chromosomes(
        best_assembly,
        output_chromosome_fasta,
        overall_output_fasta,
        ignore_list,
        stats_dict,
    )

    plasmids, circular_plasmids, total_assembly_length = sbl.write_plasmids(
        plassembler_plasmid_fasta,
        output_plasmid_fasta,
        overall_output_fasta,
        stats_dict,
        total_assembly_length,
    )

    number_of_contigs = chromosomes + plasmids - 1

    sbl.write_summary(
        hybracter_summary,
        sample,
        True,
        total_assembly_length,
        number_of_contigs,
        best_round,
        longest_contig_length,
        sbl.longest_contig_coverage(flye_info),
        int(circular_plasmids),
    )

    sbl.write_per_contig_stats(per_conting_summary, stats_dict)


if "snakemake" in globals():  # only runs under Snakemake; lets pytest import this module
    select_best_chromosome_assembly_long_complete(
        snakemake.output.hybracter_summary,
        snakemake.output.per_conting_summary,
        snakemake.output.pyrodigal_summary,
        snakemake.input.plassembler_plasmid_fasta,
        snakemake.output.chromosome_fasta,
        snakemake.output.plasmid_fasta,
        snakemake.output.total_fasta,
        snakemake.input.chrom_pre_polish_fasta,
        snakemake.wildcards.sample,
        snakemake.input.flye_info,
        snakemake.input.ignore_list,
    )
