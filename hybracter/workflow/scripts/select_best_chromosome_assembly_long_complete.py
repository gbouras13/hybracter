#!/usr/bin/env python3

import select_best_lib as sbl


def select_best_chromosome_assembly_long_complete(
    hybracter_summary,
    per_conting_summary,
    pyrodigal_summary,
    output_chromosome_fasta,
    overall_output_fasta,
    chrom_pre_polish_fasta,
    medaka_rd_1_fasta,
    medaka_rd_2_fasta,
    sample,
    flye_info,
    plassembler_fasta,
    final_plasmid_fasta,
    logic,
    ignore_list,
):
    """
    Long-read-only, complete (chromosome recovered) assembly.

    Selection is by pyrodigal mean CDS length (medaka can make things worse), so
    the default is the final medaka round and --logic best can step back to an
    earlier round. Shared IO lives in select_best_lib.
    """

    # mean CDS length per polishing round
    chrom_pre_polish_mean_cds = sbl.calculate_mean_CDS_length(chrom_pre_polish_fasta)
    medaka_rd_1_mean_cds = sbl.calculate_mean_CDS_length(medaka_rd_1_fasta)
    medaka_rd_2_mean_cds = sbl.calculate_mean_CDS_length(medaka_rd_2_fasta)

    sbl.write_single_row_tsv(
        pyrodigal_summary,
        {
            "Sample": sample,
            "pre_polish_mean_cds_length": chrom_pre_polish_mean_cds,
            "medaka_rd_1_mean_cds_length": medaka_rd_1_mean_cds,
            "medaka_rd_2_mean_cds_length": medaka_rd_2_mean_cds,
        },
    )

    # default: the final medaka round
    best_assembly = medaka_rd_2_fasta
    best_round = "medaka_rd_2"

    if logic == "best":
        if (
            medaka_rd_1_mean_cds > medaka_rd_2_mean_cds
            and chrom_pre_polish_mean_cds < medaka_rd_1_mean_cds
        ):
            best_assembly = medaka_rd_1_fasta
            best_round = "medaka_rd_1"

        if (
            medaka_rd_1_mean_cds > medaka_rd_2_mean_cds
            and medaka_rd_1_mean_cds < chrom_pre_polish_mean_cds
        ):
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
        plassembler_fasta,
        final_plasmid_fasta,
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
        snakemake.output.chromosome_fasta,
        snakemake.output.total_fasta,
        snakemake.input.chrom_pre_polish_fasta,
        snakemake.input.medaka_rd_1_fasta,
        snakemake.input.medaka_rd_2_fasta,
        snakemake.wildcards.sample,
        snakemake.input.flye_info,
        snakemake.input.plassembler_fasta,
        snakemake.output.final_plasmid_fasta,
        snakemake.params.logic,
        snakemake.input.ignore_list,
    )
