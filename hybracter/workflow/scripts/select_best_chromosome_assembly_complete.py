#!/usr/bin/env python3

import select_best_lib as sbl


def select_best_chromosome_assembly_complete(
    hybracter_summary,
    per_conting_summary,
    ale_dir,
    input_plassember_fasta,
    output_chromosome_fasta,
    output_plasmid_fasta,
    overall_output_fasta,
    ale_summary,
    chrom_pre_polish_fasta,
    medaka_rd_1_fasta,
    medaka_rd_2_fasta,
    polypolish_fasta,
    polca_fasta,
    sample,
    flye_info,
    logic,
    no_pypolca,
    ignore_list,
):
    """
    Short-read-polished, complete (chromosome recovered) assembly.

    Reads the ALE .score files, picks the round closest to zero (when --logic
    best), writes the chosen chromosome(s) and the plassembler plasmids with
    unicycler-style naming (chromosome00001, plasmid00001 ...), then writes the
    hybracter summary. Shared IO lives in select_best_lib; only the ALE-based
    best-round selection is script-specific.
    """

    filtered_score_dict, closest_to_zero_key = sbl.read_ale_scores(ale_dir, ale_summary)

    # by default the best assembly is the final round (polca, or polypolish if no pypolca)
    if no_pypolca is False:
        best_assembly = polca_fasta
        best_round = "pypolca"
    else:
        best_assembly = polypolish_fasta
        best_round = "polypolish"

    # if best is chosen then pick the round with the ALE score closest to zero
    if logic == "best":
        if not filtered_score_dict:
            print(
                f"WARNING: --logic best was set but no valid ALE scores were found; "
                f"keeping the final polishing round ({best_round})."
            )
        elif "chrom_pre_polish" in closest_to_zero_key:
            best_assembly = chrom_pre_polish_fasta
            best_round = "pre_polish"
        elif "medaka_rd_1" in closest_to_zero_key:
            best_assembly = medaka_rd_1_fasta
            best_round = "medaka_rd_1"
        elif "medaka_rd_2" in closest_to_zero_key:
            best_assembly = medaka_rd_2_fasta
            best_round = "medaka_rd_2"
        elif "polypolish" in closest_to_zero_key:
            best_assembly = polypolish_fasta
            best_round = "polypolish"
        else:  # polca
            best_assembly = polca_fasta
            best_round = "pypolca"

    stats_dict = {}

    chromosomes, longest_contig_length, total_assembly_length = sbl.write_chromosomes(
        best_assembly,
        output_chromosome_fasta,
        overall_output_fasta,
        ignore_list,
        stats_dict,
    )

    plasmids, circular_plasmids, total_assembly_length = sbl.write_plasmids(
        input_plassember_fasta,
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
    select_best_chromosome_assembly_complete(
        snakemake.output.hybracter_summary,
        snakemake.output.per_conting_summary,
        snakemake.params.ale_dir,
        snakemake.input.plassembler_fasta,
        snakemake.output.chromosome_fasta,
        snakemake.output.plasmid_fasta,
        snakemake.output.total_fasta,
        snakemake.output.ale_summary,
        snakemake.input.chrom_pre_polish_fasta,
        snakemake.params.medaka_rd_1_fasta,
        snakemake.params.medaka_rd_2_fasta,
        snakemake.params.polypolish_fasta,
        snakemake.params.polca_fasta,
        snakemake.wildcards.sample,
        snakemake.input.flye_info,
        snakemake.params.logic,
        snakemake.params.no_pypolca,
        snakemake.input.ignore_list,
    )
