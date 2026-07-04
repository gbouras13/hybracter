#!/usr/bin/env python3

import os
import sys

# snakemake only puts this script's own directory (no_medaka/) on sys.path, so
# add the parent scripts/ dir to import the shared library.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import select_best_lib as sbl


def select_best_chromosome_assembly_incomplete(
    hybracter_summary,
    per_contig_summary,
    ale_dir,
    output_fasta,
    ale_summary,
    pre_polish_fasta,
    polypolish_fasta,
    polca_fasta,
    sample,
    flye_info,
    logic,
    no_pypolca,
):
    """
    Short-read-polished, incomplete assembly, no-medaka path.

    Same as the medaka incomplete script but without the medaka round in the
    ALE-based selection. Shared IO lives in select_best_lib.
    """

    filtered_score_dict, closest_to_zero_key = sbl.read_ale_scores(ale_dir, ale_summary)

    # by default the best assembly is the final round (polca, or polypolish if no pypolca)
    if no_pypolca is False:
        best_assembly = polca_fasta
        best_round = "pypolca"
    else:
        best_assembly = polypolish_fasta
        best_round = "polypolish"

    if logic == "best":
        if not filtered_score_dict:
            print(
                f"WARNING: --logic best was set but no valid ALE scores were found; "
                f"keeping the final polishing round ({best_round})."
            )
        elif "incomp_pre_polish" in closest_to_zero_key:
            best_assembly = pre_polish_fasta
            best_round = "pre_polish"
        elif "polypolish" in closest_to_zero_key:
            best_assembly = polypolish_fasta
            best_round = "polypolish"
        else:  # polca
            best_assembly = polca_fasta
            best_round = "pypolca"

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
    select_best_chromosome_assembly_incomplete(
        snakemake.output.hybracter_summary,
        snakemake.output.per_contig_summary,
        snakemake.params.ale_dir,
        snakemake.output.fasta,
        snakemake.output.ale_summary,
        snakemake.params.pre_polish_fasta,
        snakemake.params.polypolish_fasta,
        snakemake.params.polca_fasta,
        snakemake.wildcards.sample,
        snakemake.input.flye_info,
        snakemake.params.logic,
        snakemake.params.no_pypolca,
    )
