#!/usr/bin/env python3

import select_best_lib as sbl


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
    Long-read-only, incomplete assembly.

    Selection is by pyrodigal mean CDS length; all contigs are written to a
    single FASTA (contig00001 ...). Shared IO lives in select_best_lib.
    """

    pre_polish_mean_cds = sbl.calculate_mean_CDS_length(pre_polish_fasta)
    medaka_mean_cds = sbl.calculate_mean_CDS_length(medaka_fasta)

    sbl.write_single_row_tsv(
        pyrodigal_summary,
        {
            "Sample": sample,
            "pre_polish_mean_cds_length": pre_polish_mean_cds,
            "medaka_mean_cds_length": medaka_mean_cds,
        },
    )

    # default: medaka (plausibly medaka makes it worse)
    best_assembly = medaka_fasta
    best_round = "medaka"

    if logic == "best":
        # NOTE: pre-existing latent bug, preserved verbatim by the E1 refactor.
        # This compares the FASTA *paths*, not the mean CDS lengths; it should be
        # `pre_polish_mean_cds > medaka_mean_cds`. Tracked separately (see REVIEW.md);
        # fixing it here would change behaviour, so it is intentionally untouched.
        if pre_polish_fasta > medaka_fasta:
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

    sbl.write_per_contig_stats(per_conting_summary, stats_dict)


if "snakemake" in globals():  # only runs under Snakemake; lets pytest import this module
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
