#!/usr/bin/env python3
"""
Shared helpers for the select_best_chromosome_assembly_* scripts.

Those scripts come in eight variants (complete/incomplete x short[ALE]/long[pyrodigal]
x medaka/no_medaka) that historically duplicated ~95% of their bodies. That is how
bugs like B1 (missing ALE-score guard) and B4 (plasmid circularity test) came to
exist in several copies and drift apart. The genuinely-identical IO / parsing /
writing logic lives here; each script keeps only its distinctive best-assembly
*selection* logic and calls these helpers.

Import:
  * scripts in this directory:  ``import select_best_lib as sbl``
  * scripts in ``no_medaka/``:  prepend the parent dir to sys.path via ``__file__``
    first (snakemake only puts the script's own directory on sys.path), e.g.::

        import os, sys
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        import select_best_lib as sbl
"""

import glob
import os
import sys

import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


# determines whether a file is empty
def is_file_empty(file):
    """Return True if the file is zero bytes."""
    return os.stat(file).st_size == 0


# touches an empty file
def touch_file(path):
    """Create an empty file (or update its mtime if it already exists)."""
    with open(path, "a"):
        os.utime(path, None)


def calculate_mean_CDS_length(filepath_in):
    """Calculate the mean CDS length (bp) in a FASTA file, as a float.

    Trains pyrodigal on the sequences and averages the predicted gene lengths.
    Assemblies under ~20 kb are too short for pyrodigal training and return 0.00.

    :param filepath_in: path to the input FASTA file.
    :return: mean CDS length as a float (0.00 if the assembly is too short).
    """
    import pyrodigal

    prodigal_metamode = False
    coding_table = 11

    seqs = [bytes(record.seq) for record in SeqIO.parse(filepath_in, "fasta")]

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
        total_genes = 0
        total_length = 0
        for record in SeqIO.parse(filepath_in, "fasta"):
            genes = orf_finder.find_genes(str(record.seq))
            for gene in genes:
                total_genes += 1
                total_length += len(gene.sequence())  # add the length of the gene called

        mean_cds_len = float(total_length / total_genes)
        mean_cds_len = float("{:.2f}".format(mean_cds_len))

    return mean_cds_len


def read_ale_scores(ale_dir, ale_summary):
    """Read all ``*.score`` files in ``ale_dir``, write the sorted summary TSV, and
    return ``(filtered_score_dict, closest_to_zero_key)``.

    ``filtered_score_dict`` drops files whose first line is not a valid float;
    ``closest_to_zero_key`` is the key of the best (closest-to-zero) ALE score, or
    ``None`` when there are no valid scores (B1 - callers must handle ``None``).
    """
    file_list = glob.glob(os.path.join(ale_dir, "*.score"))

    score_dict = {}
    for file_path in file_list:
        file_name = os.path.splitext(os.path.basename(file_path))[0]
        score = None
        with open(file_path, "r") as file:
            first_line = file.readline().strip()
        try:
            score = float(first_line)
        except ValueError:
            pass  # not a valid float -> score stays None
        score_dict[file_name] = score

    filtered_score_dict = {k: v for k, v in score_dict.items() if v is not None}

    closest_to_zero_key = None
    if filtered_score_dict:
        closest_to_zero_key = min(
            filtered_score_dict, key=lambda k: abs(filtered_score_dict[k] - 0)
        )

    # df with scores and files; sorted ascending (worst top, best bottom)
    scores_df = pd.DataFrame(list(score_dict.items()), columns=["Key", "Score"])
    scores_df.sort_values(by="Score", ascending=True, inplace=True)
    scores_df.to_csv(ale_summary, index=False, sep="\t")

    return filtered_score_dict, closest_to_zero_key


def longest_contig_coverage(flye_info):
    """Return the coverage of the longest contig from a Flye assembly_info.txt."""
    flye_df = pd.read_csv(flye_info, sep="\t")
    longest_contig_row = flye_df[flye_df["length"] == flye_df["length"].max()]
    return longest_contig_row["cov."].values[0]


def write_chromosomes(
    best_assembly,
    output_chromosome_fasta,
    overall_output_fasta,
    ignore_list,
    stats_dict,
):
    """Write the chromosome records of a *complete* assembly.

    Records are renamed ``chromosome00001`` etc and written to both the
    chromosome FASTA (truncated) and the overall FASTA (truncated - it is opened
    "w" here and appended to later by :func:`write_plasmids`). Circularity comes
    from ``ignore_list`` (contigs listed there are non-circular). Per-contig stats
    are recorded into ``stats_dict``.

    :return: ``(chromosomes, longest_contig_length, total_assembly_length)`` where
        ``chromosomes`` is the post-loop counter (1 + number of chromosomes),
        preserving the historical ``number_of_contigs = chromosomes + plasmids - 1``.
    """
    number_of_chromosomes = sum(1 for _ in SeqIO.parse(best_assembly, "fasta"))
    if number_of_chromosomes <= 0:
        sys.exit(f"The assembly FASTA {best_assembly} is empty.")

    # ignore list - strip '\n'
    with open(ignore_list, "r") as file:
        non_circular_chromosome_contig_list = [line.strip() for line in file.readlines()]

    chromosomes = 1
    longest_contig_length = 0
    total_assembly_length = 0

    with open(output_chromosome_fasta, "w") as output_handle:
        with open(overall_output_fasta, "w") as output_handle_overall:
            for record in SeqIO.parse(best_assembly, "fasta"):
                # keep the original id for the circularity match
                contig_id = record.id
                # rename to the 00001 form favoured for parsing
                record.id = f"chromosome{chromosomes:05}"

                sequence_length = len(record.seq)
                gc_content = round(gc_fraction(record.seq) * 100, 2)
                total_assembly_length += sequence_length

                if number_of_chromosomes == 1:
                    longest_contig_length = sequence_length
                elif sequence_length > longest_contig_length:
                    longest_contig_length = sequence_length

                if contig_id in non_circular_chromosome_contig_list:
                    record.description = f"len={sequence_length}"
                else:
                    record.description = f"len={sequence_length} circular=true"

                SeqIO.write(record, output_handle, "fasta")
                SeqIO.write(record, output_handle_overall, "fasta")

                stats_dict[record.id] = {
                    "contig_type": "chromosome",
                    "length": sequence_length,
                    "gc": gc_content,
                    "circular": (
                        "False"
                        if contig_id in non_circular_chromosome_contig_list
                        else "True"
                    ),
                }

                chromosomes += 1

    return chromosomes, longest_contig_length, total_assembly_length


def write_plasmids(
    plassembler_fasta,
    output_plasmid_fasta,
    overall_output_fasta,
    stats_dict,
    total_assembly_length,
):
    """Append plassembler plasmid records to a *complete* assembly's output.

    Records are renamed ``plasmid00001`` etc and appended to both the plasmid
    FASTA and the (already-written) overall FASTA. Circularity uses a
    case-insensitive ``circular=true`` test (B4): plassembler writes
    ``circular=true`` on the short/hybrid path but ``circular=True`` on the
    long-only path, and both must count (while ``circular=false`` must not). If
    the plassembler output is empty the plasmid FASTA is just touched.

    :return: ``(plasmids, circular_plasmids, total_assembly_length)``.
    """
    plasmids = 0
    circular_plasmids = 0

    if is_file_empty(plassembler_fasta) is False:  # plassembler output not empty
        with open(output_plasmid_fasta, "w") as output_handle:
            with open(overall_output_fasta, "a") as output_handle_overall:  # append
                for record in SeqIO.parse(plassembler_fasta, "fasta"):
                    plasmids += 1
                    record.id = f"plasmid{plasmids:05}"
                    sequence_length = len(record.seq)
                    total_assembly_length += sequence_length
                    gc_content = round(gc_fraction(record.seq) * 100, 2)

                    # drop the plassembler contig id (1, 2, 3 ...) from the header
                    record.description = record.description.split(" ", 1)[1]

                    # Case-insensitive: plassembler emits `circular=true` on the
                    # short/hybrid path but `circular=True` on the long-only path
                    # (and via canu). Match both, but still not `circular=false`.
                    completeness_flag = False
                    if "circular=true" in record.description.lower():
                        completeness_flag = True
                        circular_plasmids += 1

                    SeqIO.write(record, output_handle, "fasta")
                    SeqIO.write(record, output_handle_overall, "fasta")

                    stats_dict[record.id] = {
                        "contig_type": "plasmid",
                        "length": sequence_length,
                        "gc": gc_content,
                        "circular": str(completeness_flag),
                    }
    else:
        touch_file(output_plasmid_fasta)

    return plasmids, circular_plasmids, total_assembly_length


def write_contigs(best_assembly, output_fasta, stats_dict):
    """Write the records of an *incomplete* assembly.

    Records are renamed ``contig00001`` etc and written to a single output FASTA
    (no chromosome/plasmid split, no circularity). Per-contig stats (length, gc)
    are recorded into ``stats_dict``.

    :return: ``(number_of_contigs, longest_contig_length, total_assembly_length)``.
    """
    number_of_contigs = 0
    longest_contig_length = 0
    total_assembly_length = 0

    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(best_assembly, "fasta"):
            number_of_contigs += 1
            record.id = f"contig{number_of_contigs:05}"

            sequence_length = len(record.seq)
            total_assembly_length += sequence_length
            gc_content = round(gc_fraction(record.seq) * 100, 2)

            if number_of_contigs == 1:
                longest_contig_length = sequence_length
            elif sequence_length > longest_contig_length:
                longest_contig_length = sequence_length

            record.description = f"len={sequence_length}"
            SeqIO.write(record, output_handle, "fasta")

            stats_dict[record.id] = {
                "length": sequence_length,
                "gc": gc_content,
            }

    return number_of_contigs, longest_contig_length, total_assembly_length


def write_summary(
    hybracter_summary,
    sample,
    complete,
    total_assembly_length,
    number_of_contigs,
    best_round,
    longest_contig_length,
    longest_contig_coverage,
    number_circular_plasmids,
):
    """Write the single-row hybracter_summary TSV.

    ``complete`` is a bool (True -> "Complete" column "True"). For incomplete
    assemblies ``number_circular_plasmids`` is the string ``"Unknown"``; for
    complete ones it is an int.
    """
    summary_dict = {
        "Sample": sample,
        "Complete": "True" if complete else "False",
        "Total_assembly_length": total_assembly_length,
        "Number_of_contigs": number_of_contigs,
        "Most_accurate_polishing_round": best_round,
        "Longest_contig_length": longest_contig_length,
        "Longest_contig_coverage": longest_contig_coverage,
        "Number_circular_plasmids": number_circular_plasmids,
    }
    summary_df = pd.DataFrame([summary_dict])
    summary_df.to_csv(hybracter_summary, index=False, sep="\t")


def write_single_row_tsv(path, row_dict):
    """Write a one-row dict to a TSV (e.g. the long-read pyrodigal mean-CDS summary)."""
    pd.DataFrame([row_dict]).to_csv(path, index=False, sep="\t")


def write_per_contig_stats(per_contig_summary, stats_dict):
    """Write the per-contig stats TSV with ``contig_name`` as the first column."""
    stats_df = pd.DataFrame.from_dict(stats_dict, orient="index")
    stats_df["contig_name"] = stats_df.index
    # reorder with 'contig_name' first
    stats_df = stats_df[
        ["contig_name"] + [col for col in stats_df.columns if col != "contig_name"]
    ]
    stats_df.to_csv(per_contig_summary, index=False, sep="\t")
