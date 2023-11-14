#!/usr/bin/env python3
"""
Ryan Wick 2023 with modifications by George Bouras 2023

This script produces a human-readable output showing the differences between two alternative
assemblies of a genome built for integration with Hybracter.

The two assemblies must have the same number of contigs and corresponding contigs must have the same strand and starting position.

The script will automatically sort your contigs in order of size (to handle cases where an isolate has multiple chromosomes - e.g. Vibrios).

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. This program is distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should
have received a copy of the GNU General Public License along with this program. If not, see
<https://www.gnu.org/licenses/>.
"""

import datetime
# import edlib
import gzip
import os
import re
import sys
import textwrap

import mappy
from Bio import SeqIO


def run_compare(
    assembly_1,
    assembly_2,
    padding,
    merge,
    aligner,
    outputfile,
    reference_polishing_round,
    query_polishing_round,
):
    assembly_1, assembly_2 = load_assemblies(
        assembly_1, assembly_2, reference_polishing_round, query_polishing_round
    )
    touch_file(outputfile)
    align_sequences(
        assembly_1,
        assembly_2,
        padding,
        merge,
        aligner,
        outputfile,
        reference_polishing_round,
        query_polishing_round,
    )


def touch_file(path):
    with open(path, "a"):
        os.utime(path, None)


def load_assemblies(
    assembly_1_filename,
    assembly_2_filename,
    reference_polishing_round,
    query_polishing_round,
):
    section_header("Loading assemblies")
    log(assembly_1_filename)
    assembly_1 = load_fasta(assembly_1_filename)
    for name, seq in assembly_1:
        log(f"{reference_polishing_round}  {name}: {len(seq):,} bp")
    log()

    log(assembly_2_filename)
    assembly_2 = load_fasta(assembly_2_filename)
    for name, seq in assembly_2:
        log(f"{query_polishing_round}  {name}: {len(seq):,} bp")
    log()
    if len(assembly_1) != len(assembly_2):
        quit_with_error(
            "Error: the assembly_1 and assembly_2 assemblies need to contain the same number "
            "of sequences"
        )
    for b, a in zip(assembly_1, assembly_2):
        assembly_1_name, assembly_1_seq = b
        assembly_1_seq_len = len(assembly_1_seq)
        assembly_2_name, assembly_2_seq = a
        assembly_2_seq_len = len(assembly_2_seq)
        if assembly_1_seq_len == 0 or assembly_2_seq_len == 0:
            quit_with_error("Error: zero-length sequences are not allowed")
        ratio = assembly_1_seq_len / assembly_2_seq_len
        if ratio < 0.9 or ratio > 1.11111111:
            quit_with_error(
                f"Error: {assembly_1_name} and {assembly_2_name} are too different in length "
                f"- are the files in the same order?"
            )
    return assembly_1, assembly_2


def align_sequences(
    assembly_1,
    assembly_2,
    padding,
    merge,
    aligner,
    outputfile,
    reference_polishing_round,
    query_polishing_round,
):
    section_header("Aligning sequences")
    longest_label = get_longest_label(
        assembly_1, assembly_2, reference_polishing_round, query_polishing_round
    )
    for b, a in zip(assembly_1, assembly_2):
        assembly_1_name, assembly_1_seq = b
        assembly_2_name, assembly_2_seq = a
        output_differences(
            assembly_1_name,
            assembly_1_seq,
            assembly_2_name,
            assembly_2_seq,
            padding,
            merge,
            longest_label,
            aligner,
            outputfile,
            reference_polishing_round,
            query_polishing_round,
        )


def get_longest_label(
    assembly_1, assembly_2, reference_polishing_round, query_polishing_round
):
    longest_name, longest_seq, longest_round = 0, 0, 0
    for name, seq in assembly_1 + assembly_2:
        longest_name = max(longest_name, len(name))
        longest_seq = max(longest_seq, len(str(len(seq))))
    longest_round = max(longest_round, len(reference_polishing_round))
    # add 5 in case
    return longest_round + 5 + longest_name + 2 * longest_seq + 3


def output_differences(
    assembly_1_name,
    assembly_1_seq,
    assembly_2_name,
    assembly_2_seq,
    padding,
    merge,
    longest_label,
    aligner,
    outputfile,
    reference_polishing_round,
    query_polishing_round,
):
    log(
        f"Aligning {reference_polishing_round}:{assembly_1_name} to {query_polishing_round}:{assembly_2_name}:"
    )
    (
        assembly_1_aligned,
        assembly_2_aligned,
        differences,
        assembly_1_pos,
        assembly_2_pos,
        diff_pos,
    ) = get_aligned_seqs(assembly_1_seq, assembly_2_seq, aligner)
    if len(diff_pos) == 1:
        log(f"  1 difference")
        with open(outputfile, "a") as out:
            out.write(f"  1 difference")
            out.write("\n")
            out.close()
    else:
        log(f"  {len(diff_pos):,} differences")
        with open(outputfile, "a") as out:
            out.write(f"  {len(diff_pos):,} differences")
            out.write("\n")
            out.close()
    log()

    aligned_len = len(assembly_1_aligned)
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)

    for start, end in diff_ranges:
        try:
            # Convert positions in alignment to positions in unaligned sequences:
            assembly_1_start, assembly_1_end = (
                assembly_1_pos[start],
                assembly_1_pos[end],
            )
            assembly_2_start, assembly_2_end = (
                assembly_2_pos[start],
                assembly_2_pos[end],
            )

            # Sanity check:
            assert (
                assembly_1_aligned[start:end].replace("-", "")
                == assembly_1_seq[assembly_1_start:assembly_1_end]
            )
            assert (
                assembly_2_aligned[start:end].replace("-", "")
                == assembly_2_seq[assembly_2_start:assembly_2_end]
            )

            # Add 1 to starts to convert from 0-based exclusive ranges to 1-based inclusive ranges.
            assembly_1_label = f"{reference_polishing_round}:{assembly_1_name} {assembly_1_start+1}-{assembly_1_end}:"
            assembly_2_label = f"{query_polishing_round}:{assembly_2_name} {assembly_2_start+1}-{assembly_2_end}:"
            assert len(assembly_1_label) <= longest_label
            assert len(assembly_2_label) <= longest_label
            assembly_1_label = assembly_1_label.rjust(longest_label)
            assembly_2_label = assembly_2_label.rjust(longest_label)

            print(f"{assembly_1_label}", assembly_1_aligned[start:end], flush=True)
            print(f"{assembly_2_label}", assembly_2_aligned[start:end], flush=True)
            print(" " * longest_label, differences[start:end], flush=True)
            print(flush=True)

            # Open the output file in write mode
            with open(outputfile, "a") as out:
                out.write(f"{assembly_1_label} {assembly_1_aligned[start:end]}\n")
                out.write(f"{assembly_2_label} {assembly_2_aligned[start:end]}\n")
                out.write(" " * longest_label + differences[start:end] + "\n")
                out.write("\n")
                out.close()
        except IndexError:
            print(f"Warning: Error writing {assembly_1_label} or {assembly_2_label}")
            print("This isn't a fatal error, but the difference won't be visualised.")
            continue

    log()


def make_diff_ranges(diff_pos, padding, merge, aligned_len):
    diff_ranges = []
    last_diff_pos = None
    for p in diff_pos:
        start = max(0, p - padding)
        end = min(aligned_len, p + padding + 1)
        if not last_diff_pos:  # this is the first diff
            diff_ranges.append((start, end))
        elif p - last_diff_pos <= merge:  # this diff is close to the previous diff
            prev_start = diff_ranges[-1][0]
            diff_ranges.pop()
            diff_ranges.append((prev_start, end))
        else:  # this diff is far from the previous diff
            diff_ranges.append((start, end))
        last_diff_pos = p
    return diff_ranges


def get_aligned_seqs(assembly_1_seq, assembly_2_seq, aligner):
    cigar = get_cigar(assembly_1_seq, assembly_2_seq, aligner)
    expanded_cigar = get_expanded_cigar(cigar)
    i, j = 0, 0
    assembly_1_aligned, assembly_2_aligned, differences = [], [], []
    assembly_1_positions, assembly_2_positions, diff_positions = [], [], []
    for c in expanded_cigar:
        assembly_1_positions.append(i)
        assembly_2_positions.append(j)
        if c == "M":
            b_1 = assembly_1_seq[i]
            b_2 = assembly_2_seq[j]
            if b_1 == b_2:
                diff = " "
            else:
                diff = "*"
                diff_positions.append(len(differences))
            i += 1
            j += 1
        elif c == "=":
            b_1 = assembly_1_seq[i]
            b_2 = assembly_2_seq[j]
            diff = " "
            i += 1
            j += 1
            assert b_1 == b_2
        elif c == "X":
            b_1 = assembly_1_seq[i]
            b_2 = assembly_2_seq[j]
            diff = "*"
            diff_positions.append(len(differences))
            i += 1
            j += 1
            assert b_1 != b_2
        elif c == "I":
            b_1 = assembly_1_seq[i]
            b_2 = "-"
            diff = "*"
            diff_positions.append(len(differences))
            i += 1
        elif c == "D":
            b_1 = "-"
            b_2 = assembly_2_seq[j]
            diff = "*"
            diff_positions.append(len(differences))
            j += 1
        else:
            assert False
        assembly_1_aligned.append(b_1)
        assembly_2_aligned.append(b_2)
        differences.append(diff)
    assert i == len(assembly_1_seq)
    assert j == len(assembly_2_seq)
    assembly_1_aligned = "".join(assembly_1_aligned)
    assembly_2_aligned = "".join(assembly_2_aligned)
    differences = "".join(differences)
    for p in diff_positions:
        assert differences[p] == "*"
    assert assembly_1_aligned.replace("-", "") == assembly_1_seq
    assert assembly_2_aligned.replace("-", "") == assembly_2_seq
    return (
        assembly_1_aligned,
        assembly_2_aligned,
        differences,
        assembly_1_positions,
        assembly_2_positions,
        diff_positions,
    )


def get_cigar(assembly_1_seq, assembly_2_seq, aligner):
    if aligner == "mappy":
        return get_cigar_with_mappy(assembly_1_seq, assembly_2_seq)
    elif aligner == "edlib":
        return get_cigar_with_edlib(assembly_1_seq, assembly_2_seq)
    else:
        assert False


def get_cigar_with_mappy(assembly_1_seq, assembly_2_seq):
    a = mappy.Aligner(seq=assembly_2_seq, preset="map-ont")
    for result in a.map(assembly_1_seq):
        full_length_query = result.q_st == 0 and result.q_en == len(assembly_1_seq)
        full_length_ref = result.r_st == 0 and result.r_en == len(assembly_2_seq)
        if full_length_query and full_length_ref:
            return result.cigar_str
    quit_with_error("Error: mappy alignment failed, try using --aligner edlib")


def get_cigar_with_edlib(assembly_1_seq, assembly_2_seq):
    result = edlib.align(assembly_1_seq, assembly_2_seq, mode="NW", task="path")
    return result["cigar"]


def get_expanded_cigar(cigar):
    expanded_cigar = []
    cigar_parts = re.findall(r"\d+[IDX=M]", cigar)
    for p in cigar_parts:
        size = int(p[:-1])
        letter = p[-1]
        expanded_cigar.append(letter * size)
    return "".join(expanded_cigar)


def finished_message(start_time):
    section_header("Finished!")
    time_to_run = datetime.datetime.now() - start_time
    log(f"Time to run: {time_to_run}")
    log()


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {
        "gz": (b"\x1f", b"\x8b", b"\x08"),
        "bz2": (b"\x42", b"\x5a", b"\x68"),
        "zip": (b"\x50", b"\x4b", b"\x03", b"\x04"),
    }
    max_len = max(len(x) for x in magic_dict)
    unknown_file = open(str(filename), "rb")
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = "plain"
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == "bz2":
        sys.exit("\nError: cannot use bzip2 format - use gzip instead")
    if compression_type == "zip":
        sys.exit("\nError: cannot use zip format - use gzip instead")
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == "gz":
        return gzip.open
    else:  # plain text
        return open


# replaces Ryan's function
# reads in FASTA
# should handle the case (vibrio) with 2 chroms, but where the polisher (medaka) swapped the order for some reason
def load_fasta(fasta_filename):
    fasta_seqs = []
    records = list(SeqIO.parse(fasta_filename, "fasta"))

    # Sort the records by sequence length
    records.sort(key=lambda x: len(x.seq))

    for record in records:
        fasta_seqs.append((record.id, str(record.seq)))

    return fasta_seqs


# def load_fasta(fasta_filename):
#     fasta_seqs = []
#     with get_open_func(fasta_filename)(fasta_filename, "rt") as fasta_file:
#         name = ""
#         sequence = []
#         for line in fasta_file:
#             line = line.strip()
#             if not line:
#                 continue
#             if line[0] == ">":  # Header line = start of new contig
#                 if name:
#                     fasta_seqs.append((name.split()[0], "".join(sequence)))
#                     sequence = []
#                 name = line[1:]
#             else:
#                 sequence.append(line.upper())
#         if name:
#             fasta_seqs.append((name.split()[0], "".join(sequence)))
#     return fasta_seqs


def log(message="", end="\n"):
    print(message, file=sys.stderr, flush=True, end=end)


def section_header(text):
    log()
    time = get_timestamp()
    time_str = dim("(" + time + ")")
    header = bold_yellow_underline(text)
    print(header + " " + time_str, file=sys.stderr, flush=True)


END_FORMATTING = "\033[0m"
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
YELLOW = "\033[93m"
DIM = "\033[2m"


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def quit_with_error(text):
    terminal_width, _ = get_terminal_size_stderr()
    log()
    for line in textwrap.wrap(text, width=terminal_width - 1):
        log(line)
    log()
    sys.exit()


def get_terminal_size_stderr(fallback=(80, 24)):
    """
    Unlike shutil.get_terminal_size, which looks at stdout, this looks at stderr.
    """
    try:
        size = os.get_terminal_size(sys.__stderr__.fileno())
    except (AttributeError, ValueError, OSError):
        size = os.terminal_size(fallback)
    return size


def get_timestamp():
    return "{:%Y-%m-%d %H:%M:%S}".format(datetime.datetime.now())


def check_python_version():
    if sys.version_info.major < 3 or sys.version_info.minor < 6:
        sys.exit("\nError: compare_assemblies.py requires Python 3.6 or later")


# Unit tests for Pytest
# =====================


def test_get_longest_label():
    assembly_1 = [("seq1", "ACGT"), ("seq2", "ACGTACGTACGTACGT")]
    assembly_2 = [("seq1_polished", "ACGT"), ("seq2_polished", "ACGT")]
    reference_polishing_round = "polypolish"
    query_polishing_round = "pypolca"
    assert (
        get_longest_label(
            assembly_1, assembly_2, reference_polishing_round, query_polishing_round
        )
        == 20
    )


def test_get_expanded_cigar():
    assert get_expanded_cigar("5=") == "====="
    assert get_expanded_cigar("3=2X4=1I6=3D3=") == "===XX====I======DDD==="


def test_make_diff_ranges():
    diff_pos = [100, 110]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 121)]

    diff_pos = [100, 120]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 131)]

    diff_pos = [100, 121]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 111), (111, 132)]

    diff_pos = [100, 150]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(90, 111), (140, 161)]

    diff_pos = [2, 195]
    padding = 10
    merge = 20
    aligned_len = 200
    diff_ranges = make_diff_ranges(diff_pos, padding, merge, aligned_len)
    assert diff_ranges == [(0, 13), (185, 200)]


#### to actually run the rule
run_compare(
    snakemake.input.reference,
    snakemake.input.assembly,
    15,
    30,
    "mappy",
    snakemake.output.diffs,
    snakemake.params.reference_polishing_round,
    snakemake.params.query_polishing_round,
)
