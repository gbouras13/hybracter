"""
Unit tests for the pure helper functions in compare_assemblies.py.

These tests previously lived inside compare_assemblies.py itself, where pytest
could never collect them: importing the module ran the module-level
`run_compare(snakemake...)` and died (Q1/Q3). Now that the script is import-safe
and `mappy` is imported lazily, they live here and actually run.
"""

import compare_assemblies as ca


def test_get_longest_label():
    assembly_1 = [("seq1", "ACGT"), ("seq2", "ACGTACGTACGTACGT")]
    assembly_2 = [("seq1_polished", "ACGT"), ("seq2_polished", "ACGT")]
    reference_polishing_round = "polypolish"
    query_polishing_round = "pypolca"
    # Budget = longest_round(10) + 5 + longest_name(13) + 2*longest_seq(2) + 3 = 35,
    # which must be >= the worst-case label "polypolish:seq2_polished 1-16:" (30 chars).
    # (The original in-script assertion said 20; it predated the current formula and
    # was never collected by pytest, so it had silently gone stale - see Q3.)
    assert (
        ca.get_longest_label(
            assembly_1, assembly_2, reference_polishing_round, query_polishing_round
        )
        == 35
    )


def test_get_expanded_cigar():
    assert ca.get_expanded_cigar("5=") == "====="
    assert ca.get_expanded_cigar("3=2X4=1I6=3D3=") == "===XX====I======DDD==="


def test_make_diff_ranges():
    padding = 10
    merge = 20
    aligned_len = 200

    assert ca.make_diff_ranges([100, 110], padding, merge, aligned_len) == [(90, 121)]
    assert ca.make_diff_ranges([100, 120], padding, merge, aligned_len) == [(90, 131)]
    assert ca.make_diff_ranges([100, 121], padding, merge, aligned_len) == [
        (90, 111),
        (111, 132),
    ]
    assert ca.make_diff_ranges([100, 150], padding, merge, aligned_len) == [
        (90, 111),
        (140, 161),
    ]
    assert ca.make_diff_ranges([2, 195], padding, merge, aligned_len) == [
        (0, 13),
        (185, 200),
    ]
