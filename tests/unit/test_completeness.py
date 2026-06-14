"""
Unit/regression tests for check_completeness.get_completeness.

Pins the completeness (C/I) classification logic, including the
--circular_chromosome behaviour (length-only by default vs length-AND-circular).
"""

import check_completeness


def _write_fasta(path, records):
    """records: list of (seq_id, length)."""
    with open(path, "w") as f:
        for seq_id, length in records:
            f.write(f">{seq_id}\n" + "A" * length + "\n")


def _write_flye_info(path, rows):
    """rows: list of (seq_name, length, circ) where circ is 'Y' or 'N'."""
    with open(path, "w") as f:
        f.write("#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n")
        for seq_name, length, circ in rows:
            f.write(f"{seq_name}\t{length}\t30\t{circ}\tN\t1\t*\t*\n")


def _run(tmp_path, seq_len, circ, min_chrom_length, circular_chromosome):
    fasta = tmp_path / "assembly.fasta"
    info = tmp_path / "assembly_info.txt"
    out = tmp_path / "completeness.txt"
    _write_fasta(fasta, [("contig_1", seq_len)])
    _write_flye_info(info, [("contig_1", seq_len, circ)])
    check_completeness.get_completeness(
        str(fasta), str(out), min_chrom_length, str(info), circular_chromosome
    )
    return out.read_text().strip()


def test_long_noncircular_complete_by_default(tmp_path):
    # default: length alone -> complete
    assert _run(tmp_path, 200000, "N", 100000, False) == "C"


def test_long_circular_complete_by_default(tmp_path):
    assert _run(tmp_path, 200000, "Y", 100000, False) == "C"


def test_long_noncircular_incomplete_with_circular_flag(tmp_path):
    # --circular_chromosome: requires Flye circularity too
    assert _run(tmp_path, 200000, "N", 100000, True) == "I"


def test_long_circular_complete_with_circular_flag(tmp_path):
    assert _run(tmp_path, 200000, "Y", 100000, True) == "C"


def test_short_always_incomplete(tmp_path):
    for circ in ("Y", "N"):
        for flag in (True, False):
            assert _run(tmp_path, 50000, circ, 100000, flag) == "I"
