"""
Unit/regression tests for extract_chromosome.get_chromosome_plasmids.

Pins the chromosome-extraction + ignore_list behaviour:
- default: extract every contig above min_chrom_length; non-circular ones are
  added to the ignore_list (so dnaapler skips reorientation for them).
- --circular_chromosome: only extract contigs that are long AND Flye-circular.
"""

import extract_chromosome


def _write_fasta(path, records):
    with open(path, "w") as f:
        for seq_id, length in records:
            f.write(f">{seq_id}\n" + "A" * length + "\n")


def _write_flye_info(path, rows):
    with open(path, "w") as f:
        f.write("#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n")
        for seq_name, length, circ in rows:
            f.write(f"{seq_name}\t{length}\t30\t{circ}\tN\t1\t*\t*\n")


def _run(tmp_path, seq_len, circ, min_chrom_length, circular_chromosome):
    fasta = tmp_path / "assembly.fasta"
    info = tmp_path / "assembly_info.txt"
    chrom_out = tmp_path / "chromosome.fasta"
    ignore = tmp_path / "ignore_list.txt"
    _write_fasta(fasta, [("contig_1", seq_len)])
    _write_flye_info(info, [("contig_1", seq_len, circ)])
    extract_chromosome.get_chromosome_plasmids(
        str(fasta), str(chrom_out), str(ignore), min_chrom_length, str(info), False, circular_chromosome
    )
    return chrom_out.read_text(), ignore.read_text().strip()


def test_noncircular_long_extracted_and_ignored_default(tmp_path):
    extracted, ignored = _run(tmp_path, 200000, "N", 100000, False)
    assert "contig_1" in extracted
    assert "contig_1" in ignored


def test_circular_long_extracted_not_ignored_default(tmp_path):
    extracted, ignored = _run(tmp_path, 200000, "Y", 100000, False)
    assert "contig_1" in extracted
    assert ignored == ""


def test_noncircular_long_not_extracted_with_circular_flag(tmp_path):
    extracted, ignored = _run(tmp_path, 200000, "N", 100000, True)
    assert "contig_1" not in extracted
    assert ignored == ""


def test_circular_long_extracted_with_circular_flag(tmp_path):
    extracted, ignored = _run(tmp_path, 200000, "Y", 100000, True)
    assert "contig_1" in extracted
    assert ignored == ""


def test_short_contig_never_extracted(tmp_path):
    for circ in ("Y", "N"):
        for flag in (True, False):
            extracted, _ = _run(tmp_path, 50000, circ, 100000, flag)
            assert "contig_1" not in extracted
