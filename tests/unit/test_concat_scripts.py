"""
Round-trip unit tests for the polars-backed concat/summary scripts (§6):
add_sample_plassembler, combine_plassembler_info, create_final_hybracter_summary.

These read per-sample TSVs and concatenate/write them. The pandas->polars
migration must preserve values exactly - create_final produces the FINAL_OUTPUT
hybracter_summary the golden tests read, so ints must stay ints (e.g. "2", not
"2.0") and the "Unknown" Number_circular_plasmids of incomplete assemblies must
survive the concat.
"""

import csv

import pytest

pytest.importorskip("polars")

import add_sample_plassembler as asp
import combine_plassembler_info as cmb
import create_final_hybracter_summary as cfs

SUMMARY_HEADER = (
    "Sample\tComplete\tTotal_assembly_length\tNumber_of_contigs\t"
    "Most_accurate_polishing_round\tLongest_contig_length\t"
    "Longest_contig_coverage\tNumber_circular_plasmids\n"
)


def test_add_sample_prepends_sample_column(tmp_path):
    inp = tmp_path / "plas.tsv"
    inp.write_text("contig\tlength\tcircular\n1\t2473\tTrue\n2\t1000\tFalse\n")
    out = tmp_path / "out.tsv"
    asp.add_sample_plassembler(str(inp), str(out), "Sample1")
    assert out.read_text() == (
        "Sample\tcontig\tlength\tcircular\n"
        "Sample1\t1\t2473\tTrue\n"
        "Sample1\t2\t1000\tFalse\n"
    )


def test_add_sample_empty_input_touches_output(tmp_path):
    inp = tmp_path / "empty.tsv"
    inp.write_text("")
    out = tmp_path / "out.tsv"
    asp.add_sample_plassembler(str(inp), str(out), "Sample1")
    assert out.exists() and out.read_text() == ""


def test_combine_multiple_samples(tmp_path):
    sdir = tmp_path / "s"
    sdir.mkdir()
    (sdir / "s1.tsv").write_text("Sample\tcontig\tlength\nS1\t1\t2473\n")
    (sdir / "s2.tsv").write_text("Sample\tcontig\tlength\nS2\t1\t1500\nS2\t2\t900\n")
    out = tmp_path / "out.tsv"
    cmb.combine_sample_plassembler(str(sdir), str(out))
    lines = out.read_text().splitlines()
    assert lines[0] == "Sample\tcontig\tlength"
    assert len(lines) == 4  # header + 3 data rows (glob order is arbitrary)
    assert "S1\t1\t2473" in lines
    assert "S2\t2\t900" in lines


def test_combine_single_sample_passthrough(tmp_path):
    sdir = tmp_path / "s"
    sdir.mkdir()
    (sdir / "s1.tsv").write_text("Sample\tcontig\tlength\nS1\t1\t2473\n")
    out = tmp_path / "out.tsv"
    cmb.combine_sample_plassembler(str(sdir), str(out))
    assert out.read_text() == "Sample\tcontig\tlength\nS1\t1\t2473\n"


def test_combine_no_plasmids_touches_output(tmp_path):
    sdir = tmp_path / "s"
    sdir.mkdir()  # no *.tsv
    out = tmp_path / "out.tsv"
    cmb.combine_sample_plassembler(str(sdir), str(out))
    assert out.exists() and out.read_text() == ""


def test_create_final_preserves_ints_and_unknown(tmp_path):
    cdir = tmp_path / "c"
    cdir.mkdir()
    idir = tmp_path / "i"
    idir.mkdir()
    (cdir / "Sample1_summary.tsv").write_text(
        SUMMARY_HEADER + "Sample1\tTrue\t73079\t2\tpypolca\t70606\t61\t1\n"
    )
    (idir / "Sample2_summary.tsv").write_text(
        SUMMARY_HEADER + "Sample2\tFalse\t70606\t1\tpypolca\t70606\t61\tUnknown\n"
    )
    out = tmp_path / "final.tsv"
    cfs.make_final_summary(str(out), str(cdir), str(idir))

    rows = {r["Sample"]: r for r in csv.DictReader(open(out), delimiter="\t")}
    # exact-match golden fields must not be reformatted
    assert rows["Sample1"]["Number_of_contigs"] == "2"  # not "2.0"
    assert rows["Sample1"]["Number_circular_plasmids"] == "1"
    assert rows["Sample1"]["Complete"] == "True"
    # the incomplete sentinel must survive the concat
    assert rows["Sample2"]["Number_circular_plasmids"] == "Unknown"
    assert rows["Sample2"]["Complete"] == "False"
