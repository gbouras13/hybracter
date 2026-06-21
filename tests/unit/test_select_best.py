"""
Unit/regression tests for the select_best chromosome-assembly logic.

Pins B1: with `--logic best` but no valid ALE scores (empty score dir),
selection must fall back to the final polishing round instead of crashing with
`NameError: closest_to_zero_key`.

Pins B4: plasmid circularity must be detected with the exact `circular=true`
test, so `circular=false` headers are NOT miscounted as circular.
"""

import csv

import select_best_chromosome_assembly_complete as sbc
import select_best_chromosome_assembly_incomplete as sb


def _write_fasta(path, seq_id, length):
    path.write_text(f">{seq_id}\n" + "A" * length + "\n")


def _write_flye_info(path):
    # Flye assembly_info.txt: header + one contig row (needs length + cov. cols)
    path.write_text(
        "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
        "contig_1\t100000\t50\tN\tN\t1\t*\t*\n"
    )


def test_logic_best_with_no_ale_scores_falls_back(tmp_path):
    """--logic best + empty ALE dir must not raise (B1) and keeps the final round."""
    ale_dir = tmp_path / "ale"
    ale_dir.mkdir()  # empty -> no *.score files -> filtered_score_dict is empty

    polca = tmp_path / "polca.fasta"
    _write_fasta(polca, "contig_1", 100000)
    flye_info = tmp_path / "assembly_info.txt"
    _write_flye_info(flye_info)

    summary = tmp_path / "hybracter_summary.tsv"
    per_contig = tmp_path / "per_contig.tsv"
    ale_summary = tmp_path / "ale_summary.tsv"
    out_fasta = tmp_path / "out.fasta"

    # Pre-fix this raised NameError: closest_to_zero_key
    sb.select_best_chromosome_assembly_incomplete(
        str(summary),
        str(per_contig),
        str(ale_dir),
        str(out_fasta),
        str(ale_summary),
        str(polca),  # pre_polish_fasta (unused on this path)
        str(polca),  # medaka_fasta (unused)
        str(polca),  # polypolish_fasta (unused)
        str(polca),  # polca_fasta (the fallback assembly)
        "Sample1",
        str(flye_info),
        "best",
        False,  # no_pypolca -> default best_round is "pypolca"
    )

    # the script renames contigs to contig00001 etc; just confirm output exists
    assert "contig00001" in out_fasta.read_text()
    row = next(csv.DictReader(open(summary), delimiter="\t"))
    assert row["Complete"] == "False"
    assert row["Most_accurate_polishing_round"] == "pypolca"


def test_plasmid_circularity_only_counts_circular_true(tmp_path):
    """B4: only `circular=true` plasmids count; `circular=false`/absent must not.

    plassembler headers carry `circular=true` (or, in principle, `circular=false`).
    The old substring test `"circular" in ...` also matched `circular=false` and
    inflated Number_circular_plasmids. Pin the exact `circular=true` behaviour
    across all three cases: true, false, and absent.
    """
    ale_dir = tmp_path / "ale"
    ale_dir.mkdir()  # empty -> logic="last" just takes the polca round

    # chromosome (the selected best assembly)
    polca = tmp_path / "polca.fasta"
    _write_fasta(polca, "contig_1", 100000)

    # three plasmids: circular=true, circular=false, and no circular tag.
    # the first whitespace-delimited token is stripped (plassembler contig id),
    # so each header needs an id token followed by the description.
    plassembler = tmp_path / "plassembler.fasta"
    plassembler.write_text(
        ">1 len=1000 circular=true copy_number=2\n" + "A" * 1000 + "\n"
        ">2 len=2000 circular=false copy_number=1\n" + "C" * 2000 + "\n"
        ">3 len=3000 copy_number=1\n" + "G" * 3000 + "\n"
    )

    ignore_list = tmp_path / "ignore.txt"
    ignore_list.write_text("")  # nothing ignored
    flye_info = tmp_path / "assembly_info.txt"
    _write_flye_info(flye_info)

    summary = tmp_path / "hybracter_summary.tsv"
    per_contig = tmp_path / "per_contig.tsv"
    ale_summary = tmp_path / "ale_summary.tsv"
    chrom_fasta = tmp_path / "chromosome.fasta"
    plasmid_fasta = tmp_path / "plasmid.fasta"
    overall_fasta = tmp_path / "overall.fasta"

    sbc.select_best_chromosome_assembly_complete(
        str(summary),
        str(per_contig),
        str(ale_dir),
        str(plassembler),
        str(chrom_fasta),
        str(plasmid_fasta),
        str(overall_fasta),
        str(ale_summary),
        str(polca),  # chrom_pre_polish_fasta (unused on this path)
        str(polca),  # medaka_rd_1_fasta (unused)
        str(polca),  # medaka_rd_2_fasta (unused)
        str(polca),  # polypolish_fasta (unused)
        str(polca),  # polca_fasta (the selected assembly)
        "Sample1",
        str(flye_info),
        "last",  # not "best" -> just take the polca round
        False,  # no_pypolca -> best_round "pypolca"
        str(ignore_list),
    )

    row = next(csv.DictReader(open(summary), delimiter="\t"))
    # only the circular=true plasmid counts (false + absent must not)
    assert row["Number_circular_plasmids"] == "1"

    per = {r["contig_name"]: r for r in csv.DictReader(open(per_contig), delimiter="\t")}
    plasmids = {k: v for k, v in per.items() if v["contig_type"] == "plasmid"}
    assert sorted(v["circular"] for v in plasmids.values()) == ["False", "False", "True"]
