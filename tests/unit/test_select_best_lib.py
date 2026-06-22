"""
Unit tests for the shared select_best_lib helpers (E1).

These pin the IO/parsing/writing behaviour that the eight
select_best_chromosome_assembly_* scripts delegate to, so the de-duplication
refactor is verified to be behaviour-preserving.
"""

import csv

import select_best_lib as sbl


def _write_fasta(path, records):
    """records: list of (header_without_gt, sequence)."""
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records))


def _read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# --------------------------------------------------------------------------- #
# tiny helpers
# --------------------------------------------------------------------------- #
def test_is_file_empty(tmp_path):
    empty = tmp_path / "empty"
    empty.write_text("")
    full = tmp_path / "full"
    full.write_text("x")
    assert sbl.is_file_empty(str(empty)) is True
    assert sbl.is_file_empty(str(full)) is False


def test_touch_file_creates(tmp_path):
    p = tmp_path / "touched"
    assert not p.exists()
    sbl.touch_file(str(p))
    assert p.exists() and p.read_text() == ""


# --------------------------------------------------------------------------- #
# ALE scores
# --------------------------------------------------------------------------- #
def test_read_ale_scores_picks_closest_to_zero(tmp_path):
    ale = tmp_path / "ale"
    ale.mkdir()
    (ale / "chrom_pre_polish.score").write_text("-50.0\n")
    (ale / "polypolish.score").write_text("-12.0\n")  # closest to zero
    (ale / "polca.score").write_text("-30.0\n")
    (ale / "broken.score").write_text("not_a_float\n")  # filtered out

    summary = tmp_path / "ale_summary.tsv"
    filtered, closest = sbl.read_ale_scores(str(ale), str(summary))

    assert closest == "polypolish"
    assert "broken" not in filtered
    assert set(filtered) == {"chrom_pre_polish", "polypolish", "polca"}
    # summary written and sorted ascending (worst first)
    rows = _read_tsv(summary)
    assert rows[0]["Key"] == "chrom_pre_polish"  # -50 is most negative


def test_read_ale_scores_empty_returns_none(tmp_path):
    ale = tmp_path / "ale"
    ale.mkdir()  # no .score files
    summary = tmp_path / "ale_summary.tsv"
    filtered, closest = sbl.read_ale_scores(str(ale), str(summary))
    assert filtered == {} and closest is None  # B1: callers must handle None


# --------------------------------------------------------------------------- #
# flye coverage
# --------------------------------------------------------------------------- #
def test_longest_contig_coverage(tmp_path):
    flye = tmp_path / "assembly_info.txt"
    flye.write_text(
        "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
        "contig_1\t100000\t37\tN\tN\t1\t*\t*\n"
        "contig_2\t5000\t99\tN\tN\t1\t*\t*\n"  # higher cov but shorter
    )
    # returns the longest contig's cov. as the raw string from the file (it is only
    # ever written into a summary TSV), picking contig_1 (37) over the shorter contig_2
    assert sbl.longest_contig_coverage(str(flye)) == "37"


# --------------------------------------------------------------------------- #
# chromosome writer (complete)
# --------------------------------------------------------------------------- #
def test_write_chromosomes(tmp_path):
    best = tmp_path / "best.fasta"
    _write_fasta(best, [("contigA", "A" * 1000), ("contigB", "C" * 500)])
    ignore = tmp_path / "ignore.txt"
    ignore.write_text("contigB\n")  # contigB is non-circular

    chrom_fa = tmp_path / "chrom.fasta"
    overall_fa = tmp_path / "overall.fasta"
    stats = {}
    chromosomes, longest, total = sbl.write_chromosomes(
        str(best), str(chrom_fa), str(overall_fa), str(ignore), stats
    )

    assert chromosomes == 3  # post-loop counter (1 + 2 chromosomes)
    assert longest == 1000 and total == 1500
    # renamed + circularity
    text = chrom_fa.read_text()
    assert ">chromosome00001 len=1000 circular=true" in text
    assert ">chromosome00002 len=500" in text and "chromosome00002 len=500 circular" not in text
    assert stats["chromosome00001"]["circular"] == "True"
    assert stats["chromosome00002"]["circular"] == "False"
    assert stats["chromosome00001"]["contig_type"] == "chromosome"
    # overall fasta got both records too
    assert overall_fa.read_text().count(">chromosome") == 2


# --------------------------------------------------------------------------- #
# plasmid writer (complete)
# --------------------------------------------------------------------------- #
def test_write_plasmids_counts_only_circular_true(tmp_path):
    plas = tmp_path / "plassembler.fasta"
    _write_fasta(
        plas,
        [
            ("1 len=1000 circular=true", "A" * 1000),
            ("2 len=2000 circular=false", "C" * 2000),
            ("3 len=3000", "G" * 3000),
        ],
    )
    out_plas = tmp_path / "plasmid.fasta"
    overall = tmp_path / "overall.fasta"
    overall.write_text(">chromosome00001 len=5 circular=true\nACGTA\n")  # pre-existing
    stats = {}
    plasmids, circular, total = sbl.write_plasmids(
        str(plas), str(out_plas), str(overall), stats, 5
    )
    assert plasmids == 3
    assert circular == 1  # B4: only circular=true counts
    assert total == 5 + 1000 + 2000 + 3000
    assert stats["plasmid00001"]["circular"] == "True"
    assert stats["plasmid00002"]["circular"] == "False"
    # overall fasta was appended to, not truncated
    assert ">chromosome00001" in overall.read_text()
    assert overall.read_text().count(">plasmid") == 3


def test_write_plasmids_empty_touches_file(tmp_path):
    plas = tmp_path / "plassembler.fasta"
    plas.write_text("")  # empty
    out_plas = tmp_path / "plasmid.fasta"
    overall = tmp_path / "overall.fasta"
    overall.write_text(">chromosome00001\nACGT\n")
    stats = {}
    plasmids, circular, total = sbl.write_plasmids(
        str(plas), str(out_plas), str(overall), stats, 4
    )
    assert plasmids == 0 and circular == 0 and total == 4
    assert out_plas.exists() and out_plas.read_text() == ""


# --------------------------------------------------------------------------- #
# contig writer (incomplete)
# --------------------------------------------------------------------------- #
def test_write_contigs(tmp_path):
    best = tmp_path / "best.fasta"
    _write_fasta(best, [("x", "A" * 300), ("y", "C" * 900)])
    out = tmp_path / "out.fasta"
    stats = {}
    n, longest, total = sbl.write_contigs(str(best), str(out), stats)
    assert n == 2 and longest == 900 and total == 1200
    text = out.read_text()
    assert ">contig00001 len=300" in text and ">contig00002 len=900" in text
    assert "circular" not in text  # incomplete: no circularity
    assert set(stats["contig00001"]) == {"length", "gc"}


# --------------------------------------------------------------------------- #
# summary / per-contig writers
# --------------------------------------------------------------------------- #
def test_write_summary_complete_vs_incomplete(tmp_path):
    comp = tmp_path / "comp.tsv"
    sbl.write_summary(str(comp), "S1", True, 1500, 2, "pypolca", 1000, 37, 1)
    row = _read_tsv(comp)[0]
    assert row["Complete"] == "True" and row["Number_circular_plasmids"] == "1"

    inc = tmp_path / "inc.tsv"
    sbl.write_summary(str(inc), "S1", False, 1500, 2, "medaka", 1000, 37, "Unknown")
    row = _read_tsv(inc)[0]
    assert row["Complete"] == "False" and row["Number_circular_plasmids"] == "Unknown"


def test_write_single_row_tsv(tmp_path):
    out = tmp_path / "pyrodigal.tsv"
    sbl.write_single_row_tsv(
        str(out), {"Sample": "S1", "pre_polish_mean_cds_length": 950.5}
    )
    row = _read_tsv(out)[0]
    assert row["Sample"] == "S1" and row["pre_polish_mean_cds_length"] == "950.5"


def test_write_per_contig_stats_contig_name_first(tmp_path):
    out = tmp_path / "per_contig.tsv"
    stats = {
        "chromosome00001": {"contig_type": "chromosome", "length": 1000, "gc": 50.0, "circular": "True"},
    }
    sbl.write_per_contig_stats(str(out), stats)
    header = out.read_text().splitlines()[0].split("\t")
    assert header[0] == "contig_name"
    row = _read_tsv(out)[0]
    assert row["contig_name"] == "chromosome00001" and row["contig_type"] == "chromosome"
