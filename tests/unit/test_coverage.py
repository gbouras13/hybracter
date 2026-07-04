"""
Unit/regression tests for the coverage-estimation scripts.

Guards the #142 class of bug (seqkit TSV parsing must not depend on a fragile
pandas/numpy combination) for both the long-read and short-read paths.
"""

import pytest

import lr_coverage
import sr_coverage


def _write_seqkit(path, sum_len, extra_cols=False):
    """Write a minimal `seqkit stats -T` (tab-separated) table with one data row.

    `extra_cols=True` appends N50/N90 columns (as `seqkit stats -N 50 -N 90`
    produces for long reads) to prove parsing is by column *name*, not position.
    """
    cols = ["file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len", "max_len"]
    vals = ["reads.fq.gz", "FASTQ", "DNA", "100", str(sum_len), "50", "150.0", "300"]
    if extra_cols:
        cols += ["N50", "N90"]
        vals += ["180", "70"]
    path.write_text("\t".join(cols) + "\n" + "\t".join(vals) + "\n")


def test_lr_get_sum_len(tmp_path):
    f = tmp_path / "long.tsv"
    _write_seqkit(f, 1234567)
    assert lr_coverage.get_sum_len(str(f)) == 1234567


def test_lr_get_sum_len_with_extra_n_columns(tmp_path):
    # long-read seqkit output carries N50/N90 columns; parse must stay by-name
    f = tmp_path / "long.tsv"
    _write_seqkit(f, 9_000_000, extra_cols=True)
    assert lr_coverage.get_sum_len(str(f)) == 9_000_000


def test_sr_get_sum_len(tmp_path):
    f = tmp_path / "r1.tsv"
    _write_seqkit(f, 555)
    assert sr_coverage.get_sum_len(str(f)) == 555


def test_lr_coverage_value_written(tmp_path):
    seqkit = tmp_path / "l.tsv"
    _write_seqkit(seqkit, 5_000_000)
    out = tmp_path / "cov.txt"
    # min_bases low enough to proceed; 5,000,000 / 1,000,000 = 5x
    lr_coverage.get_quick_coverage_estimate_long(str(seqkit), "1000000", "100", str(out))
    assert out.read_text().strip() == "5"


def test_lr_coverage_below_min_bases_exits(tmp_path):
    seqkit = tmp_path / "l.tsv"
    _write_seqkit(seqkit, 100)
    out = tmp_path / "cov.txt"
    # long_bases (100) < min_bases (100000) -> hybracter should stop
    with pytest.raises(SystemExit):
        lr_coverage.get_quick_coverage_estimate_long(str(seqkit), "1000", "100000", str(out))


def test_sr_coverage_value_written(tmp_path):
    r1 = tmp_path / "r1.tsv"
    r2 = tmp_path / "r2.tsv"
    _write_seqkit(r1, 3_000_000)
    _write_seqkit(r2, 2_000_000)
    out = tmp_path / "cov.txt"
    # (3,000,000 + 2,000,000) / 1,000,000 = 5x
    sr_coverage.get_quick_coverage_estimate(str(r1), str(r2), "1000000", str(out))
    assert out.read_text().strip() == "5"
