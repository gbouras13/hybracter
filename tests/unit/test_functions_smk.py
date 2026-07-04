"""
Unit tests for the helpers in workflow/rules/preflight/functions.smk.

functions.smk is plain Python that Snakemake `include:`s; its module-level code is
only imports + defs (the Snakemake globals ``config`` / ``dictReads`` are referenced
*inside* functions), so it can be loaded and the pure helpers tested directly.

Pins B3 (is_fasta_file: plain FASTA, gzipped FASTA, non-FASTA, missing file, gzipped
non-FASTA) and the E4 lrge genome-size helpers (auto mode).
"""

import gzip
import importlib.machinery
import importlib.util
from pathlib import Path

SMK = (
    Path(__file__).resolve().parents[2]
    / "hybracter"
    / "workflow"
    / "rules"
    / "preflight"
    / "functions.smk"
)


def _load_functions_smk():
    loader = importlib.machinery.SourceFileLoader("functions_smk", str(SMK))
    spec = importlib.util.spec_from_loader("functions_smk", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod


fns = _load_functions_smk()


# ---------- B3: is_fasta_file ----------


def test_is_fasta_plain(tmp_path):
    p = tmp_path / "a.fasta"
    p.write_text(">c1\nACGT\n")
    assert fns.is_fasta_file(str(p)) is True


def test_is_fasta_gzipped(tmp_path):
    p = tmp_path / "a.fasta.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">c1\nACGT\n")
    assert fns.is_fasta_file(str(p)) is True


def test_non_fasta_is_false(tmp_path):
    p = tmp_path / "notfasta.txt"
    p.write_text("this is not a fasta\njust some text\n")
    assert fns.is_fasta_file(str(p)) is False


def test_missing_file_is_false(tmp_path):
    assert fns.is_fasta_file(str(tmp_path / "does_not_exist.fasta")) is False


def test_gzipped_non_fasta_is_false(tmp_path):
    p = tmp_path / "notfasta.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("not a fasta at all\n")
    assert fns.is_fasta_file(str(p)) is False


# ---------- E4: lrge genome-size helpers (auto mode) ----------


def test_get_min_chrom_length_auto(tmp_path):
    lrge = tmp_path / "size.txt"
    lrge.write_text("1000000\n")
    assert fns.getMinChromLength(str(lrge), "s", True) == int(1000000 * 0.8)


def test_get_min_bases_auto(tmp_path):
    lrge = tmp_path / "size2.txt"
    lrge.write_text("1000000\n")
    assert fns.getMinBases(str(lrge), "s", True, 10) == int(10 * 1000000 * 0.8)


def test_get_target_bases_auto(tmp_path):
    lrge = tmp_path / "size3.txt"
    lrge.write_text("1000000\n")
    assert fns.getTargetBases(str(lrge), "s", True, 5) == int(5 * 1000000 * 0.8)
