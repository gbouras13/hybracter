"""
# to run on mac
pytest --run_mac .
"""

import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pytest

__author__ = "George Bouras"
__copyright__ = "Copyright 2023, Hybracter"
__license__ = "MIT"
__type__ = "Test Script"
__maintainer__ = "George Bouras"
__email__ = "george.bouras@adelaide.edu.au"

# Every test in this module shells out to a full hybracter run (needs conda +
# the assembler tool stack), so the whole file is the slow "integration" lane.
# The fast unit tests live in tests/unit/. Run these with: pytest -m integration
pytestmark = pytest.mark.integration

# test data
test_data_path = Path("hybracter/test_data")
db_dir: Path = "plassembler_db"


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")

# to run on macs
@pytest.fixture(scope="session")
def run_mac(pytestconfig):
    return pytestconfig.getoption("run_mac")

threads = 2


def remove_directory(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_hybracter_t_hybrid(run_mac):
    """hybracter test-hybrid"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-hybrid --threads {threads} --output {outdir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_hybrid_last(run_mac):
    """hybracter test-hybrid"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --logic last"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_t_hybrid_no_medaka(run_mac):
    """hybracter test-hybrid no medaka"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --no_medaka"
    if run_mac:
        cmd = f"hybracter test-hybrid -h"
        #cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_t_subsample_depth(run_mac):
    """hybracter test-hybrid no medaka"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --no_medaka --subsample_depth 50"
    if run_mac:
        cmd = f"hybracter test-hybrid -h"
        #cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_t_depth_filter(run_mac):
    """hybracter test-hybrid depth filter 2"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --no_medaka --depth_filter 2"
    if run_mac:
        cmd = f"hybracter test-hybrid -h"
        #cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_hybrid_last_no_medaka(run_mac):
    """hybracter test-hybrid - don't worry about MAC as no medaka"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --logic last --no_medaka --skip_qc"
    if run_mac:
        cmd = f"hybracter test-hybrid -h"
        #cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_t_hybrid_no_pypolca(run_mac):
    """hybracter test-hybrid no pypolca"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --no_pypolca --skip_qc"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_t_hybrid_no_pypolca_last_logic(run_mac):
    """hybracter test-hybrid no pypolca"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --no_pypolca --logic last --skip_qc"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_t_long(run_mac):
    """hybracter test-long"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-long --threads {threads} --output {outdir} --skip_qc"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_t_long_last(run_mac):
    """hybracter test-long"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-long --threads {threads} --output {outdir} --logic last --skip_qc"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_t_long_depth_filter(run_mac):
    """hybracter test-long with depth_filter - also methylation model for linux"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-long --threads {threads} --output {outdir} --skip_qc  --depth_filter 2"
    if run_mac:
        cmd += " --mac"
    else:
        cmd += " --medakaModel r1041_e82_400bps_bacterial_methylation"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_t_long_no_medaka(run_mac):
    """hybracter test-long no medaka - don't worry about mac"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-long --threads {threads} --output {outdir} --no_medaka"
    if run_mac:
        cmd = f"hybracter test-hybrid -h"
        #cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_install(run_mac):
    """test hybracter install and install dirs for later tests"""
    cmd = f"hybracter install  --databases {db_dir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)


"""
test with csv
"""


def test_hybracter_hybrid_csv(run_mac):
    """test hybracter hybrid default"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input.csv"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_hybrid_csv_auto(run_mac):
    """test hybracter hybrid auto

    Note that in 'auto' mode, it is going to be very hard for incomplete samples to happen
    This is because the required No. unique counter k-mers to be count is set to 10 
    Likely to get longer assemblies than 10kmers (in general), done some tests on subsampled data
    So warn people to be careful using kmc

    Sample 1
    kmc
    No. of unique counted k-mers       :        72697
    which should give a chromosome size of 58157  ( assembly is 70ish)

    Sample 2
    No. of unique counted k-mers       : 4542  
    which should give a chromosome size of 3633( assembly is 70ish)

    """
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input_auto.csv"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir} --auto"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_hybrid_csv_skip_qc_not_gzipped(run_mac):
    """test hybracter hybrid with --skip_qc and a long read fastq that isn't gzipped """
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input_not_gzipped.csv"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_hybrid_csv_ultra_low_cov(run_mac):
    """test hybracter hybrid default ultra low illumina coverage (4x) - to make sure pypolca is skipped"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input_ultra_low_illumina_coverage.csv"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_hybrid_csv_low_cov(run_mac):
    """test hybracter hybrid default low illumina coverage"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input_low_illumina_coverage.csv"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_no_pypolca_no_medaka(run_mac):
    """test hybracter hybrid no pypolca and no medaka"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input.csv"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --no_pypolca --databases {db_dir} --no_medaka "
    if run_mac:
        cmd = f"hybracter test-hybrid -h"
        #cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_no_medaka(run_mac):
    """test hybracter hybrid no medaka"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input.csv"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir} --no_medaka"
    if run_mac:
        cmd = f"hybracter test-hybrid -h"
        #cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_no_pypolca(run_mac):
    """test hybracter hybrid no pypolca"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input.csv"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --no_pypolca --databases {db_dir} "
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_subsample_depth(run_mac):
    """test hybracter hybrid subsample_Depth"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input.csv"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --no_pypolca --databases {db_dir} --subsample_depth 50"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_hybrid_datadir_1dir(run_mac):
    """test hybracter hybrid with --datadir - 1 dir"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input_datadir.csv"
    datadir: Path = "hybracter/test_data/Fastqs/"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir} --datadir {datadir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_hybrid_datadir_2dirs(run_mac):
    """test hybracter hybrid with --datadir - 2 dirs"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_hybrid_input_datadir.csv"
    # 2 dirs
    datadir: Path = "hybracter/test_data/Fastqs/,hybracter/test_data/Fastqs/"
    cmd = f"hybracter hybrid --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir} --datadir {datadir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

"""
long
"""

def test_hybracter_long(run_mac):
    """test hybracter long"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_long_input.csv"
    cmd = f"hybracter long --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

# test override
# bit useless for mac (as --bacteria won't work) but won't break anything so just leave

def test_hybracter_long_medaka_override(run_mac):
    """test hybracter long"""
    outdir: Path = "test_hybracter_output_medaka_override"
    input_csv: Path = test_data_path / "test_long_input.csv"
    cmd = f"hybracter long --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir} --medaka_override"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

# test --extra_params_flye

def test_hybracter_long_medaka_override(run_mac):
    """test hybracter long"""
    outdir: Path = "test_hybracter_output_medaka_override"
    input_csv: Path = test_data_path / "test_long_input.csv"
    cmd = f"hybracter long --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir} --extra_params_flye \"--meta\" "
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_long_auto(run_mac):
    """

    Note that in 'auto' mode, it is going to be very hard for incomplete samples to happen
    This is because the required No. unique counter k-mers to be count is set to 10 
    Likely to get longer assemblies than 10kmers (in general), done some tests on subsampled data
    If the chrom size is really low --> suggests that you need more sequencing data :)
    So warn people to be careful using kmc

    Sample 1
    kmc
    No. of unique counted k-mers       :        72697
    which should give a chromosome size of 58157  ( assembly is 70ish)

    Sample 2
    No. of unique counted k-mers       : 4542  
    which should give a chromosome size of 3633( assembly is 70ish)

    """
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_long_input_auto.csv"
    cmd = f"hybracter long --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir} --auto"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)

def test_hybracter_long_datadir(run_mac):
    """test hybracter long with --datadir"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_long_input_datadir.csv"
    datadir: Path = "hybracter/test_data/Fastqs/"
    cmd = f"hybracter long --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir} --datadir {datadir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)
    

def test_hybracter_long_no_medaka(run_mac):
    """test hybracter long"""
    outdir: Path = "test_hybracter_output"
    input_csv: Path = test_data_path / "test_long_input.csv"
    cmd = f"hybracter long --input {input_csv} --threads {threads} --output {outdir} --databases {db_dir} --no_medaka"
    if run_mac:
        cmd = f"hybracter test-hybrid -h"
        #cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


"""
single
"""


def test_hybracter_hybrid_single(run_mac):
    """test hybracter hybrid single"""
    outdir: Path = "test_hybracter_output_hybr_single"
    longreads: Path = test_data_path / "Fastqs/test_long_reads.fastq.gz"
    reads1: Path = test_data_path / "Fastqs/test_short_reads_R1.fastq.gz"
    reads2: Path = test_data_path / "Fastqs/test_short_reads_R2.fastq.gz"
    cmd = f"hybracter hybrid-single -l {longreads} -1 {reads1} -2 {reads2} -c 50000 -s Sample1 --threads {threads} --output {outdir} --databases {db_dir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_hybrid_single_no_medaka(run_mac):
    """test hybracter hybrid single no medaka"""
    outdir: Path = "test_hybracter_output_hybr_single"
    longreads: Path = test_data_path / "Fastqs/test_long_reads.fastq.gz"
    reads1: Path = test_data_path / "Fastqs/test_short_reads_R1.fastq.gz"
    reads2: Path = test_data_path / "Fastqs/test_short_reads_R2.fastq.gz"
    cmd = f"hybracter hybrid-single -l {longreads} -1 {reads1} -2 {reads2} -c 50000 -s Sample1 --threads {threads} --output {outdir} --databases {db_dir} --no_medaka"
    if run_mac:
        cmd = f"hybracter test-hybrid -h"
        #cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_long_single(run_mac):
    """test hybracter long"""
    outdir: Path = "test_hybracter_output_long_single"
    longreads: Path = test_data_path / "Fastqs/test_long_reads.fastq.gz"
    cmd = f"hybracter long-single -l {longreads} -c 50000 -s Sample1  --threads {threads} --output {outdir} --databases {db_dir}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_long_single_no_medaka_subsample_depth(run_mac):
    """test hybracter long no medaka and with --subsample_depth"""
    outdir: Path = "test_hybracter_output_long_single"
    longreads: Path = test_data_path / "Fastqs/test_long_reads.fastq.gz"
    cmd = f"hybracter long-single -l {longreads} -c 50000 -s Sample1  --threads {threads} --output {outdir} --databases {db_dir} --no_medaka --subsample_depth 50"
    if run_mac:
        cmd = f"hybracter test-hybrid -h"
        #cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


"""
custom db
"""


def test_hybracter_dnaapler_custom_db(run_mac):
    """test hybracter hybrid-single with custom_db"""
    outdir: Path = "test_hybracter_output_custom_db"
    longreads: Path = test_data_path / "Fastqs/test_long_reads.fastq.gz"
    reads1: Path = test_data_path / "Fastqs/test_short_reads_R1.fastq.gz"
    reads2: Path = test_data_path / "Fastqs/test_short_reads_R2.fastq.gz"
    dnaapler_custom_db: Path = test_data_path / "dnaapler_custom_db/all.faa"
    cmd = f"hybracter hybrid-single -l {longreads} -1 {reads1} -2 {reads2} -c 50000 -s Sample1 --threads {threads} --output {outdir} --databases {db_dir} --dnaapler_custom_db {dnaapler_custom_db}"
    if run_mac:
        cmd += " --mac"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_t_hybrid_circular_chromosome(run_mac):
    """hybracter test-hybrid with --circular_chromosome"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --no_medaka --skip_qc --circular_chromosome"
    if run_mac:
        cmd = "hybracter test-hybrid -h"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_t_long_circular_chromosome(run_mac):
    """hybracter test-long with --circular_chromosome"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-long --threads {threads} --output {outdir} --no_medaka --skip_qc --circular_chromosome"
    if run_mac:
        cmd = "hybracter test-long -h"
    exec_command(cmd)
    remove_directory(outdir)


# ---------------------------------------------------------------------------
# Unit tests for check_completeness and extract_chromosome logic
# These test the Python functions directly without going through Snakemake.
#
# Snakemake normally injects a `snakemake` builtin into scripts at runtime.
# We replicate that here so the module-level calls in each script don't fail
# on import. The mock points at real (but empty/throwaway) temp files so
# pandas / BioPython don't raise FileNotFoundError.
# ---------------------------------------------------------------------------

def _make_snakemake_mock(tmpdir, circ="Y"):
    """Return a SimpleNamespace that looks enough like snakemake for module import."""
    import types
    fasta = os.path.join(tmpdir, "_mock.fa")
    info = os.path.join(tmpdir, "_mock_info.txt")
    with open(fasta, "w") as f:
        f.write(">c0\nAAAA\n")
    with open(info, "w") as f:
        f.write("# col1\n")
        f.write(f"c0\t4\t1\t{circ}\tN\t1\t*\t*\n")
    return types.SimpleNamespace(
        input=types.SimpleNamespace(fasta=fasta, info=info),
        output=types.SimpleNamespace(
            completeness_check=os.path.join(tmpdir, "_mock_comp.txt"),
            fasta=os.path.join(tmpdir, "_mock_chrom.fa"),
            ignore_list=os.path.join(tmpdir, "_mock_ignore.txt"),
        ),
        params=types.SimpleNamespace(
            min_chrom_length="10",
            circular_chromosome=False,
            polypolish_flag=False,
        ),
    )


def _write_fasta(path, seq_id, length):
    with open(path, "w") as f:
        f.write(f">{seq_id}\n" + "A" * length + "\n")


def _write_flye_info(path, seq_id, circ):
    with open(path, "w") as f:
        f.write("# col1\n")
        f.write(f"{seq_id}\t100\t50\t{circ}\tN\t1\t*\t*\n")


def _import_scripts():
    """Import check_completeness and extract_chromosome with a mocked snakemake builtin.
    Safe to call multiple times — returns cached modules after the first import."""
    import builtins
    import importlib
    import types

    scripts_dir = str(Path("hybracter/workflow/scripts"))
    if scripts_dir not in sys.path:
        sys.path.insert(0, scripts_dir)

    # We need a real temp dir for the mock so module-level calls don't blow up
    _tmpdir = tempfile.mkdtemp()
    builtins.snakemake = _make_snakemake_mock(_tmpdir)

    # Force fresh import (in case a previous test already imported them)
    for mod in ("check_completeness", "extract_chromosome"):
        if mod in sys.modules:
            del sys.modules[mod]

    import check_completeness as cc
    import extract_chromosome as ec

    # Clean up the builtin — we don't need it after import
    if hasattr(builtins, "snakemake"):
        del builtins.snakemake
    shutil.rmtree(_tmpdir, ignore_errors=True)

    return cc, ec


class TestCheckCompleteness:
    """Unit tests for check_completeness.get_completeness."""

    @pytest.fixture(autouse=True, scope="class")
    def load_module(self, request):
        cc, _ = _import_scripts()
        request.cls.cc = cc

    def _run(self, seq_length, circ, min_chrom_length, circular_chromosome):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = os.path.join(tmpdir, "assembly.fasta")
            info = os.path.join(tmpdir, "assembly_info.txt")
            out = os.path.join(tmpdir, "completeness.txt")
            _write_fasta(fasta, "contig_1", seq_length)
            _write_flye_info(info, "contig_1", circ)
            self.cc.get_completeness(fasta, out, min_chrom_length, info, circular_chromosome)
            return open(out).read().strip()

    def test_noncircular_long_default_complete(self):
        """Non-circular contig above threshold → Complete with default settings (length only)."""
        assert self._run(200000, "N", 100000, False) == "C"

    def test_circular_long_default_complete(self):
        """Circular contig above threshold → Complete with default settings."""
        assert self._run(200000, "Y", 100000, False) == "C"

    def test_noncircular_long_circular_flag_incomplete(self):
        """Non-circular contig above threshold → Incomplete with --circular_chromosome."""
        assert self._run(200000, "N", 100000, True) == "I"

    def test_circular_long_circular_flag_complete(self):
        """Circular contig above threshold → Complete with --circular_chromosome."""
        assert self._run(200000, "Y", 100000, True) == "C"

    def test_short_always_incomplete(self):
        """Contig below threshold → Incomplete regardless of circularity or flag."""
        assert self._run(50000, "Y", 100000, False) == "I"
        assert self._run(50000, "Y", 100000, True) == "I"
        assert self._run(50000, "N", 100000, False) == "I"
        assert self._run(50000, "N", 100000, True) == "I"


class TestExtractChromosome:
    """Unit tests for extract_chromosome.get_chromosome_plasmids."""

    @pytest.fixture(autouse=True, scope="class")
    def load_module(self, request):
        _, ec = _import_scripts()
        request.cls.ec = ec

    def _run(self, seq_length, circ, min_chrom_length, circular_chromosome):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = os.path.join(tmpdir, "assembly.fasta")
            info = os.path.join(tmpdir, "assembly_info.txt")
            chrom_out = os.path.join(tmpdir, "chromosome.fasta")
            ignore_list = os.path.join(tmpdir, "ignore_list.txt")
            _write_fasta(fasta, "contig_1", seq_length)
            _write_flye_info(info, "contig_1", circ)
            self.ec.get_chromosome_plasmids(
                fasta, chrom_out, ignore_list, min_chrom_length, info, False, circular_chromosome
            )
            extracted = open(chrom_out).read()
            ignored = open(ignore_list).read().strip()
            return extracted, ignored

    def test_noncircular_long_extracted_default(self):
        """Non-circular contig above threshold IS extracted by default and added to ignore_list."""
        extracted, ignored = self._run(200000, "N", 100000, False)
        assert "contig_1" in extracted
        assert "contig_1" in ignored

    def test_circular_long_extracted_default_not_ignored(self):
        """Circular contig above threshold is extracted by default and NOT in ignore_list."""
        extracted, ignored = self._run(200000, "Y", 100000, False)
        assert "contig_1" in extracted
        assert ignored == ""

    def test_noncircular_long_not_extracted_circular_flag(self):
        """Non-circular contig above threshold is NOT extracted with --circular_chromosome."""
        extracted, ignored = self._run(200000, "N", 100000, True)
        assert "contig_1" not in extracted
        assert ignored == ""

    def test_circular_long_extracted_circular_flag_not_ignored(self):
        """Circular contig above threshold is extracted with --circular_chromosome and NOT in ignore_list."""
        extracted, ignored = self._run(200000, "Y", 100000, True)
        assert "contig_1" in extracted
        assert ignored == ""

    def test_short_contig_never_extracted(self):
        """Contig below threshold is never extracted regardless of circularity or flag."""
        for circ in ("Y", "N"):
            for flag in (True, False):
                extracted, _ = self._run(50000, circ, 100000, flag)
                assert "contig_1" not in extracted


def test_citation():
    """test hybracter citation"""
    cmd = "hybracter citation"
    exec_command(cmd)


def test_version():
    """test hybracter version"""
    cmd = "hybracter version"
    exec_command(cmd)



class TestErrors:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.threads = threads  

    def test_no_db(self, run_mac):
        """test hybracter with no db"""
        with pytest.raises(RuntimeError):
            outdir: Path = "test_hybracter_output"
            empty_db: Path = "tests/empty_db"
            cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --database {empty_db} --no_pypolca"
            if run_mac:
                cmd += " --mac"
            exec_command(cmd)
            remove_directory(outdir)

    def test_hybracter_t_hybrid_md(run_mac):
        """hybracter test-hybrid with --min_depth 1000"""
        with pytest.raises(RuntimeError):
            outdir: Path = "test_hybracter_output"
            cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --min_depth 1000"
            if run_mac:
                cmd += " --mac"
            exec_command(cmd)
            remove_directory(outdir)

    def test_hybracter_t_long_md(run_mac):
        """hybracter test-long --min_depth 1000 also test """
        with pytest.raises(RuntimeError):
            outdir: Path = "test_hybracter_output"
            cmd = f"hybracter test-long --threads {threads} --output {outdir} --min_depth 1000"
            if run_mac:
                cmd += " --mac"
            exec_command(cmd)
            remove_directory(outdir)


# cleanup
remove_directory(db_dir)
