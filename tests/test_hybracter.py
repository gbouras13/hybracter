"""
# to run on mac
pytest --run_mac .  
"""

import os
import shutil
import subprocess
import unittest
from pathlib import Path

import pytest

__author__ = "George Bouras"
__copyright__ = "Copyright 2023, Hybracter"
__license__ = "MIT"
__type__ = "Test Script"
__maintainer__ = "George Bouras"
__email__ = "george.bouras@adelaide.edu.au"

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
    """hybracter test-long with depth_filter"""
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-long --threads {threads} --output {outdir}--skip_qc --depth_filter 2"
    if run_mac:
        cmd += " --mac"
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
