import os
import shutil
import subprocess
from pathlib import Path
import unittest
import pytest

__author__ = "George Bouras"
__copyright__ = "Copyright 2023, Hybracter"
__license__ = "MIT"
__type__ = "Test Script"
__maintainer__ = "George Bouras"
__email__ = "george.bouras@adelaide.edu.au"



@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    """set the working directory for all tests"""
    monkeypatch.chdir(tmp_dir)


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


def test_hybracter_hybrid():
    """test hybracter hybrid"""
    threads = 4
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --no_polca"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_long():
    """test hybracter long"""
    threads = 4
    outdir: Path = "test_hybracter_output"
    cmd = f"hybracter test-long --threads {threads} --output {outdir}"
    exec_command(cmd)
    remove_directory(outdir)


def test_hybracter_install():
    """test hybracter install"""
    outdir: Path = "plassembler_db"
    cmd = f"hybracter install  --databases {outdir}"
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


class errors(unittest.TestCase):

    def test_no_db(self):
        """test hybracter with no db"""
        with self.assertRaises(RuntimeError):
            threads = 4
            outdir: Path = "test_hybracter_output"
            empty_db: Path = "tests/empty_db"
            cmd = f"hybracter test-hybrid --threads {threads} --output {outdir} --database {empty_db} --no_polca"
            exec_command(cmd)
            remove_directory(outdir)