"""
Golden (assembly-output regression) tests.

Runs the bundled `hybracter test-long` end-to-end and compares the FINAL_OUTPUT
summary against a committed reference, so structural regressions surface in CI:
lost/extra contigs (cf. #133), completeness misclassification, plasmid-count
drift (cf. #154), or large assembly-length changes.

Assertions are deliberately tolerant of single-base polishing differences
(which legitimately change with medaka/tool versions):
  * exact: Complete, Number_of_contigs, Number_circular_plasmids
  * within tolerance: Total_assembly_length, Longest_contig_length (+/- 2%)

The references are generated on Linux. canu's long-read plasmid assembly is
machine-dependent on the tiny simulated plasmid (a given machine is
self-consistent, but different machines can recover or miss it), so OFF the
reference platform (e.g. the macOS CI runner) only the chromosome-level fields
(Complete, Longest_contig_length) are asserted - plasmid-recovery-dependent
fields are checked strictly on Linux only.

Bootstrapping the reference (run once, then commit the generated file):
    HYBRACTER_UPDATE_GOLDEN=1 pytest -m integration tests/test_golden.py
"""

import csv
import os
import platform
import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).resolve().parents[1]
GOLDEN_DIR = REPO_ROOT / "tests" / "golden"

# The golden references are generated on Linux. On the reference platform we
# assert everything strictly. Off it (e.g. macOS / Apple Silicon CI) we assert
# only the rock-stable chromosome-level fields: the long-read plasmid assembly
# (canu) is *machine-dependent* on the tiny simulated plasmid - a given machine
# is self-consistent, but different machines can recover or miss the plasmid
# (it sits far above the depth filter, so this is the assembler, not filtering).
# Plasmid-recovery-dependent fields (Number_of_contigs, Number_circular_plasmids,
# Total_assembly_length) therefore can't be matched against a Linux-generated
# golden cross-machine. Plasmid-loss regressions are still caught on Linux.
REFERENCE_PLATFORM = "Linux"

EXACT_FIELDS = ["Complete", "Number_of_contigs", "Number_circular_plasmids"]
TOLERANT_FIELDS = {"Total_assembly_length": 0.02, "Longest_contig_length": 0.02}

# chromosome-only subsets used off the reference platform
EXACT_FIELDS_OFF_REF = ["Complete"]
TOLERANT_FIELDS_OFF_REF = {"Longest_contig_length": 0.02}

THREADS = 4


def _summary_path(outdir):
    return Path(outdir) / "FINAL_OUTPUT" / "hybracter_summary.tsv"


def _read_summary(path):
    with open(path) as fh:
        return {row["Sample"]: row for row in csv.DictReader(fh, delimiter="\t")}


def _run_and_compare(subcommand, golden_name, tmp_path):
    outdir = tmp_path / subcommand
    cmd = f"hybracter {subcommand} --threads {THREADS} --output {outdir} --skip_qc"
    # Ensure TERM is set: without it, unicycler emits `tput: No value for $TERM`
    # on stderr, which plassembler's dependency check merges into stdout and then
    # mis-parses ("Unicycler not found"). CI sets TERM=linux; do the same here so
    # the test is self-contained when run outside CI.
    env = {**os.environ, "TERM": os.environ.get("TERM") or "xterm"}
    subprocess.run(cmd, shell=True, check=True, env=env)

    produced = _summary_path(outdir)
    assert produced.is_file(), f"no summary produced at {produced}"

    golden = GOLDEN_DIR / golden_name

    # Bootstrap mode: write the observed summary as the new reference and skip.
    if os.environ.get("HYBRACTER_UPDATE_GOLDEN"):
        GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
        shutil.copy(produced, golden)
        pytest.skip(f"HYBRACTER_UPDATE_GOLDEN set: wrote reference -> {golden}")

    if not golden.is_file():
        pytest.skip(
            f"No golden reference at {golden}. Generate it once with "
            f"HYBRACTER_UPDATE_GOLDEN=1 pytest -m integration tests/test_golden.py and commit it."
        )

    observed = _read_summary(produced)
    expected = _read_summary(golden)

    assert set(observed) == set(expected), (
        f"sample set changed: produced {sorted(observed)} vs golden {sorted(expected)}"
    )

    if platform.system() == REFERENCE_PLATFORM:
        exact_fields, tolerant_fields = EXACT_FIELDS, TOLERANT_FIELDS
    else:
        # off the reference platform: chromosome-level fidelity only (canu
        # long-read plasmid recovery is machine-dependent on the test data)
        exact_fields, tolerant_fields = EXACT_FIELDS_OFF_REF, TOLERANT_FIELDS_OFF_REF

    for sample, exp in expected.items():
        obs = observed[sample]
        for field in exact_fields:
            assert obs[field] == exp[field], (
                f"{sample}: {field} regressed ({obs[field]!r} != golden {exp[field]!r})"
            )
        for field, tol in tolerant_fields.items():
            o, e = float(obs[field]), float(exp[field])
            assert abs(o - e) <= tol * e, (
                f"{sample}: {field} drifted {o} vs golden {e} (> {tol * 100:.0f}%)"
            )


def test_golden_long(tmp_path):
    _run_and_compare("test-long", "test_long_summary.tsv", tmp_path)


def test_golden_hybrid(tmp_path):
    _run_and_compare("test-hybrid", "test_hybrid_summary.tsv", tmp_path)
