"""
CLI smoke tests for the hybracter click group.

`--help` for the group and every subcommand must exit 0 - this catches broken
option definitions and import-time errors in __main__.py without running any
pipeline. Also pins the #119 sample-name sanitisation helper.
"""

import pytest
from click.testing import CliRunner

from hybracter.__main__ import cli
from hybracter.util import sanitise_sample_name


def test_group_help_exits_zero():
    result = CliRunner().invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "hybracter" in result.output.lower()


@pytest.mark.parametrize("command", sorted(cli.commands))
def test_subcommand_help_exits_zero(command):
    result = CliRunner().invoke(cli, [command, "--help"])
    assert result.exit_code == 0, f"{command} --help failed:\n{result.output}"


def test_version_command_exits_zero():
    result = CliRunner().invoke(cli, ["version"])
    assert result.exit_code == 0


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("my sample ", "my_sample"),  # #119: trailing space + internal space
        ("  leading", "leading"),
        ("trailing  ", "trailing"),
        ("a b c", "a_b_c"),
        ("clean", "clean"),
    ],
)
def test_sanitise_sample_name(raw, expected):
    assert sanitise_sample_name(raw) == expected
