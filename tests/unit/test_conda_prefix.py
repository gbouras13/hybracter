"""
Pins the #166 conda-prefix fallback: a read-only default prefix (shared conda
env) falls back to a user cache dir, unless it already holds built envs
(the container) or the user passed --conda-prefix explicitly.
"""

import os

import pytest

import hybracter.util as util


@pytest.fixture
def default_prefix(tmp_path, monkeypatch):
    pkg = tmp_path / "pkg"
    (pkg / "workflow").mkdir(parents=True)
    monkeypatch.setattr(util, "snake_base", lambda rel: str(pkg / rel))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    yield pkg / "workflow" / "conda"
    # restore permissions so pytest can clean up tmp_path
    for p in (pkg / "workflow", pkg / "workflow" / "conda"):
        if p.exists():
            os.chmod(p, 0o755)


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores permissions")
def test_writable_default_is_kept(default_prefix):
    assert util.resolve_conda_prefix(str(default_prefix)) == str(default_prefix)
    assert default_prefix.is_dir()


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores permissions")
def test_readonly_empty_default_falls_back(default_prefix, tmp_path):
    default_prefix.mkdir()
    os.chmod(default_prefix, 0o555)
    os.chmod(default_prefix.parent, 0o555)
    resolved = util.resolve_conda_prefix(str(default_prefix))
    assert resolved == str(tmp_path / "cache" / "hybracter" / "conda")
    assert os.path.isdir(resolved)


@pytest.mark.skipif(os.geteuid() == 0, reason="root ignores permissions")
def test_readonly_default_with_built_envs_is_kept(default_prefix):
    (default_prefix / "0123abcd_").mkdir(parents=True)
    os.chmod(default_prefix, 0o555)
    assert util.resolve_conda_prefix(str(default_prefix)) == str(default_prefix)


def test_explicit_prefix_is_honoured(default_prefix, tmp_path):
    custom = str(tmp_path / "does" / "not" / "matter")
    assert util.resolve_conda_prefix(custom) == custom
