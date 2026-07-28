"""Data-directory resolution.

Before the package layout landed, ``base_data_dir`` was derived from the
location of ``bin/gkmqc.py``, so it was correct wherever you ran from.
Deriving it from the CWD instead made it order-dependent, and made
"resolved to the wrong place" a silent failure.  These tests pin the
resolution order down and keep the failure loud.
"""

import os

import pytest

from gkmqc import _paths


def test_cwd_named_data_resolves_to_itself(tmp_path):
    """``cd data && gkmqc buildidx ...`` -- the documented README flow.

    The naive ``$CWD/data`` rule looks for ``<repo>/data/data`` here,
    misses, and silently falls through to a package-adjacent guess.
    """
    d = tmp_path / "data"
    d.mkdir()
    os.chdir(d)
    assert _paths.set_data_dir() == str(d)


def test_data_under_cwd_resolves(tmp_path):
    """Running from the repository root."""
    d = tmp_path / "data"
    d.mkdir()
    os.chdir(tmp_path)
    assert _paths.set_data_dir() == str(d)


def test_env_var_beats_cwd(tmp_path):
    cwd_data = tmp_path / "data"
    cwd_data.mkdir()
    other = tmp_path / "elsewhere"
    other.mkdir()
    os.chdir(tmp_path)
    os.environ["GKMQC_DATA_DIR"] = str(other)
    assert _paths.set_data_dir() == str(other)


def test_explicit_beats_env_var(tmp_path):
    env_dir = tmp_path / "from_env"
    env_dir.mkdir()
    cli_dir = tmp_path / "from_cli"
    cli_dir.mkdir()
    os.environ["GKMQC_DATA_DIR"] = str(env_dir)
    os.chdir(tmp_path)
    assert _paths.set_data_dir(str(cli_dir)) == str(cli_dir)


def test_unresolvable_raises_rather_than_guessing(tmp_path, monkeypatch):
    """The regression that motivated this module.

    The old resolver returned a package-adjacent path unconditionally.
    Under a non-editable install that is ``<site-packages>/data``, so the
    run died much later inside ``os.mkdir`` with a bare FileNotFoundError
    instead of saying what was actually wrong.
    """
    work = tmp_path / "no_data_here"
    work.mkdir()
    os.chdir(work)
    # Simulate an installed package with no adjacent data/ directory.
    monkeypatch.setattr(_paths, "__file__", str(tmp_path / "fake_pkg" / "_paths.py"))
    with pytest.raises(RuntimeError, match="could not locate the gkmQC data directory"):
        _paths.set_data_dir()


def test_missing_explicit_dir_rejected_unless_creating(tmp_path):
    """``evaluate`` needs an existing index; ``buildidx`` may create one.

    Without this split, a typo in ``-D`` silently builds an empty tree and
    fails much later with a confusing "index not found".
    """
    missing = tmp_path / "typo"
    with pytest.raises(RuntimeError, match="does not exist"):
        _paths.set_data_dir(str(missing))
    assert _paths.set_data_dir(str(missing), create_ok=True) == str(missing)


def test_resolution_is_exported_for_child_processes(tmp_path):
    """Regression: the resolved directory must reach Pool workers.

    A module global does not survive the ``spawn`` and ``forkserver``
    start methods -- the child re-imports and would resolve again from
    its own CWD, which by then is the gkmQC output directory.  Exporting
    into the environment is what makes the value inheritable.
    """
    d = tmp_path / "data"
    d.mkdir()
    _paths.set_data_dir(str(d))
    assert os.environ["GKMQC_DATA_DIR"] == str(d)


def test_get_data_dir_is_lazy(tmp_path):
    """Nothing resolves at import time.

    ``cli.main()`` chdir()s into the output directory before the null-seq
    pool runs, so a value computed at import time would be computed
    against the wrong directory depending on import order.
    """
    d = tmp_path / "data"
    d.mkdir()
    os.chdir(tmp_path)
    assert _paths._DATA_DIR is None
    assert _paths.get_data_dir() == str(d)
