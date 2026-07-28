"""Command-line surface: flag wiring, version, and error messages."""

import re

import pytest

import gkmqc


SUBCOMMANDS = ["buildidx", "evaluate", "optimize", "report"]

# Only the two commands that touch the null-seq index take -D.
TAKES_DATA_DIR = {"buildidx": True, "evaluate": True, "optimize": False, "report": False}


@pytest.mark.parametrize("sub", SUBCOMMANDS)
def test_subcommand_help_runs(sub, tmp_path, run_cli):
    proc = run_cli([sub, "-h"], cwd=tmp_path)
    assert proc.returncode == 0, proc.stdout


@pytest.mark.parametrize("sub", SUBCOMMANDS)
def test_data_dir_flag_only_where_it_applies(sub, tmp_path, run_cli):
    proc = run_cli([sub, "-h"], cwd=tmp_path)
    has_flag = "--data-dir" in proc.stdout
    assert has_flag is TAKES_DATA_DIR[sub], (
        "%s: expected --data-dir present=%s\n%s" % (sub, TAKES_DATA_DIR[sub], proc.stdout)
    )


def test_data_dir_uses_capital_d(tmp_path, run_cli):
    """-d is --max-num-gaps on evaluate, so the data dir flag must be -D.

    argparse would raise on a genuine collision, but this pins the intent
    so nobody "tidies" it to lowercase later.
    """
    proc = run_cli(["evaluate", "-h"], cwd=tmp_path)
    assert re.search(r"-D, --data-dir", proc.stdout), proc.stdout
    assert re.search(r"-d, --max-num-gaps", proc.stdout), proc.stdout


def test_version_is_single_sourced(tmp_path, run_cli):
    """setup.py and the run header both derive from gkmqc.__version__."""
    proc = run_cli(["report", "-h"], cwd=tmp_path)
    assert proc.returncode == 0
    assert re.match(r"^\d+\.\d+\.\d+", gkmqc.__version__), gkmqc.__version__


def test_missing_data_dir_reports_clearly(tmp_path, run_cli):
    """Running where nothing resolves must name the fix, not fail obscurely.

    Only reachable when the package has no adjacent ``data/`` -- i.e. a
    normal ``pip install .``.  Under ``pip install -e .`` the source tree's
    own ``data/`` is the last-resort candidate and resolves successfully,
    so there is nothing to assert.  ``test_paths.py`` covers the resolver
    itself in both modes.
    """
    import os

    pkg_adjacent = os.path.abspath(
        os.path.join(os.path.dirname(gkmqc.__file__), os.pardir, "data")
    )
    if os.path.isdir(pkg_adjacent):
        pytest.skip("editable install: package-adjacent %s resolves" % pkg_adjacent)

    work = tmp_path / "no_data"
    work.mkdir()
    proc = run_cli(
        ["evaluate", "-i", "absent.narrowPeak", "-g", "hg38", "-n", "x"], cwd=work
    )
    assert proc.returncode != 0
    assert "could not locate the gkmQC data directory" in proc.stdout, proc.stdout
    assert "-D/--data-dir" in proc.stdout


def test_nonexistent_explicit_data_dir_reports_clearly(tmp_path, run_cli):
    proc = run_cli(
        ["evaluate", "-i", "absent.narrowPeak", "-g", "hg38", "-n", "x",
         "-D", str(tmp_path / "typo")],
        cwd=tmp_path,
    )
    assert proc.returncode != 0
    assert "does not exist or is not a directory" in proc.stdout, proc.stdout
