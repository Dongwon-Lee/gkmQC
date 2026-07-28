"""Shared fixtures for the gkmQC test suite.

The suite is deliberately self-contained: every test builds the data it
needs under pytest's ``tmp_path``.  Nothing here reads a real genome or
a prebuilt null-seq index, so the whole suite runs in seconds on a
laptop and needs no reference downloads.
"""

import os
import random
import subprocess
import sys
import tarfile

import pytest


# ``gkmqc.cli.main()`` writes to module-level state (it declares
# ``global HEADER`` and appends to it), so driving it more than once in a
# single interpreter accumulates output.  CLI-level tests therefore run
# it in a subprocess.
_CLI_RUNNER = "import sys; from gkmqc.cli import main; main()"


@pytest.fixture(autouse=True)
def reset_data_dir_state():
    """Undo the global state ``_paths.set_data_dir()`` establishes.

    ``set_data_dir()`` caches into a module global *and* exports
    ``GKMQC_DATA_DIR`` so multiprocessing workers inherit it.  Without
    this fixture, whichever test ran first would decide the answer for
    every test after it.
    """
    from gkmqc import _paths

    saved = os.environ.get("GKMQC_DATA_DIR")
    _paths._DATA_DIR = None
    os.environ.pop("GKMQC_DATA_DIR", None)
    yield
    _paths._DATA_DIR = None
    if saved is None:
        os.environ.pop("GKMQC_DATA_DIR", None)
    else:
        os.environ["GKMQC_DATA_DIR"] = saved


@pytest.fixture
def run_cli():
    """Run the gkmQC CLI in a subprocess; return the CompletedProcess.

    ``python -c CODE a b c`` leaves ``sys.argv == ['-c', 'a', 'b', 'c']``,
    which is exactly what argparse expects, so the argv list is passed
    through untouched.
    """

    def _run(argv, cwd, env=None, check=False):
        full_env = dict(os.environ)
        full_env.pop("GKMQC_DATA_DIR", None)
        if env:
            full_env.update(env)
        return subprocess.run(
            [sys.executable, "-c", _CLI_RUNNER] + [str(a) for a in argv],
            cwd=str(cwd), env=full_env, check=check,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        )

    return _run


@pytest.fixture
def tiny_genome(tmp_path):
    """A 2-chromosome synthetic genome archive, as buildidx expects.

    Sequences are short (3 kb) and the tests drive buildidx with a small
    ``-w`` so the null-seq index stays tiny.  The alphabet includes ``N``
    and lowercase (repeat-masked) bases so the N/GC/repeat bit arrays all
    get exercised.
    """
    rng = random.Random(7)
    src = tmp_path / "src_fa"
    src.mkdir()
    names = ["chrA", "chrB"]
    for chrom in names:
        seq = "".join(rng.choice("ACGTacgtN") for _ in range(3000))
        with open(src / ("%s.fa" % chrom), "w") as fh:
            fh.write(">%s\n" % chrom)
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    archive = tmp_path / "tiny.chromFa.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for chrom in names:
            tar.add(src / ("%s.fa" % chrom), arcname="%s.fa" % chrom)

    return {"archive": archive, "chroms": names, "genome": "tiny"}


@pytest.fixture
def data_dir(tmp_path):
    """An empty, writable gkmQC data directory."""
    d = tmp_path / "data"
    d.mkdir()
    return d
