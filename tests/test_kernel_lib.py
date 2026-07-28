"""Locating the compiled gkm-kernel library.

``gkmkern_pylib.so`` is built by ``cd src && make && make install`` and
shipped as package data by ``pip install .``.  Running those two the
other way round installs a package with no C library, which used to
surface as a bare ``OSError`` from ctypes.
"""

import numpy as np
import pytest

from gkmqc import gkmsvm


# C, L, k, d, M, H, gamma, posfile, negfile, nthreads, verbosity
GKM_ARGS = [4, 10, 6, 3, 50, 0.5, 2.0, "pos.fa", "neg.fa", 1, 1]


def test_library_ships_next_to_the_package():
    """bin_dir must point into the package, not a repo-relative ../bin.

    The pre-package layout looked for the .so beside the checkout, which
    does not survive installation into site-packages.
    """
    import os
    from gkmqc import __file__ as pkg_file
    assert gkmsvm.bin_dir == os.path.dirname(os.path.realpath(pkg_file))


def test_missing_library_explains_the_build_order(tmp_path, monkeypatch):
    monkeypatch.setattr(gkmsvm, "bin_dir", str(tmp_path))  # empty directory
    with pytest.raises(RuntimeError) as excinfo:
        gkmsvm.computeGkmKernel(list(GKM_ARGS))
    message = str(excinfo.value)
    assert "gkmkern_pylib.so not found" in message
    assert "make install" in message
    assert "pip install" in message


@pytest.mark.skipif(
    not __import__("os").path.exists(
        __import__("os").path.join(gkmsvm.bin_dir, "gkmkern_pylib.so")
    ),
    reason="C library not built; run 'cd src && make && make install'",
)
def test_library_loads_when_present():
    np.ctypeslib.load_library("gkmkern_pylib.so", gkmsvm.bin_dir)
