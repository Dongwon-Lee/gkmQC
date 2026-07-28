"""Multiprocessing start-method regressions.

Python 3.14 changed the default start method on Linux from ``fork`` to
``forkserver``.  Under ``fork`` a Pool worker inherits the parent's
memory; under ``forkserver`` and ``spawn`` it re-imports the module and
sees module globals unset.  gkmQC shares two things with its workers
that way, and both broke.

Each case runs in a subprocess because the start method is process-wide
and cannot be changed back once pools exist.
"""

import subprocess
import sys

import pytest

START_METHODS = ["fork", "forkserver", "spawn"]


CROSSVALIDATE_UNDER = """
import sys, multiprocessing as mp
def main(method):
    mp.set_start_method(method, force=True)
    import logging, numpy as np
    logging.disable(logging.INFO)
    from gkmqc import gkmsvm
    n = 40
    rs = np.random.RandomState(0)
    kmat = np.eye(n * 2) + rs.rand(n * 2, n * 2) * 0.01
    kmat = np.maximum(kmat, kmat.T)
    # C, tol, shrinking, cache, ncv, repeats, fast_estimation, seed, procs
    args_svm = (1.0, 0.001, 1, 100, 5, 1, 0, 1, 2)
    auc, _ = gkmsvm.crossValidate(args_svm, kmat, n, n)
    print("AUC=%.6f" % auc)
if __name__ == "__main__":
    main(sys.argv[1])
"""


@pytest.mark.parametrize("method", START_METHODS)
def test_crossvalidate_survives_start_method(method, tmp_path):
    """Regression: ``NameError: name 'kmat' is not defined``.

    ``crossValidate`` hands the kernel matrix to its workers through a
    module global rather than as an argument -- pickling a multi-GB
    matrix per worker would cost far more than the fork.  That sharing
    only works under ``fork``, so the pool must request it explicitly
    rather than inherit whatever the platform defaults to.

    On 3.14 this failed *after* the full kernel computation, so a user
    paid the entire runtime before seeing the error.
    """
    script = tmp_path / "cv.py"
    script.write_text(CROSSVALIDATE_UNDER)
    proc = subprocess.run(
        [sys.executable, str(script), method],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=300,
    )
    assert proc.returncode == 0, "crossValidate failed under %r:\n%s" % (method, proc.stdout)
    assert "AUC=" in proc.stdout


def test_crossvalidate_result_is_start_method_independent(tmp_path):
    """Requesting fork must not change the numbers it produces."""
    script = tmp_path / "cv.py"
    script.write_text(CROSSVALIDATE_UNDER)
    seen = {}
    for method in START_METHODS:
        proc = subprocess.run(
            [sys.executable, str(script), method],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=300,
        )
        assert proc.returncode == 0, proc.stdout
        line = [l for l in proc.stdout.splitlines() if l.startswith("AUC=")][-1]
        seen[method] = line
    assert len(set(seen.values())) == 1, "AUC differs by start method: %r" % seen


DATA_DIR_IN_WORKER = """
import sys, os, multiprocessing as mp
def worker(_):
    from gkmqc import _paths
    return _paths.get_data_dir()
def main(method, data_dir, elsewhere):
    mp.set_start_method(method, force=True)
    from gkmqc import _paths
    _paths.set_data_dir(data_dir)
    os.chdir(elsewhere)          # cli.main() chdir()s to the output dir
    with mp.Pool(2) as pool:
        seen = set(pool.map(worker, [0, 1]))
    print("WORKER_SAW=%s" % ("|".join(sorted(seen))))
if __name__ == "__main__":
    main(*sys.argv[1:])
"""


@pytest.mark.parametrize("method", START_METHODS)
def test_data_dir_reaches_workers(method, tmp_path):
    """Regression: ``could not locate the gkmQC data directory``.

    The null-seq sampler calls ``get_data_dir()`` inside the worker.
    ``cli.main()`` has already chdir()'d into the output directory by
    then, so a worker that re-resolves gets it wrong -- or, once the
    resolver stopped guessing, raises outright.
    """
    script = tmp_path / "dd.py"
    script.write_text(DATA_DIR_IN_WORKER)
    data = tmp_path / "data"
    data.mkdir()
    out = tmp_path / "output_dir"
    out.mkdir()
    proc = subprocess.run(
        [sys.executable, str(script), method, str(data), str(out)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=300,
    )
    assert proc.returncode == 0, "worker lost the data dir under %r:\n%s" % (method, proc.stdout)
    assert "WORKER_SAW=%s" % data in proc.stdout, proc.stdout
