# gkmQC test suite

Self-contained: every test builds what it needs in a temporary directory.
No genome download, no prebuilt null-seq index, no network. The whole
suite runs in about two minutes.

## Running

```bash
pip install -e ".[test]"      # or: pip install pytest
pytest tests/                 # or: ./tests/run_tests.sh
```

Build the C library first if you want `test_kernel_lib.py` to check that
the shipped `gkmkern_pylib.so` actually loads — that one test skips when
the library is absent, everything else runs either way:

```bash
cd src && make && make install && cd ..
pip install -e ".[test]"
```

## What is covered

| File | Covers |
|---|---|
| `test_paths.py` | Data-directory resolution order, and that an unresolvable directory raises instead of guessing |
| `test_workers.py` | Both multiprocessing start-method regressions, under `fork`, `forkserver` and `spawn` |
| `test_cli.py` | `-D/--data-dir` wiring per subcommand, flag letters, error messages |
| `test_buildidx.py` | End-to-end index build on a synthetic 2-chromosome genome, including the `.fai` |
| `test_kernel_lib.py` | Locating `gkmkern_pylib.so`, and the build-order error message |

## Why the start-method tests exist

Python 3.14 changed the default `multiprocessing` start method on Linux
from `fork` to `forkserver`. Under `fork` a pool worker inherits the
parent's memory; under `forkserver` and `spawn` it re-imports the module
and sees module globals unset.

gkmQC shared two things with workers that way, and both broke on 3.14:

- the resolved data directory, which the null-seq sampler reads inside
  the worker — now exported through the environment so children inherit
  it under any start method;
- the kernel matrix in `crossValidate`, passed via a module global to
  avoid pickling a multi-GB array to every worker — the pool now requests
  the `fork` context explicitly rather than taking the platform default.

The second failed *after* the full kernel computation, so the whole
runtime was spent before the error appeared. `test_workers.py`
parametrises over all three start methods so a future default change
cannot reintroduce either silently.

## Not covered here

A full `evaluate` run needs a real genome index and takes 1–2 hours. It
is not part of this suite. To run one by hand:

```bash
export GKMQC_DATA_DIR=/path/to/gkmQC/data
cd test
gkmqc evaluate -i foo.narrowPeak -g hg38 -n foo -@ 10
```
