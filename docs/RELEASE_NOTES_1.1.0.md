# gkmQC 1.1.0

Python 3 packaging overhaul. gkmQC now installs as a proper Python package with a
`gkmqc` console script, and no longer depends on the unmaintained `pyfasta`.

Thanks to @chlee-tabin for #5, which this release is built on.

## Breaking changes

Existing scripts and pipelines need updating:

- **The command is now `gkmqc`, not `gkmqc.py`.** Replace `../bin/gkmqc.py <cmd>` with
  `gkmqc <cmd>`. There is no `gkmqc.py` alias — shipping both would let the plain script
  shadow the package on import.
- **Installation now requires `pip install .`**, and the order matters:
  ```bash
  cd src && make && make install   # builds the C library into the package
  cd .. && pip install .
  ```
  Running `pip install` first produces a package with no C library. If that happens, gkmQC
  now tells you so explicitly instead of raising a bare `OSError`.
- **Python >= 3.10** (was 3.7, EOL June 2023). Nothing in gkmQC's own code forced this;
  3.7 simply cannot install the current numpy/scipy/scikit-learn/matplotlib stack.
- **`pyfasta` is replaced by `pyfaidx`.** Update your environment
  (`conda env create -f environment.yml` handles it). Sequence output is unchanged.
- **Module imports moved into the package**: `import seqs_nullgen` becomes
  `from gkmqc import seqs_nullgen`. Only affects code importing gkmQC internals directly.
- **The data directory is now resolved explicitly.** Order: `-D/--data-dir`, then
  `$GKMQC_DATA_DIR`, then `./data`, then the current directory if it is itself named `data`.
  If none resolve, gkmQC stops with an error instead of guessing.

  This is the one change most likely to bite an existing workflow: running from a directory
  with no `data/` in it — the README's `cd test && gkmqc evaluate ...` flow — now needs
  `-D` or `$GKMQC_DATA_DIR`:
  ```bash
  export GKMQC_DATA_DIR=/path/to/gkmQC/data
  ```
  `buildidx` and `evaluate` accept `-D`; `optimize` and `report` do not use the data
  directory.

## Fixes

- `buildidx` run as `cd data && gkmqc buildidx ...` wrote to the wrong location under a
  non-editable install — it resolved the data directory into `site-packages` and failed.
- `make install` failed on a fresh clone when no conda environment was active, because the
  `bin/` directory no longer exists in a clean checkout.
- A missing `gkmkern_pylib.so` now produces an actionable message naming the build steps.
- `buildidx` writes the `.fai` FASTA index as it builds, rather than leaving the first
  `evaluate` run to create it.

## Known issues

- The precomputed index tarballs still ship `pyfasta`'s `.flat`/`.gdx` files and no `.fai`;
  `pyfaidx` builds the `.fai` on first use, which requires the reference directory to be
  writable. Regenerating the tarballs is planned.

- On Python older than 3.13, `pytest tests/` reports one failure in
  `test_cli.py::test_data_dir_uses_capital_d`. This is a defect in the test, not in gkmQC:
  the test matched `argparse`'s help layout, which changed in 3.13 (`-D DATA_DIR,
  --data-dir DATA_DIR` before, `-D, --data-dir DATA_DIR` after). The flags themselves are
  correct on every supported version. Fixed on `main` after this tag.
