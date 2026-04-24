"""Path resolution helpers for gkmQC.

Pre-package layout, ``base_data_dir`` was derived from the location of
``bin/gkmqc.py`` via ``os.path.dirname(os.path.dirname(__file__))``.  Once
the scripts are installed as a package under site-packages, that trick
points into the Python install — not where the user stored their null-seq
indexes.

Resolution order (first non-empty wins):

1. ``$GKMQC_DATA_DIR`` environment variable (preferred — set once and
   reuse across commands).
2. ``./data`` under the current working directory (preserves the README
   flow that does ``cd data && ../bin/gkmqc.py ...``; if the user is one
   level above ``data/`` that still resolves correctly).
3. Package-adjacent ``../data`` — only useful when running out of an
   editable install from a source checkout.
"""

import os


def _resolve_base_data_dir() -> str:
    env = os.environ.get("GKMQC_DATA_DIR")
    if env:
        return env

    cwd_data = os.path.join(os.getcwd(), "data")
    if os.path.isdir(cwd_data):
        return cwd_data

    # fall back to a package-adjacent ``data/`` so editable installs
    # (``pip install -e .``) keep working out of the source tree.
    pkg_adjacent = os.path.abspath(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, "data")
    )
    return pkg_adjacent


base_data_dir = _resolve_base_data_dir()
