"""Path resolution helpers for gkmQC.

Pre-package layout, ``base_data_dir`` was derived from the location of
``bin/gkmqc.py`` via ``os.path.dirname(os.path.dirname(__file__))``.  Once
the scripts are installed as a package under site-packages, that trick
points into the Python install — not where the user stored their null-seq
indexes.

Resolution order (highest precedence first):

1. ``--data-dir`` / ``-D`` on the command line, via :func:`set_data_dir`.
2. ``$GKMQC_DATA_DIR`` environment variable (set once, reuse across
   commands).
3. ``./data`` under the current working directory — the flow that runs
   gkmQC from the repository root.
4. The current working directory itself when it is named ``data`` — the
   README flow that does ``cd data && gkmqc buildidx ...``.
5. Package-adjacent ``../data`` — only useful when running out of an
   editable install from a source checkout.

If none of these resolve, :func:`get_data_dir` raises rather than handing
back a path that does not exist; the previous implementation returned the
package-adjacent candidate unconditionally, which under a non-editable
install pointed into site-packages and failed much later with a confusing
``FileNotFoundError``.

Resolution is deliberately *lazy*: nothing is computed at import time, so
a module that imports this one after an ``os.chdir()`` still sees the
directory that was resolved up front by :func:`set_data_dir`.
"""

import os

_DATA_DIR = None


def _resolve() -> str:
    env = os.environ.get("GKMQC_DATA_DIR")
    if env:
        return os.path.abspath(env)

    cwd = os.getcwd()

    cwd_data = os.path.join(cwd, "data")
    if os.path.isdir(cwd_data):
        return cwd_data

    # ``cd data && gkmqc buildidx ...`` — the CWD *is* the data directory.
    if os.path.basename(cwd) == "data":
        return cwd

    # fall back to a package-adjacent ``data/`` so editable installs
    # (``pip install -e .``) keep working out of the source tree.
    pkg_adjacent = os.path.abspath(
        os.path.join(os.path.dirname(os.path.realpath(__file__)), os.pardir, "data")
    )
    if os.path.isdir(pkg_adjacent):
        return pkg_adjacent

    raise RuntimeError(
        "could not locate the gkmQC data directory. Pass -D/--data-dir, or set "
        "GKMQC_DATA_DIR, or run from a directory containing 'data/'."
    )


def set_data_dir(path=None, create_ok=False) -> str:
    """Resolve the data directory once, eagerly.

    Call this from ``main()`` *before* any ``os.chdir()``, so that the
    CWD-relative candidates above are evaluated against the directory the
    user actually invoked gkmQC from.

    ``path`` is the explicit ``--data-dir`` value, or ``None`` to fall
    through to the resolution order documented above.  When ``path`` is
    given but does not exist, ``create_ok`` decides: ``buildidx`` may
    create its data directory, everything else must be pointed at one
    that already exists (a typo would otherwise silently produce an empty
    tree and a much later "index not found" error).
    """
    global _DATA_DIR

    if path:
        resolved = os.path.abspath(path)
        if not create_ok and not os.path.isdir(resolved):
            raise RuntimeError(
                "--data-dir does not exist or is not a directory: %s" % resolved
            )
        _DATA_DIR = resolved
    else:
        _DATA_DIR = _resolve()

    return _DATA_DIR


def get_data_dir() -> str:
    """Return the resolved data directory, resolving on first use."""
    global _DATA_DIR
    if _DATA_DIR is None:
        _DATA_DIR = _resolve()
    return _DATA_DIR
