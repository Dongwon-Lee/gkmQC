"""Setup script for gkmQC.

Install the Python side of gkmQC.  Before running ``pip install .`` (or
``pip install -e .``), build the C library so ``gkmkern_pylib.so`` is
copied next to the package:

    cd src && make && make install

The Makefile's ``install`` target drops ``gkmkern_pylib.so`` into
``gkmqc/`` so this setup.py can ship it as package data, and drops the
standalone ``gkmkern`` CLI into ``$CONDA_PREFIX/bin`` (or ``bin/`` for
non-conda flows).
"""

import os
import re

from setuptools import setup, find_packages


def _version():
    """Single-source the version from gkmqc/__init__.py.

    Read rather than imported: importing the package at build time would
    pull in numpy/pyfaidx before install_requires has been satisfied.
    """
    init = os.path.join(os.path.dirname(__file__), "gkmqc", "__init__.py")
    with open(init) as fh:
        m = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', fh.read(), re.M)
    if not m:
        raise RuntimeError("could not find __version__ in %s" % init)
    return m.group(1)


setup(
    name="gkmqc",
    version=_version(),
    description="gkmQC: gapped k-mer-SVM quality check for chromatin accessibility data",
    author="Seong Kyu Han, Dongwon Lee",
    author_email="dongwon.lee@childrens.harvard.edu",
    url="https://github.com/Dongwon-Lee/gkmQC",
    license="GPLv3",
    python_requires=">=3.10",
    packages=find_packages(include=["gkmqc", "gkmqc.*"]),
    package_data={
        # Ship the compiled C kernel and the SLURM sbatch wrapper
        # alongside the Python package so gkmsvm.computeGkmKernel and
        # the SLURM launch path can locate them via
        # os.path.dirname(__file__) — works for editable installs,
        # wheel installs, and source-tree execution.
        "gkmqc": ["gkmkern_pylib.so", "gkmsvm_slurm.sh"],
    },
    install_requires=[
        "numpy",
        "scipy",
        "scikit-learn",
        "bitarray",
        "matplotlib",
        # pyfasta (the previous dep) was abandoned in 2015 and is
        # Python-2-only; pyfaidx is a maintained drop-in for the one
        # call site this project uses.
        "pyfaidx",
    ],
    entry_points={
        "console_scripts": [
            "gkmqc=gkmqc.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
