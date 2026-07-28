#!/bin/bash
# Run the gkmQC test suite. Any extra arguments are passed to pytest,
# e.g. ./tests/run_tests.sh -k buildidx -v
set -o errexit
set -o nounset

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! python -c "import pytest" 2>/dev/null; then
    echo "pytest not found. Install it with:  pip install -e \".[test]\"" >&2
    exit 1
fi

if ! python -c "import gkmqc" 2>/dev/null; then
    echo "gkmqc not importable. Install it with:  pip install -e ." >&2
    exit 1
fi

exec python -m pytest "$here" "$@"
