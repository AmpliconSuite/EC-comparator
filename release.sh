#!/usr/bin/env bash
set -euo pipefail

PACKAGE="EC-comparator"

python -m pip install --upgrade build twine pytest setuptools wheel

echo "========================================"
echo "1. Running tests"
echo "========================================"

python -m pytest -q tests/test_cli_examples.py

echo
echo "Tests passed."

echo "========================================"
echo "2. Cleaning previous builds"
echo "========================================"

rm -rf build dist *.egg-info eccomparator.egg-info

echo "========================================"
echo "3. Building package"
echo "========================================"

python -m build

echo "========================================"
echo "4. Checking package"
echo "========================================"

python -m twine check dist/*

echo "========================================"
echo "5. Uploading to TestPyPI"
echo "========================================"

python -m twine upload \
    --repository testpypi \
    dist/*

echo "========================================"
echo "6. Installing from TestPyPI"
echo "========================================"

python -m pip uninstall -y "${PACKAGE}" || true

python -m pip install \
    --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ \
    "${PACKAGE}"

echo "========================================"
echo "7. Testing TestPyPI installation"
echo "========================================"

EC-comparator --version
EC-comparator --help

echo "TestPyPI installation works."

echo "========================================"
echo "8. Uploading to PyPI"
echo "========================================"

python -m twine upload dist/*

echo "========================================"
echo "9. Installing from PyPI"
echo "========================================"

python -m pip uninstall -y "${PACKAGE}" || true

python -m pip install "${PACKAGE}"

echo "========================================"
echo "10. Testing PyPI installation"
echo "========================================"

EC-comparator --version
EC-comparator --help

echo
echo "========================================"
echo "Release completed successfully."
echo "========================================"