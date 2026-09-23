### Build and install (for developers)

Install the packaging tools:

```bash
python -m pip install  --upgrade installer toml build twine
```

Go to the repository root:

```bash
cd EC-comparator
```

Remove previous builds:

```bash
rm -rf build dist *.egg-info eccomparator.egg-info
```

Build the package:

```bash
python -m build
```

The `dist/` directory should contain:

```text
dist/
├── ec_comparator-0.1.0-py3-none-any.whl
└── ec_comparator-0.1.0.tar.gz
```

Check the built packages:

```bash
python -m twine check dist/*
```

Both the wheel and source distribution should report:

```text
PASSED
```

### Test the package locally

Create a clean test environment:

```bash
conda create -n eccomparator-test -c conda-forge -c bioconda python=3.10 bedtools=2.31.1
conda activate eccomparator-test
```

Install the wheel:

```bash
python -m pip install dist/ec_comparator-0.1.0-py3-none-any.whl
```

Check the installed CLI:

```bash
EC-comparator --version
```

It should output:

```text
EC-comparator 0.1.0
```

### Upload to TestPyPI

Upload the package:

```bash
python -m twine upload --repository testpypi dist/*
```

Use your TestPyPI credentials/API token when prompted:

https://test.pypi.org/

To test installation from TestPyPI, first remove the locally installed package:

```bash
python -m pip uninstall EC-comparator -y
```

Install from TestPyPI while using the regular PyPI index for dependencies:

```bash
python -m pip install \
    --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ \
    EC-comparator
```

Verify the installation:

```bash
EC-comparator --version
```

## Upload to PyPI

After successfully testing TestPyPI:

```bash
python -m twine upload dist/*
```

Test the official PyPI installation:

```bash
python -m pip uninstall EC-comparator -y
python -m pip install EC-comparator

EC-comparator --version
```



### Running tests (for developers)

The CLI example tests generate PDF reports and require `xhtml2pdf`.

If the tests fail with:

```text
ModuleNotFoundError: No module named 'xhtml2pdf'
```

install it in the active environment:

```bash
python -m pip install xhtml2pdf
```

Run the CLI example tests:

```bash
pytest -q tests/test_cli_examples.py
```

If a test fails while generating reports, confirm that `xhtml2pdf` is installed in the same Python environment used to run `pytest`.

