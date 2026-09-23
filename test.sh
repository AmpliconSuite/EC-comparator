for py in 3.10 3.11 3.12; do
    echo "======================================"
    echo "Testing Python $py"
    echo "======================================"

    conda create -n "eccomp-py${py//./}" python="$py" -y

    conda run -n "eccomp-py${py//./}" \
        python -m pip install -r requirements.txt

    conda run -n "eccomp-py${py//./}" \
        python -m pip install -e .

    conda run -n "eccomp-py${py//./}" \
        python -m pip install pytest xhtml2pdf

    conda run -n "eccomp-py${py//./}" \
        pytest -q tests/test_cli_examples.py
done
