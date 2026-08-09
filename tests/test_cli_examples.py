import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLES = [
    ("examples/test4_1/true.bed", "examples/test4_1/reconstruct.bed", "examples/test4_1/output"),
    ("examples/ecdna1/true.bed", "examples/ecdna1/reconstructed.bed", "examples/ecdna1/output"),

    ("examples/ecdna2/true_format.bed", "examples/ecdna2/reconstructed_format.bed", "examples/ecdna2/output"),
    ("examples/ecdna2/true_format.bed", "examples/ecdna2/reconstructed_format1.bed", "examples/ecdna2/output"),
    ("examples/ecdna2/true_format.bed", "examples/ecdna2/reconstructed_format2.bed", "examples/ecdna2/output"),

    ("examples/ecdna3/true_format.bed", "examples/ecdna3/reconstructed_format.bed", "examples/ecdna3/output"),
    ("examples/ecdna3/true_format.bed", "examples/ecdna3/reconstructed_format1.bed", "examples/ecdna3/output"),
    ("examples/ecdna3/true_format.bed", "examples/ecdna3/reconstructed_format2.bed", "examples/ecdna3/output"),

    ("examples/ecdna4/true_format.bed", "examples/ecdna4/reconstructed_format.bed", "examples/ecdna4/output"),
    ("examples/ecdna5/true_format.bed", "examples/ecdna5/reconstructed_format.bed", "examples/ecdna5/output"),
    ("examples/ecdna6/true_format.bed", "examples/ecdna6/reconstructed_format.bed", "examples/ecdna6/output"),
    ("examples/ecdna7/true_format.bed", "examples/ecdna7/reconstructed_format.bed", "examples/ecdna7/output"),
    ("examples/ecdna8/true_format.bed", "examples/ecdna8/reconstructed_format.bed", "examples/ecdna8/output"),
    ("examples/ecdna9/true_format.bed", "examples/ecdna9/reconstructed_format.bed", "examples/ecdna9/output"),
    ("examples/ecdna10/true_format.bed", "examples/ecdna10/reconstructed_format.bed", "examples/ecdna10/output"),

    ("examples/ecdna11/true_format.bed", "examples/ecdna11/reconstructed_format.bed", "examples/ecdna11/output"),
    ("examples/ecdna11/true_format.bed", "examples/ecdna11/reconstructed_format2.bed", "examples/ecdna11/output"),

    ("examples/ecdna12/true_format.bed", "examples/ecdna12/reconstructed_format.bed", "examples/ecdna12/output"),
    ("examples/ecdna13/true_format.bed", "examples/ecdna13/reconstructed_format.bed", "examples/ecdna13/output"),

    ("examples/ecdna13_reshuffled/true_format.bed", "examples/ecdna13_reshuffled/reconstructed_format.bed", "examples/ecdna13_reshuffled/output"),
    ("examples/ecdna13_reshuffled/true_format_reverse.bed", "examples/ecdna13_reshuffled/reconstructed_format.bed", "examples/ecdna13_reshuffled/output"),
    ("examples/ecdna13_reshuffled/true_format_acyclic.bed", "examples/ecdna13_reshuffled/reconstructed_format.bed", "examples/ecdna13_reshuffled/output"),
    ("examples/ecdna13_reshuffled/true_format_rotate.bed", "examples/ecdna13_reshuffled/reconstructed_format.bed", "examples/ecdna13_reshuffled/output"),
    ("examples/ecdna13_reshuffled/true_format_reverse_wrong.bed", "examples/ecdna13_reshuffled/reconstructed_format.bed", "examples/ecdna13_reshuffled/output"),

    ("examples/ecdna15/true.bed", "examples/ecdna15/reconstruct.bed", "examples/ecdna15/output"),
]


@pytest.mark.parametrize("a_b_c", EXAMPLES)
def test_cli_generates_report(tmp_path, a_b_c):
    a, b, c = a_b_c
    a = Path(a)
    b = Path(b)
    outdir = Path(c)
    
    assert a.exists(), f"example file {a} missing"
    assert b.exists(), f"example file {b} missing"

    outdir.mkdir(parents=True, exist_ok=True)
    print(sys.executable)
    cmd = [sys.executable, "-m", "eccomparator.main", "-a", str(a), "-b", str(b), "-d", str(outdir), "--report","--plot"]
    
    res = subprocess.run(cmd, check=False, capture_output=True, text=True)

    # CLI should exit without crashing
    assert res.returncode == 0, f"CLI failed: {res.stderr}\n{res.stdout}"

    report = outdir / "report.html"
    assert report.exists(), f"report not created at {report}"
