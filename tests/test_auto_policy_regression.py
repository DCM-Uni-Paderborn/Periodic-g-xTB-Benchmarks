import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = (
    ROOT
    / "validation"
    / "accelerated_exchange"
    / "auto_policy_regression_20260722"
)


def test_exact_build_automatic_policy_archive():
    result = subprocess.run(
        [sys.executable, str(ARCHIVE / "verify.py")],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip() == (
        "automatic-policy exact-build regression: PASS (196/196)"
    )
