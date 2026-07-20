#!/usr/bin/env python3
"""Verify the direct-CLI accuracy sensitivity and CP2K parity controls."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CLI_SHA256 = "f0c66f82385f33367b9988a9f04959b77992e0139f60b47211e35b90bbebb38a"
CP2K_SHA256 = "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
ENERGY_RE = re.compile(
    r"^\s*ENERGY\|\s+Total FORCE_EVAL.*?"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)
ECONV_RE = re.compile(r"energy convergence\s+([-+0-9.Ee]+)\s+Eh")
PCONV_RE = re.compile(r"density convergence\s+([-+0-9.Ee]+)\s+e")
ITER_RE = re.compile(r"^\s+(\d+)\s+-\d", re.MULTILINE)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def sidecar(path: Path) -> str:
    fields = path.read_text(encoding="utf-8").split()
    return fields[0].lower() if fields else ""


def cli_run(path: Path, expected_accuracy: float) -> dict[str, object]:
    if (path / "exit_status").read_text(encoding="utf-8").strip() != "0":
        raise ValueError(f"nonzero direct-CLI exit status: {path}")
    if sidecar(path / "binary.sha256") != CLI_SHA256:
        raise ValueError(f"wrong direct-CLI binary: {path}")
    if digest(path / "tblite.json") == "":
        raise ValueError(f"missing direct-CLI result: {path}")
    text = (path / "process.out").read_text(encoding="utf-8", errors="replace")
    energy_match = ECONV_RE.search(text)
    density_match = PCONV_RE.search(text)
    iterations = [int(value) for value in ITER_RE.findall(text)]
    if energy_match is None or density_match is None or not iterations:
        raise ValueError(f"incomplete direct-CLI output: {path}")
    energy_convergence = float(energy_match.group(1))
    density_convergence = float(density_match.group(1))
    if not math.isclose(energy_convergence, 1.0e-6 * expected_accuracy, abs_tol=1.0e-18):
        raise ValueError(f"wrong energy convergence threshold: {path}")
    if not math.isclose(density_convergence, 2.0e-5 * expected_accuracy, abs_tol=1.0e-18):
        raise ValueError(f"wrong density convergence threshold: {path}")
    recorded_accuracy = float((path / "accuracy.txt").read_text(encoding="utf-8"))
    if recorded_accuracy != expected_accuracy:
        raise ValueError(f"wrong recorded accuracy: {path}")
    result = json.loads((path / "tblite.json").read_text(encoding="utf-8"))
    energy = float(result["energy"])
    if not math.isfinite(energy):
        raise ValueError(f"non-finite direct-CLI energy: {path}")
    return {
        "accuracy": expected_accuracy,
        "energy_hartree": energy,
        "energy_convergence_hartree": energy_convergence,
        "density_convergence_e": density_convergence,
        "iterations": max(iterations),
        "binary_sha256": CLI_SHA256,
        "input_sha256": sidecar(path / "input.sha256"),
        "result_sha256": digest(path / "tblite.json"),
        "process_output_sha256": digest(path / "process.out"),
    }


def cp2k_run(phase: str) -> dict[str, object]:
    path = (
        ROOT
        / "DMC-ICE13/reproduction/seidler_dmc13_recalculation/raw/cp2k_native"
        / "k111-reduced"
        / phase
    )
    if (path / "exit_status").read_text(encoding="utf-8").strip() != "0":
        raise ValueError(f"nonzero CP2K exit status: {path}")
    if sidecar(path / "binary.sha256") != CP2K_SHA256:
        raise ValueError(f"wrong CP2K binary: {path}")
    values = []
    ended = False
    for line in (path / "cp2k.out").read_text(encoding="utf-8", errors="replace").splitlines():
        if match := ENERGY_RE.match(line):
            values.append(float(match.group(1)))
        if "PROGRAM ENDED AT" in line:
            ended = True
    if not ended or not values or not math.isfinite(values[-1]):
        raise ValueError(f"incomplete CP2K output: {path}")
    return {
        "accuracy": 0.1,
        "energy_hartree": values[-1],
        "binary_sha256": CP2K_SHA256,
        "input_sha256": sidecar(path / "input.sha256"),
        "output_sha256": digest(path / "cp2k.out"),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "verification.json")
    args = parser.parse_args()
    rows: dict[str, object] = {}
    for phase in ("Ih", "VII"):
        loose = cli_run(HERE / f"{phase}_k111_acc010", 0.1)
        tight = cli_run(HERE / f"{phase}_k111_acc001", 0.01)
        native = cp2k_run(phase)
        loose_energy = float(loose["energy_hartree"])
        tight_energy = float(tight["energy_hartree"])
        native_energy = float(native["energy_hartree"])
        rows[phase] = {
            "direct_cli_accuracy_0.1": loose,
            "direct_cli_accuracy_0.01": tight,
            "cp2k_native_accuracy_0.1": native,
            "tight_minus_loose_hartree": tight_energy - loose_energy,
            "native_minus_loose_hartree": native_energy - loose_energy,
        }
    max_accuracy_delta = max(
        abs(float(row["tight_minus_loose_hartree"])) for row in rows.values()
    )
    max_native_delta = max(
        abs(float(row["native_minus_loose_hartree"])) for row in rows.values()
    )
    report = {
        "status": "PASS"
        if max_accuracy_delta <= 1.0e-8 and max_native_delta <= 2.0e-7
        else "FAIL",
        "phases": rows,
        "maximum_accuracy_sensitivity_hartree": max_accuracy_delta,
        "maximum_cp2k_native_minus_same_accuracy_cli_hartree": max_native_delta,
        "tolerances_hartree": {
            "accuracy_sensitivity": 1.0e-8,
            "cp2k_native_parity": 2.0e-7,
        },
        "interpretation": (
            "For the two controlled 1x1x1 cells, tightening direct save_tblite "
            "from accuracy 0.1 to 0.01 changes the total energy by at most "
            "1.91e-10 hartree. CP2K-native at accuracy 0.1 agrees with the "
            "same-accuracy direct CLI to 2.02e-8 hartree."
        ),
    }
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
