#!/usr/bin/env python3
"""Verify the qualified-build DMC energy sentinels and provenance."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "raw"
ENERGY_RE = re.compile(r"ENERGY\| Total FORCE_EVAL.*?([-+0-9.Ee]+)\s*$")
REFERENCE_BINARY = "b0dacc7dea4035ea5fb817eb1054f2b288016bfb63c9a96bceca878a44524c2f"
QUALIFIED_BINARY = "ab69000f995fbec138fd988f38b3cfe6ebde2dc58b7102e56f4bc52b3687e536"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def energy(run: Path) -> float:
    output = run / "cp2k.out"
    text = output.read_text()
    assert "PROGRAM ENDED AT" in text, f"incomplete output: {output}"
    assert (run / "exit_status").read_text().strip() == "0", run
    values = [float(match.group(1)) for line in text.splitlines()
              if (match := ENERGY_RE.search(line))]
    assert values, f"missing energy: {output}"
    return values[-1]


def recorded_digest(path: Path) -> str:
    return path.read_text().split()[0]


def check_affinity(run: Path) -> None:
    text = (run / "affinity_preexec.txt").read_text()
    header = text.splitlines()[0]
    fields = dict(item.split("=", 1) for item in header.split()
                  if "=" in item)
    assert fields["expected_cpu"] == fields["allowed"], (run, header)
    assert "," not in fields["allowed"] and "-" not in fields["allowed"], header
    assert f"Cpus_allowed_list:\t{fields['allowed']}" in text, (run, text)


def main() -> None:
    reference = {
        "k222_Ih": energy(RAW / "reference-k222-Ih"),
        "k222_VII": energy(RAW / "reference-k222-VII"),
        "k333_VII": energy(RAW / "baseline-k333-VII"),
    }
    candidate = {
        "k222_Ih": energy(RAW / "qualified-k222-Ih"),
        "k222_VII": energy(RAW / "qualified-k222-VII"),
        "k333_VII": energy(RAW / "qualified-k333-VII"),
    }
    deltas = {key: candidate[key] - reference[key] for key in reference}
    maximum = max(map(abs, deltas.values()))
    assert maximum <= 1.0e-12, maximum

    archived = json.loads((RAW / "comparison.json").read_text())
    assert archived["status"] == "PASS"
    assert archived["reference_hartree"] == reference
    assert archived["qualified_hartree"] == candidate
    assert archived["qualified_minus_reference_hartree"] == deltas
    assert archived["maximum_absolute_delta_hartree"] == maximum

    run = RAW / "baseline-k333-VII"
    assert recorded_digest(run / "binary.sha256") == REFERENCE_BINARY
    check_affinity(run)
    for name in ("qualified-k222-Ih", "qualified-k222-VII", "qualified-k333-VII"):
        run = RAW / name
        assert recorded_digest(run / "binary.sha256") == QUALIFIED_BINARY
        check_affinity(run)

    input_map = {
        "qualified-k222-Ih": ROOT / "inputs/k222-Ih.inp",
        "qualified-k222-VII": ROOT / "inputs/k222-VII.inp",
        "baseline-k333-VII": ROOT / "inputs/k333-VII.inp",
        "qualified-k333-VII": ROOT / "inputs/k333-VII.inp",
    }
    for name, input_path in input_map.items():
        assert recorded_digest(RAW / name / "input.sha256") == digest(input_path)

    print(json.dumps({
        "qualified_minus_reference_hartree": deltas,
        "maximum_absolute_delta_hartree": maximum,
        "status": "PASS",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
